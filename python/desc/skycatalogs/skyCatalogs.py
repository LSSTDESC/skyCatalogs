import os
import sys
import re
import yaml
from collections import namedtuple
import healpy
import numpy as np
import numpy.ma as ma
import pyarrow.parquet as pq
from astropy import units
from desc.skycatalogs.objects import *
from desc.skycatalogs.readers import *
from desc.skycatalogs.readers import ParquetReader
from desc.skycatalogs.utils.sed_utils import MagNorm, create_cosmology

__all__ = ['SkyCatalog', 'open_catalog', 'Box', 'Disk']


Box = namedtuple('Box', ['ra_min', 'ra_max', 'dec_min', 'dec_max'])

# radius is measured in arcseconds
Disk = namedtuple('Disk', ['ra', 'dec', 'radius_as'])

_aDisk = Disk(1.0, 1.0, 1.0)
_aBox = Box(-1.0, 1.0, -2.0, 2.0)

# This function should maybe be moved to utils
def _get_intersecting_hps(hp_ordering, nside, region):
    '''
    Given healpixel structure defined by hp_ordering and nside, find
    all healpixels which instersect region, where region may be either
    a box or a disk.
    Return as some kind of iterable
    Note it's possible extra hps which don't actually intersect the region
    will be returned
    '''
    # First convert region description to an array (4,3) with (x,y,z) coords
    # for each vertex
    if type(region) == type(_aBox):
        vec = healpy.pixelfunc.ang2vec([region.ra_min, region.ra_max,
                                        region.ra_max, region.ra_min],
                                       [region.dec_min, region.dec_min,
                                        region.dec_max, region.dec_max],
                                       lonlat=True)

        return healpy.query_polygon(nside, vec, inclusive=True, nest=False)
    if type(region) == type(_aDisk):
        # Convert inputs to the types query_disk expects
        center = healpy.pixelfunc.ang2vec(region.ra, region.dec,
                                          lonlat=True)
        radius_rad = (region.radius_as * units.arcsec).to_value('radian')

        pixels = healpy.query_disc(nside, center, radius_rad, inclusive=True,
                                   nest=False)
        return pixels
def _compute_mask(region, ra, dec):
    '''
    Compute mask according to region for provided data
    Parameters
    ----------
    region         Supported shape or None
    ra,dec         Coordinates for data to be masked
    Returns
    -------
    mask of elements to be omitted

    '''
    mask = None
    if type(region) == type(_aBox):
        mask = np.logical_or((ra < region.ra_min),
                             (ra > region.ra_max))
        mask = np.logical_or(mask, (dec < region.dec_min))
        mask = np.logical_or(mask, (dec > region.dec_max))
    if type(region) == type(_aDisk):
        # Change positions to 3d vectors to measure distance
        p_vec = healpy.pixelfunc.ang2vec(ra, dec, lonlat=True)

        c_vec = healpy.pixelfunc.ang2vec(region.ra,
                                         region.dec,
                                         lonlat=True)
        radius_rad = (region.radius_as * units.arcsec).to_value('radian')


        # Rather than comparing arcs, it is equivalent to compare chords
        # (or square of chord length)
        diff = p_vec - c_vec
        obj_chord_sq = np.sum(np.square(p_vec - c_vec),1)

        # This is to be compared to square of chord for angle a corresponding
        # to disk radius.  That's 4(sin(a/2)^2)
        rad_chord_sq = 4 * np.square(np.sin(0.5 * radius_rad) )
        mask = obj_chord_sq > rad_chord_sq

    return mask

class SkyCatalog(object):
    '''
    A base class with derived classes for galaxies, static (w.r.t. coordinates)
    point sources, SSOs

    '''
    def __init__(self, config, mp=False, verbose=False):
        '''
        Parameters
        ----------
        config:  dict.  Typically the result of loading a yaml file
        mp:      boolean  Default False. Set True to enable multiprocessing.
        '''
        self._config = config
        self._mp = mp
        # There may be more to do at this point but not too much.
        # In particular, don't read in anything from data files

        self.verbose = verbose
        self._validate_config()

        # Outer dict: hpid for key. Value is another dict
        #    with keys 'files', 'object_types', each with value another dict
        #    for 'files', map filepath to handle (initially None)
        #    for 'object_types', map object type to filepath
        self._hp_info = dict()
        hps = self._find_all_hps()

        cosmology = create_cosmology(config['Cosmology'])
        self._magnorm_f = MagNorm(cosmology)

    @property
    def mag_norm_f(self):
        '''
        Return function object used to calculate mag_norm
        '''
        return self._magnorm_f

    @property
    def raw_config(self):
        '''
        Return config, typically uploaded from yaml.
        '''
        return self._config

    # One might check that the config is complete and well-formed
    #  - for example require certain keys, such as catalog_name,
    #  data_file_type, area_parition, root_directory, object_types -
    # to exist, and check that that the data directory (value of
    # root_directory) exists.
    def _validate_config(self):
        pass

    def _find_all_hps(self):
        '''
        For each healpix with files matching pattern in the directory,
        update _hp_info as needed to keep track of all files for that healpix
        and the object types included in those files.   Also..

        Returns
        -------
        Set of healpix pixels with at least one file in the directory

        '''

        cat_dir = self._config['root_directory']

        # If major organization is by healpix, healpix # will be in
        # subdirectory name.  Otherwise there may be no subdirectories
        # or data may be organized by component type.
        # Here only handle case where data files are directly in root dir
        files = os.listdir(cat_dir)
        o_types = self._config['object_types']

        hp_set = set()
        for f in files:
            for ot in o_types:
                if 'file_template' in o_types[ot]:
                    m = re.match(o_types[ot]['file_template'], f)
                    if m:
                        hp = int(m['healpix'])
                        hp_set.add(hp)

                        if hp not in self._hp_info:
                            self._hp_info[hp] = {'files' : {f : None},
                                                 'object_types' : {ot : f}}
                        else:
                            this_hp = self._hp_info[hp]
                            if f not in this_hp['files'] :
                                this_hp['files'][f] = None
                            this_hp['object_types'][ot] = f
        return hp_set

    def get_hps_by_region(self, region):
        '''
        Parameters
        ----------
        Region can be a box (named 4-tuple (min-ra, max-ra, min-dec, max-dec))
        or a circle (named 3-tuple (ra, dec, radius))
        Catalog area partition must be by healpix

        Returns
        -------
        Set of healpixels intersecting the region
        '''
        # If area_partition doesn't use healpix, raise exception

        return _get_intersecting_hps(
            self._config['area_partition']['ordering'],
            self._config['area_partition']['nside'],
            region)
            #.intersection(self._hps.keys())

    def get_object_type_names(self):
        '''
        Returns
        -------
        All object type names in the catalog's config
        '''
        return set(self._config['object_types'].keys())

    # Add more functions to return parts of config of possible interest
    # to user

    def get_objects_by_region(self, region, obj_type_set=None, datetime=None):
        '''
        Parameters
        ----------
        region         region is a named tuple.  May be box or circle
        obj_type_set   Return only these objects. Defaults to all available
        datetime       Python datetime object. Ignored except for SSO

        Returns
        -------
        ObjectList containing sky objects visible in the region
        [at the specified time]
        '''
        # Take intersection of obj_type_list and available object types
        # Determine healpix intersecting the region

        if self.verbose:
            print("Region ", region)
            print("obj_type_set ", obj_type_set)
        if self._config['area_partition']['type'] == 'healpix':
            hps = self.get_hps_by_region(region)

        # otherwise raise a not-supported exception

        object_list = ObjectList()
        if obj_type_set is None:
            obj_types = self.get_object_type_names()
        else:
            obj_types = self.get_object_type_names().intersection(obj_type_set)
        for hp in hps:
            # Maybe have a multiprocessing switch? Run-time option when
            # catalog is opened?
            c = self.get_objects_by_hp(hp, region, obj_types, datetime)
            if (len(c)) > 0:
                object_list.append_object_list(c)

        return object_list

    def get_object_iterator_by_region(self, region=None, obj_type_set=None,
                                      max_chunk=None, datetime=None):
        '''
        Parameters
        ----------
        region         Either a box or a circle (each represented as
                       named tuple)
        obj_type_set   Return only these objects. Defaults to all available
        max_chunk      If specified, iterator will return no more than this
                       number of objections per iteration
        datetime       Python datetime object. Ignored for all but SSO

        Returns
        -------
        An iterator
        '''
        raise NotImplementedError('get_object_iterator_by_region not implemented yet. See get_objects_by_region instead for now')

    def get_objects_by_hp(self, hp, region=None, obj_type_set=None,
                          datetime=None):
        '''
        Find all object
        # Find the right Sky Catalog file or files (depends on obj_type_set)
        # Get file handle(s) and store if we don't already have it (them)

        Returns
        -------
        ObjectList containing sky objects in the region and the hp
        '''
        object_list = ObjectList()

        G_COLUMNS = ['galaxy_id', 'ra', 'dec']
        PS_COLUMNS = ['object_type', 'id', 'ra', 'dec']
        if self.verbose:
            print('Working on healpix pixel ', hp)

        obj_types = obj_type_set
        if obj_types is None:
            obj_types = self._config['object_types'].keys()
        else:
            parents = set()
            for ot in obj_types:
                if 'parent' in self._config['object_types'][ot]:
                    parents.add(self._config['object_types'][ot]['parent'])
            obj_types = obj_types.union(parents)

        # Associate object types with readers.  May be > one type per reader
        rdr_ot = dict()
        root_dir = self._config['root_directory']
        for ot in obj_types:
            if 'file_template' in self._config['object_types'][ot]:
                f = self._hp_info[hp]['object_types'][ot]
            elif 'parent' in self._config['object_types'][ot]:
                f = self._hp_info[hp]['object_types'][self._config['object_types'][ot]['parent']]
            if f not in self._hp_info[hp]:
                the_reader = parquet_reader.ParquetReader(os.path.join(root_dir,f), mask=None)
                self._hp_info[hp][f] = the_reader
            the_reader = self._hp_info[hp][f]
            if the_reader in rdr_ot:
                rdr_ot[the_reader].add(ot)
            else:
                rdr_ot[the_reader] = set([ot])

        # Now get minimal columns for objects using the readers
        for rdr in rdr_ot:
            if 'galaxy' in rdr_ot[rdr]:
                arrow_t = rdr.read_columns(G_COLUMNS, None) # or read_row_group
                # Make a boolean array, value set to 1 for objects
                # outside the region
                if region is not None:
                    # This belongs in a separate routine
                    mask = _compute_mask(region, arrow_t['ra'], arrow_t['dec'])

                else:
                    mask = None

                if mask is not None:
                    masked_ra = ma.array(arrow_t['ra'], mask=mask)
                    if self.verbose:
                        print("Masked array size: ", masked_ra.size)
                        print("Masked array compressed size: ", masked_ra.compressed().size)
                    ra_compress = masked_ra.compressed()
                    if ra_compress.size > 0:
                        dec_compress = ma.array(arrow_t['dec'], mask=mask).compressed()
                        id_compress = ma.array(arrow_t['galaxy_id'], mask=mask).compressed()
                    else:
                        continue
                else:
                    ra_compress = arrow_t['ra']
                    dec_compress = arrow_t['dec']
                    id_compress = arrow_t['galaxy_id']

                new_collection = ObjectCollection(ra_compress,
                                                  dec_compress,
                                                  id_compress,
                                                  'galaxy',
                                                  hp,
                                                  self,
                                                  region=region,
                                                  mask=mask,
                                                  reader=rdr)
                object_list.append_collection(new_collection)

                # Now do the same for point sources
            if 'star' in rdr_ot[rdr]:
                arrow_t = rdr.read_columns(PS_COLUMNS, None) # or read_row_group
                # Make a boolean array, value set to 1 for objects
                # outside the region
                if region is not None:
                    # This belongs in a separate routine
                    mask = _compute_mask(region, arrow_t['ra'], arrow_t['dec'])
                else:
                    mask = None

                if mask is not None:
                    masked_ra = ma.array(arrow_t['ra'], mask=mask)
                    ra_compress = masked_ra.compressed()
                    if ra_compress.size > 0:
                        dec_compress = ma.array(arrow_t['dec'], mask=mask).compressed()
                        id_compress = ma.array(arrow_t['id'], mask=mask).compressed()
                    else:
                        continue
                else:
                    ra_compress = arrow_t['ra']
                    dec_compress = arrow_t['dec']
                    id_compress = arrow_t['id']

                new_collection = ObjectCollection(ra_compress,
                                                  dec_compress,
                                                  id_compress,
                                                  'star',
                                                  hp,
                                                  self,
                                                  region=region,
                                                  mask=mask,
                                                  reader=rdr)
                object_list.append_collection(new_collection)



        return object_list

    # For generator version, do this a row group at a time
    #    but if region cut leaves too small a list, read more rowgroups
    #    to achieve a reasonable size list (or exhaust the file)
    def get_object_iterator_by_hp(self, hp, obj_type_set=None,
                                  max_chunk=None, datetime=None):
        '''
        Parameters
        ----------
        datetime       Python datetime object.
        hp             A healpix id
        obj_type_set   Return only these objects. Defaults to all available
        max_chunk      If specified, iterator will return no more than this
                       number of objections per iteration
        Returns
        -------
        An iterator
        '''
        pass


def open_catalog(config_file, mp=False):
    '''
    Parameters
    ----------
    yaml file containing config

    Returns
    -------
    SkyCatalog
    '''
    with open(config_file) as f:
        return SkyCatalog(yaml.safe_load(f), mp)

if __name__ == '__main__':
    ##cfg_file_name = 'to_translate.yaml'
    cfg_file_name = 'latest.yaml'

    if len(sys.argv) > 1:
        cfg_file_name = sys.argv[1]
    cfg_file = os.path.join('/global/homes/j/jrbogart/desc_git/skyCatalogs/cfg',
                            cfg_file_name)
    ##cfg_file = '/global/homes/j/jrbogart/desc_git/skyCatalogs/cfg/to_translate.yaml'

    # For tract 3828
    #   55.73604 < ra < 57.563452
    #  -37.19001 < dec < -35.702481


    cat = open_catalog(cfg_file)
    hps = cat._find_all_hps()
    print('Found {} healpix pixels '.format(len(hps)))
    for h in hps: print(h)

    ra_min_tract = 55.736
    ra_max_tract = 57.564
    dec_min_tract = -37.190
    dec_max_tract = -35.702
    ##ra_min_small = 56.0
    ##ra_max_small = 56.2
    ra_min_small = 55.9
    ra_max_small = 56.1
    dec_min_small = -36.2
    dec_max_small = -36.0

    sed_fmt = 'lambda: {:.1f}  f_lambda: {:g}'

    rgn = Box(ra_min_small, ra_max_small, dec_min_small, dec_max_small)

    intersect_hps = _get_intersecting_hps('ring', 32, rgn)

    print("For region ", rgn)
    print("intersecting pixels are ", intersect_hps)

    print('Invoke get_objects_by_region with box region')
    object_list = cat.get_objects_by_region(rgn)
                                            ##### temporary obj_type_set={'galaxy', 'star'} )
    #                                        obj_type_set=set(['galaxy']) )
    # Try out get_objects_by_hp with no region
    #colls = cat.get_objects_by_hp(9812, None, set(['galaxy']) )

    print('Number of collections returned:  ', object_list.collection_count)

    colls = object_list.get_collections()
    for c in colls:
        n_obj = len(c)
        print("For hpid ", c.get_partition_id(), "found ", n_obj, " objects")
        print("First object: ")
        print(c[0], '\nid=', c[0].id, ' ra=', c[0].ra, ' dec=', c[0].dec,
              ' belongs_index=', c[0]._belongs_index,
              ' object_type: ', c[0].object_type )

        print("Slice [1:3]")
        slice13 = c[1:3]
        for o in slice13:
            print('id=',o.id, ' ra=',o.ra, ' dec=',o.dec, ' belongs_index=',
                  o._belongs_index,  ' object_type: ', o.object_type)
            print(o.object_type)
            if o.object_type == 'star':
                print(o.get_instcat_entry())
                (lmbda, f_lambda, magnorm) = o.get_sed(resolution=1.0)
                print('For star magnorm: ', magnorm)
                if magnorm < 1000:
                    print('Length of sed: ', len(lmbda))
                    for i in range(10):
                        print(sed_fmt.format(lmbda[i], f_lambda[i]))
                    mid = int(len(lmbda)/2)
                    for i in range(mid-5, mid+5):
                        print("ix=",i,"  ", sed_fmt.format(lmbda[i], f_lambda[i]))

            else:
                for cmp in ['disk', 'bulge', 'knots']:
                    print(cmp)
                    if cmp in o.subcomponents:
                        print(o.get_instcat_entry(component=cmp))
                        (lmbda, f_lambda, magnorm) = o.get_sed(cmp)
                        print('magnorm: ', magnorm)
                        if magnorm < 1000:
                            for i in range(10):
                                print(sed_fmt.format(lmbda[i], f_lambda[i]))

        if n_obj > 200:
            print("Object 200")
            print(c[200], '\nid=', c[200].id, ' ra=', c[200].ra, ' dec=',
                  c[200].dec,
                  ' belongs_index=', c[200]._belongs_index)
        if n_obj > 163997:
            slice_late = c[163994:163997]
            print('\nobjects indexed 163994 through 163996')
            for o in slice_late:
                print('id=',o.id, ' ra=',o.ra, ' dec=',o.dec, ' belongs_index=',
                      o._belongs_index)

    print('Total object count: ', len(object_list))

    obj = object_list[0]
    print("Type of element in object_list:", type(obj))

    if object_list[0].object_type == 'galaxy':
        redshift0 = object_list[0].get_native_attribute('redshift')
        print('First redshift: ', redshift0)


    sum = 0
    for obj in object_list:
        sum = sum + 1
        if sum > 7216:
            print("Sum is now ", sum)
            print("obj id: ", obj.id)
        if sum > 7220:
            break

    print(f'Object list len:  {len(object_list)}')
    print(f'Objects found with "in":  {sum}')

    segment = object_list[2:5]
    print('Information for slice 2:5 ')
    for o in segment:
        print(f'object {o.id} of type {o.object_type} belongs to collection {o._belongs_to}')

    print('\nInformation for slice 285:300')
    segment = object_list[285:300]
    for o in segment:
        print(f'object {o.id} of type {o.object_type} belongs to collection {o._belongs_to}')

    #ixes = ([3,5,8],)
    ixes = (np.array([3,5,8, 300, 303]),)
    print(f'\nObjects with indexes {ixes[0]}')
    for o in object_list[ixes]:
        print(o.id)
    print(f'\nObjects in slice [3:9]')
    for o in object_list[3:9]:
        print(o.id)
    print(f'\nObjects in slice [300:304]')
    for o in object_list[300:304]:
        print(o.id)
