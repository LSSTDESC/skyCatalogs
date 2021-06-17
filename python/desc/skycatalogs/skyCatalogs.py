import os
import re
import yaml
from collections import namedtuple
import healpy
import numpy as np
import numpy.ma as ma
import pyarrow.parquet as pq
from astropy import units
from objects import *
from readers import *
from readers import ParquetReader

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
    all healpixels which instersect region, defined by min/max ra and dec
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
        center = healpy.pixelfunc.ang2vec([region.ra], [region.dec],
                                          lonlat=True)
        radius_rad = (region.radius_as * units.arcsec).to_value('radian')

        return healpy.query_disc(nside, center, radius_rad, inclusive=True,
                                 nest=False)

class SkyCatalog(object):
    '''
    A base class with derived classes for galaxies, static (w.r.t. coordinates)
    point sources, SSOs

    '''
    def __init__(self, config, mp=False):
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
        # One might check that the config is complete and well-formed
        #  - for example require certain keys, such as catalog_name,
        #  data_file_type, area_parition, root_directory, object_types -
        # to exist, and check that that the data directory (value of
        # root_directory) exists.
        # create an empty dict for
        # per-HEALpix pixel information and so forth.

        self._validate_config()

        # Outer dict: hpid for key. Value is another dict
        #    with keys 'files', 'object_types', each with value another dict
        #    for 'files', map filepath to handle (initially None)
        #    for 'object_types', map object type to filepath
        self._hp_info = dict()
        self._find_all_hps()

    def _validate_config(self):
        pass

    def _find_all_hps(self):
        # for each healpix with files matching pattern in the directory,
        # update _hp_info
        #   self._hp_info[hpid]['files'] =
        #     {'relpath' : the_filename, 'handle' : None}    and,
        #  for each object_type represented in the file,
        #   self._hp_info[hpid]['object_types'][ot] = the_filename
        # When file is open, set handle to the Parquet file object
        # (or perhaps something else if underlying format is not Parquet)

        cat_dir = self._config['root_directory']

        # If major organization is by healpix, healpix # will be in
        # subdirectory name.  Otherwise there may be no subdirectories
        # or data may be organized by component type.
        # Here only handle case where data files are directly in root dir
        files = os.listdir(cat_dir)
        o_types = self._config['object_types']
        #####o_types_enum = set([OBJECT_TYPES[t] for t in o_types])
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
        Region can be a box (named 4-tuple (min-ra, max-ra, min-dec, max-dec))
        or a circle (named 3-tuple (ra, dec, radius))
        Catalog area partition must be by healpix
        '''
        # If area_partition doesn't use healpix, raise exception

        return _get_intersecting_hps(
            self._config['area_partition']['ordering'],
            self._config['area_partition']['nside'],
            region)
            #.intersection(self._hps.keys())

    def get_object_type_names(self):
        return set(self._config['object_types'].keys())

    # Add more functions to return parts of config of possible interest
    # to user

    def get_objects_by_region(self, datetime, region, obj_type_set=None):
        '''
        Parameters
        ----------
        datetime       Python datetime object.
        region         region is a named tuple.  May be box or circle
        obj_type_set   Return only these objects. Defaults to all available

        Returns
        -------
        Collection of SkyObjects visible in the region at the specified time
        '''
        # Take intersection of obj_type_list and available object types
        # Determine healpix intersecting the region

        print("Region ", region)
        print("obj_type_set ", obj_type_set)
        if self._config['area_partition']['type'] == 'healpix':
            hps = self.get_hps_by_region(region)

        # otherwise raise a not-supported exception

        obj_colls = []
        if obj_type_set is None:
            obj_types = self.get_object_type_names()
        else:
            obj_types = self.get_object_type_names().intersection(obj_type_set)
        for hp in hps:
            # Maybe have a multiprocessing switch? Run-time option when
            # catalog is opened?
            c = self.get_objects_by_hp(datetime, hp, region, obj_types)
            if (len(c)) > 0:
                obj_colls = obj_colls + c

        return obj_colls

    def get_object_iterator_by_region(self, datetime, region=None,
                                      obj_type_set=None, max_chunk=None):
        '''
        Parameters
        ----------
        datetime       Python datetime object.
        region         min & max for ra and dec
        obj_type_set   Return only these objects. Defaults to all available
        max_chunk      If specified, iterator will return no more than this
                       number of objections per iteration
        Returns
        -------
        An iterator
        '''
        pass

    def get_objects_by_hp(self, datetime, hp, region=None, obj_type_set=None):
        # Find the right Sky Catalog file or files (depends on obj_type_set)
        # Get file handle(s) and store if we don't already have it (them)

        # Return list of object collections

        G_COLUMNS = ['galaxy_id', 'ra', 'dec']
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
                f = self._hp_info[hp]['object_types'][obj_types[ot]['parent']]
            if f not in self._hp_info[hp]:
                ##self._hp_info[hp][f] = pq.ParquetFile(os.path.join(root_dir,f))
                # or maybe something like
                the_reader = parquet_reader.ParquetReader(os.path.join(root_dir,f))
                self._hp_info[hp][f] = the_reader
            the_reader = self._hp_info[hp][f]
            if the_reader in rdr_ot:
                rdr_ot[the_reader].add(ot)
                #print("added object type for reader ", rdr_ot)
            else:
                rdr_ot[the_reader] = set([ot])
                #print("added reader ", the_reader, "and object type ", ot)

        # Now get minimal columns for objects using the readers
        obj_collect = []
        for rdr in rdr_ot:
            if 'galaxy' in rdr_ot[rdr]:
                arrow_t = rdr.read_columns(G_COLUMNS) # or read_row_group
                # Make a boolean array, value set to 1 for objects
                # outside the region
                if region is not None:
                    if type(region) == type(_aBox):
                        mask = np.logical_or((arrow_t['ra'] < region.ra_min),
                                             (arrow_t['ra'] > region.ra_max))
                        mask = np.logical_or(mask, (arrow_t['dec'] < region.dec_min))
                        mask = np.logical_or(mask, (arrow_t['dec'] > region.dec_max))
                    if type(region) == type(_aDisk):
                        # Change positions to 3d vectors to measure distance
                        p_vec = healpy.pixelfunc.ang2vec(arrow_t['ra'],
                                                         arrow_t['dec'],
                                                         lonlat=True)

                        c_vec = healpy.pixelfunc.ang1vec([region.ra],
                                                         [region.dec],
                                                         lonlat=True)[0]
                        # change disk radius to radians
                        radius_rad = (region.radius_as * units.arcsec).to_value('radian')
                        inners = [np.dot(pos, c_vec) for pos in p_vec]
                        mask = np.arccos(inners) > radius_rad

                else:
                    mask = None

                # Any future reads should compress output using the mask
                rdr.set_mask(mask)

                if mask is not None:
                    masked_ra = ma.array(arrow_t['ra'], mask=mask)
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

                obj_collect.append(BaseObjectCollection(ra_compress,
                                                        dec_compress,
                                                        id_compress,
                                                        'galaxy',
                                                        include_mask=mask,
                                                        hp_id=hp,
                                                        region=region,
                                                        reader=rdr))
        return obj_collect


        # For generator version, do this a row group at a time
        #    but if region cut leaves too small a list, read more rowgroups
        #    to achieve a reasonable size list (or exhaust the file)

    def get_object_iterator_by_hp(self, datetime, hp, obj_type_set=None,
                                  max_chunk=None):
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
    cfg_file = '/global/homes/j/jrbogart/Joanne_git/skyCatalogs/cfg/galaxy.yaml'

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
    ra_min_small = 56.0
    ra_max_small = 56.2
    dec_min_small = -36.2
    dec_max_small = -36.0

    rgn = Box(ra_min_small, ra_max_small, dec_min_small, dec_max_small)

    intersect_hps = _get_intersecting_hps('ring', 32, rgn)

    print("For region ", rgn)
    print("intersecting pixels are ", intersect_hps)

    print('Invoke get_objects_by_hp with box region')
    colls = cat.get_objects_by_region(0, rgn, obj_type_set=set(['galaxy']) )
    # Try out get_objects_by_hp with no region
    #print('Invoke get_objects_by_hp with region=None')
    #colls = cat.get_objects_by_hp(0, 9812, None, obj_type_set=set(['galaxy']) )
    print('Number of collections returned:  ', len(colls))

    for c in colls:
        print("For hpid ", c.get_hpid(), "found ", len(c), " objects")
        print("First object: ")
        print(c[0], '\nid=', c[0].id, ' ra=', c[0].ra, ' dec=', c[0].dec,
              ' belongs_index=', c[0]._belongs_index)

        print("Slice [1:3]")
        slice13 = c[1:3]
        for o in slice13:
            print('id=',o.id, ' ra=',o.ra, ' dec=',o.dec, ' belongs_index=',
                  o._belongs_index)
        print("Object 1000")
        print(c[1000], '\nid=', c[1000].id, ' ra=', c[1000].ra, ' dec=',
              c[1000].dec,
              ' belongs_index=', c[1000]._belongs_index)
        slice_late = c[163994:163997]
        print('\nobjects indexed 163994 through 163996')
        for o in slice_late:
            print('id=',o.id, ' ra=',o.ra, ' dec=',o.dec, ' belongs_index=',
                  o._belongs_index)

    coll = colls[0]
    print(type(coll))
    #reds = coll.redshifts()
    #print('redshift length: ', len(reds))
    #print('first few: ', reds[0], ' ', reds[1], ' ', reds[2])
    redshift0 = coll[0].redshift
    print('First redshift: ', redshift0)
