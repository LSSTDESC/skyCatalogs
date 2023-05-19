import os
import sys
import re
import yaml
import logging
from collections import namedtuple
import healpy
import numpy as np
import numpy.ma as ma
import pyarrow.parquet as pq
from astropy import units as u
from desc.skycatalogs.objects.base_object import load_lsst_bandpasses
from desc.skycatalogs.objects.base_object import  CatalogContext
from desc.skycatalogs.objects import *
from desc.skycatalogs.readers import *
from desc.skycatalogs.readers import ParquetReader
from desc.skycatalogs.utils.sed_tools import ObservedSedFactory
from desc.skycatalogs.utils.sed_tools import MilkyWayExtinction
from desc.skycatalogs.utils.config_utils import Config

__all__ = ['SkyCatalog', 'open_catalog', 'Box', 'Disk', 'PolygonalRegion']


Box = namedtuple('Box', ['ra_min', 'ra_max', 'dec_min', 'dec_max'])

# radius is measured in arcseconds
Disk = namedtuple('Disk', ['ra', 'dec', 'radius_as'])

from lsst.sphgeom import ConvexPolygon, UnitVector3d, LonLat
class PolygonalRegion:

    def __init__(self, vertices_radec=None, convex_polygon=None):
        '''
        Supply either an object of type lsst.sphgeom.ConvexPolygon
        or a list of vertices, each a tuple (ra, dec) in degrees,
        which describe a convex polygon
        '''
        if convex_polygon:
            if isinstance(convex_polygon, ConvexPolygon):
                self._convex_polygon = convex_polygon
                return
        if vertices_radec:
            if not isinstance(vertices_radec, list):
                raise TypeError(f'PolygonalRegion: Argument {vertices_radec} is not a list')
            vertices = [UnitVector3d(LonLat.fromDegrees(v_rd[0], v_rd[1])) for v_rd in vertices_radec]
            self._convex_polygon = ConvexPolygon(vertices)
            return
        raise ValueError('PolygonalRegion: Either vertices_radec or convex_polygon must have an acceptable value')

    def get_vertices(self):
        '''
        Return vertices as list of 3d vectors
        '''
        return self._convex_polygon.getVertices()

    def get_vertices_radec(self):
        v3d = self.get_vertices()
        vertices_radec = []
        for v in v3d:
            vertices_radec.append((LonLat.longitudeOf(v).asDegrees(),
                                   LonLat.latitudeOf(v).asDegrees()))
        return vertices_radec

    def get_containment_mask(self, ra, dec, included=True):
        '''
        Parameters
        ----------
        ra, dec      parallel float arrays, units are degrees. Together
                     they define the list of points to be checked for
                     containment
        included     boolean   If true, mask bit will be set to True for
                               contained points, else False.    Reverse
                               the settings if included is False
        '''
        # convert to radians
        ra = [(r * u.degree).to_value(u.radian) for r in ra]
        dec = [(d * u.degree).to_value(u.radian) for d in dec]

        mask = self._convex_polygon.contains(ra, dec)
        if included:
            return mask
        else:
            return np.logical_not(mask)

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
    if isinstance(region, Box):
        vec = healpy.pixelfunc.ang2vec([region.ra_min, region.ra_max,
                                        region.ra_max, region.ra_min],
                                       [region.dec_min, region.dec_min,
                                        region.dec_max, region.dec_max],
                                       lonlat=True)

        return healpy.query_polygon(nside, vec, inclusive=True, nest=False)
    if isinstance(region, Disk):
        # Convert inputs to the types query_disk expects
        center = healpy.pixelfunc.ang2vec(region.ra, region.dec,
                                          lonlat=True)
        radius_rad = (region.radius_as * u.arcsec).to_value('radian')

        pixels = healpy.query_disc(nside, center, radius_rad, inclusive=True,
                                   nest=False)
        return pixels

    if isinstance(region, PolygonalRegion):
        return healpy.query_polygon(nside, region.get_vertices(),
                                    inclusive=True, nest=False)

def _compress_via_mask(tbl, id_column, region, galaxy=True):
    '''
    Parameters
    ----------
    tbl          data table including columns named "ra", "dec", and id_column
    id_column    string
    region       mask should restrict to this region (or not at all if None)

    Returns
    -------
    4 values for galaxies: ra, dec, id, mask
    5 values for pointsources: ra, dec, id, object_type, mask
    If objects are in the region, ra, dec, id correspond to those objects.
    mask will mask off unused objects
    If there are no objects in the region, all return values are None

    '''
    if region is not None:
        if isinstance(region, PolygonalRegion):        # special case
            # Find bounding box
            vertices = region.get_vertices_radec()  # list of (ra, dec)
            ra_max = max([v[0] for v in vertices])
            ra_min = min([v[0] for v in vertices])
            dec_max = max([v[1] for v in vertices])
            dec_min = min([v[1] for v in vertices])
            bnd_box = Box(ra_min, ra_max, dec_min, dec_max)
            # Compute mask for that box
            mask = _compute_mask(bnd_box, tbl['ra'], tbl['dec'])
            if all(mask): # even bounding box doesn't intersect table rows
                if galaxy:
                    return None, None, None, None
                else:
                    return None, None, None, None, None

            # Get compressed ra, dec
            ra_compress = ma.array(tbl['ra'], mask=mask).compressed()
            dec_compress = ma.array(tbl['dec'], mask=mask).compressed()
            poly_mask = _compute_mask(region, ra_compress, dec_compress)

            # Get indices of unmasked
            ixes = np.where(~mask)

            # Update bounding box mask with polygonal mask
            mask[ixes] |= poly_mask
        else:
            mask = _compute_mask(region, tbl['ra'], tbl['dec'])

        if all(mask):
            if galaxy:
                return None, None, None, None
            else:
                return None, None, None, None, None
        else:
            ra_compress = ma.array(tbl['ra'], mask=mask).compressed()
            dec_compress = ma.array(tbl['dec'], mask=mask).compressed()
            id_compress = ma.array(tbl[id_column], mask=mask).compressed()
            if galaxy:
                return ra_compress, dec_compress, id_compress, mask
            else:
                object_type_compress = ma.array(tbl['object_type'],
                                                mask=mask).compressed()
                return ra_compress, dec_compress, id_compress, object_type_compress, mask
    else:
        if galaxy:
            return tbl['ra'], tbl['dec'], tbl[id_column],None
        else:
            return tbl['ra'], tbl['dec'], tbl[id_column],tbl['object_type'],None

def _compute_mask(region, ra, dec):
    '''
    Compute mask according to region for provided data
    Parameters
    ----------
    region         Supported shape (box, disk)  or None
    ra,dec         Coordinates for data to be masked, in degrees
    Returns
    -------
    mask of elements in ra, dec arrays to be omitted

    '''
    mask = None
    if isinstance(region, Box):
        mask = np.logical_or((ra < region.ra_min),
                             (ra > region.ra_max))
        mask = np.logical_or(mask, (dec < region.dec_min))
        mask = np.logical_or(mask, (dec > region.dec_max))
    if isinstance(region, Disk):
        # Change positions to 3d vectors to measure distance
        p_vec = healpy.pixelfunc.ang2vec(ra, dec, lonlat=True)

        c_vec = healpy.pixelfunc.ang2vec(region.ra,
                                         region.dec,
                                         lonlat=True)
        radius_rad = (region.radius_as * u.arcsec).to_value('radian')


        # Rather than comparing arcs, it is equivalent to compare chords
        # (or square of chord length)
        diff = p_vec - c_vec
        obj_chord_sq = np.sum(np.square(p_vec - c_vec),axis=1)

        # This is to be compared to square of chord for angle a corresponding
        # to disk radius.  That's 4(sin(a/2)^2)
        rad_chord_sq = 4 * np.square(np.sin(0.5 * radius_rad) )
        mask = obj_chord_sq > rad_chord_sq
    if isinstance(region, PolygonalRegion):
        mask = region.get_containment_mask(ra, dec, included=False)
    return mask

class SkyCatalog(object):
    '''
    A base class with derived classes for galaxies, static (w.r.t. coordinates)
    point sources, SSOs

    '''
    def __init__(self, config, mp=False, skycatalog_root=None, verbose=False,
                 loglevel='INFO'):
        '''
        Parameters
        ----------
        config:  dict.  Typically the result of loading a yaml file
        mp:      boolean  Default False. Set True to enable multiprocessing.
        skycatalog_root: If not None, overrides value in config or
                         in environment variable SKYCATALOG_ROOT
        '''
        self._config = Config(config)
        self._logger = logging.getLogger('skyCatalogs.client')

        if not self._logger.hasHandlers():
            self._logger.setLevel(loglevel)
            ch = logging.StreamHandler()
            ch.setLevel(loglevel)
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            ch.setFormatter(formatter)
            self._logger.addHandler(ch)

        self._mp = mp
        if 'schema_version' not in config and 'schema_version' not in config['provenance']['versioning']:
            self._cat_dir = config['root_directory']
        else:
            sky_root = config['skycatalog_root']        # default
            if skycatalog_root:
                sky_root = skycatalog_root
            else:
                sky_root_env = os.getenv('SKYCATALOG_ROOT', None)
                if sky_root_env:
                    sky_root = sky_root_env

            self._cat_dir = os.path.join(sky_root, config['catalog_dir'])

        self._logger.info(f'Catalog data will be read from {self._cat_dir}')
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

        self._observed_sed_factory =\
            ObservedSedFactory(self._config.get_tophat_parameters(),
                               config['Cosmology'])
        self._extinguisher = MilkyWayExtinction(self.observed_sed_factory)

        # Make our properties accessible to BaseObject, etc.
        self.catalog_context = CatalogContext(self)

    @property
    def observed_sed_factory(self):
        return self._observed_sed_factory

    @property
    def extinguisher(self):
        return self._extinguisher

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
        and the object types included in those files.

        Returns
        -------
        Set of healpix pixels with at least one file in the directory

        '''

        # If major organization is by healpix, healpix # will be in
        # subdirectory name.  Otherwise there may be no subdirectories
        # or data may be organized by component type.
        # Here only handle case where data files are directly in root dir
        files = os.listdir(self._cat_dir)
        o_types = self._config['object_types']

        hp_set = set()
        for f in files:
            for ot in o_types:
                # find all keys containing the string 'file_template'
                template_keys = [k for k in o_types[ot] if 'file_template' in k]
                for k in template_keys:
                    m = re.fullmatch(o_types[ot][k], f)
                    if m:
                        hp = int(m['healpix'])
                        hp_set.add(hp)

                        if hp not in self._hp_info:
                            self._hp_info[hp] = {'files' : {f : None},
                                                 'object_types' : {ot : [f]}}
                        else:
                            this_hp = self._hp_info[hp]
                            # Value of 'object_types' is now a list
                            if f not in this_hp['files'] :
                                this_hp['files'][f] = None
                            if ot in this_hp['object_types']:
                                this_hp['object_types'][ot].append(f)
                            else:
                                this_hp['object_types'][ot] = [f]
        return hp_set

    def get_hps_by_region(self, region):
        '''
        Parameters
        ----------
        Region can be a box (named 4-tuple (min-ra, max-ra, min-dec, max-dec)),
        a circle (named 3-tuple (ra, dec, radius)) or of type
        PolygonalRegion.
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
        region         region is a named tuple(may be box or circle)
                       or object of type PolygonalRegion
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
        Not yet implemented.

        Parameters
        ----------
        region         Either a box or a circle (each represented as
                       named tuple) or object of type PolygonalRegion
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
        Find all objects in the healpixl (and region if specified)
        of requested types

        Parameters
        ----------
        hp     The healpix pixel.
        region If supplied defines a disk, box or convex polygon.
               Return only objects contained in it
        obj_type_set   Type of objects to fetch.  Defaults to all
        datetime       Ignored for now; no support for moving objects

        Returns
        -------
        ObjectList containing sky objects in the region and the hp
        '''
        object_list = ObjectList()

        if hp not in self._hp_info:
            print(f'WARNING: In SkyCatalog.get_objects_by_hp healpixel {hp} intersects region but has no catalog file')
            return object_list

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

        # Find the right Sky Catalog file or files (depends on obj_type_set)
        # Get file handle(s) and store if we don't already have it (them)

        # Associate object types with files and files with readers.
        # May be > one type per reader (e.g. different kinds of point sources
        # stored in a single file)
        # Also may be more than one file per type (galaxy has two)
        rdr_ot = dict()   # maps readers to set of object types it reads

        for ot in obj_types:
            # FInd files associated with this type
            if 'file_template' in self._config['object_types'][ot]:
                f_list = self._hp_info[hp]['object_types'][ot]
            elif 'parent' in self._config['object_types'][ot]:
                f_list = self._hp_info[hp]['object_types'][self._config['object_types'][ot]['parent']]

            for f in f_list:
                if self._hp_info[hp]['files'][f] is None:            # no reader yet
                    full_path = os.path.join(self._cat_dir, f)
                    the_reader = parquet_reader.ParquetReader(full_path, mask=None)
                    self._hp_info[hp]['files'][f] = the_reader
                else:
                    the_reader = self._hp_info[hp]['files'][f]
                # associate object type with this reader
                if the_reader in rdr_ot:
                    rdr_ot[the_reader].add(ot)
                else:
                    rdr_ot[the_reader] = set([ot])

        # Now get minimal columns for objects using the readers
        galaxy_readers = []
        pointsource_readers = []

        # First make a pass to separate out the readers
        for rdr in rdr_ot:
            if 'galaxy' in rdr_ot[rdr]:
                galaxy_readers.append(rdr)
            elif ot in rdr_ot[rdr]:
                pointsource_readers.append(rdr)

        # Make galaxy collections
        for rdr in galaxy_readers:
            if 'ra' not in rdr.columns:
                continue

            # Make a collection for each row group

            for rg in range(rdr.n_row_groups):
                arrow_t = rdr.read_columns(G_COLUMNS, None, rg)

                ra_c, dec_c, id_c, mask = _compress_via_mask(arrow_t,
                                                             'galaxy_id',
                                                             region)

                if ra_c is not None:
                    new_collection = ObjectCollection(ra_c, dec_c, id_c,
                                                      'galaxy', hp, self,
                                                      region=region, mask=mask,
                                                      readers=galaxy_readers,
                                                      row_group=rg)
                    object_list.append_collection(new_collection)

        # Make point source collection
        for rdr in pointsource_readers:
            if 'ra' not in rdr.columns:
                continue

            # Make a collection for each row group
            for rg in range(rdr.n_row_groups):
                arrow_t = rdr.read_columns(PS_COLUMNS, None, rg)

                ra_c, dec_c, id_c, object_type,mask = _compress_via_mask(arrow_t, 'id', region, galaxy=False)

                if ra_c is not None and object_type[0] in obj_type_set:
                    new_collection = ObjectCollection(ra_c, dec_c, id_c,
                                                      object_type[0], hp, self,
                                                      region=region, mask=mask,
                                                      readers=pointsource_readers, row_group=rg)
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


def open_catalog(config_file, mp=False, skycatalog_root=None, verbose=False):
    '''
    Parameters
    ----------
    config_file   yaml file containing config
    skycatalog_root optional override of skycatalog_root. This value may also
                    be supplied by setting the environment variable
                    SKYCATALOG_ROOT.  Precedence is 1) argument
                    2) environment variable 3) value in config file for
                    key skycatalog_root.  However set, this value
                    joined to value of the key catalog_dir will be used
                    to find the catalog data.

    Returns
    -------
    SkyCatalog
    '''
    # Get LSST bandpasses in case we need to compute fluxes
    load_lsst_bandpasses()
    with open(config_file) as f:
        return SkyCatalog(yaml.safe_load(f), skycatalog_root=skycatalog_root,
                          mp=mp, verbose=verbose)

if __name__ == '__main__':
    import time
    ###cfg_file_name = 'latest.yaml'
    cfg_file_name = 'future_latest.yaml'
    skycatalog_root = os.path.join(os.getenv('SCRATCH'),'desc/skycatalogs')

    if len(sys.argv) > 1:
        cfg_file_name = sys.argv[1]
    cfg_file = os.path.join('/global/homes/j/jrbogart/desc_git/skyCatalogs/cfg',
                            cfg_file_name)

    write_sed = False
    if len(sys.argv) > 2:
        write_sed = True

    # For tract 3828
    #   55.73604 < ra < 57.563452
    #  -37.19001 < dec < -35.702481


    cat = open_catalog(cfg_file, skycatalog_root=skycatalog_root)
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
    vertices = [(ra_min_small, dec_min_small), (ra_min_small, dec_max_small),
                (ra_max_small, dec_max_small), (ra_max_small, dec_min_small)]
    rgn_poly = PolygonalRegion(vertices_radec=vertices)

    intersect_hps = _get_intersecting_hps('ring', 32, rgn)

    print("For region ", rgn)
    print("intersecting pixels are ", intersect_hps)

    intersect_poly_hps = _get_intersecting_hps('ring', 32, rgn_poly)
    print("For region ", rgn_poly)
    print("intersecting pixels are ", intersect_poly_hps)


    print('Invoke get_objects_by_region with box region')
    t0 = time.time()
    object_list = cat.get_objects_by_region(rgn)
    t_done = time.time()
    print('Took ', t_done - t0)
                                            ##### temporary obj_type_set={'galaxy', 'star'} )
    #                                        obj_type_set=set(['galaxy']) )
    # Try out get_objects_by_hp with no region
    #colls = cat.get_objects_by_hp(9812, None, set(['galaxy']) )

    print('Number of collections returned for box:  ', object_list.collection_count)
    print('Object count for box: ', len(object_list))

    print('Invoke get_objects_by_region with polygonal region')
    t0 = time.time()
    object_list_poly = cat.get_objects_by_region(rgn_poly)
    t_done = time.time()
    print('Took ', t_done - t0)

    print('Number of collections returned for polygon:  ',
          object_list_poly.collection_count)
    assert(object_list.collection_count == object_list_poly.collection_count)
    print('Object count for polygon: ', len(object_list_poly))

    fudge = 5
    assert(len(object_list_poly) > len(object_list) - fudge)
    assert(len(object_list_poly) < len(object_list) + fudge)

    print('Now try inscribed polygon which is not a box')
    ra_avg = (ra_min_small + ra_max_small) * 0.5
    dec_avg = (dec_min_small + dec_max_small) * 0.5
    diamond_vertices = [(ra_min_small, dec_avg), (ra_avg, dec_max_small),
                        (ra_max_small, dec_avg), (ra_avg, dec_min_small)]
    rgn_diamond = PolygonalRegion(vertices_radec=diamond_vertices)

    intersect_diamond_hps = _get_intersecting_hps('ring', 32, rgn_diamond)
    print("for diamond region ", rgn_diamond)
    print("intersecting pixels are ", intersect_diamond_hps)

    print('Invoke get_objects_by_region with diamond region')
    t0 = time.time()
    object_list_diamond = cat.get_objects_by_region(rgn_diamond)
    t_done = time.time()
    print('Took ', t_done - t0)

    print('Number of collections returned for diamond:  ',
          object_list_diamond.collection_count)
    print('Object count for diamond: ', len(object_list_diamond))


    #### TEMP FOR DEBUGGING
    ### exit(0)


    colls = object_list.get_collections()
    got_a_sed = False
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
                #(lmbda, f_lambda, magnorm) = o.get_sed(resolution=1.0)
                sed, magnorm = o.get_sed(resolution=1.0)
                print('For star magnorm: ', magnorm)
                if magnorm < 1000:
                    print('Length of sed: ', len(sed.wave_list))
            else:
                for cmp in ['disk', 'bulge', 'knots']:
                    print(cmp)
                    if cmp in o.subcomponents:
                        print(o.get_instcat_entry(component=cmp))
                        sed, _ = o.get_sed(cmp)
                        if sed:
                            print('Length of sed table: ', len(sed.wave_list))
                            if not got_a_sed:
                                got_a_sed = True
                                th = o.get_native_attribute(f'sed_val_{cmp}')
                                print('Tophat values: ', th)
                                sed, _ = o.get_sed(component=cmp)
                                print('Simple sed wavelengths:')
                                print(sed.wave_list)
                                print('Simple sed values:')
                                print([sed(w) for w in sed.wave_list])
                                if write_sed:
                                    o.write_sed('simple_sed.txt', component=cmp)
                                sed_fine, _ = o.get_sed(component=cmp,
                                                        resolution=1.0)
                                print('Bin width = 1 nm')
                                print('Initial wl values', sed_fine.wave_list[:20])
                                print('Start at bin 100', sed_fine.wave_list[100:120])
                                print('Initial values')
                                print([sed_fine(w) for w in sed_fine.wave_list[:20]])
                                print('Start at bin 100')
                                print([sed_fine(w) for w in sed_fine.wave_list[100:120]])
                        else:
                            print('All-zero sed')

                # Try out old wrapper functions
                print("\nget_dust:")
                i_av, i_rv, g_av, g_rv = o.get_dust()
                print(f'i_av={i_av} i_rv={i_rv} g_av={g_av} g_rv={g_rv}')
                print("\nget_wl_params")
                g1, g2, mu = o.get_wl_params()
                print(f'g1={g1} g2={g2} mu={mu}')
                print("\nget_gsobject_components. Keys of returned dict:")
                gs_dict = o.get_gsobject_components()
                print(gs_dict.keys())
                print("\nget_observer_sed_components.  Keys of returned dict:")
                o_seds = o.get_observer_sed_components()
                print(o_seds.keys())

                f = o.get_LSST_flux('i')
                print(f'Flux for i bandpass: {f}')
                fluxes = o.get_LSST_fluxes()
                for k,v in fluxes.items():
                    print(f'Bandpass {k} has flux {v}')

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
