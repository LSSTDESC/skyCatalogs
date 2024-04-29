import os
import re
import logging
import healpy
import numpy as np
import numpy.ma as ma
from astropy import units as u
import lsst.sphgeom
from skycatalogs.objects.base_object import load_lsst_bandpasses, load_roman_bandpasses
from skycatalogs.utils.catalog_utils import CatalogContext
from skycatalogs.objects.base_object import ObjectList, ObjectCollection
from skycatalogs.objects.gaia_object import GaiaObject, GaiaCollection
from skycatalogs.objects.sso_object import SsoObject, SsoCollection
from skycatalogs.objects.sso_object import EXPOSURE_DEFAULT
# from skycatalogs.objects.sso_object import find_sso_files
from skycatalogs.readers import ParquetReader
from skycatalogs.utils.sed_tools import TophatSedFactory, DiffskySedFactory
from skycatalogs.utils.sed_tools import SsoSedFactory
from skycatalogs.utils.sed_tools import MilkyWayExtinction
from skycatalogs.utils.config_utils import Config
from skycatalogs.utils.shapes import Box, Disk, PolygonalRegion
from skycatalogs.utils.shapes import compute_region_mask
from skycatalogs.objects.sncosmo_object import SncosmoObject, SncosmoCollection
from skycatalogs.objects.star_object import StarObject
from skycatalogs.objects.galaxy_object import GalaxyObject
from skycatalogs.objects.diffsky_object import DiffskyObject
from skycatalogs.objects.snana_object import SnanaObject, SnanaCollection

__all__ = ['SkyCatalog', 'open_catalog']


# This function should maybe be moved to utils
def _get_intersecting_hps(hp_ordering, nside, region):
    '''
    Given healpixel structure defined by hp_ordering and nside, find
    all healpixels which intersect region, where region may be a box,
    a disk or a polygonal region.
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

        pixels = healpy.query_polygon(nside, vec, inclusive=True, nest=False)
    elif isinstance(region, Disk):
        # Convert inputs to the types query_disk expects
        center = healpy.pixelfunc.ang2vec(region.ra, region.dec,
                                          lonlat=True)
        radius_rad = (region.radius_as * u.arcsec).to_value('radian')

        pixels = healpy.query_disc(nside, center, radius_rad, inclusive=True,
                                   nest=False)

    elif isinstance(region, PolygonalRegion):
        pixels = healpy.query_polygon(nside, region.get_vertices(),
                                      inclusive=True, nest=False)
    else:
        raise ValueError('Unsupported region type')

    # ensure pixels are always presented in the same order
    pixels = list(pixels)
    pixels.sort()
    return pixels


def _compress_via_mask(tbl, id_column, region, source_type='galaxy',
                       mjd=None, exposure=EXPOSURE_DEFAULT):
    '''
    Parameters
    ----------
    tbl          data table including columns named "ra", "dec", and id_column
                 (and also "object_type" colum if source_type is "star"
                 or "Gaia_star", possibly
                 start_mjd and end_mjd if galaxy is False and mjd is not
                 None)
    id_column    string
    region       mask should restrict to this region (or not at all if None)
    source_type  string of expected object type
    mjd          if not none, may be used to filter transient or variable
                 objects
    exposure     length of exposure if mjd is not None

    Returns
    -------
    4 values for galaxies and snana: ra, dec, id, mask
    5 values for pointsources: ra, dec, id, object_type, mask
    5 values for SSO: ra, dec, id, mjd, mask
    If objects are in the region, ra, dec, id correspond to those objects.
    mask will mask off unused objects
    If there are no objects in the region, all return values are None

    '''
    if isinstance(tbl[id_column][0], (int, np.int64)):
        tbl[id_column] = [str(an_id) for an_id in tbl[id_column]]

    no_obj_type_return = (source_type in {'galaxy', 'diffsky_galaxy',
                                          'snana', 'sso'})
    no_mjd_return = (source_type != 'sso')   # for now
    transient_filter = ('start_mjd' in tbl) and ('end_mjd' in tbl) and mjd is not None
    variable_filter = ('mjd' in tbl)

    if region is not None:
        if isinstance(region, PolygonalRegion):
            # Using native PolygonalRegion selection is slow, so
            # pre-mask the RA, Dec values using an acceptance
            # cone that encloses all of the vertices.

            # Compute the mean direction from the region vertices.
            vertices = region.get_vertices()
            mean_dir = vertices[0]
            for vertex in vertices[1:]:
                mean_dir += vertex
            mean_dir = lsst.sphgeom.UnitVector3d(mean_dir)

            # Find largest vertex offset from the mean direction in arcsec.
            max_offset = np.degrees(np.arccos(min([mean_dir.dot(_)
                                                   for _ in vertices])))*3600.

            # Construct a "Disk" enclosing the polygonal region.
            lon_lat = lsst.sphgeom.LonLat(mean_dir)
            ra = lon_lat.getLon().asDegrees()
            dec = lon_lat.getLat().asDegrees()
            disk_region = Disk(ra, dec, max_offset)

            mask = compute_region_mask(disk_region, tbl['ra'], tbl['dec'])
            if all(mask):
                if no_obj_type_return and no_mjd_return:
                    return None, None, None, None
                else:    # currently if object type is returned, mjd is not
                    return None, None, None, None, None

            # Get compressed ra, dec
            ra_compress = ma.array(tbl['ra'], mask=mask).compressed()
            dec_compress = ma.array(tbl['dec'], mask=mask).compressed()
            poly_mask = compute_region_mask(region, ra_compress, dec_compress)

            # Get indices of unmasked
            ixes = np.where(~mask)

            # Update bounding box mask with polygonal mask
            mask[ixes] |= poly_mask
        else:
            mask = compute_region_mask(region, tbl['ra'], tbl['dec'])

        if transient_filter:
            time_mask = _compute_transient_mask(mjd, tbl['start_mjd'],
                                                tbl['end_mjd'])
            mask = np.logical_or(mask, time_mask)
        elif variable_filter:
            time_mask = _compute_variable_mask(mjd, tbl['mjd'], exposure)
            mask = np.logical_or(mask, time_mask)

        if all(mask):
            if no_obj_type_return and no_mjd_return:
                return None, None, None, None
            else:
                return None, None, None, None, None
        else:
            ra_compress = ma.array(tbl['ra'], mask=mask).compressed()
            dec_compress = ma.array(tbl['dec'], mask=mask).compressed()
            id_compress = ma.array(tbl[id_column], mask=mask).compressed()
            if no_obj_type_return:
                if no_mjd_return:
                    return ra_compress, dec_compress, id_compress, mask
                else:
                    mjd_compress = ma.array(tbl['mjd'], mask=mask).compressed()
                    return ra_compress, dec_compress, id_compress, mjd_compress, mask
            else:

                object_type_compress = ma.array(tbl['object_type'],
                                                mask=mask).compressed()
                return ra_compress, dec_compress, id_compress, object_type_compress, mask
    else:
        if no_obj_type_return:
            if transient_filter:
                time_mask = _compute_transient_mask(mjd, tbl['start_mjd'],
                                                    tbl['end_mjd'])
                ra_compress = ma.array(tbl['ra'], mask=time_mask).compressed()
                dec_compress = ma.array(tbl['dec'],
                                        mask=time_mask).compressed()
                id_compress = ma.array(tbl[id_column],
                                       mask=time_mask).compressed()
                return ra_compress, dec_compress, id_compress, time_mask
            elif variable_filter:
                time_mask = _compute_variable_mask(mjd, tbl['mjd'], exposure)
                if time_mask is not None:
                    ra_compress = ma.array(tbl['ra'], mask=time_mask).compressed()
                    dec_compress = ma.array(tbl['dec'],
                                            mask=time_mask).compressed()
                    id_compress = ma.array(tbl[id_column],
                                           mask=time_mask).compressed()
                    mjd_compress = ma.array(tbl['mjd'], mask=time_mask).compressed()
                    return ra_compress, dec_compress, id_compress, mjd_compress, time_mask
                else:
                    return tbl['ra'], tbl['dec'], tbl[id_column], tbl['mjd'], None
            else:
                return tbl['ra'], tbl['dec'], tbl[id_column], None
        else:
            return tbl['ra'], tbl['dec'], tbl[id_column], tbl['object_type'], None


def _compute_transient_mask(current_mjd, start_mjd, end_mjd):
    '''
    Starting with an existing mask of excluded objects, exclude additional
    objects not visible at time current_mjd
    Parameters
    ----------
    current_mjd    Float  Time for which mask will be computed
    start_mjd      Array of float. Bound on when object starts being
                   detectable
    end_mjd        Array of float. Bound on when object ceases being
                   visible
    Returns
    -------
    mask of objects detectable at specified time
    '''
    mask = np.logical_or((current_mjd > end_mjd), (current_mjd < start_mjd))

    return mask


SECONDS_PER_DAY = 24.0*3600.0
MJD_EPS = 0.1/SECONDS_PER_DAY
JD_SEC = 1.0/SECONDS_PER_DAY


def _compute_variable_mask(current_mjd, mjd_column, exposure, epsilon=MJD_EPS):
    '''
    Compute mask to exclude all entries with
    mjd_column  > current_mjd - epsilon and mjd < current_mjd + exposure

    Parameters
    ----------
    current_mjd  float            mjd of interest
    mjd_column   array of float   mjd for each entry
    exposure     float            exposure in seconds
    epsilon      float            tolerance for matching mjd entry

    Returns
    -------
    mask
    '''
    if not current_mjd:
        return None

    exposure_jd = exposure * JD_SEC
    mask = np.logical_or(mjd_column < (current_mjd - epsilon),
                         mjd_column > (current_mjd + exposure_jd))
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
        self._global_partition = None
        if 'area_partition' in self._config.keys():       # old-style config
            self._global_partition = self._config['area_partition']

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

        self._sky_root = os.path.abspath(sky_root)

        self._logger.info(f'Catalog data will be read from {self._cat_dir}')
        # There may be more to do at this point but not too much.
        # In particular, don't read in anything from data files

        self.verbose = verbose
        self._validate_config()

        # Outer dict: hpid for key. Value is another dict
        #    with keys 'files', 'object_types', each with value another dict
        #       for 'files', map filepath to handle (initially None)
        #       for 'object_types', map object type to filepath
        self._hp_info = dict()
        _ = self._find_all_hps()

        # NOTE: the use of TophatSedFactory is appropriate *only* for an
        # input galaxy catalog with format like cosmoDC2, which includes
        # definitions of tophat SEDs. A different implementation will
        # be needed for newer galaxy catalogs
        th_parameters = self._config.get_tophat_parameters()
        if th_parameters:
            self._observed_sed_factory =\
                TophatSedFactory(th_parameters, config['Cosmology'])
        elif 'diffsky_galaxy' in config['object_types']:
            self._observed_sed_factory =\
                DiffskySedFactory(self._cat_dir,
                                  config['object_types']['diffsky_galaxy']
                                  ['sed_file_template'],
                                  config['Cosmology'])
        if 'sso' in config['object_types']:
            self._sso_sed_path = config['provenance']['inputs'].get('sso_sed',
                                                                    'sso_sed.db')

            self._sso_sed_factory = SsoSedFactory(self._sso_sed_path)
        self._extinguisher = MilkyWayExtinction()

        # Make our properties accessible to BaseObject, etc.
        self.cat_cxt = CatalogContext(self)

        # register object types which are in the config
        if 'gaia_star' in config['object_types']:
            self.cat_cxt.register_source_type('gaia_star',
                                              object_class=GaiaObject,
                                              collection_class=GaiaCollection,
                                              custom_load=True)
        if 'sncosmo' in config['object_types']:
            self.cat_cxt.register_source_type(
                'sncosmo',
                object_class=SncosmoObject,
                collection_class=SncosmoCollection)
        if 'star' in config['object_types']:
            self.cat_cxt.register_source_type('star',
                                              object_class=StarObject)
        if 'galaxy' in config['object_types']:
            self.cat_cxt.register_source_type('galaxy',
                                              object_class=GalaxyObject)
        if 'snana' in config['object_types']:
            self.cat_cxt.register_source_type('snana',
                                              object_class=SnanaObject,
                                              collection_class=SnanaCollection)
        if 'diffsky_galaxy' in config['object_types']:
            self.cat_cxt.register_source_type('diffsky_galaxy',
                                              object_class=DiffskyObject)
        if 'sso' in config['object_types']:
            self.cat_cxt.register_source_type('sso',
                                              object_class=SsoObject,
                                              collection_class=SsoCollection)

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

    def _find_hps_by_type(self, name, type_config):
        '''
        Parameters
        ----------
        name        string    object type name
        type_config      dict      config information about name

        Returns
        -------
        set of healpixels with data for specified type

        Side-effects: Update self._hp_info
        '''
        # Assume flux files, if they exist, are in same dir. as main files
        if 'file_template' not in type_config:
            return set()
        tmpl = type_config['file_template']
        rel_dir, tmpl = os.path.split(tmpl)
        tmpl_keys = [tmpl]
        if 'flux_file_template' in type_config:
            flux_tmpl = type_config['flux_file_template']
            flux_dir, flux_tmpl = os.path.split(flux_tmpl)
            if flux_dir != rel_dir:
                raise ValueError(f'Flux and main files for {name} in different directories')
            tmpl_keys.append(flux_tmpl)
        to_search = self._cat_dir
        if rel_dir != '':
            to_search = os.path.join(self._cat_dir, rel_dir)

        files = os.listdir(to_search)
        hp_set = set()
        for f in files:
            if rel_dir == '':
                fpath = f
            else:
                fpath = os.path.join(rel_dir, f)
            for k in tmpl_keys:
                if k.find('?P<healpix>') != -1:
                    m = re.fullmatch(k, f)
                    if m:
                        hp = int(m['healpix'])
                        hp_set.add(hp)

                        if hp not in self._hp_info:
                            self._hp_info[hp] = {'files': {fpath: None},
                                                 'object_types': {name: [fpath]}}
                        else:
                            this_hp = self._hp_info[hp]
                            # Value of 'object_types' is now a list
                            if fpath not in this_hp['files']:
                                this_hp['files'][fpath] = None
                            if name in this_hp['object_types']:
                                this_hp['object_types'][name].append(fpath)
                            else:
                                this_hp['object_types'][name] = [fpath]
        return hp_set

    def _find_all_hps(self):
        '''
        For each healpix with files matching pattern in or under the
        directory containing the config file, update _hp_info as needed to keep
        track of all files for that healpix and the object types included in
        those files.

        Returns
        -------
        Sorted list of healpix pixels with at least one file in the directory

        '''
        if len(self._hp_info) > 0:
            hp_list = list(self._hp_info.keys())
            hp_list.sort()
            return hp_list

        o_types = list(self.toplevel_only((self._config['object_types'])))
        hp_set = set()
        for (k, v) in [(o_t, self._config['object_types'][o_t]) for o_t in o_types]:
            new_hps = self._find_hps_by_type(k, v)
            hp_set.update(new_hps)

        hp_list = list(hp_set)
        hp_list.sort()
        return hp_list

    def hps_by_type(self, object_type):
        '''
        Parameters
        ----------
        object_type   string

        Returns
        -------
        list of healpixels (int) having data for object type object_type
        (or its parent type if it has one)
        '''
        if not self._hp_info:
            return []
        hps = []
        if 'parent' in self._config['object_types'][object_type]:
            object_type = self._config['object_types'][object_type]['parent']
        for hp, val in self._hp_info.items():
            if object_type in val['object_types']:
                hps.append(hp)
        hps.sort()
        return hps

    def get_hps_by_region(self, region, object_type='galaxy'):
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
        if self._global_partition is None:
            partition = self._config['object_types'][object_type]['area_partition']
        else:
            partition = self._global_partition
        return _get_intersecting_hps(partition['ordering'], partition['nside'],
                                     region)

    def get_object_type_names(self):
        '''
        Returns
        -------
        All object type names in the catalog's config
        '''
        names = list(set(self._config['object_types'].keys()))
        names.sort()
        return names

    def default_object_type_set(self):
        if 'default_object_types' in self._config.keys():
            return set(self._config['default_object_types'])
        else:
            return self.get_object_type_names()

    def toplevel_only(self, object_types):
        '''
        Parameters
        ----------
        object_types     Set of object type names
        Remove object types with a parent.  Add in the parent.

        Return the resulting types (as list) and their values
        '''
        objs_copy = set(object_types)
        for obj in object_types:
            parent = self._config.get_object_parent(obj)
            if parent is not None:
                objs_copy.remove(obj)
                objs_copy.add(parent)
        return objs_copy

    def get_objects_by_region(self, region, obj_type_set=None, mjd=None,
                              exposure=EXPOSURE_DEFAULT):
        '''
        Parameters
        ----------
        region         region is a named tuple(may be box or circle)
                       or object of type PolygonalRegion
        obj_type_set   Return only these objects. Defaults to value in config if
                       specified; else default to all defined in config
        mjd            MJD of observation epoch.
        exposure       exposure length (seconds)

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

        object_list = ObjectList()
        if obj_type_set is None:
            obj_types = self.default_object_type_set()
        else:
            obj_types = set(self.get_object_type_names()).intersection(obj_type_set)
        obj_types = self.toplevel_only(obj_types)

        # Ensure they're always ordered the same way
        obj_types = list(obj_types)
        obj_types.sort()

        for ot in obj_types:
            new_list = self.get_object_type_by_region(region, ot, mjd=mjd,
                                                      exposure=exposure)
            object_list.append_object_list(new_list)

        return object_list

    def get_object_type_by_region(self, region, object_type, mjd=None,
                                  exposure=EXPOSURE_DEFAULT):
        '''
        Parameters
        ----------
        region        box, circle or PolygonalRegion. Supported region
                      types made depend on object_type
        object_type   known object type without parent
        mjd           MJD of observation epoch.
        exposure      exposure (seconds)

        Returns
        -------
        an ObjectList containing all objects found.

        '''

        out_list = ObjectList()
        if self._global_partition is not None:
            partition = self._global_partition
        else:
            partition = self._config['object_types'][object_type]['area_partition']

        coll_type = self.cat_cxt.lookup_collection_type(object_type)
        if coll_type is not None:
            if self.cat_cxt.use_custom_load(object_type):
                coll = coll_type.load_collection(region, self, mjd=mjd,
                                                 exposure=EXPOSURE_DEFAULT)
                if isinstance(coll, ObjectCollection):
                    out_list.append_collection(coll)
                else:  # ObjectList
                    out_list.append_object_list(coll)
                return out_list

        if partition != 'None':
            if partition['type'] == 'healpix':
                hps = self.get_hps_by_region(region, object_type)
                for hp in hps:
                    c = self.get_object_type_by_hp(hp, object_type, region,
                                                   mjd,
                                                   exposure=EXPOSURE_DEFAULT)
                    if len(c) > 0:
                        out_list.append_object_list(c)
                return out_list
        else:
            raise NotImplementedError(f'Unsupported object type {object_type}')

    def get_object_type_by_hp(self, hp, object_type, region=None, mjd=None,
                              exposure=EXPOSURE_DEFAULT):
        object_list = ObjectList()

        #  Do we need to check more specifically by object type?
        # if hp not in self._hp_info:
        if hp not in self.hps_by_type(object_type):
            msg = f'In SkyCatalog.get_object_type_by_hp, healpix {hp}  '
            msg += f'intersects region but has no data file for {object_type}'
            self._logger.warning(msg)
            return object_list

        if object_type in ['galaxy', 'diffsky_galaxy']:
            columns = ['galaxy_id', 'ra', 'dec']
            id_name = 'galaxy_id'
        elif object_type in ['snana']:
            columns = ['id', 'ra', 'dec', 'start_mjd', 'end_mjd']
            id_name = 'id'
        elif object_type in ['star', 'sncosmo']:
            columns = ['object_type', 'id', 'ra', 'dec']
            id_name = 'id'
        elif object_type in ['sso']:
            id_name = 'id'
            columns = ['id', 'ra', 'dec', 'mjd']
        else:
            raise NotImplementedError(f'Unsupported object type {object_type}')

        coll_class = self.cat_cxt.lookup_collection_type(object_type)

        if self.verbose:
            print('Working on healpix pixel ', hp)
        rdr_ot = dict()   # maps readers to set of object types it reads

        if 'file_template' in self._config['object_types'][object_type]:
            f_list = self._hp_info[hp]['object_types'][object_type] \
                if object_type in self._hp_info[hp]['object_types'] else []
        elif 'parent' in self._config['object_types'][object_type]:
            f_list = self._hp_info[hp]['object_types'][self._config['object_types'][object_type]['parent']]

        for f in f_list:
            if self._hp_info[hp]['files'][f] is None:        # no reader yet
                full_path = os.path.join(self._cat_dir, f)
                the_reader = ParquetReader(full_path, mask=None)
                self._hp_info[hp]['files'][f] = the_reader
            else:
                the_reader = self._hp_info[hp]['files'][f]
                # associate object type with this reader
            if the_reader in rdr_ot:
                rdr_ot[the_reader].add(object_type)
            else:
                rdr_ot[the_reader] = set([object_type])

        # Find readers needed to get minimal columns for our object type
        the_readers = []
        for rdr in rdr_ot:
            if object_type in rdr_ot[rdr]:
                the_readers.append(rdr)

        # Unfortunately galaxies and point sources can't be handled quite
        # the same way since we need to read an extra column for pointsources
        # There will have to be some "if galaxy... else ... code"
        for rdr in the_readers:
            if 'ra' not in rdr.columns:
                continue

            # Make a collection for each row group
            for rg in range(rdr.n_row_groups):
                arrow_t = rdr.read_columns(columns, None, rg)
                if object_type in {'galaxy', 'diffsky_galaxy', 'snana'}:
                    ra_c, dec_c, id_c, mask =\
                        _compress_via_mask(arrow_t,
                                           id_name,
                                           region,
                                           source_type=object_type,
                                           mjd=mjd)
                    if ra_c is not None:
                        new_collection = coll_class(ra_c, dec_c, id_c,
                                                    object_type, hp, self,
                                                    region=region,
                                                    mjd=mjd,
                                                    mask=mask,
                                                    readers=the_readers,
                                                    row_group=rg)
                        if object_type == 'snana':
                            # file pattern really should come from cfg
                            base = f'snana_{hp}.hdf5'
                            SED_file = os.path.join(self._cat_dir, base)
                            new_collection.set_SED_file(SED_file)
                        object_list.append_collection(new_collection)
                elif object_type in {'sso'}:
                    ra_c, dec_c, id_c, mjd_c, mask =\
                        _compress_via_mask(arrow_t,
                                           id_name,
                                           region,
                                           source_type=object_type,
                                           mjd=mjd, exposure=exposure)
                    if ra_c is not None:
                        new_collection = SsoCollection(ra_c, dec_c, id_c, hp,
                                                       self,
                                                       mjd_individual=mjd_c,
                                                       region=region,
                                                       mjd=mjd, mask=mask,
                                                       readers=the_readers,
                                                       row_group=rg,
                                                       exposure=exposure)
                        object_list.append_collection(new_collection)
                else:
                    ra_c, dec_c, id_c, object_type_c, mask =\
                        _compress_via_mask(arrow_t, id_name, region,
                                           source_type={object_type})
                    if ra_c is not None and object_type_c[0] == object_type:
                        new_collection = coll_class(ra_c, dec_c, id_c,
                                                    object_type_c[0], hp,
                                                    self, region=region,
                                                    mjd=mjd,
                                                    mask=mask,
                                                    readers=the_readers,
                                                    row_group=rg)
                        object_list.append_collection(new_collection)

        return object_list

    # For generator version, do this a row group at a time
    #    but if region cut leaves too small a list, read more rowgroups
    #    to achieve a reasonable size list (or exhaust the file)

    def get_object_iterator_by_hp(self, hp, obj_type_set=None,
                                  max_chunk=None, mjd=None):
        '''
        Parameters
        ----------
        mjd            MJD of observation epoch
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
    # Get bandpasses in case we need to compute fluxes
    _ = load_lsst_bandpasses()
    _ = load_roman_bandpasses()

    from skycatalogs.utils.config_utils import open_config_file

    config_dict = open_config_file(config_file)
    return SkyCatalog(config_dict, skycatalog_root=skycatalog_root, mp=mp,
                      verbose=verbose)