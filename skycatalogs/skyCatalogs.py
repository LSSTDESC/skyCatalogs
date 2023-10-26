import os
import sys
import re
import yaml
from yamlinclude import YamlIncludeConstructor
import logging
import healpy
import numpy as np
import numpy.ma as ma
from astropy import units as u
from skycatalogs.objects.base_object import load_lsst_bandpasses
from skycatalogs.utils.catalog_utils import CatalogContext
from skycatalogs.objects.base_object import ObjectList
from skycatalogs.objects.gaia_object import GaiaObject, GaiaCollection
from skycatalogs.readers import ParquetReader
from skycatalogs.utils.sed_tools import TophatSedFactory
from skycatalogs.utils.sed_tools import MilkyWayExtinction
from skycatalogs.utils.config_utils import Config
from skycatalogs.utils.shapes import Box, Disk, PolygonalRegion
from skycatalogs.objects.sncosmo_object import SncosmoObject
from skycatalogs.objects.star_object import StarObject
from skycatalogs.objects.galaxy_object import GalaxyObject
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


def _compress_via_mask(tbl, id_column, region, source_type={'galaxy'},
                       mjd=None):
    '''
    Parameters
    ----------
    tbl          data table including columns named "ra", "dec", and id_column
                 (and also "object_type" colum if galaxy is False, possibly
                 start_mjd and end_mjd if galaxy is False and mjd is not
                 None)
    id_column    string
    region       mask should restrict to this region (or not at all if None)
    source_type  string or set (or other container) of expected object type(s).
                 For now, must be a singleton
    mjd          if not none, may be used to filter transient objects

    Returns
    -------
    4 values for galaxies and snana: ra, dec, id, mask
    5 values for pointsources: ra, dec, id, object_type, mask
    If objects are in the region, ra, dec, id correspond to those objects.
    mask will mask off unused objects
    If there are no objects in the region, all return values are None

    '''
    if isinstance(tbl[id_column][0], (int, np.int64)):
        tbl[id_column] = [str(an_id) for an_id in tbl[id_column]]

    if not isinstance(source_type, str):
        source_type = set(source_type)
        if len(source_type) != 1:
            raise NotImplementedError('_compress_via_mask only accepts singleton object_type')
        source_type = source_type.pop()
    no_obj_type_return = (source_type in {'galaxy', 'snana'})
    time_filter = ('start_mjd' in tbl) and ('end_mjd' in tbl) and mjd is not None

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
            mask = _compute_region_mask(bnd_box, tbl['ra'], tbl['dec'])
            if all(mask):    # even bounding box doesn't intersect table rows
                if no_obj_type_return:
                    return None, None, None, None
                else:
                    return None, None, None, None, None

            # Get compressed ra, dec
            ra_compress = ma.array(tbl['ra'], mask=mask).compressed()
            dec_compress = ma.array(tbl['dec'], mask=mask).compressed()
            poly_mask = _compute_region_mask(region, ra_compress, dec_compress)

            # Get indices of unmasked
            ixes = np.where(~mask)

            # Update bounding box mask with polygonal mask
            mask[ixes] |= poly_mask
        else:
            mask = _compute_region_mask(region, tbl['ra'], tbl['dec'])

        if time_filter:
            time_mask = _compute_time_mask(mjd, tbl['start_mjd'],
                                           tbl['end_mjd'])
            mask = np.logical_or(mask, time_mask)

        if all(mask):
            if no_obj_type_return:
                return None, None, None, None
            else:
                return None, None, None, None, None
        else:
            ra_compress = ma.array(tbl['ra'], mask=mask).compressed()
            dec_compress = ma.array(tbl['dec'], mask=mask).compressed()
            id_compress = ma.array(tbl[id_column], mask=mask).compressed()
            if no_obj_type_return:
                return ra_compress, dec_compress, id_compress, mask
            else:

                object_type_compress = ma.array(tbl['object_type'],
                                                mask=mask).compressed()
                return ra_compress, dec_compress, id_compress, object_type_compress, mask
    else:
        if no_obj_type_return:
            if time_filter:
                time_mask = _compute_time_mask(mjd, tbl['start_mjd'],
                                               tbl['end_mjd'])
                ra_compress = ma.array(tbl['ra'], mask=time_mask).compressed()
                dec_compress = ma.array(tbl['dec'],
                                        mask=time_mask).compressed()
                id_compress = ma.array(tbl[id_column],
                                       mask=time_mask).compressed()
                return ra_compress, dec_compress, id_compress, time_mask
            else:
                return tbl['ra'], tbl['dec'], tbl[id_column], None
        else:
            return tbl['ra'], tbl['dec'], tbl[id_column], tbl['object_type'], None


def _compute_region_mask(region, ra, dec):
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
        obj_chord_sq = np.sum(np.square(p_vec - c_vec), axis=1)

        # This is to be compared to square of chord for angle a corresponding
        # to disk radius.  That's 4(sin(a/2)^2)
        rad_chord_sq = 4 * np.square(np.sin(0.5 * radius_rad))
        mask = obj_chord_sq > rad_chord_sq
    if isinstance(region, PolygonalRegion):
        mask = region.get_containment_mask(ra, dec, included=False)
    return mask


def _compute_time_mask(current_mjd, start_mjd, end_mjd):
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
        self._observed_sed_factory =\
            TophatSedFactory(th_parameters, config['Cosmology'])

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
            self.cat_cxt.register_source_type('sncosmo',
                                              object_class=SncosmoObject)
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

    def get_objects_by_region(self, region, obj_type_set=None, mjd=None):
        '''
        Parameters
        ----------
        region         region is a named tuple(may be box or circle)
                       or object of type PolygonalRegion
        obj_type_set   Return only these objects. Defaults to value in config if
                       specified; else default to all defined in config
        mjd            MJD of observation epoch.

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
            new_list = self.get_object_type_by_region(region, ot, mjd=mjd)
            object_list.append_object_list(new_list)

        return object_list

    def get_object_type_by_region(self, region, object_type, mjd=None):
        '''
        Parameters
        ----------
        region        box, circle or PolygonalRegion. Supported region
                      types made depend on object_type
        object_type   known object type without parent
        mjd           MJD of observation epoch.

        Returns all objects found
        '''

        out_list = ObjectList()
        if self._global_partition is not None:
            partition = self._global_partition
        else:
            partition = self._config['object_types'][object_type]['area_partition']

        coll_type = self.cat_cxt.lookup_collection_type(object_type)
        if coll_type is not None:
            if self.cat_cxt.use_custom_load(object_type):
                out_list.append_collection(coll_type.load_collection(region,
                                                                     self,
                                                                     mjd=mjd))
                return out_list

        if partition != 'None':
            if partition['type'] == 'healpix':
                hps = self.get_hps_by_region(region, object_type)
                for hp in hps:
                    c = self.get_object_type_by_hp(hp, object_type,
                                                   region, mjd)
                    if len(c) > 0:
                        out_list.append_object_list(c)
                return out_list
        else:
            raise NotImplementedError(f'Unsupported object type {object_type}')

    def get_object_type_by_hp(self, hp, object_type, region=None, mjd=None):
        object_list = ObjectList()

        #  Do we need to check more specifically by object type?
        # if hp not in self._hp_info:
        if hp not in self.hps_by_type(object_type):
            msg = f'In SkyCatalog.get_object_type_by_hp, healpix {hp}  '
            msg += f'intersects region but has no data file for {object_type}'
            self._logger.warning(msg)
            return object_list

        if object_type in ['galaxy']:
            COLUMNS = ['galaxy_id', 'ra', 'dec']
            id_name = 'galaxy_id'
        elif object_type in ['snana']:
            COLUMNS = ['id', 'ra', 'dec', 'start_mjd', 'end_mjd']
            id_name = 'id'
        elif object_type in ['star', 'sncosmo']:
            COLUMNS = ['object_type', 'id', 'ra', 'dec']
            id_name = 'id'
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
            # ##f_list = self._hp_info[hp]['object_types'][self._config['object_types'][ot]['parent']]
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
                arrow_t = rdr.read_columns(COLUMNS, None, rg)
                if object_type in {'galaxy', 'snana'}:
                    ra_c, dec_c, id_c, mask =\
                        _compress_via_mask(arrow_t,
                                           id_name,
                                           region,
                                           source_type={object_type},
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
    # Get LSST bandpasses in case we need to compute fluxes
    _ = load_lsst_bandpasses()
    base_dir = os.path.dirname(config_file)
    YamlIncludeConstructor.add_to_loader_class(loader_class=yaml.SafeLoader,
                                               base_dir=base_dir)
    with open(config_file) as f:
        return SkyCatalog(yaml.safe_load(f), skycatalog_root=skycatalog_root,
                          mp=mp, verbose=verbose)


if __name__ == '__main__':
    import time
    cfg_file_name = 'skyCatalog.yaml'
    skycatalog_root = os.getenv('SKYCATALOG_ROOT')
    catalog_dir = 'reorg'
    if len(sys.argv) > 1:
        catalog_dir = sys.argv[1]
    if len(sys.argv) > 2:
        cfg_file_name = sys.argv[2]

    cfg_file = os.path.join(skycatalog_root, catalog_dir, cfg_file_name)

    write_sed = False
    if len(sys.argv) > 3:
        write_sed = True

    # For tract 3828
    #   55.73604 < ra < 57.563452
    #  -37.19001 < dec < -35.702481

    cat = open_catalog(cfg_file, skycatalog_root=skycatalog_root)
    hps = cat._find_all_hps()
    print('Found {} healpix pixels '.format(len(hps)))
    for h in hps:
        print(h)
    ra_min_tract = 55.736
    ra_max_tract = 57.564
    dec_min_tract = -37.190
    dec_max_tract = -35.702
    # ra_min_small = 56.0
    # ra_max_small = 56.2
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

    at_slac = os.getenv('HOME').startswith('/sdf/home/')
    if not at_slac:
        obj_types = {'star', 'galaxy', 'sncosmo', 'snana'}
    else:
        obj_types = {'star', 'galaxy', 'sncosmo', 'snana', 'gaia_star'}

    print('Invoke get_objects_by_region with box region, no gaia')
    t0 = time.time()
    object_list = cat.get_objects_by_region(rgn,
                                            obj_type_set={'star', 'galaxy',
                                                          'sncosmo'})
    t_done = time.time()
    print('Took ', t_done - t0)
    # #### temporary obj_type_set={'galaxy', 'star'} )
    #                                        obj_type_set=set(['galaxy']) )
    # Try out get_objects_by_hp with no region
    # colls = cat.get_objects_by_hp(9812, None, set(['galaxy']) )

    print('Number of collections returned for box:  ',
          object_list.collection_count)
    print('Object count for box: ', len(object_list))

    print('Invoke get_objects_by_region with polygonal region')
    t0 = time.time()
    object_list_poly = cat.get_objects_by_region(rgn_poly,
                                                 obj_type_set=obj_types)
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
    object_list_diamond = cat.get_objects_by_region(rgn_diamond,
                                                    obj_type_set=obj_types)
    t_done = time.time()
    print('Took ', t_done - t0)

    print('Number of collections returned for diamond:  ',
          object_list_diamond.collection_count)
    print('Object count for diamond: ', len(object_list_diamond))

    # ### TEMP FOR DEBUGGING
    # ## exit(0)

    # For now SIMS_SED_LIBRARY_DIR is undefined at SLAC, making it impossible
    # to get SEDs for stars. So (crudely) determine whether or not
    # we're running at SLAC

    colls = object_list.get_collections()
    got_a_sed = False
    for c in colls:
        n_obj = len(c)
        print("For hpid ", c.get_partition_id(), "found ", n_obj, " objects")
        if (n_obj) < 1:
            continue
        print("First object: ")
        print(c[0], '\nid=', c[0].id, ' ra=', c[0].ra, ' dec=', c[0].dec,
              ' belongs_index=', c[0]._belongs_index,
              ' object_type: ', c[0].object_type)

        if (n_obj < 3):
            continue
        print("Slice [1:3]")
        slice13 = c[1:3]
        for o in slice13:
            print('id=', o.id, ' ra=', o.ra, ' dec=', o.dec, ' belongs_index=',
                  o._belongs_index,  ' object_type: ', o.object_type)
            print(o.object_type)
            if o.object_type == 'star':
                if not at_slac:
                    print(o.get_instcat_entry())
                    sed, magnorm = o._get_sed()
                    print('For star magnorm: ', magnorm)
                    if magnorm < 1000:
                        print('Length of sed: ', len(sed.wave_list))
            elif o.object_type == 'sncosmo':
                print(o.get_instcat_entry())
            elif o.object_type == 'galaxy':
                for cmp in ['disk', 'bulge', 'knots']:
                    print(cmp)
                    if cmp in o.subcomponents:
                        # broken for galaxies currently
                        # print(o.get_instcat_entry(component=cmp))
                        sed, _ = o._get_sed(cmp)
                        if sed:
                            print('Length of sed table: ', len(sed.wave_list))
                            if not got_a_sed:
                                got_a_sed = True
                                th = o.get_native_attribute(f'sed_val_{cmp}')
                                print('Tophat values: ', th)
                                sed, _ = o._get_sed(component=cmp)
                                print('Simple sed wavelengths:')
                                print(sed.wave_list)
                                print('Simple sed values:')
                                print([sed(w) for w in sed.wave_list])
                                if write_sed:
                                    o.write_sed('simple_sed.txt',
                                                component=cmp)
                                sed_fine, _ = o._get_sed(component=cmp,
                                                         resolution=1.0)
                                print('Bin width = 1 nm')
                                print('Initial wl values',
                                      sed_fine.wave_list[:20])
                                print('Start at bin 100',
                                      sed_fine.wave_list[100:120])
                                print('Initial values')
                                print([sed_fine(w) for w in sed_fine.wave_list[:20]])
                                print('Start at bin 100')
                                print([sed_fine(w) for w in sed_fine.wave_list[100:120]])
                        else:
                            print('All-zero sed')

                # Try out old wrapper functions
                print("\nget_dust:")
                i_av, i_rv, g_av, g_rv = o._get_dust()
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
                for k, v in fluxes.items():
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
                print('id=', o.id, ' ra=', o.ra, ' dec=', o.dec,
                      ' belongs_index=', o._belongs_index)

    print('Total object count: ', len(object_list))

    if len(object_list) == 0:
        print('Empty object list. All done')
        exit(0)

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

    if len(object_list) < 5:
        print('Very short object list (< 5 elements).  done')
        exit(0)

    segment = object_list[2:5]
    print('Information for slice 2:5 ')
    for o in segment:
        print(f'object {o.id} of type {o.object_type} belongs to collection {o._belongs_to}')

    if len(object_list) < 304:
        print('Object list len < 304.  All done')
        exit(0)
    print('\nInformation for slice 285:300')
    segment = object_list[285:300]
    for o in segment:
        print(f'object {o.id} of type {o.object_type} belongs to collection {o._belongs_to}')

    # ixes = ([3,5,8],)
    ixes = (np.array([3, 5, 8, 300, 303]),)
    print(f'\nObjects with indexes {ixes[0]}')
    for o in object_list[ixes]:
        print(o.id)
    print('\nObjects in slice [3:9]')
    for o in object_list[3:9]:
        print(o.id)
    print('\nObjects in slice [300:304]')
    for o in object_list[300:304]:
        print(o.id)

    print('all done')
