import yaml
from collections import namedtuple
import healpy

__all__ = ['SkyCatalogs', 'open_catalog', 'Region']

Region = namedtuple('Region', ['ra_min', 'ra_max', 'dec_min', 'dec_max'])

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
    vec = healpy.pixelfunc.ang2vec([region.ra_min, region.ra_max,
                                    region.ra_max, region.ra_min],
                                   [region.dec_min, region.dec_min,
                                    region.dec_max, region.dec_max],
                                   lonlat=True)

    return healpy.query_polygon(nside, vec, inclusive=True)

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

        self._hps = self._find_all_hps()

    def _find_all_hps(self):
        # for each file matching pattern in the directory, make an
        # entry in self._hps
        #   self._hps[hpid] = {'relpath' : the_filename, gal_handle : None}
        # When file is open, set gal_handle to the Parquet file object
        # (or perhaps something else if underlying format is not Parquet)

        pass         # for now

    def get_hps_by_region(self, region):
        '''
        Region is a named 4-tuple (min-ra, max-ra, min-dec, max-dec)
        Catalog area partition must be by healpix
        '''
        # If area_partition doesn't use healpix, raise exception

        return _get_intersecting_hps(
            self._config['area_partition']['ordering'],
            self._config['area_partition']['nside'],
            region).intersetion(self._hps.keys())

    def get_object_type_names(self):
        return set(self._config['object_types'].keys())

    # Add more functions to return parts of config of possible interest
    # to user

    def get_objects_by_region(self, datetime, region, obj_type_set=None):
        '''
        Parameters
        ----------
        datetime       Python datetime object.
        region         min & max for ra and dec
        obj_type_set   Return only these objects. Defaults to all available

        Returns
        -------
        Collection of SkyObjects visible in the region at the specified time
        '''
        # Take intersection of obj_type_list and available object types
        # Determine healpix intersecting the region

        if self._config['area_partition']['type'] == 'healpix':
            hps = self.get_hps_by_region(region)

        # otherwise raise a not-supported exception

        objects = []
        if obj_type_set is None:
            obj_types = self.get_object_type_names()
        else:
            obj_types = self.get_object_type_names().intersection(obj_types)
        for hp in hps:
            # Maybe have a multiprocessing switch? Run-time option when
            # catalog is opened?
            objects = objects + self.get_objects_by_healpix(datetime, hp,
                                                            obj_types)

        return objects

    def get_object_iterator_by_region(self, datetime, region,
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

    def get_objects_by_hp(self, datetime, region, hp, obj_type_set=None):
        # Determine if healpix pixel is entirely contained within region.
        # If so we don't have to check whether individual object are
        # in the region or not.

        # Get file handle(s) and store if we don't already have it (them)
        # Read into dataframe
        # transpose;
        # create object for each row, plus galaxy object with subcomponents
        # for each duple or triple with same object id
        # return collection of objects (and cache under self._hps[hp])




        pass     # implementation to be filled in later

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
