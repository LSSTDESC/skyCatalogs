import yaml

__all__ = ['SkyCatalogs', open_catalog]


class SkyCatalogs(object):
    '''

    '''
    def __init__(self, config):
        '''
        Parameters
        ----------
        config:  dict.  Typically the result of loading a yaml file
        '''
        self._config = config
        # There may be more to do at this point but not too much.
        # In particular, don't read in anything from data files
        # One might check that the config is complete and well-formed
        #  - for example require certain keys, such as catalog_name,
        #  data_file_type, area_parition, root_directory, object_types -
        # to exist, and check that that the data directory (value of
        # root_directory) exists.
        # create an empty dict for
        # per-HEALpix pixel information and so forth.
        self._hps = {}

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
            hps = self.get_intersecting_hps(region)   # returns an iterable
        # otherwise raise a not-supported exception

        objects = []
        if obj_type_set is None:
            obj_types = self.get_object_type_names()
        else:
            obj_types = self.get_object_type_names().intersection(obj_types)
        for hp in hps:
            objects = objects + self.get_objects_by_healpix(datetime, hp,
                                                            obj_types)

        return objects

    def get_objects_by_healpix(self, datetime, hp, obj_type_set):
        # Determine if healpix pixel is entirely contained within region.
        # If so we don't have to check whether individual object are
        # in the region or not.






        pass     # implementation to be filled in later


def open_catalog(config_file):
    '''
    Parameters
    ----------
    yaml file containing config

    Returns
    -------
    SkyCatalog
    '''
    with open(config_file) as f:
        return SkyCatalog(yaml.safe_load(f))
