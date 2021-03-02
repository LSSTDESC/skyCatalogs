import yaml

__all__ = ['SkyCatalogs', open_catalog]


class SkyCatalogs(object):
    '''

    '''
    def __init__(self, config):
        self._config = config
        # There may be more to do at this point but not too much.
        # In particular, don't read in anything from data files
        # One might check that the config is complete and well-formed
        # and that the data directory exists, create an empty dict for
        # per-HEALpix pixel information and so forth.


    def get_objects_by_region(datetime, region, obj_type_list=None):
        '''
        Parameters
        ----------
        datetime       Python datetime object.
        region         min & max for ra and dec
        obj_type_list  Return only these objects. Defaults to all available

        Returns
        -------
        Collection of SkyObjects visible in the region at the specified time
        '''
        # Take intersection of obj_type_list and available object types
        # Determine healpix intersecting the region

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
