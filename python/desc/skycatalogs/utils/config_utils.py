import yaml
from collections import namedtuple

__all__ = ['Config', 'open_config_file', 'Tophat']

def open_config_file(config_file):
    '''
    Given path to config file, return a Config object
    '''
    with open(config_file) as f:
        return Config(yaml.safe_load(f))

Tophat = namedtuple('Tophat', ['start', 'width'])


class Config(object):
    '''
    A wrapper around the dict which is the contents of a Sky Catalog
    config file which understands some of the semantics
    '''
    def __init__(self, cfg):
        self._cfg = cfg

    def list_sed_models(self):
        return self._cfg['SED_models'].keys()

    def list_object_types(self):
        return self._cfg['object_types'].keys()

    def get_sed_model(self, modelname):
        return self._cfg['SED_models'][modelname]

    def object_is_composite(self, objectname):
        return 'composite' in self._cfg['object_types'][objectname]

    def get_object_parent(self, objectname):
        if 'parent' in self._cfg['object_types'][objectname]:
            return self._cfg['object_type'][objectname]['parent']
        else:
            return None

    def get_object_sedmodel(self, objectname):
        if 'sed_model' in self._cfg['object_types'][objectname]:
            return self._cfg['object_type'][objectname]['sed_model']
        else:
            return None

    def get_tophat_parameters(self):
        '''
        Return list of named tuples
        Should maybe be part of Sky Catalogs API
        '''
        raw_bins = self._cfg['SED_models'][0]['tophat']['bins']

        return [ Tophat(b[0], b[1]) for b in raw_bins]
