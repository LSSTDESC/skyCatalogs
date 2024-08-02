import yaml
import os
from pathlib import PurePath

'''
Per-object type classes which know how to write a suitable config
fragment for a particular object type
'''

__all__ = ['BaseConfigFragment']
_FILE_PATH = str(PurePath(__file__))
_SKYCATALOGS_DIR = _FILE_PATH[:_FILE_PATH.rindex('/skycatalogs')]

_TEMPLATE_DIR = os.path.join(_SKYCATALOGS_DIR, 'skycatalogs', 'data',
                             'cfg_templates')


class BaseConfigFragment():
    def __init__(self, prov, object_type_name=None, template_name=None,
                 area_partition=None, data_file_type=None):
        self._object_type_name = object_type_name
        self._prov = prov
        self._opt_dict = dict()
        if area_partition:
            self._opt_dict['area_partition'] = area_partition
        if data_file_type:
            self._opt_dict['data_file_type'] = data_file_type
        self._template_name = template_name
        if not template_name:
            if object_type_name:
                self._template_name = object_type_name + '_template.yaml'
        if object_type_name:
            self._frag_name = self._object_type_name + '.yaml'
        else:
            self._frag_name = None

    @property
    def object_type(self):
        return self._object_type_name

    @property
    def fragment_name(self):
        return self._frag_name

    def make_fragment(self):
        '''
        Returns a dict which, when dumped to yaml, will be the config fragment
        for the object type.
        Must be implemented by subclass
        '''
        return self.generic_create()

    def generic_create(self):
        template_path = os.path.join(_TEMPLATE_DIR, self._template_name)
        with open(template_path, 'r') as f:
            data = yaml.load(f, Loader=yaml.SafeLoader)

        opt = self._opt_dict
        other = dict()
        for key in opt:
            if opt[key] is not None:
                other[key] = opt[key]
        data.update(other)

        data['provenance'] = self._prov
        return data
