import yaml
import jsonschema.validate
import os

from desc.skycatalogs.utils.exceptions import NoSchemaVersionError, ConfigDuplicateKeyError

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
        raw_bins = self._cfg['SED_models']['tophat']['bins']

        return [ Tophat(b[0], b[1]) for b in raw_bins]

    def get_config_value(self, key_path):
        '''
        Find value belonging to key in a config.
        Parameters
        ----------
        key_path    string   of form 'a/b/c/../q' where all but last
                             component must be a dict

        Returns
        -------
        Value associated with key_path if config contains such a path
        '''
        path_items = key_path.split('/')
        d = self._cfg
        for i in path_items[:-1]:
            if not i in d:
                raise ValueError(f'Item {i} not found')
            d = d[i]
            if not isinstance(d, dict):
                raise ValueError(f'intermediate {d} is not a dict')
        return d[path_items[-1]]

    def validate(self, schema=None):
        '''
        Parameters
        ----------
        schema    Identifies schema to validate against.
                  If string, open file with this file path
                  If dict, use "as is"
                  If None, attempt to deduce schema version (and hence
                  file path) from schema_version keyword
        Returns:  None
                  Raise exception if validation fails for any reason:
                  desc.skycatalogs.exception.NoSchema if schema specification
                  can't be found
                  OSError if schema file can't be read
                  yaml.YAMLerror if schema can't be loaded
                  jsonschema.exceptions.SchemaError if schema is invalid
                  jsonschema.exceptions.ValidationError otherwise
        '''
        fpath = None
        if schema is None:
            if 'schema_version' not in self._cfg:
                raise NoSchemaVersionError
            fpath = _find_schema_path(self._cfg["schema_version"])
            if not fpath:
                raise NoSchemaVersionError('Schema file not found')
        elif isinstance(schema, string):
            fpath = schema
        if fpath:
            try:
                f = open(fpath)
                sch = yaml.safe_load(f)
            except OSError as e:
                raise NoSchemaVersionError('Schema file not found or unreadable')
            except yaml.YAMLError as ye:
                raise ye
        if isinstance(schema, dict):
            sch = schema

        jsonschema.validate(self._cfg, schema.dict)

    def add_key(self, k, d):
        '''
        Parameters
        ----------
        k        Top-level key belonging to the schema
        v        value for the key
        '''

        if k in self._cfg[k]:
            raise ConfigDuplicateKeyError(k)
        self._cfg[k] = v

    def write_config(self, dirpath, filename=None):
        '''
        Export self to yaml document and write to specified directory.
        If filename is None the file will be named after catalog_name
        '''
        self.validate()

        if not filename:
            filename = self._cfg['catalog_name'] + '.yaml'




def _find_schema_path(schema_spec):
    '''
    Given a schema version specification, return the file path
    where the file describing it belongs
    '''
    fname = f'skycatalogs_config_{self._cfg["schema_spec"]}'
    here = os.path.dirname(__file__)
    return = os.path.join(here, '../../../../cfg', fname)

def create_config(catalog_name, schema_version):
    return Config({'catalog_name' : catalog_name,
                   'schema_version' : schema_version})
