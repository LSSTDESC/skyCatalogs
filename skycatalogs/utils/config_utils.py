import os
import yaml
import git
import logging
from jsonschema import validate

from .exceptions import NoSchemaVersionError, ConfigDuplicateKeyError

from collections import namedtuple

__all__ = ['Config', 'open_config_file', 'Tophat', 'create_config',
           'assemble_SED_models', 'assemble_MW_extinction',
           'assemble_cosmology', 'assemble_object_types', 'assemble_provenance',
           'assemble_variability_models', 'write_yaml', 'CURRENT_SCHEMA_VERSION']

CURRENT_SCHEMA_VERSION='1.2.0'

def open_config_file(config_file):
    '''
    Given path to config file, return a Config object
    '''
    with open(config_file) as f:
        return Config(yaml.safe_load(f))

Tophat = namedtuple('Tophat', ['start', 'width'])


def _custom_dir(c, add):
    '''
    dir should return
       * functions defined for c's class
       * instance variables belonging to c
       * public stuff from the contained object (some table-like thing)
    '''
    return dir(type(c)) + list(c.__dict__.keys()) + add


class DelegatorBase:
    '''
    This class allows derived classes to delegate operations to a
    member.   Here it is used to delegate dict operations to
    the member self._cfg
    '''

    @property
    def _delegate(self):
        pub = [o for o in dir(self.default) if not o.startswith('_')]
        return pub
    def __getattr__(self, k):
        if k in self._delegate:
            return getattr(self.default, k)
        raise AttributeError(k)
    def __dir__(self):
        return _custom_dir(self, self._delegate)

    def __init__(self):
        pass

class Config(DelegatorBase):
    '''
    A wrapper around the dict which is the contents of a Sky Catalog
    config file (a dict) which understands some of the semantics

    '''
    def __init__(self, cfg, logname=None):
        '''
        Parameters
        ----------
        cfg  is normally a dict, but it can also itself be a Config.
        In this case, store its dict to delegate to.
        '''
        if isinstance(cfg, Config):
            cfg = cfg._cfg

        self._cfg = cfg
        self._logname = logname

        super().__init__()

        self.default = cfg

    def __getitem__(self, k):
        '''
        Specific override of __getitem__, delegating to dict member
        so that [ ] syntax will be handled properly
        '''
        return self._cfg.__getitem__(k)

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
            return self._cfg['object_types'][objectname]['parent']
        else:
            return None

    def get_object_sedmodel(self, objectname):
        if 'sed_model' in self._cfg['object_types'][objectname]:
            return self._cfg['object_types'][objectname]['sed_model']
        else:
            return None

    def get_tophat_parameters(self):
        '''
        Return list of named tuples
        Should maybe be part of Sky Catalogs API
        '''
        if not self.get_config_value('SED_models/tophat', silent=True):
            return None
        raw_bins = self._cfg['SED_models']['tophat']['bins']

        return [ Tophat(b[0], b[1]) for b in raw_bins]

    def get_config_value(self, key_path, silent=False):
        '''
        Find value belonging to key in a config.
        Parameters
        ----------
        key_path    string   of form 'a/b/c/../q' where all but last
                             component must be a dict
        silent      boolean  If False (default) complain when key is not
                             found.  Else return None

        Returns
        -------
        Value associated with key_path if config contains such a path
        '''
        path_items = key_path.split('/')
        d = self._cfg
        for i in path_items[:-1]:
            if not i in d.keys():
                if silent:
                    return None
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
                  skycatalogs.exception.NoSchema if schema specification
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

    def add_key(self, k, v):
        '''
        Parameters
        ----------
        k        Top-level key belonging to the schema
        v        value for the key
        '''

        if k in self._cfg:
            raise ConfigDuplicateKeyError(k)
        self._cfg[k] = v

    def write_config(self, dirpath, filename=None, overwrite=False):
        '''
        Export self to yaml document and write to specified directory.
        If filename is None the file will be named after catalog_name

        Parameters
        ----------
        dirpath        Directory to which file will be written
        filename       If supplied, use for filename
        overwrite      By default do not overwrite existing file with same path

        Return
        ------
        Full path of output config
        '''
        ###self.validate()   skip for now

        if not filename:
            filename = self._cfg['catalog_name'] + '.yaml'

        return write_yaml(self._cfg, os.path.join(dirpath, filename),
                          overwrite=overwrite, logname=self._logname)


def write_yaml(input_dict, outpath, overwrite=False, logname=None):
        if not overwrite:
            try:
                with open(outpath, mode='x') as f:
                    yaml.dump(input_dict, f)
            except FileExistsError:
                if logname:
                    logger = logging.getLogger(logname)
                    logger.warning('Config.write_yaml: Will not overwrite pre-existing config file ' + outpath)
                else:
                    print('Config.write_yaml: Will not overwrite pre-existing config file ' + outpath)
                return

        else:
            with open(outpath, mode='w') as f:
                yaml.dump(input_dict, f)

        return outpath


def _find_schema_path(schema_spec):
    '''
    Given a schema version specification, return the file path
    where the file describing it belongs
    '''
    fname = f'skycatalogs_schema_{self._cfg["schema_spec"]}'
    here = os.path.dirname(__file__)
    return os.path.join(here, '../../cfg', fname)

def create_config(catalog_name, logname=None):
    return Config({'catalog_name' : catalog_name}, logname)
#                  'schema_version' : schema_version,
#                  'code_version' : desc.skycatalogs.__version__}, logname)

def assemble_cosmology(cosmology):
    d = {k : cosmology.__getattribute__(k) for k in ('Om0', 'Ob0', 'sigma8',
                                                     'n_s')}
    d['H0'] = float(cosmology.H0.value)
    return d

def assemble_MW_extinction():
    av = {'mode' : 'data'}
    rv = {'mode' : 'constant', 'value' : 3.1}
    return {'r_v' : rv, 'a_v' : av}

def assemble_object_types(pkg_root):
    '''
    Include all supported object types even though a particular catalog
    might not use them all
    '''
    t_path = os.path.join(pkg_root, 'cfg', 'object_types.yaml')
    with open(t_path) as f:
        d = yaml.safe_load(f)
        return d['object_types']

def assemble_SED_models(bins):
    tophat_d = { 'units' : 'angstrom', 'bin_parameters' : ['start', 'width']}
    tophat_d['bins'] = bins
    file_nm_d = {'units' : 'nm'}
    return {'tophat' : tophat_d, 'file_nm' : file_nm_d}

def assemble_provenance(pkg_root, inputs={}, schema_version=None):

    if not schema_version:
        schema_version = CURRENT_SCHEMA_VERSION
    import skycatalogs
    version_d = {'schema_version' : schema_version}
    if '__version__' in dir(skycatalogs):
        code_version = skycatalogs.__version__
    else:
        code_version = 'unknown'
    version_d['code_version'] = code_version

    repo = git.Repo(pkg_root)
    has_uncommited = repo.is_dirty()
    has_untracked = (len(repo.untracked_files) > 0)


    git_d = {}
    git_d['git_hash'] = repo.commit().hexsha
    git_d['git_branch'] = repo.active_branch.name
    status = []
    if has_uncommited:
        status.append('UNCOMMITTED_FILES')
    if has_untracked:
        status.append('UNTRACKED_FILES')
    if len(status) == 0:
        status.append('CLEAN')
    git_d['git_status'] = status

    if inputs:
        return {'versioning' : version_d,'skyCatalogs_repo' : git_d,
                'inputs' : inputs}
    else:
        return{'versioning' : version_d, 'skyCatalogs_repo' : git_d}

# In config just keep track of models by object type. Information
# about the parameters they require is internal to the code.
_AGN_MODELS = ['agn_random_walk']
_SNCOSMO_MODELS = ['sn_salt2_extended']

def assemble_variability_models(object_types):
    '''
    Add information about all known variability models for supplied object
    types.
    Parameters
    ----------
    object_types: iterable of (pointsource) object type names

    Returns
    -------
    A dict with object type names for keys.  Values are also dicts with
    keys = model name and values defining struct of parameters for that
    model, e.g. ordered dict with parameter names for keys and parameter
    data types for values
    '''
    models = dict()
    if 'agn' in object_types:
        models['agn'] = _AGN_MODELS
    if 'sncosmo' in object_types:
        models['sncosmo'] = _SNCOSMO_MODELS

    return models
