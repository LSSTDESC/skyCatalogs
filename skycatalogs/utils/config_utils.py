import os
import sys
import yaml
import logging
from typing import Any
from .exceptions import ConfigDuplicateKeyError
from collections import namedtuple

__all__ = ['Config', 'open_config_file', 'Tophat', 'create_config',
           'assemble_SED_models', 'assemble_MW_extinction',
           'assemble_cosmology', 'assemble_object_types',
           'assemble_provenance', 'assemble_variability_models', 'write_yaml',
           'CURRENT_SCHEMA_VERSION']

CURRENT_SCHEMA_VERSION = '1.2.0'


class YamlIncludeLoader(yaml.SafeLoader):

    """YAML Loader that supports file include directives.

    Uses ``!include`` directive in a YAML file to point to another
    YAML file to be included. The path in the include directive is relative
    to the top-level file

        storageClasses: !include storageClasses.yaml

    Examples
    --------
    >>> with open("document.yaml", "r") as f:
           data = yaml.load(f, Loader=YamlIncludeLoader)

    Parameters
    ----------
    stream :  text io stream
        The stream to parse.

    This code was adapted from the LSST Science Pipelines Butler.
    See in particular the Loader class in
    daf_butler/python/lsst/daf/butler/_config.py in the daf_butler repo
    https://github.com/lsst/daf_butler
    """
    def __init__(self, filestream):
        super().__init__(filestream)
        self._logger =  logging.getLogger('YamlIncludeLoader')
        self._current_dir = os.path.dirname(filestream.name)
        self.add_constructor("!include", YamlIncludeLoader.include)

    def include(self, node: yaml.Node) -> list[Any] | dict[str, Any]:
        result: list[Any] | dict[str, Any]
        if isinstance(node, yaml.ScalarNode):
            return self.extractFile(self.construct_scalar(node))  # type: ignore[arg-type]

        elif isinstance(node, yaml.SequenceNode):
            result = []
            for filename in self.construct_sequence(node):
                result.append(self.extractFile(filename))
            return result

        elif isinstance(node, yaml.MappingNode):
            result = {}
            for k, v in self.construct_mapping(node).items():
                if not isinstance(k, str):
                    raise TypeError(f"Expected only strings in YAML mapping; got {k!r} of type {type(k)}.")
                result[k] = self.extractFile(v)
            return result

        else:
            self._logger.error("Unrecognised node type in !include statement",
                  file=sys.stderr)
            raise yaml.constructor.ConstructorError

    def extractFile(self, filepath: str) -> Any:
        if os.path.isabs(filepath):
            actual_path = filepath
        else:
            actual_path = os.path.join(self._current_dir, filepath)
        self._logger.info("Opening YAML file via !include: %s", actual_path)

        # Read all the data from the resource
        with open(actual_path) as f:
            content = yaml.load(f, YamlIncludeLoader)
        return content


def open_config_file(config_file):
    '''
    Given path to config file, return a Config object
    '''
    with open(config_file) as f:
        content = yaml.load(f, Loader=YamlIncludeLoader)
        return Config(content)


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

    def __contains__(self, k):
        return k in self._cfg

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

        return [Tophat(b[0], b[1]) for b in raw_bins]

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
            if i not in d.keys():
                if silent:
                    return None
                raise ValueError(f'Item {i} not found')
            d = d[i]
            if not isinstance(d, dict):
                raise ValueError(f'intermediate {d} is not a dict')

        if path_items[-1] in d:
            return d[path_items[-1]]
        else:
            if silent:
                return None
            raise ValueError(f'Item {i} not found')

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
        # ##self.validate()   skip for now

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
            txt = 'Config.write_yaml: Will not overwrite pre-existing config '
            if logname:
                logger = logging.getLogger(logname)
                logger.warning(txt + outpath)
            else:
                print(txt + outpath)
            return
    else:
        with open(outpath, mode='w') as f:
            yaml.dump(input_dict, f)

    return outpath


def create_config(catalog_name, logname=None):
    return Config({'catalog_name': catalog_name}, logname)


def assemble_cosmology(cosmology):
    d = {k: cosmology.__getattribute__(k) for k in ('Om0', 'Ob0', 'sigma8',
                                                    'n_s')}
    d['H0'] = float(cosmology.H0.value)
    return d


def assemble_MW_extinction():
    av = {'mode': 'data'}
    rv = {'mode': 'constant', 'value': 3.1}
    return {'r_v': rv, 'a_v': av}


def assemble_object_types(pkg_root, galaxy_nside=32):
    '''
    Include all supported object types even though a particular catalog
    might not use them all
    '''
    t_path = os.path.join(pkg_root, 'cfg', 'object_types.yaml')
    with open(t_path) as f:
        d = yaml.safe_load(f)
        d['object_types']['galaxy']['area_partition']['nside'] = galaxy_nside
        return d['object_types']


def assemble_SED_models(bins):
    tophat_d = {'units': 'angstrom', 'bin_parameters': ['start', 'width']}
    tophat_d['bins'] = bins
    file_nm_d = {'units': 'nm'}
    return {'tophat': tophat_d, 'file_nm': file_nm_d}


def assemble_provenance(pkg_root, inputs={}, schema_version=None):
    import git

    if not schema_version:
        schema_version = CURRENT_SCHEMA_VERSION
    import skycatalogs
    version_d = {'schema_version': schema_version}
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
        return {'versioning': version_d, 'skyCatalogs_repo': git_d,
                'inputs': inputs}
    else:
        return{'versioning': version_d, 'skyCatalogs_repo': git_d}


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
