import os
import sys
import yaml
import json
import pyarrow.parquet as pq
import logging
from typing import Any
from .exceptions import ConfigDuplicateKeyError
from collections import namedtuple

__all__ = ['Config', 'open_config_file', 'Tophat',
           'get_file_metadata', 'CURRENT_SCHEMA_VERSION',
           'YamlPassthruIncludeLoader']

CURRENT_SCHEMA_VERSION = '1.3.0'


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
        self._logger = logging.getLogger('YamlIncludeLoader')
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


class YamlPassthruIncludeLoader(yaml.SafeLoader):

    """YAML Loader that just returns scalar nodes like
           !include the_path
       as is, without attempting to include the contents of the_path


    Examples
    --------
    >>> with open("document.yaml", "r") as f:
           data = yaml.load(f, Loader=YamlPassthruIncludeLoader)

    Parameters
    ----------
    stream :  text io stream
        The stream to parse.

    Modified version of YamlIncludeLoader above.
    """
    def __init__(self, filestream):
        super().__init__(filestream)
        self._logger = logging.getLogger('YamlPassthruIncludeLoader')
        self._current_dir = os.path.dirname(filestream.name)
        self.add_constructor("!include", YamlPassthruIncludeLoader.include)

    def include(self, node: yaml.Node) -> list[Any] | dict[str, Any]:
        result: list[Any] | dict[str, Any]
        if isinstance(node, yaml.ScalarNode):
            # The following returns only what we need.
            return node.tag + ' ' + node.value

        elif isinstance(node, yaml.SequenceNode):
            result = []
            for filename in self.construct_sequence(node):
                result.append(filename)
            return result

        elif isinstance(node, yaml.MappingNode):
            result = {}
            for k, v in self.construct_mapping(node).items():
                if not isinstance(k, str):
                    raise TypeError(f"Expected only strings in YAML mapping; got {k!r} of type {type(k)}.")
                result[k] = str(v)              # self.extractFile(v)
            return result

        else:
            self._logger.error("Unrecognised node type in !include statement",
                               file=sys.stderr)
            raise yaml.constructor.ConstructorError


def open_config_file(config_file):
    '''
    Given path to config file, return a Config object. See Config class below
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
    config file (a dict) which understands some of the dict semantics

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

        # cache schema version
        self._schema_version = None
        if 'schema_version' in cfg:
            self._schema_version = cfg['schema_version']
        else:
            self._schema_version = self.get_config_value('provenance/versioning/schema_version', silent=True)

        self._schema_pre_130 = True
        if self._schema_version:
            self._cmps = [int(c) for c in self._schema_version.split('.')]
            if (self._cmps[0] > 1) or (self._cmps[0] == 1 and self._cmps[1] > 2):
                self._schema_pre_130 = False

    def __getitem__(self, k):
        '''
        Specific override of __getitem__, delegating to dict member
        so that [ ] syntax will be handled properly
        '''
        return self._cfg.__getitem__(k)

    def __contains__(self, k):
        return k in self._cfg

    @property
    def schema_version(self):
        return self._schema_version

    def schema_pre_130(self):
        ''' If schema version is older than 1.3.0, structure is different'''
        return self._schema_pre_130

    def list_object_types(self):
        return self._cfg['object_types'].keys()

    def object_is_composite(self, objectname):
        return 'composite' in self._cfg['object_types'][objectname]

    def get_tophat_parameters(self, object_type='galaxy'):
        '''
        Return list of named tuples
        Should maybe be part of Sky Catalogs API
        '''
        if self._schema_pre_130:
            tophat_path = 'SED_models/tophat'
        else:
            tophat_path = f'object_types/{object_type}/tophat'

        tophat = self.get_config_value(tophat_path, silent=True)

        if not tophat:
            return None
        raw_bins = tophat['bins']

        return [Tophat(b[0], b[1]) for b in raw_bins]

    def get_cosmology(self, object_type=None):
        '''
        Return cosmology parameters.  Location in config will depend
        on schema version and object type

        Parameters
        ----------
        schema_version   string or None.   Of form x.y.z
        object_type      string or None    If object type specified, use it.
                                           Else presume only one galaxy-like
                                           object type present
        Return          dict or None
        '''
        if self._schema_pre_130:
            return self._cfg['Cosmology']

        if not object_type:
            if 'galaxy' in self._cfg['object_types']:
                object_type = 'galaxy'
            elif 'diffsky_galaxy' in self._cfg['object_types']:
                object_type = 'diffsky_galaxy'
            elif 'skysim5000' in self._cfg['object_types']:
                object_type = 'skysim5000'
            else:
                return None
        return self._cfg['object_types'][object_type]['Cosmology']

    def get_sso_sed_path(self):
        if self._schema_pre_130:
            return self._cfg['provenance']['inputs'].get('sso_sed',
                                                         'sso_sed.db')
        else:
            return self._cfg['object_types']['sso']['provenance']['inputs'].get('sso_sed', 'sso_sed.db')

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


def get_file_metadata(fpath, key='provenance'):
    '''
    Extract metadata from a parquet catalog file. Entries in the file
    metadata are key-value pairs.  Return a dict with the specified
    key as its only key
    '''
    meta = pq.read_table(fpath, columns=[]).schema.metadata
    to_return = dict()

    if key.encode() in meta:
        to_return[key] = json.loads(meta[key.encode()])
    return to_return
