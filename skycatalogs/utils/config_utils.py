import os
import sys
import yaml
import json
import pyarrow.parquet as pq
import logging
from typing import Any
from .exceptions import ConfigDuplicateKeyError
from collections import namedtuple

__all__ = ['Config', 'open_config_file', 'Tophat', 'create_config',
           'assemble_MW_extinction', 'assemble_cosmology',
           'assemble_provenance',
           'get_file_metadata', 'CURRENT_SCHEMA_VERSION',
           'YamlPassthruIncludeLoader', 'ConfigWriter']

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

    def get_tophat_parameters(self):
        '''
        Return list of named tuples
        Should maybe be part of Sky Catalogs API
        '''
        if self._schema_pre_130:
            tophat_path = 'SED_models/tophat'
        else:
            tophat_path = 'object_types/galaxy/tophat'

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


def create_config(catalog_name, logname=None):
    return Config({'catalog_name': catalog_name}, logname)


# A collection of utilities used by CatalogCreator to assemble and write
# configs
def assemble_cosmology(cosmology):
    d = {k: cosmology.__getattribute__(k) for k in ('Om0', 'Ob0', 'sigma8',
                                                    'n_s') if k in dir(cosmology)}
    d['H0'] = float(cosmology.H0.value)
    return d


def assemble_MW_extinction():
    av = {'mode': 'data'}
    rv = {'mode': 'constant', 'value': 3.1}
    return {'r_v': rv, 'a_v': av}


def assemble_provenance(pkg_root, inputs={}, run_options=None,
                        schema_version=None):
    '''
    Assemble provenance information, usually pertaining to a sinlge
    object type

    Parameters
    ----------
    pkg_root  string
        top directory of the git package
    inputs    dict
        For a name like 'star_truth' the path to the corresponding input file.
    run_options dict or None
         Options the script create_sc.py was called with
    schema_verion string or None
         If None (usual case) use current schema version

    Return
    ------
    dict
    '''
    import skycatalogs
    try:
        import git
        have_git = True
    except ImportError:
        have_git = False

    if not schema_version:
        schema_version = CURRENT_SCHEMA_VERSION
    version_d = {'schema_version': schema_version}
    if '__version__' in dir(skycatalogs):
        code_version = skycatalogs.__version__
    else:
        code_version = 'unknown'
    version_d['code_version'] = code_version

    repo = git.Repo(pkg_root)
    has_uncommited = repo.is_dirty()
    has_untracked = (len(repo.untracked_files) > 0)

    if have_git:
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

        to_return = dict()
        to_return['versioning'] = version_d
        to_return['skyCatalogs_repo'] = git_d

    if inputs:
        to_return['inputs'] = inputs
    if run_options:
        to_return['run_options'] = run_options

    return to_return


def assemble_file_metadata(pkg_root, inputs=None, run_options=None,
                           flux_file=False, throughputs_versions=None):
    '''
    Assemble the metadata to be included in a skyCatalogs binary data file
    '''
    to_return = assemble_provenance(pkg_root, inputs=inputs,
                                    run_options=run_options)
    if flux_file:
        # add a section 'flux_dependencies' containing at least
        # galsim version and, if possible, throughputs version
        to_return['flux_dependencies'] = dict()
        from galsim import version as galsim_version
        to_return['flux_dependencies']['galsim_version'] = galsim_version
        if throughputs_versions:
            for k, v in throughputs_versions.items():
                to_return['flux_dependencies'][k] = v

    return to_return


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


def _read_yaml(inpath, silent=True, resolve_include=True):
    '''
    Parameters
    ----------
    inpath     string            path to file
    silent     boolean           if file not found and silent is true,
                                 return None.  Else raise exception
    resolve_include boolean      if False, return values like
                                 !include star.yaml
                                 literally. If True, replace with content
                                 of references file

    Returns
    -------
    dict representing contents of file, or None if silent and file not found
    '''
    if resolve_include:
        ldr = YamlIncludeLoader
    else:
        ldr = YamlPassthruIncludeLoader
    try:
        with open(inpath, mode='r') as f:
            data = yaml.load(f, Loader=ldr)
            return data
    except FileNotFoundError as ex:
        if silent:
            return None
        else:
            raise ex


class ConfigWriter:
    '''
    Saves a little context needed by utilities
    '''
    def __init__(self, skycatalog_root, catalog_dir, top_name,
                 overwrite=False, logname=None):
        self._skycatalog_root = skycatalog_root
        self._catalog_dir = catalog_dir
        self._out_dir = os.path.join(skycatalog_root, catalog_dir)
        self._top_name = top_name
        self._overwrite = overwrite
        if not logname:
            logname = 'skyCatalogs:ConfigUtils'
        self._logger = logging.getLogger(logname)

    def write_yaml(self, input_dict, outpath):
        '''
        Write yaml file if
          * it doesn't already exist     or
          * we're allowed to overwrite

        Parameters
        ----------
        input_dict    dict   Contents to be output to yaml
        outpath       string Where to write the output

        Returns
        -------
        output path (same as argument) if a file is written, else None
        '''
        if self._overwrite:
            return self.update_yaml(input_dict, outpath)

        try:
            with open(outpath, mode='x') as f:
                yaml.dump(input_dict, f)
        except FileExistsError:
            txt = 'write_yaml: Will not overwrite pre-existing config'
            self._logger.warning(txt + outpath)
            return None

        return outpath

    def update_yaml(self, input_dict, outpath):
        '''
        Write yaml regardless of whether file is present or not
        '''
        with open(outpath, mode='w') as f:
            yaml.dump(input_dict, f)
        return outpath

    def find_schema_version(self, top):
        if 'schema_version' not in top.keys():
            return None
        return top['schema_version']

    def schema_compatible(self, schema_version):
        # For now just require that major versions match
        # Schema designations are of the form M.m.p where M, m  and p are
        # integers denoting major, minor and patch version
        current = CURRENT_SCHEMA_VERSION.split('.')
        incoming = schema_version.split('.')
        return current[0] == incoming[0]

    def write_configs(self, config_fragment):
        '''
        Write yaml fragment for specified object type and write or update
        top yaml file if
        * fragment doesn't already exist
        * or overwrite is true

        Parameters
        ----------
        config_fragment   BaseConfigFragment   Knows how to create
                                               config_fragment for a particular
                                               object type
        '''

        # Need the following machinery (class IncludeValue; routine
        # include_representer) in order to output values like
        #       !include  some_path.yaml
        # properly
        class IncludeValue(str):
            def __new__(cls, a):
                return str.__new__(cls, a)

            def __repr__(self):
                return "IncludeValue(%s)" % self

        def include_representer(dumper, value):
            # To avoid outputting any quotes, use style='|'
            return dumper.represent_scalar(u'!include', u'%s' % value,
                                           style='|')

        overwrite = self._overwrite
        top_path = os.path.join(self._out_dir, self._top_name + '.yaml')

        top = _read_yaml(top_path, silent=True, resolve_include=False)
        if top:
            top_exists = True
            object_type_exists = config_fragment.object_type in top['object_types']
            schema_version = self.find_schema_version(top)
        else:
            top_exists = False
            object_type_exists = False
            schema_version = None

        if top_exists:
            if not self.schema_compatible(schema_version):
                if not overwrite:
                    raise RuntimeError('Incompatible skyCatalogs config versions')
                else:
                    self._logger.warning('Overwriting config with incompatible schema version')
            if object_type_exists and not overwrite:
                return

        # Generate and write fragment for the object type
        fragment_name = config_fragment.fragment_name
        frag = config_fragment.make_fragment()
        frag_path = os.path.join(self._out_dir, fragment_name)
        self.write_yaml(frag, frag_path)

        # Write or update top file if necessary
        object_type = config_fragment.object_type
        value = IncludeValue(fragment_name)
        yaml.add_representer(IncludeValue, include_representer)
        if top_exists and not overwrite:
            if object_type_exists and top['object_types'][object_type] == value:
                # No change necessary
                return

            # Otherwise need to add or modify value for our object type
            # First have to fix values for any other object types already
            # mentions.  Value read in looks like "!include an_obj_type.yaml"
            for k, v in top['object_types'].items():
                cmps = v.split(' ')
                new_value = IncludeValue(cmps[1])
                top['object_types'][k] = new_value

            top['object_types'][object_type] = value
            self.update_yaml(top, top_path)
            return

        # Write out top file from scratch, ignoring other object types
        # which may have been referenced in older version
        top_dict = {'skycatalog_root': self._skycatalog_root,
                    'catalog_dir': self._catalog_dir,
                    'catalog_name': self._top_name,
                    'schema_version': CURRENT_SCHEMA_VERSION}
        objs = {object_type: value}
        top_dict['object_types'] = objs
        self.write_yaml(top_dict, top_path)
