import os
import sys
import yaml
import json
import pyarrow.parquet as pq
import logging
from typing import Any
from .exceptions import ConfigDuplicateKeyError
from collections import namedtuple
from .source_config import create_galaxy_config
from skycatalogs.utils.source_config import create_diffsky_galaxy_config
from .source_config import create_star_config
from .source_config import create_sso_config
from .source_config import create_snana_config
from .source_config import create_gaia_star_butler_config
from .source_config import create_gaia_star_direct_config

__all__ = ['Config', 'open_config_file', 'Tophat', 'create_config',
           'assemble_SED_models', 'assemble_MW_extinction',
           'assemble_cosmology',       # 'assemble_object_types',
           'assemble_provenance',   # 'write_yaml',
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

    """YAML Loader that just returns !include value as is as a string


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
        self._logger = logging.getLogger('YamlPassthruIncludeLoader')
        self._current_dir = os.path.dirname(filestream.name)
        self.add_constructor("!include", YamlPassthruIncludeLoader.include)

    def include(self, node: yaml.Node) -> list[Any] | dict[str, Any]:
        result: list[Any] | dict[str, Any]
        if isinstance(node, yaml.ScalarNode):
            return str(node)
            # return self.extractFile(self.construct_scalar(node))  # type: ignore[arg-type]

        elif isinstance(node, yaml.SequenceNode):
            result = []
            for filename in self.construct_sequence(node):
                # result.append(self.extractFile(filename))
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

    # def write_config(self, dirpath, filename=None, overwrite=False):
    #     '''
    #     Export self to yaml document and write to specified directory.
    #     If filename is None the file will be named after catalog_name

    #     Parameters
    #     ----------
    #     dirpath        Directory to which file will be written
    #     filename       If supplied, use for filename
    #     overwrite      By default do not overwrite existing file with same path

    #     Return
    #     ------
    #     Full path of output config
    #     '''
    #     # ##self.validate()   skip for now

    #     if not filename:
    #         filename = self._cfg['catalog_name'] + '.yaml'

    #     return write_yaml(self._cfg, os.path.join(dirpath, filename),
    #                       overwrite=overwrite, logname=self._logname)


def create_config(catalog_name, logname=None):
    return Config({'catalog_name': catalog_name}, logname)


def assemble_cosmology(cosmology):
    d = {k: cosmology.__getattribute__(k) for k in ('Om0', 'Ob0', 'sigma8',
                                                    'n_s') if k in dir(cosmology)}
    d['H0'] = float(cosmology.H0.value)
    return d


def assemble_MW_extinction():
    av = {'mode': 'data'}
    rv = {'mode': 'constant', 'value': 3.1}
    return {'r_v': rv, 'a_v': av}


# # won't need this
# def assemble_object_types(pkg_root, galaxy_nside=32):
#     '''
#     Include all supported object types even though a particular catalog
#     might not use them all
#     '''
#     t_path = os.path.join(pkg_root, 'cfg', 'object_types.yaml')
#     with open(t_path) as f:
#         d = yaml.safe_load(f)
#         d['object_types']['galaxy']['area_partition']['nside'] = galaxy_nside
#         return d['object_types']

# this will change or be eliminated
def assemble_SED_models(bins):
    to_return = dict()
    file_nm_d = {'units': 'nm'}
    to_return['file_nm'] = file_nm_d
    if bins:
        tophat_d = {'units': 'angstrom', 'bin_parameters': ['start', 'width']}
        tophat_d['bins'] = bins
        to_return['tophat'] = tophat_d
    return to_return


def assemble_provenance(pkg_root, inputs={}, run_options=None,
                        schema_version=None):
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
                           flux_file=False):
    to_return = assemble_provenance(pkg_root, inputs=inputs,
                                    run_options=run_options)
    if flux_file:
        # add a section 'flux_dependencies' containing at least
        # galsim version and, if possible, throughputs version
        to_return['flux_dependencies'] = dict()
        from galsim import version as galsim_version
        to_return['flux_dependencies']['galsim_version'] = galsim_version
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


class ConfigWriter():
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

        # associate fragment name with config builder
        self._make_frag = {'star': create_star_config,
                           'galaxy': create_galaxy_config,
                           'diffsky_galaxy': create_diffsky_galaxy_config,
                           'sso': create_sso_config,
                           'snana': create_snana_config,
                           'gaia_star_butler': create_gaia_star_butler_config,
                           'gaia_star_direct': create_gaia_star_direct_config}

    def write_yaml(self, input_dict, outpath):
        if not self._overwrite:
            try:
                with open(outpath, mode='x') as f:
                    yaml.dump(input_dict, f)
            except FileExistsError:
                txt = 'write_yaml: Will not overwrite pre-existing config'
                if self._logname:
                    logger = logging.getLogger(self._logname)
                    logger.warning(txt + outpath)
                else:
                    print(txt + outpath)
                return
        else:
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

    def write_configs(self, object_type, provenance, cosmology=None,
                      tophat_bins=None, fragment_name=None, id_prefix=None,
                      butler_parameters=None, basename_template=None):
        '''
        Write yaml fragment for specified object type and write or update
        top yaml file if
        * fragment doesn't already exist
        * or overwrite is true

        Parameters
        ----------
        object_type   string       An object type
        provenance    dict
        cosmology     dict         Ignored for all but galaxy, diffsky_galaxy
        tophat_bins   list         Ignored for all but galaxy
        fragment_name string       Name (not include '.yaml') of file to be
                                   written.  Defaults to object_type
        id_prefix     string       Ignored for all but gaia_star
        butler_parameters dict     Ignored for all but gaia_star (via butler)
        basename_template string   Ignored for all but giaa_star (direct FITS)
        '''
        overwrite = self._overwrite
        top_path = os.path.join(self._out_dir, self._top_name + '.yaml')

        top = _read_yaml(top_path, silent=True, resolve_include=False)
        if top:
            top_exists = True
            object_type_exists = object_type in top['object_types']
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
        if not fragment_name:
            fragment_name = object_type
        if fragment_name in {'star', 'sso', 'snana'}:
            frag = self._make_frag[fragment_name](provenance)
        elif fragment_name == 'diffsky_galaxy':
            frag = self._make_frag[fragment_name](provenance, cosmology)
        elif fragment_name == 'galaxy':
            frag = self._make_frag[fragment_name](provenance, cosmology,
                                                  tophat_bins)
        elif fragment_name == 'gaia_star_butler':
            frag = self._make_frag[fragment_name](provenance,
                                                  id_prefix=id_prefix,
                                                  butler_parameters=butler_parameters)
        elif fragment_name == 'gaia_star_direct':
            frag = self._make_frag[fragment_name](provenance,
                                                  id_prefix=id_prefix,
                                                  basename_template=basename_template)
        else:
            raise ValueError(f'ConfigWriter.write_configs: unknown fragment {fragment_name}')
        basen = fragment_name + '.yaml'
        frag_path = os.path.join(self._out_dir, basen)
        self.write_yaml(frag, frag_path)

        # Write or update top file if necessary
        value = '!include ' + basen
        if top_exists and not overwrite:
            if object_type_exists and top['object_types'][object_type] == value:
                # No change necessary
                return

            # Otherwise need to add or modify value for our object type
            top['object_types'][object_type] = value
            self.write_yaml(top, top_path)
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
