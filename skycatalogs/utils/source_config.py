import yaml
import os
from pathlib import PurePath

'''
Utilities to create suitable config fragments for various source types
'''

__all__ = ['create_galaxy_config', 'create_star_config',
           'create_diffsky_galaxy_config', 'create_gaia_star_butler_config',
           'create_gaia_star_direct_config', 'create_snana_config',
           'create_sso_config']
_FILE_PATH = str(PurePath(__file__))
_SKYCATALOGS_DIR = _FILE_PATH[:_FILE_PATH.rindex('/skycatalogs')]

_TEMPLATE_DIR = os.path.join(_SKYCATALOGS_DIR, 'skycatalogs', 'data',
                             'cfg_templates')


def create_galaxy_config(provenance, cosmology, tophat_bins,
                         area_partition=None, data_file_type=None):
    '''
    Make a dict of DC2 galaxy parameters

    Parameters
    ----------
    provenance:   a dict including keys for inputs, code version, etc.
    cosmology:    a dict of cosmology parameters
    tophat_bins: A list of integer duples (start, width) describing the bins
    area_partition: Dict describing partition.  Defaults to standard scheme
    data_file_type: string   Defaults to 'parquet'

    Returns
    -------
    A dict

The resulting dict looks much like the file galaxy_template.yaml under
skycatalogs/cfg_templates with added tophat bins, provenance and cosmology

    '''
    opt_dict = {'area_partition': area_partition,
                'data_file_type': data_file_type}
    data = _generic_create('galaxy_template.yaml', provenance,
                           opt_dict)

    # add tophat bins, cosmology
    data['tophat']['bins'] = tophat_bins
    data['cosmology'] = cosmology
    return data


def create_star_config(provenance, area_partition=None, data_file_type=None):
    '''
    Make a dict of DC2 star parameters.

    Parameters
    ----------
    provenance:   a dict including keys for inputs, code version, etc.
    area_partition: Dict describing partition.  Defaults to standard scheme
    data_file_type: string   Defaults to 'parquet'

    Returns
    -------
    A dict

    '''
    opt_dict = {'area_partition': area_partition,
                'data_file_type': data_file_type}
    return _generic_create('star_template.yaml', provenance,
                           opt_dict)


def create_diffsky_galaxy_config(provenance, cosmology, area_partition=None,
                                 data_file_type=None):
    '''
    Make a dict of diffsky galaxy parameters

    Parameters
    ----------
    provenance:   a dict including keys for inputs, code version, etc. Required
    cosmology:    a dict of cosmology parameters.  Required
    area_partition: Dict describing partition.  Defaults to standard scheme
    data_file_type: string   Defaults to 'parquet'

    Returns
    -------
    A dict

    The resulting dict looks much like the file galaxy_diffsky template.yaml
    under skycatalogs/cfg_templates with added sections for cosmology and
    provenance

    '''
    opt_dict = {'area_partition': area_partition,
                'data_file_type': data_file_type}
    data = _generic_create('diffsky_galaxy_template.yaml', provenance,
                           opt_dict)

    data['cosmology'] = cosmology
    return data


def create_gaia_star_butler_config(provenance, id_prefix=None,
                                   butler_parameters=None):
    '''
    Create config for access to gaia stars via butler

    Parameters
    ----------
    provenance dict  (required)
    id_prefix  string   String to prepend to ids in the data  (optional)
    butler_parameters  dict    (optional)

    Defaults for optional arguments can be found in the template file
    gaia_star_butler_template.yaml

    Returns
    -------
    A dict
    '''
    opt_dict = {'id_prefix': id_prefix, 'butler_parameters': butler_parameters}

    return _generic_create('gaia_star_butler_template.yaml',
                           provenance, opt_dict)


def create_gaia_star_direct_config(provenance, id_prefix=None,
                                   area_partition=None, data_dir=None,
                                   basename_template=None):
    '''
    Create config for access to gaia stars via butler

    Parameters
    ----------
    provenance dict  (required)
    id_prefix  string   String to prepend to ids in the data  (optional)
    area_partition  dict    (optional)
    data_dir        string   (optional)
    basename_template string  (optional)

    Defaults for optional arguments can be found in the template file
    gaia_star_direct_template.yaml

    Returns
    -------
    A dict
    '''
    opt_dict = {'id_prefix': id_prefix, 'area_partition': area_partition,
                'data_dir': data_dir, 'basename_template': basename_template}

    return _generic_create('gaia_star_direct_template.yaml',
                           provenance, opt_dict)


def create_snana_config(provenance, area_partition=None, data_file_type=None):
    '''
    Make a dict of DC2 star parameters.

    Parameters
    ----------
    provenance:  dict including keys for inputs, code version, etc. Required
    area_partition: dict describing partition.  Optional
    data_file_type: string   Optional

    Defaults for optional arguments can be found in the template file
    gaia_star_direct_template.yaml

    Returns
    -------
    A dict

    '''
    opt_dict = {'area_partition': area_partition,
                'data_file_type': data_file_type}
    return _generic_create('snana_template.yaml', provenance, opt_dict)


def create_sso_config(provenance, area_partition=None, data_file_type=None):
    opt_dict = {'area_partition': area_partition,
                'data_file_type': data_file_type}
    return _generic_create('sso_template_yaml', provenance, opt_dict)


def _generic_create(template_filename, provenance, opt):
    '''
    Parameters
    ----------
    template_filename  string   filename of template config
    provenance       dict      Contents of the provenance node.  Cannot be
                               None
    opt              dict      keywords whose values may need to be replaced
                               in the template config

    '''
    template_path = os.path.join(_TEMPLATE_DIR, template_filename)
    with open(template_path, 'r') as f:
        data = yaml.load(f, Loader=yaml.SafeLoader)
    other = dict()
    for key in opt:
        if opt[key] is not None:
            other[key] = opt[key]
    if len(other.keys()) > 0:
        data.update(other)

    data['provenance'] = provenance
    return data
