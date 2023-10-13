from collections import namedtuple
from enum import Enum
import numpy as np

__all__ = ['column_finder', 'check_file', 'write_to_instance', 'SourceType',
           'STAR_FMT', 'CMP_FMT',
           'form_star_instance_columns', 'form_cmp_instance_columns']


STAR_FMT = '''{:s} {:s} {:.14f} {:.14f} {:.8f} {:s} {:d} {:d} {:d} {:d} {:d}
{:d} {:s} {:s} {:s} {:.8f} {:f}'''

CMP_FMT = '''{:s} {:s} {:.14f} {:.14f} {:.8f} {:s} {:.9g} {:.9g} {:.9g} {:.9g}
 {:d} {:d} {:s} {:.9g} {:.9g} {:f} {:.0f} {:s} {:s} {:.8f} {:f}'''


def form_star_instance_columns(band):
    star_instance = [column_finder('prefix', SourceType.FIXED,
                                   ('object', np.dtype('U6'))),
                     # column_finder('uniqueId', SourceType.DATA, 'id'),
                     column_finder('uniquePsId', SourceType.COMPUTE, ['id']),
                     column_finder('raPhoSim', SourceType.DATA, 'ra'),
                     column_finder('decPhoSim', SourceType.DATA, 'dec'),
                     column_finder('maskedMagNorm', SourceType.DATA,
                                   'magnorm'),
                     column_finder('sedFilepath', SourceType.DATA,
                                   'sed_filepath'),
                     column_finder('redshift', SourceType.FIXED, (0, int)),
                     column_finder('gamma1', SourceType.FIXED, (0, int)),
                     column_finder('gamma2', SourceType.FIXED, (0, int)),
                     column_finder('kappa', SourceType.FIXED, (0, int)),
                     column_finder('raOffset', SourceType.FIXED, (0, int)),
                     column_finder('decOffset', SourceType.FIXED, (0, int)),
                     column_finder('spatialmodel', SourceType.FIXED,
                                   ('point', np.dtype('U5'))),
                     column_finder('internalExtinctionModel', SourceType.FIXED,
                                   ('none', np.dtype('U4'))),
                     column_finder('galacticExtinctionModel',
                                   SourceType.CONFIG,
                                   ('object_types/star/MW_extinction', str)),
                     column_finder('galactivAv', SourceType.DATA,
                                   'MW_av'),
                     column_finder('galacticRv', SourceType.CONFIG,
                                   ('MW_extinction_values/r_v/value', float))]
    return star_instance


def _form_knots_instance_columns(cmp, band):
    cmp_instance = [column_finder('prefix', SourceType.FIXED, ('object',
                                                               np.dtype('U6'))),
                    column_finder('uniqueId', SourceType.COMPUTE,
                                  ['galaxy_id', f'{cmp}']),
                    column_finder('raPhoSim', SourceType.DATA, 'ra'),
                    column_finder('decPhoSim', SourceType.DATA, 'dec'),
                    column_finder('phoSimMagNorm', SourceType.DATA,
                                  'knots_magnorm'),
                    column_finder('sedFilepath', SourceType.COMPUTE,
                                  [f'sed_val_{cmp}', 'redshift_hubble']),
                    column_finder('redshift', SourceType.DATA, 'redshift'),
                    column_finder('gamma1', SourceType.DATA, 'shear_1'),
                    column_finder('gamma2', SourceType.DATA, 'shear_2'),
                    column_finder('kappa', SourceType.DATA, 'convergence'),
                    column_finder('raOffset', SourceType.FIXED, (0, int)),
                    column_finder('decOffset', SourceType.FIXED, (0, int)),
                    column_finder('spatialmodel', SourceType.CONFIG,
                                  (f'object_types/{cmp}_basic/spatial_model',
                                   'str')),
                    column_finder('majorAxis', SourceType.DATA,
                                  'size_disk_true'),
                    column_finder('minorAxis', SourceType.DATA,
                                  'size_minor_disk_true'),
                    column_finder('positionAngle', SourceType.DATA,
                                  'position_angle_unlensed'),
                    column_finder('sindex', SourceType.DATA, 'n_knots'),
                    column_finder('internalExtinctionModel', SourceType.FIXED,
                                  ('none', np.dtype('U4'))),
                    column_finder('galacticExtinctionModel', SourceType.CONFIG,
                                  (f'object_types/{cmp}_basic/MW_extinction',
                                   'str')),
                    column_finder('galactivAv', SourceType.DATA,
                                  'MW_av'),
                    column_finder('galacticRv', SourceType.CONFIG,
                                  ('MW_extinction_values/r_v/value', float))]
    return cmp_instance


def form_cmp_instance_columns(cmp, band):
    if cmp == 'knots':
        return _form_knots_instance_columns(cmp, band)

    cmp_instance = [column_finder('prefix', SourceType.FIXED,
                                  ('object', np.dtype('U6'))),
                    column_finder('uniqueId', SourceType.COMPUTE,
                                  ['galaxy_id', f'{cmp}']),
                    column_finder('raPhoSim', SourceType.DATA, 'ra'),
                    column_finder('decPhoSim', SourceType.DATA, 'dec'),
                    column_finder('phoSimMagNorm', SourceType.DATA,
                                  f'{cmp}_magnorm'),
                    column_finder('sedFilepath', SourceType.COMPUTE,
                                  [f'sed_val_{cmp}', 'redshift_hubble']),
                    column_finder('redshift', SourceType.DATA, 'redshift'),
                    column_finder('gamma1', SourceType.DATA, 'shear_1'),
                    column_finder('gamma2', SourceType.DATA, 'shear_2'),
                    column_finder('kappa', SourceType.DATA, 'convergence'),
                    column_finder('raOffset', SourceType.FIXED, (0, int)),
                    column_finder('decOffset', SourceType.FIXED, (0, int)),
                    column_finder('spatialmodel', SourceType.CONFIG,
                                  (f'object_types/{cmp}_basic/spatial_model',
                                   'str')),
                    column_finder('majorAxis', SourceType.DATA,
                                  f'size_{cmp}_true'),
                    column_finder('minorAxis', SourceType.DATA,
                                  f'size_minor_{cmp}_true'),
                    column_finder('positionAngle', SourceType.DATA,
                                  'position_angle_unlensed'),
                    column_finder('sindex', SourceType.DATA, f'sersic_{cmp}'),
                    column_finder('internalExtinctionModel', SourceType.FIXED,
                                  ('none', np.dtype('U4'))),
                    column_finder('galacticExtinctionModel', SourceType.CONFIG,
                                  (f'object_types/{cmp}_basic/MW_extinction',
                                   'str')),
                    column_finder('galactivAv', SourceType.DATA,
                                  'MW_av'),
                    column_finder('galacticRv', SourceType.CONFIG,
                                  ('MW_extinction_values/r_v/value', float))]
    return cmp_instance


column_finder = namedtuple('ColumnFinder', ['instance_name', 'source_type',
                                            'source_parm'])

SourceType = Enum('SourceType', 'DATA CONFIG FIXED COMPUTE')
'''
Used in source_type field of column_finder to describe source of each
instance catalog value

DATA    Sky catalog column
CONFIG  Sky catalog config field
FIXED   Fixed value supplied inline
COMPUTE Arbitrary computation which may involve any of the above
'''


def check_file(path):
    '''Look for a file that should not exist'''
    try:
        f = open(path, mode='r')
    except FileNotFoundError:
        return

    raise ValueError(f'File for {path} already exists')


def write_to_instance(handle, ordered_data_dict, fmt):
    '''
    Parameters
    ----------
    handle             file handle to which text will be written
    ordered_data_dict  columns in the order they should be written
    '''
    col_list = [ordered_data_dict[k] for k in ordered_data_dict]
    for row in zip(*col_list):
        handle.write(fmt.format(*row) + '\n')


def form_object_string(obj, band, component):
    '''
    parse columns for this object/component type
    fetch data and config values as needed
    form values into a list.
    '''
    row = []

    if obj.object_type in ['star', 'sncosmo']:
        fmt = STAR_FMT
        columns = form_star_instance_columns(band)
    elif obj.object_type in ['galaxy', 'disk', 'bulge']:
        fmt = CMP_FMT
        if obj.object_type == 'galaxy':
            cmp = component
        else:
            cmp = obj.object_type
        columns = form_cmp_instance_columns(cmp, band)
    else:
        raise NotImplementedError(f'translate_utils.form_object_string: Unsupported object type {obj.object_type}')

    for c in columns:
        if c.source_type == SourceType.DATA:
            row.append(obj.get_native_attribute(c.source_parm))
        elif c.source_type == SourceType.CONFIG:
            v = obj.belongs_to.config.get_config_value(c.source_parm[0])
            t = c.source_parm[1]
            if str(t) in ['float', 'int']:
                q = eval(f'{t}({v})')
            else:
                q = v
            row.append(q)
        elif c.source_type == SourceType.FIXED:
            v = c.source_parm[0]
            t = c.source_parm[1]
            if str(t) in ['float', 'int']:
                q = eval(f'{t}({v})')
            else:
                q = v
            row.append(q)
        elif c.source_type == SourceType.COMPUTE:
            if c.instance_name not in ['sedFilepath', 'uniqueId',
                                       'uniquePsId']:
                raise ValueError(f'translate_utils.form_object_string: Bad COMPUTE entry {c.instance_name}')
            if c.instance_name == 'sedFilepath':
                row.append(f'{obj.get_native_attribute("galaxy_id")}_{cmp}_sed.txt')
            elif c.instance_name == 'uniqueId':
                if cmp not in ['disk', 'bulge', 'knots']:
                    raise ValueError(f'translate_utils.form_object_string: Bad COMPUTE entry {c.instance_name} for component {cmp}')
                row.append(f'{str(obj.get_native_attribute("galaxy_id"))}_{cmp}')
            else:    # uniquePsId.  Output must be string but input may be int
                row.append(f'{str(obj.get_native_attribute("id"))}')
        else:
            raise ValueError(f'translate_utils.form_object_string: Unknown source type {c.source_type}')

    return write_to_string(row, fmt)


def write_to_string(row, fmt):
    '''
    Output single string, composed from input values in row.
    Types of values in row must match those expected by fmt.
y    Parameters
    ----------
    row     list of values
    fmt     fmt string.

    Returns
    -------
    Formatted string

    '''
    return(fmt.format(*row))
