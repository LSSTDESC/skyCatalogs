from collections import namedtuple, OrderedDict
from enum import Enum
import numpy as np

__all__ = ['column_finder', 'check_file', 'write_to_instance', 'SourceType', 'STAR_FMT', 'CMP_FMT',
           'form_star_instance_column', 'form_cmp_instance_columns']


STAR_FMT = '{:s} {:d} {:.14f} {:.14f} {:.8f} {:s} {:d} {:d} {:d} {:d} {:d} {:d} {:s} {:s} {:s} {:.8f} {:f}'

CMP_FMT = '{:s} {:d} {:.14f} {:.14f}, {:.8f}, {:s} {:.9g} {:.9g} {:.9g} {:.9g} {:d} {:d} {:s} {:.9g} {:.9g} {:f} {:s} {:s} {:.8f} {:f}'

def form_star_instance_columns(band):
    star_instance = [column_finder('prefix', SourceType.FIXED, ('object', np.dtype('U6'))),
                     column_finder('uniqueId', SourceType.DATA, 'id'),
                     column_finder('raPhoSim', SourceType.DATA, 'ra'),
                     column_finder('decPhoSim', SourceType.DATA, 'dec'),
                     column_finder('maskedMagNorm', SourceType.DATA, 'magnorm'),
                     column_finder('sedFilepath',SourceType.DATA, 'sed_filepath'),
                     column_finder('redshift', SourceType.FIXED, (0, int)),
                     column_finder('gamma1', SourceType.FIXED, (0, int)),
                     column_finder('gamma2', SourceType.FIXED, (0, int)),
                     column_finder('kappa', SourceType.FIXED, (0, int)),
                     column_finder('raOffset', SourceType.FIXED, (0, int)),
                     column_finder('decOffset', SourceType.FIXED, (0, int)),
                     column_finder('spatialmodel', SourceType.FIXED, ('point', np.dtype('U5'))),
                     column_finder('internalExtinctionModel', SourceType.FIXED, ('none', np.dtype('U4'))),
                     column_finder('galacticExtinctionModel', SourceType.CONFIG,
                                   'object_types/star/MW_extinction'),
                     column_finder('galactivAv', SourceType.DATA,
                                   f'MW_av_lsst_{band}'),
                     column_finder('galacticRv', SourceType.CONFIG,
                                   'MW_extinction_values/r_v/value')]
    return star_instance

def form_cmp_instance_columns(cmp, band):
    cmp_instance = [column_finder('prefix', SourceType.FIXED, ('object', np.dtype('U6'))),
                    column_finder('uniqueId', SourceType.DATA, 'galaxy_id'),
                    column_finder('raPhoSim', SourceType.DATA, 'ra'),
                    column_finder('decPhoSim', SourceType.DATA, 'dec'),
                    column_finder('phoSimMagNorm', SourceType.DATA, f'{cmp}_magnorm'),
                    column_finder('sedFilepath',SourceType.COMPUTE, [f'sed_val_{cmp}','redshift_hubble']),
                    column_finder('redshift', SourceType.DATA, 'redshift'),
                    column_finder('gamma1', SourceType.DATA, 'shear_1'),
                    column_finder('gamma2', SourceType.DATA, 'shear_2'),
                    column_finder('kappa', SourceType.DATA, 'convergence'),
                    column_finder('raOffset', SourceType.FIXED, (0, int)),
                    column_finder('decOffset', SourceType.FIXED, (0, int)),
                    column_finder('spatialmodel', SourceType.CONFIG,
                                  f'object_types/{cmp}_basic/spatial_model'),
                    column_finder('majorAxis', SourceType.DATA, f'size_{cmp}_true'),
                    column_finder('minorAxis', SourceType.DATA, f'size_minor_{cmp}_true'),
                    column_finder('positionAngle', SourceType.DATA, 'position_angle_unlensed'),
                    column_finder('internalExtinctionModel', SourceType.FIXED, ('none', np.dtype('U4'))),
                    column_finder('galacticExtinctionModel', SourceType.CONFIG,
                                  f'object_types/{cmp}_basic/MW_extinction'),
                    column_finder('galactivAv', SourceType.DATA,
                                  f'MW_av_lsst_{band}'),
                    column_finder('galacticRv', SourceType.CONFIG,
                                  'MW_extinction_values/r_v/value')]
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
    except FileNotFoundError as e:
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

def write_to_string(row, fmt):
    '''
    Output single string, composed from input values in row.
    Types of values in row must match those expected by fmt.
    Parameters
    ----------
    row     list of values
    fmt     fmt string.

    Returns
    -------
    Formatted string

    '''
    return(fmt.format(*row))
