from collections import namedtuple, OrderedDict
from enum import Enum

__all__ = ['column_finder', 'check_file', 'write_to_instance', 'SourceType' ]

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
        handle.write(fmt.format(*row))
