from collections import namedtuple, OrderedDict

__all__ = ['column_finder', 'check_file', 'write_to_instance' ]

column_finder = namedtuple('ColumnFinder', ['instance_name', 'source_type',
                                            'source_parm'])

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
