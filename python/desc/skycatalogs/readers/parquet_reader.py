from collections import OrderedDict
import pyarrow as pa
import pyarrow.parquet as pq
import numpy as np
import numpy.ma as ma
import warnings

__all__ = ['ParquetReader']

class ParquetReader:
    '''
    Handle reads for a particular file
    Create a separate ParquetReader instance for each file, region
    (as specified by optional mask argument)
    '''
    def __init__(self, filepath, mask=None):
        self._filepath = filepath
        # meta data includes num_columns, num_rows, num_row_groups
        self._meta = pq.read_metadata(filepath)
        self._schema = pq.read_schema(filepath)
        self._columns = set(self._schema.names)
        #  To support iterator also save # row groups, length of each.
        #  If mask, must be able to get slices corresponding to the different
        #  row groups
        self._open()

    def _open(self):
        self._pqfile = pq.ParquetFile(self._filepath)

    def close(self):
        '''
        There doesn't seem to be any explicit way to close the file
        via pyarrow.parquet.  Maybe delete the object?   Until this is
        understood do nothing.
        '''

        if self._pqfile:
            pass


    def read_columns(self, cols, mask):
        '''
        Parameters
        -----------
        cols       list of column names belonging to the file
        mask       if not None, use it and return compressed array

        Returns
        -------
        OrdereDict where keys are column names, values are numpy arrays
        '''

        if not set(cols).issubset(self._columns):
            # raise exception?   For now, just
            return None

        d = OrderedDict()
        if not self._pqfile:
            self._open()
        tbl = self._pqfile.read(columns = cols)
        for c in cols:
            c_data = np.array(tbl[c])
            if mask is not None:
                d[c] = ma.array(c_data, mask=mask).compressed()
            else:
                d[c] = c_data

        return d
