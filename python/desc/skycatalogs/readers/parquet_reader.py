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
        self._pqfile = self._open()
        self._mask = mask

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


    def read_columns(self, cols, apply_mask=True):
        '''
        Parameters
        -----------
        cols       list of column names belonging to the file
        apply_mask if True and mask has been set, return compressed array

        Returns
        -------
        dict where keys are column names, values are numpy arrays
        '''

        d = dict()
        use_mask = apply_mask and (self._mask is not None)
        if not self._pqfile:
            self._open()
        tbl = self._pqfile.read(columns = cols)
        for c in cols:
            c_data = np.array(tbl[c])
            if use_mask:
                d[c] = ma.array(c_data, mask=self._mask).compressed()
            else:
                d[c] = c_data

        return d

    def set_mask(self, mask):
        if not self._mask:
            self._mask = mask
        else:
            warnings.warn('ParquetReader: override of mask existing setting not allowed. Mask is unchanged')
