import pyarrow as pa
import pyarrow.parquet as pq

__all__ = ['ParquetReader']

class ParquetReader:
    '''
    Handle reads for a particular file
    '''
    def __init__(self, filepath):
        self._filepath = filepath
        self._pqfile = self._open()

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


    def read_columns(self, cols):
        '''
        Parameters
        -----------
        cols      list of column names belonging to the file

        Returns
        -------
        pyarrow Table
        '''

        if not self._pqfile:
            self._open()
        return self._pqfile.read(columns = cols)
