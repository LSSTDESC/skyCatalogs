import os
import pyarrow.parquet as pq
import pyarrow as pa
import pandas as pd

from skycatalogs.utils.creator_utils import make_MW_extinction_av
from skycatalogs.utils.creator_utils import make_MW_extinction_rv


class AddExtinction():
    def __init__(self, in_dir, out_dir, starts_with):
        '''
        Rewrite sky catalog-like parquet files, adding columns MW_ra, MW_rv

        Parameters
        ----------
        in_dir       string      directory where input files may be found
        out_dir      string      directory where outputs are to be written
        starts_with  string      Form of input an output file name is always
                                 <starts_with><pixel>.parquet
        '''
        self._in_dir = in_dir
        self._out_dir = out_dir
        self._starts_with = starts_with

    def write(self, pixel):
        fname = f'{self._starts_with}{str(pixel)}.parquet'
        infile = pq.ParquetFile(os.path.join(self._in_dir, fname))
        arrow_schema = (infile.schema).to_arrow_schema()
        out_schema = arrow_schema.append(pa.field('MW_av', pa.float32()))
        out_schema = out_schema.append(pa.field('MW_rv', pa.float32()))

        writer = pq.ParquetWriter(os.path.join(self._out_dir, fname),
                                  out_schema)
        n_row_group = infile.metadata.num_row_groups

        for g in range(n_row_group):
            out_dict = dict()
            tbl = infile.read_row_group(g)
            for c in arrow_schema.names:
                out_dict[c] = [i.as_py() for i in tbl[c]]
            out_dict['MW_av'] = make_MW_extinction_av(out_dict['ra'],
                                                      out_dict['dec'])
            out_dict['MW_rv'] = make_MW_extinction_rv(out_dict['ra'],
                                                      out_dict['dec'])
            out_df = pd.DataFrame.from_dict(out_dict)
            out_table = pa.Table.from_pandas(out_df, schema=out_schema)
            writer.write_table(out_table)

        writer.close()
