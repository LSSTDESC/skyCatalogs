import sys
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pandas as pd


def make_parquet(input_path):
    '''
    Given a text file in table format where columns are wavelength (nanometers)
    and flux, convert to parquet
    '''

    wv, flux = np.genfromtxt(input_path, unpack=True)
    # wv32 = np.array(wv, np.float32)

    df = pd.DataFrame({'wavelength': wv, 'flux': flux})
    # df = pd.DataFrame({'wavelength': wv32, 'flux': flux})
    # This method produces a file somewhat larger than if we left wv alone

    table = pa.Table.from_pandas(df)

    # as does this
    # table = pa.Table.from_arrays([wv32, flux], names=['wavelength', 'flux'])

    out_path = input_path + '.parquet'
    pq.write_table(table, out_path)


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print('Requires filepath argument')
        exit(0)

    fname = sys.argv[1]

    print('Called with filepath argument ', fname)

    mode = 'parquet'

    if len(sys.argv) > 2:
        mode = sys.argv[2]

    if mode == 'parquet':
        make_parquet(fname)
