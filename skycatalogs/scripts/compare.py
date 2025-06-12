import os
import numpy as np
import pyarrow.parquet as pq
import argparse

'''
Compare contents of two parquet files which are supposed to be the same
(except for metadata).
'''

parser = argparse.ArgumentParser(
    description='''
    Compare contents of two skyCatalog files''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('file1', help='full path to the first file')
parser.add_argument('file2',
                    help='''
full path to the second file OR just containing directory (including final
slash) if basename is the same as for file1''')
parser.add_argument('--object-type', default=None,
                    choices=['galaxy', 'star', 'sso', 'trilegal', None],
                    help='If None deduce from file name')
parser.add_argument('--catalog-type', default=None,
                    choices=['main', 'flux', None],
                    help='flux or main.  If None, deduce from file name')
parser.add_argument('--columns', nargs='*', default=[],
                    help='Columns to compare.  If None, deduce from object and catalog types')
parser.add_argument('--debug', action='store_true',
                    help='print output when operations succeed')

args = parser.parse_args()

file1 = args.file1
dir1, base1 = os.path.split(file1)
dir2, base2 = os.path.split(args.file2)
if not base2:
    file2 = os.path.join(dir2, base1)
else:
    file2 = args.file2

if args.debug:
    print(f'Comparing {file1} and {file2}')

# Basenames are of the form {object-type}_{hp}.parquet or
# {object-type}_flux_{hp}.parquet
cmps = base1.split('_')
obj = args.object_type
if not obj:
    obj = cmps[0]
if obj == 'pointsource':
    obj = 'star'
cat_type = args.catalog_type
if not cat_type:
    if len(cmps) == 2:
        cat_type = 'main'
    else:
        cat_type = 'flux'

cols = args.columns
if len(cols) == 0:
    if cat_type == 'flux':
        cols.append('lsst_flux_r')
        if obj == 'galaxy':
            cols.append('galaxy_id')
        else:
            cols.append('id')
    else:      # main
        cols.extend(['ra', 'dec'])
        if obj == 'galaxy':
            cols.extend(['galaxy_id', 'MW_av', 'redshift'])
        elif obj == 'star':
            cols.extend(['id', 'MW_av', 'magnorm'])
        elif obj == 'sso':
            cols.extend(['id', 'mjd', 'trailed_source_mag'])
        elif obj == 'trilegal':
            cols.extend(['id', 'logT', 'imag'])

pq_file1 = pq.ParquetFile(file1)
pq_file2 = pq.ParquetFile(file2)

assert pq_file1.metadata.num_rows == pq_file2.metadata.num_rows, 'Same number of rows'

if args.debug:
    print("Number of rows is identical")

tbl1 = pq_file1.read_row_group(0, columns=cols)
tbl2 = pq_file2.read_row_group(0, columns=cols)

for c in cols:
    assert tbl1[c] == tbl2[c], f'Column {c} is identical'
    if args.debug:
        print(f"Column {c} is identical")
