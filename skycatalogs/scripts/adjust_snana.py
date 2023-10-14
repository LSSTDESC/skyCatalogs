'''
Read in original SNANA parquet file; output file identical except for
added columns MW_av, MW_rv
'''
import os
import argparse
from skycatalogs.utils.add_extinction import AddExtinction

parser = argparse.ArgumentParser(description='''
Make a new skyCatalogs parquet file from an old, adding columns for MW_av, MW_rv
''')

parser.add_argument('indir',
                    help='Directory containing input file(s). Required')
parser.add_argument('outdir', help='Directory for output. Required')
parser.add_argument('--pixels', type=int, nargs='*', default=[],
                    help='healpix pixels for which new files will be created. Required')
parser.add_argument('--starts-with',
                    help='That part of the filename preceding healpixel',
                    default='snana_')

args = parser.parse_args()

if os.path.abspath(args.indir) == os.path.abspath(args.outdir):
    raise ValueError('Input and output directories must be different')

writer = AddExtinction(args.indir, args.outdir, args.starts_with)

for p in args.pixels:
    writer.write(p)
