import argparse
from skycatalogs.trilegal_catalog_creator import TrilegalSEDGenerator

parser = argparse.ArgumentParser(
    description='''
    Create Trilegal SEDs for specified healpixels. Assumes coresponding
    main SkyCatalogs files exist for these healpixels''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--pixels', type=int, nargs='*', default=[9556],
                    help='healpix pixels for which catalogs will be created')
parser.add_argument('--input-dir', help='where to find main files')
parser.add_argument('--output-dir',
                    help='where to write output. Defaults to input-dir',
                    default=None)

args = parser.parse_args()
output_dir = args.output_dir
if not output_dir:
    output_dir = args.input_dir

generator = TrilegalSEDGenerator(args.input_dir, output_dir)

generator.generate(args.pixels)
