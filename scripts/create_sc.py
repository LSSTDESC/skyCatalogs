# Try out CatalogCreator
# Note: root dir for SED files is $SIMS_SED_LIBRARY_DIR, defined by
#       LSST Science Pipelines setup
#       This applies for galaxy, star and AGN SEDs.
'''
Create sky catalogs for one or more healpixels. Invoke with --help
for details
'''
# For now partitioning is fixed

import os
import argparse
from desc.skycatalogs.catalog_creator import CatalogCreator
from desc.skycatalogs.utils.common_utils import print_date, print_callinfo


area_partition = {'type' : 'healpix', 'ordering' : 'ring', 'nside' : 32}

parser = argparse.ArgumentParser(description='''
Create Sky Catalogs. By default create a galaxy catalog for a
single healpix pixel''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--pointsource', action='store_true',
                    help='if used, create point source catalog(s)')
parser.add_argument('--no-galaxies', action='store_true',
                    help='if used galaxy catalogs will NOT be created')
parser.add_argument('--pixels', type=int, nargs='*', default=[9556],
                    help='healpix pixels for which catalogs will be created')
out_dir = os.path.join(os.getenv('SCRATCH'), 'desc', 'skycatalogs', 'test')
parser.add_argument('--output-dir', help='directory for output files',
                    default=out_dir)
parser.add_argument('--sed-subdir', help='subdirectory to prepend to paths of galaxy SEDs as written to the sky catalog', default='galaxyTopHatSED')
parser.add_argument('--verbose', help='print more output if true',
                    action='store_true')
parser.add_argument('--galaxy-magnitude-cut', default=None,
                    help='If supplied exclude galaxies with magnitude above')

args = parser.parse_args()

print_callinfo('create_sc', args)

output_dir = args.output_dir


parts = args.pixels

creator = CatalogCreator(parts, area_partition, output_dir=output_dir,
                         verbose=args.verbose,
                         mag_cut=args.galaxy_magnitude_cut,
                         sed_subdir=args.sed_subdir)
print('Starting with healpix pixel ', parts[0])
if not args.no_galaxies:
    print("Creating galaxy catalogs")
    creator.create('galaxy')

if args.pointsource:
    print("Creating point source catalogs")
    creator.create('pointsource')

print('All done')
print_date()
