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
parser.add_argument('--no-pointsources', action='store_true',
                    help='if used, point source catalogs will NOT be created')
parser.add_argument('--no-galaxies', action='store_true',
                    help='if used galaxy catalogs will NOT be created')
parser.add_argument('--pixels', type=int, nargs='*', default=[9556],
                    help='healpix pixels for which catalogs will be created')
parser.add_argument('--skycatalog_root',
                    help='Root directory for sky catalogs, typically site-dependent. If not specified, use value environment variable SKYCATALOG_ROOT')
parser.add_argument('--catalog-dir', help='directory for output files relative to skycatalog_root',
                    default='.')
parser.add_argument('--sed-subdir', help='subdirectory to prepend to paths of galaxy SEDs as written to the sky catalog', default='galaxyTopHatSED')
parser.add_argument('--verbose', help='print more output if true',
                    action='store_true')
parser.add_argument('--galaxy-magnitude-cut', default=29.0, type=float,
                    help='Exclude galaxies with r-magnitude above this value')
parser.add_argument('--knots-magnitude-cut', default=27.0, type=float,
                    help='Galaxies with i-magnitude above this cut get no knots')
parser.add_argument('--no-knots', action='store_true', help='If supplied omit knots component. Default is False')

parser.add_argument('--write-config', action='store_true',
                    help='''If present output config for the new catalog.
                    Path is determined by --config-path.''')
parser.add_argument('--config-path', default=None, help='''
                    Output path for config file. If no value,
                    config will be written to same location as data,
                    with filenmame config.yaml''')

args = parser.parse_args()

print_callinfo('create_sc', args)

skycatalog_root = args.skycatalog_root
if not skycatalog_root:
    skycatalog_root = os.getenv('SKYCATALOG_ROOT')

parts = args.pixels

creator = CatalogCreator(parts, area_partition, skycatalog_root=skycatalog_root,
                         catalog_dir=args.catalog_dir,
                         write_config=args.write_config,
                         config_path=args.config_path,
                         verbose=args.verbose,
                         mag_cut=args.galaxy_magnitude_cut,
                         sed_subdir=args.sed_subdir,
                         knots_mag_cut=args.knots_magnitude_cut,
                         knots=(not args.no_knots))
print('Starting with healpix pixel ', parts[0])
if not args.no_galaxies:
    print("Creating galaxy catalogs")
    creator.create('galaxy')

if not args.no_pointsources:
    print("Creating point source catalogs")
    creator.create('pointsource')

print('All done')
print_date()
