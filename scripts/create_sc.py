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
import logging
from desc.skycatalogs.catalog_creator import CatalogCreator
from desc.skycatalogs.utils.common_utils import print_date, log_callinfo


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
                    help='Root directory for sky catalogs, typically site-dependent. If not specified, use value of environment variable SKYCATALOG_ROOT')
parser.add_argument('--catalog-dir', help='directory for output files relative to skycatalog_root',
                    default='.')
parser.add_argument('--sed-subdir', help='subdirectory to prepend to paths of galaxy SEDs as written to the sky catalog', default='galaxyTopHatSED')
parser.add_argument('--log-level', help='controls logging output',
                    default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'])
parser.add_argument('--galaxy-magnitude-cut', default=29.0, type=float,
                    help='Exclude galaxies with r-magnitude above this value')
parser.add_argument('--knots-magnitude-cut', default=27.0, type=float,
                    help='Galaxies with i-magnitude above this cut get no knots')
parser.add_argument('--no-knots', action='store_true', help='If supplied omit knots component. Default is False')

parser.add_argument('--config-path', default=None, help='''
                    Output path for config file. If no value,
                    config will be written to same location as data,
                    with filenmame taken from catalog-name.
                    A config file is written iff galaxies are requested''')
parser.add_argument('--catalog-name', default='skyCatalog',
                    help='''If write-config option is used, this value
                    will appear in the config and will be part of
                    the filename''')

parser.add_argument('--overwrite-data', action='store_true',
                    help='If supplied overwrite existing data files; else skip with message')
parser.add_argument('--flux-only', action='store_true',
                    help='If supplied only do flux files. Main files must already exist')
parser.add_argument('--main-only', action='store_true',
                    help='If supplied only do main files, not flux files')
parser.add_argument('--flux-parallel', default=16, type=int,
                    help='Number of processes to run in parallel when computing fluxes')


args = parser.parse_args()
logname = 'skyCatalogs.creator'
logger = logging.getLogger(logname)
logger.setLevel(args.log_level)

ch = logging.StreamHandler()
ch.setLevel(args.log_level)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)

logger.addHandler(ch)

log_callinfo('create_sc', args, logname)

skycatalog_root = args.skycatalog_root
if not skycatalog_root:
    skycatalog_root = os.getenv('SKYCATALOG_ROOT')

parts = args.pixels

creator = CatalogCreator(parts, area_partition, skycatalog_root=skycatalog_root,
                         catalog_dir=args.catalog_dir,
                         config_path=args.config_path,
                         catalog_name=args.catalog_name,
                         mag_cut=args.galaxy_magnitude_cut,
                         sed_subdir=args.sed_subdir,
                         knots_mag_cut=args.knots_magnitude_cut,
                         knots=(not args.no_knots), logname=logname,
                         overwrite_data=args.overwrite_data,
                         flux_only=args.flux_only, main_only=args.main_only,
                         flux_parallel=args.flux_parallel)

logger.info(f'Starting with healpix pixel {parts[0]}')
if not args.no_galaxies:
    logger.info("Creating galaxy catalogs")
    creator.create('galaxy')

if not args.no_pointsources:
    logger.info("Creating point source catalogs")
    creator.create('pointsource')


logger.info('All done')
print_date()
