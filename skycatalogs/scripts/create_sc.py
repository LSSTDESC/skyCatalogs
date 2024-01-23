#!/usr/bin/env python
# Try out CatalogCreator
# Note: root dir for SED files is $SIMS_SED_LIBRARY_DIR, defined by
#       LSST Science Pipelines setup
#       This applies for galaxy, star and AGN SEDs.
'''
Create sky catalogs for one or more healpixels. Invoke with --help
for details
'''

import os
import numpy as np
import argparse
import logging
import yaml
from skycatalogs.catalog_creator import CatalogCreator
from skycatalogs.utils.common_utils import print_date, log_callinfo

parser = argparse.ArgumentParser(description='''
Create Sky Catalogs. By default create a galaxy catalog for a
single healpix pixel''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--no-pointsources', '--no-point', action='store_true',
                    help='if used, point source catalogs will NOT be created')
parser.add_argument('--no-galaxies', '--no-gal', action='store_true',
                    help='if used galaxy catalogs will NOT be created')
parser.add_argument('--pixels', type=int, nargs='*', default=[9556],
                    help='healpix pixels for which catalogs will be created')
parser.add_argument('--skycatalog-root',
                    help='''Root directory for sky catalogs, typically
                    site-dependent. If not specified, use value of
                    environment variable SKYCATALOG_ROOT''')
parser.add_argument('--catalog-dir', '--cat-dir',
                    help='output file directory relative to skycatalog_root',
                    default='.')
parser.add_argument('--sed-subdir',
                    help='''subdirectory to prepend to paths of galaxy SEDs
                    as written to the sky catalog''',
                    default='galaxyTopHatSED')
parser.add_argument('--log-level', help='controls logging output',
                    default='INFO', choices=['DEBUG', 'INFO', 'WARNING',
                                             'ERROR'])
parser.add_argument('--galaxy-magnitude-cut', '--gal-mag-cut',
                    default=29.0, type=float,
                    help='Exclude galaxies with r-magnitude above this value')
parser.add_argument('--knots-magnitude-cut', default=27.0, type=float,
                    help='Galaxies with i-magnitude above this cut get no knots')
parser.add_argument('--no-knots', action='store_true',
                    help='If supplied omit knots component. Default is False')

parser.add_argument('--config-path', default=None, help='''
                    Output path for config file. If no value,
                    config will be written to same location as data,
                    with filenmame taken from catalog-name.
                    A config file is written iff galaxies are requested''')
parser.add_argument('--catalog-name', '--cat-name', default='skyCatalog',
                    help='''If write-config option is used, this value
                    will appear in the config and will be part of
                    the filename''')

parser.add_argument('--skip-done', action='store_true',
                    help='If supplied skip existing data files; else overwrite with message')
parser.add_argument('--flux-only', action='store_true',
                    help='If supplied only do flux files. Main files must already exist')
parser.add_argument('--main-only', action='store_true',
                    help='If supplied only do main files, not flux files')
parser.add_argument('--flux-parallel', default=16, type=int,
                    help='Number of processes to run in parallel when computing fluxes')
parser.add_argument('--provenance', '--prov', choices=['yaml'], help='''
                     Persist git provenance information for each file
                     written. Only supported format currently is as a
                     small yaml file, written to the data directory.''')
parser.add_argument('--options-file', default=None, help='''
                    path to yaml file associating option names with values.
                    Values for any options included will take precedence.''')
parser.add_argument('--dc2', action='store_true',
                    help='''If supplied provide values comparable to those
                            used for DC2 run. Default is False''')
parser.add_argument('--galaxy-nside', default=32, type=int,
                    choices=2**np.arange(15),
                    help='''Pixel nsides for galaxy output. Must be power of 2''')
parser.add_argument('--galaxy-stride', default=1_000_000, type=int,
                    help='''max # rows in a galaxy row group; default 1 million''')
parser.add_argument('--galaxy-type', default='cosmodc2',
                    choices=['cosmodc2', 'diffsky'],
                    help='''Galaxies may be based on either a cosmodc2 or
                    diffsky catalog''')
parser.add_argument('--galaxy-truth', default=None, help='''Truth catalog
                    on which skyCatalogs galaxies will be based. Default
                    depends on value of --galaxy-type option''')
parser.add_argument('--include-roman-flux', action='store_true',
                    help='If supplied calculate & store Roman as well as LSST  fluxes')
parser.add_argument('--star-input-fmt', default='sqlite',
                    choices=['sqlite', 'parquet'], help='''
                    star truth may come from either sqlite db or collection
                    of parquet files''')

args = parser.parse_args()

if args.options_file:
    with open(args.options_file) as f:
        opt_dict = yaml.safe_load(f)
        for k in opt_dict:
            if k in args:
                args.__setattr__(k, opt_dict[k])
            else:
                raise ValueError(f'Unknown attribute "{k}" in options file {args.options_file}')
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
if args.provenance:
    provenance = args.provenance
else:
    provenance = None

creator = CatalogCreator(parts, area_partition=None,
                         skycatalog_root=skycatalog_root,
                         catalog_dir=args.catalog_dir,
                         config_path=args.config_path,
                         catalog_name=args.catalog_name,
                         mag_cut=args.galaxy_magnitude_cut,
                         sed_subdir=args.sed_subdir,
                         knots_mag_cut=args.knots_magnitude_cut,
                         knots=(not args.no_knots), logname=logname,
                         skip_done=args.skip_done,
                         flux_only=args.flux_only, main_only=args.main_only,
                         flux_parallel=args.flux_parallel,
                         galaxy_nside=args.galaxy_nside,
                         galaxy_stride=args.galaxy_stride,
                         provenance=provenance,
                         dc2=args.dc2, galaxy_type=args.galaxy_type,
                         galaxy_truth=args.galaxy_truth,
                         include_roman_flux=args.include_roman_flux,
                         star_input_fmt=args.star_input_fmt)
logger.info(f'Starting with healpix pixel {parts[0]}')
if not args.no_galaxies:
    logger.info("Creating galaxy catalogs")
    creator.create('galaxy')

if not args.no_pointsources:
    logger.info("Creating point source catalogs")
    creator.create('pointsource')


logger.info('All done')
print_date()
