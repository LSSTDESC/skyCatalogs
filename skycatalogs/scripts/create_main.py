'''
Create "main" sky catalogs for one or more healpixels, single object type.
Invoke with --help for details
'''

import os
import numpy as np
import argparse
import logging
import yaml
from skycatalogs.main_catalog_creator import MainCatalogCreator
from skycatalogs.utils.common_utils import print_date, log_callinfo
from skycatalogs.utils.common_utils import callinfo_to_dict

parser = argparse.ArgumentParser(
    description='''
    Create main Sky Catalogs for specified object type.''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('object_type',
                    choices=['star', 'cosmodc2_galaxy', 'diffsky_galaxy',
                             'sso'],
                    help='Object type for which catalog is to be created')
parser.add_argument('--pixels', type=int, nargs='*', default=[9556],
                    help='healpix pixels for which catalogs will be created')
parser.add_argument('--skycatalog-root',
                    help='''Root directory for sky catalogs, typically
                    site-dependent. If not specified, use value of
                    environment variable SKYCATALOG_ROOT''')
parser.add_argument('--catalog-dir', '--cat-dir',
                    help='output file directory relative to skycatalog_root',
                    default='.')
parser.add_argument('--log-level', help='controls logging output',
                    default='INFO', choices=['DEBUG', 'INFO', 'WARNING',
                                             'ERROR'])
parser.add_argument('--galaxy-magnitude-cut', '--gal-mag-cut',
                    default=29.0, type=float,
                    help='''Exclude galaxies with r-magnitude above this value.
                            ignored for non-galaxy object types''')
parser.add_argument('--knots-magnitude-cut', default=27.0, type=float,
                    help='''Galaxies with i-magnitude above this cut get
                            no knots. Ignored for non-galaxy object types''')
parser.add_argument('--no-knots', action='store_true',
                    help='''If supplied omit knots component. Default is False.
                            Ignored for non-galaxy object types''')

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
parser.add_argument('--options-file', default=None, help='''
                    path to yaml file associating option names with values.
                    Values for any options included will take precedence.''')
parser.add_argument('--dc2', action='store_true',
                    help='''If supplied provide values comparable to those
                            used for DC2 run. Default is False. Applies
                            only to object type "cosmodc2_galaxy"''')
parser.add_argument('--nside', default=32, type=int,
                    choices=2**np.arange(15),
                    help='''Pixel nsides for output. Must be power of 2''')
parser.add_argument('--stride', default=1_000_000, type=int,
                    help='''max # rows in a row group; default 1 million''')
parser.add_argument('--truth', default=None, help='''Truth catalog
                    used as input. Default
                    depends on value of object_type''')
parser.add_argument('--star-input-fmt', default='sqlite',
                    choices=['sqlite', 'parquet'], help='''
                    star truth may come from either sqlite db or collection
                    of parquet files. Applies only if object_type is "star"''')
parser.add_argument('--sso-sed', default=None, help='''
                    path to sqlite file containing SED to be used
                    for all SSOs. Ignored of object_type is not sso''')

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

log_callinfo('create_main', args, logname)

skycatalog_root = args.skycatalog_root
if not skycatalog_root:
    skycatalog_root = os.getenv('SKYCATALOG_ROOT')

parts = args.pixels

opt_dict = callinfo_to_dict(args)

creator = MainCatalogCreator(args.object_type, parts,
                             skycatalog_root=skycatalog_root,
                             catalog_dir=args.catalog_dir,
                             truth=args.truth,
                             config_path=args.config_path,
                             catalog_name=args.catalog_name,
                             mag_cut=args.galaxy_magnitude_cut,
                             knots_mag_cut=args.knots_magnitude_cut,
                             knots=(not args.no_knots), logname=logname,
                             skip_done=args.skip_done,
                             nside=args.nside,
                             stride=args.stride,
                             dc2=args.dc2,
                             star_input_fmt=args.star_input_fmt,
                             sso_sed=args.sso_sed, # probably not needed
                             run_options=opt_dict)
if len(parts) > 0:
    logger.info(f'Starting with healpix pixel {parts[0]}')
elif args.object_type == "sso":
    logger.info('Creating catalogs for all available healpixels')

creator.create()

logger.info('All done')
print_date()
