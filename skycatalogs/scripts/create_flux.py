'''
Create flux sky catalogs for one or more healpixels, single object type.
Invoke with --help for details
'''
import os
import argparse
import logging
import yaml
from skycatalogs.flux_catalog_creator import FluxCatalogCreator
from skycatalogs.utils.common_utils import print_date, log_callinfo
from skycatalogs.utils.common_utils import callinfo_to_dict

parser = argparse.ArgumentParser(
    description='''
    Create flux Sky Catalogs for specified object type''',
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
parser.add_argument('--config-path', default=None, help='''
                    Directory containing config file. If no value,
                    config is assumed to be in same location as data,
                    with filenmame taken from catalog-name.''')
parser.add_argument('--catalog-name', '--cat-name', default='skyCatalog',
                    help='''Name of config file describing the catalog''')
parser.add_argument('--log-level', help='controls logging output',
                    default='INFO', choices=['DEBUG', 'INFO', 'WARNING',
                                             'ERROR'])

parser.add_argument(
    '--skip-done', action='store_true',
    help='If supplied skip existing data files; else overwrite with message')
parser.add_argument(
    '--flux-parallel', default=16, type=int,
    help='Number of processes to run in parallel when computing fluxes')
parser.add_argument('--options-file', default=None, help='''
                    path to yaml file associating option names with values.
                    Values for any options included will take precedence.''')
parser.add_argument('--dc2', action='store_true',
                    help='''If supplied provide values comparable to those
                            used for DC2 run.''')
parser.add_argument(
    '--include-roman-flux', action='store_true',
    help='If supplied calculate & store Roman as well as LSST  fluxes')
parser.add_argument('--sso-sed', default=None, help='''
                    path to sqlite file containing SED to be used
                    for all SSOs''')

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

log_callinfo('create_flux', args, logname)

skycatalog_root = args.skycatalog_root
if not skycatalog_root:
    skycatalog_root = os.getenv('SKYCATALOG_ROOT')

parts = args.pixels

opt_dict = callinfo_to_dict(args)

creator = FluxCatalogCreator(args.object_type, parts,
                             skycatalog_root=skycatalog_root,
                             catalog_dir=args.catalog_dir,
                             config_path=args.config_path,
                             catalog_name=args.catalog_name,
                             logname=logname,
                             skip_done=args.skip_done,
                             flux_parallel=args.flux_parallel,
                             dc2=args.dc2,
                             include_roman_flux=args.include_roman_flux,
                             sso_sed=args.sso_sed,
                             run_options=opt_dict)
if len(parts) > 0:
    logger.info(f'Starting with healpix pixel {parts[0]}')

creator.create()

logger.info('All done')
print_date()
