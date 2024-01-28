# Use DiffskySedGenerator class
# For each healpixel for which SEDs are to be created, main skyCatalogs
# parquet file must already exist
'''
Create sky catalogs for one or more healpixels. Invoke with --help
for details
'''

import os
import argparse
import logging
import yaml
from skycatalogs.utils.common_utils import print_date, log_callinfo

parser = argparse.ArgumentParser(description='''
Create SEDs for diffsky galaxies.''',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--pixels', type=int, nargs='*', default=[9556],
                    help='healpix pixels for which catalogs will be created')
parser.add_argument('--config-path', help='path to a skyCatalogs config file')
parser.add_argument('--skycatalog-root',
                    help='''Root directory for sky catalogs, typically
                    site- or user-dependent. If not specified, use value of
                    environment variable SKYCATALOG_ROOT or value from config
                    file''')
parser.add_argument('--catalog-dir', '--cat-dir',
                    help='output file directory relative to skycatalog_root',
                    default='.')
parser.add_argument('--output-dir',
                    help='''If specified output SEDs here. Else write to
                    dir used for diffsky catalogs''')
parser.add_argument('--log-level', help='controls logging output',
                    default='INFO', choices=['DEBUG', 'INFO', 'WARNING',
                                             'ERROR'])
parser.add_argument('--galaxy-truth', help='''GCRCatalogs name for galaxy
                    truth catalog.  If not specified, default will be used''')
parser.add_argument('--rel-err', type=float, default=0.03, help='''
                    target tolerance for flux error''')
parser.add_argument('--wave-ang-min', type=int, default=500,
                    help='''Min wavelength to keep in SEDs''')
parser.add_argument('--wave-ang-max', type=int, default=100000,
                    help='''Max wavelength to keep in SEDs''')
parser.add_argument('--n-per', type=int, default=100000, help='''
                    number of galaxies to be processed together''')
parser.add_argument('--overwrite', action='store_true',
                    help='''If supplied overwrite existing data files;
                    else skip with message''')
parser.add_argument('--options-file', default=None, help='''
                    path to yaml file associating option names with values.
                    Values for any options included will take precedence.''')

args = parser.parse_args()

if args.options_file:
    with open(args.options_file) as f:
        opt_dict = yaml.safe_load(f)
        for k in opt_dict:
            if k in args:
                args.__setattr__(k, opt_dict[k])
            else:
                raise ValueError(f'Unknown attribute "{k}" in options file {args.options_file}')
logname = 'diffsky_sed.creator'
logger = logging.getLogger(logname)
logger.setLevel(args.log_level)

ch = logging.StreamHandler()
ch.setLevel(args.log_level)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)

logger.addHandler(ch)

log_callinfo('create_diffsky_sed', args, logname)

skycatalog_root = args.skycatalog_root
if not skycatalog_root:
    skycatalog_root = os.getenv('SKYCATALOG_ROOT')

from skycatalogs.skyCatalogs import open_catalog

sky_cat = open_catalog(args.config_path, skycatalog_root=skycatalog_root)
from skycatalogs.diffsky_sedgen import DiffskySedGenerator

# hard-code for now.  Should be able to retrieve galaxy truth from sky_cat
galaxy_truth = args.galaxy_truth
if galaxy_truth is None:
    galaxy_truth = 'roman_rubin_2023_v1.1.2_elais'

creator = DiffskySedGenerator(logname=logname, galaxy_truth=galaxy_truth,
                              output_dir=args.output_dir, sky_cat=sky_cat,
                              skip_done=(not args.overwrite),
                              rel_err=args.rel_err,
                              wave_ang_min=args.wave_ang_min,
                              wave_ang_max=args.wave_ang_max,
                              n_per=args.n_per, sed_out=args.output_dir)

for p in args.pixels:
    creator.generate_pixel(p)
    logger.info(f'Done with pixel {p}')

logger.info('All done')
print_date()
