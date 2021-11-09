
'''
Create instance catalog(s) from a sky catalog.
'''

import os
import argparse
from desc.skycatalogs.catalog_creator import CatalogCreator
from desc.skycatalogs.utils.common_utils import print_date, print_callinfo


area_partition = {'type' : 'healpix', 'ordering' : 'ring', 'nside' : 32}

parser = argparse.ArgumentParser(description='''
Translate sky catalogs to equivalent instance catalogs
single healpix pixel''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input-dir', help='directory containing sky catalogs')

parser.add_argument('--object-types', choices=['star', 'galaxy', 'disk',
                                               'bulge'], nargs='*',
                    help='''Objects for which instance catalogs should be
                           created.  Note galaxy is shorthand for disk & bulge''', default=['star', 'galaxy'])
                    help='if used, create point source catalog(s)')
parser.add_argument('--pixels', type=int, nargs='*', default=[9556],
                    help='healpix pixels for which catalogs will be created')
parser.add_argument('output-dir', help='directory for output files',
                    default='./')
parser.add_argument('--band', choices=['u', 'g', 'r', 'i', 'z', 'y']
                    default=['r'], help='lsst band for the visit')
parser.add_argument('--verbose', help='print more output if true',
                    action='store_true')

args = parser.parse_args()

print_callinfo('translate_sc', args)

output_dir = args.output_dir

#rsdir = None
#if args.random_sed_dir:
#    rsdir = args.random_sed_dir
parts = args.pixels

#creator = CatalogCreator(parts, area_partition)
creator = CatalogCreator(parts, area_partition, output_dir=output_dir,
                         verbose=args.verbose,
                         mag_cut=args.galaxy_magnitude_cut,
                         random_sed_dir=args.random_sed_dir)
print('Starting with healpix pixel ', parts[0])
if not args.no_galaxies:
    print("Creating galaxy catalogs")
    creator.create('galaxy')

if args.pointsource:
    print("Creating point source catalogs")
    creator.create('pointsource')

print('All done')
print_date()
