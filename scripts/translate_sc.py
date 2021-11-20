
'''
Create instance catalog(s) from a sky catalog.
'''

import os
import argparse
#from desc.skycatalogs.catalog_creator import CatalogCreator
from desc.skycatalogs.utils.common_utils import print_date, print_callinfo
from desc.skycatalogs.translate import Translator
from desc.skycatalogs.skyCatalogs import Disk

area_partition = {'type' : 'healpix', 'ordering' : 'ring', 'nside' : 32}

parser = argparse.ArgumentParser(description='''
Translate sky catalogs to equivalent instance catalogs for a particular
visit. Currently the only observation parameter of interest is band.
''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('visit', type=int,
                    help='visit to which generated files belong')
parser.add_argument('config_path', help='config for source sky catalogs')

parser.add_argument('outputdir', help='directory for output files')

parser.add_argument('--object-types', choices=['star', 'galaxy', 'disk',
                                               'bulge'], nargs='*',
                    metavar='OBJ_TYPE',
                    help='''Objects for which instance catalogs should be
                           created.  Note galaxy is shorthand for disk & bulge''',
                    default=['star', 'galaxy'])
parser.add_argument('--pixels', type=int, nargs='*', default=[9556],
                    metavar='PIXEL',
                    help='healpix pixels for which catalogs will be created')
parser.add_argument('--region', type=float, nargs=3,
                    metavar=('RA_deg', 'DEC_deg', 'RADIUS_as'),
                    help='''Value must be ra dec radius for a circular region
 where ra and dec define the center. ra and dec are in units of degrees;
 radius is in arcseconds.
 If supplied, this option overrides the pixels option.''')
parser.add_argument('--band', choices=['u', 'g', 'r', 'i', 'z', 'y'],
                    default='r', help='lsst band for the visit')
parser.add_argument('--verbose', help='print more output if true',
                    action='store_true')
parser.add_argument('--clear', help='remove any old files for this visit',
                    action='store_true')

args = parser.parse_args()

print_callinfo('translate_sc', args)

translator = Translator(args.visit, args.config_path, args.outputdir,
                        object_types=args.object_types,
                        band=args.band, verbose=args.verbose,
                        clear=args.clear)

if args.region:
    disk = Disk(args.region[0], args.region[1], args.region[2])
else:
    disk = None

translator.translate_visit(args.pixels, disk)

# print('Starting with healpix pixel ', parts[0])
# if not args.no_galaxies:
#     print("Creating galaxy catalogs")
#     creator.create('galaxy')

# if args.pointsource:
#     print("Creating point source catalogs")
#     creator.create('pointsource')

print('All done')
print_date()
