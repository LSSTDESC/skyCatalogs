import os
import sys
import argparse
from desc.skycatalogs.skyCatalogs import *
from desc.skycatalogs.utils.sed_utils import Cmp, LookupInfo
from desc.skycatalogs.utils.common_utils import *
from desc.skycatalogs.utils.config_utils import *

if __name__ == "__main__":
    '''
    Write SED files (bulge & disk) derived from tophat model in cosmoDC2
    The output file should look like the ones ImSim is expecting, e.g.

    # Wavelength (nm)   F_lamA (normalized erg/cm2/s/A)
    9.100 0.003413
    9.400 0.005055
        ...

    but with only 30 lines of data.
    SED files should express wavelength in nanometers, but cosmoDC2
    uses Angstroms (10 angstroms = 1 nanometer)

    F_lamA is  "spectral flux density": energy per unit time per unit area per
    unit freq.

    Also write out tables (one for bulge, one for disk) with columns
    row_number   skycat_index   galaxy_id     tophat_values   magnorm

    Here magnorm is 1.0 unless we've rescaled values before writing SED file.

    '''
    parser = argparse.ArgumentParser(description='Write text files in format expected for SEDs where input SEDs come from a Sky Catalog in parquet format')
    parser.add_argument('skycatalog_config', type=str,
                        help='path to config for Sky Catalog to be read')
    parser.add_argument('--healpix', type=int, default='9556',
                        help='healpix from which SEDs will be read')
    parser.add_argument('--n-seds', type=int, default='1',
                        help='# of seds (each for bulge and disk) to be written to individual files')
    parser.add_argument('--output-dir', type=str, default='.',
                        help='path to directory where files are written')
    parser.add_argument('--count-start', type=int, default='0',
                        help='initial N to use when writing files named fake_bulge_sed_HP_N.txt, fake_disk_sed_HP_N')

    args = parser.parse_args()
    print_callinfo('write_sed_files', args)

    # Get a LookupInfo object
    lookup = LookupInfo(os.getenv('SIMS_SED_LIBRARY_DIR'), args.healpix)

    # open sky catalog,  get object list for our hp
    cat = open_catalog(args.skycatalog_config)
    cfg = open_config_file(args.skycatalog_config)
    bins = cfg.get_tophat_parameters()

    obj_list = cat.get_objects_by_hp(0, args.healpix, None, set(['galaxy']))
    collect = obj_list.get_collections()[0]

    cmp_bulge = Cmp('bulge', collect, args.output_dir, args.healpix,
                    args.n_seds, bins, lookup)
    cmp_bulge.create(args.count_start)

    cmp_disk = Cmp('disk', collect, args.output_dir, args.healpix, args.n_seds,
                   bins, lookup)
    cmp_disk.create(args.count_start)
