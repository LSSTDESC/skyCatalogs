import os
import sys
import argparse
import numpy as np
import numpy.ma as ma
from numpy.random import default_rng
import pandas as pd
import yaml
from desc.skycatalogs.skyCatalogs import *
from desc.skycatalogs.utils.sed_utils import *
from desc.skycatalogs.utils.common_utils import *

class Cmp(object):
    '''
    Handle writing of SED files and booking for either disk or bulge
    '''
    def __init__(self, cmp_name, obj_coll, output_dir, hp, n_seds, bins):
        '''
        Parameters
        ----------
        cmp_name     string    one of 'bulge', 'disk'
        obj_coll               object collection coming from sky catalog, typically all
                           galaxies belonging to a particular pixel
        output_dir   string    where to write output SED files
        hp           int       in case we decide to embed in output filename
        n_seds       int       how many SED files to write
        bins         list      list of (start, width) tuples describing bins.
        '''
        self.cmp_name = cmp_name
        self.output_dir = output_dir
        self.hp = hp
        self.coll = obj_coll
        self.n_seds = n_seds
        self.n_seds_done = 0
        self.bins = bins

    def _write_sed(self, outpath, sed_list, bins):
        '''
        Convert cosmoDC2-style tophat SEDs to a file of the form expected by
        ImSim.
        Parameters
        ----------
        outpath   string              full path of output file
        sed_list  list of floats      list of values as they appear in cosmoDC2 catalog
        bins      list((start,width)) bin definitions

        Return
        ------
        magnorm  (usually 1.0)
        '''

        (lmbda,f_lambda,magnorm) = convert_tophat_sed(bins, sed_list)
        write_sed_file(outpath, lmbda, f_lambda, wv_unit='nm')
        return magnorm     # for now

    def _write_summary(self, ix, gal, sed, orig_magnorm, our_magnorm):
        # Filepath.  Use same output dir.
        print('Entered _srite_summary for component ', self.cmp_name)
        print_date()
        basename = f'{self.cmp_name}_random_sed_hp{self.hp}_summary.csv'
        outpath = os.path.join(self.output_dir, basename)

        # transpose SEDs
        sed_transposed = [[row[i] for row in sed] for i in range(len(sed[0]))]

        # form column names
        sed_col_names = [f'bin_{b.start}_{b.width}' for b in self.bins]
        out_dict = {}
        out_dict['chosen_ix'] = ix
        out_dict['gal_id'] = gal
        out_dict['orig_magnorm'] = orig_magnorm
        out_dict['our_magnorm'] = our_magnorm
        for ib in range(len(sed_col_names)):
            out_dict[sed_col_names[ib]] = sed_transposed[ib]

        df = pd.DataFrame(data=out_dict)
        df.to_csv(path_or_buf=outpath)

    def create(self, count_start=0):
        '''
        Create SED files as specified at init time and also table describing which
        tophat SEDs were used.
        count_start may be > 0 in case some of the required files have already been
        created and we just want to pick up where we left off.  [but initial draft
        won't support this since there are complications]
        '''
        # For debugging predictability
        seed_dict = {}
        seed_dict['bulge'] = 135711
        seed_dict['disk'] = 890123

        print('create called for component', self.cmp_name)
        print_date()
        sed_col = 'sed_val_' + self.cmp_name
        sed = np.array(self.coll.get_attribute(sed_col))
        magnorm_col = self.cmp_name + '_magnorm'
        magnorm = np.array(self.coll.get_attribute(magnorm_col))
        gal_id = np.array(self.coll.get_attribute('galaxy_id'))

        # mask off anything with magnorm infinite
        mask_inf = np.isinf(magnorm)
        good_sed = ma.array(sed, mask=mask_inf).compressed()
        good_gal_id = ma.array(gal_id, mask=mask_inf).compressed()
        good_magnorm = ma.array(magnorm, mask=mask_inf).compressed()

        # Choose entries at random
        rng = default_rng(seed_dict[self.cmp_name])
        ix_list = rng.integers(low=0, high=len(good_magnorm), size=self.n_seds)
        print("random index list: ")
        for i in ix_list:
            print(i)
        gal_chosen = [good_gal_id[i] for i in ix_list]
        # magnorm will either always be 1 or will be calculated per galaxy
        # magnorm_chosen = [magnorm[i] for i in ix_list]
        sed_chosen = [good_sed[i] for i in ix_list]
        orig_magnorm_chosen = [good_magnorm[i] for i in ix_list]
        our_magnorm = []
        for i in range(len(sed_chosen)):
            # Form output path
            filename = f'{self.cmp_name}_random_sed_{self.hp}_{i}.txt'
            outpath = os.path.join(self.output_dir, filename)
            our_magnorm.append(self._write_sed(outpath, sed_chosen[i], self.bins))
            print('Wrote file ', i)
            print_date()

        # Make summary table and write to a file
        self._write_summary(ix_list, gal_chosen, sed_chosen, orig_magnorm_chosen,
                            our_magnorm)

def _get_tophat_parameters(config):
    '''
    Return list of named tuples
    Should maybe be part of Sky Catalogs API
    '''
    with open(config) as f:
        d = yaml.safe_load(f)
        raw_bins = d['SED_models'][0]['tophat']['bins']

    return [ Tophat(b[0], b[1]) for b in raw_bins]

if __name__ == "__main__":
    '''
    Write one or more SED files (bulge & disk) derived from the tophat model in cosmoDC2
    The output file should look like the ones ImSim is expecting, e.g.

    # Wavelength (nm)   F_lamA (normalized erg/cm2/s/A)
    9.100 0.003413
    9.400 0.005055
        ...

    but with only 30 lines of data.
    SED files should express wavelength in nanometers, but cosmoDC2
    uses Angstroms (10 angstroms = 1 nanometer)

    F_lamA is  "spectral flux density": energy per unit time per unit area per unit freq.

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

    # open sky catalog,  get object list for our hp
    cat = open_catalog(args.skycatalog_config)

    bins = _get_tophat_parameters(args.skycatalog_config)

    obj_list = cat.get_objects_by_hp(0, args.healpix, None, set(['galaxy']))
    collect = obj_list.get_collections()[0]

    cmp_bulge = Cmp('bulge', collect, args.output_dir, args.healpix, args.n_seds,
                    bins)
    cmp_bulge.create(args.count_start)

    cmp_disk = Cmp('disk', collect, args.output_dir, args.healpix, args.n_seds, bins)
    cmp_disk.create(args.count_start)
