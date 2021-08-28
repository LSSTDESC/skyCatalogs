from collections import namedtuple
import os
from astropy import units as u
import h5py
import pandas as pd

import numpy as np
import numpy.ma as ma
from numpy.random import default_rng
from desc.skycatalogs.utils.config_utils import Tophat
from desc.skycatalogs.utils.common_utils import print_date

__all__ = ['SCALE_FACTOR', 'LookupInfo', 'Cmp']

# Can more or less use second equation under
# https://en.wikipedia.org/wiki/AB_magnitude#Expression_in_terms_of_f%CE%BB
# This is almost our case except for us f_nu and f_lam are not spectral flux densities:
# they're just "spectral flux".  So the scale factor is something like
#    (1 / (units for tophat value)) * c    (c in cm/sec I guess)
#    (1/4.4659) * 10^(-13) * 2.99792 * 10^10  = 6.7129 * 10^(-4)

SCALE_FACTOR = 6.7129e-4


def _convert_tophat_sed(a_bins, f_nu, granularity=None):
    '''
    Given a tophat SED, produce an equivalent SED as list of (wavelength, f_lambda)
    Parameters
    ----------
    a_bins: list of Tophat [tuples (start, width)] in Angstroms
    values: list of  f_nu
    granularity: spacing between SED values.  If None, use input start values/end values
                 plus midpoints of each bin, translated from Angstrom to nm

    What to do about extrapolation?   Tophat lambdas (in nm) range from 100 to 1740.2 + 259.6,
    so about 2000.
    A SED file has min lambda = 9.1, max = 160000 (but starting at 20000 bins are 20000 wide.
    But LSST filters only cover a range comfortably within [300 nm, 1100 nm] so this shouldn't be
    an issue.   Can just start the SED file at 300 nm and run to 1100. In this case,
    don't even use the first 5 tophats or the last  4.

    return: arrays lambda, f_lambda where lambda is in nm and f_lambda
            is in erg / (cm**2 * s * ang)
            Or return lambda, f_lambda needing to be scaled, and also the scale factor
    '''


    # For each tophat value
    # Choose lambda to use.  For first draft, just use midpoint
    lam_a = np.array([ b.start + 0.5*b.width for b in a_bins])     # Angstroms
    f_nu = np.array(f_nu)

    #        convert from f_nu to f_lambda:
    # Up to a constant - universal for all SEDs - all I need to do is divide by lambda^2
    f_lam = f_nu / (lam_a * lam_a)

    if (lam_a[0] > lam_a[1]):     # reverse
        lam_a[:] = lam_a[::-1]
        f_lam[:] = f_lam[::-1]

    magnorm = 1.0         # for now
    return 0.1 * lam_a, f_lam, magnorm


def _write_sed_file(path, wv, f_lambda, wv_unit=None, f_lambda_unit=None):
    '''
    Write a two-column text file.  First column is wavelength, second is luminosity value
    If units are supplied, write a comment line at the top
    Parameters
    ----------
    path           Where to write the file and what to call it
    wv             List or array of wavelength values
    f_lambda       List or array of luminosities.  Must be the same length as wv
    wv_unit        String describing units for first column
    f_lambda_unit  String describing units for second column
    '''

    header = '#  '
    if wv_unit:
        header += wv_unit + ' '
    else:
        header += ' lambda unit unknown '
    if f_lambda_unit:
        header += f_lambda_unit
    else:
        header += ' f_lambda unit unknown'
    header += '\n'
    with open(path, mode="w") as f:
        f.write(header)
        for i in range(len(wv)):
            line = '{:8.2f}  {:g}\n'.format(wv[i], f_lambda[i])
            #line = f'{wv[i]}   {f_lambda[i]}\n'
            f.write(line)
    f.close()

class LookupInfo(object):
    '''
    Stash information from the lookup file for a particular hp which
    will be useful for Cmp class
    '''
    def __init__(self, sed_library_dir, hp):
        self.sed_lookup_file = os.path.join(sed_library_dir,
                                            f'sed_fit_{hp}.h5')
        self.sed_names = None

    def cache_info(self):
        if self.sed_names:   return

        with h5py.File(self.sed_lookup_file) as f:
            self.sed_names = f['sed_names']
            self.disk_sed = f['disk_sed']
            self.bulge_sed = f['bulge_sed']
            self.galaxy_id = f['galaxy_id']

    def get_component_sed_file(self, cmp, galaxy_id, min_ix=0):
        # Start searching for galaxy_id starting with min_ix
        the_ix = -1
        if cmp not in ['bulge', 'disk']:
            raise ValueError(f'Unknown component type "{cmp}" ')
        for i in range(min_ix, len(self.galaxy_id)):
            if self.galaxy_id[i] == galaxy_id:
                the_ix = i
                break
        if the_ix == -1:
            raise ValueError(f'Galaxy {galaxy_id} not found')

        if cmp == 'bulge':
            return self.sed_names[self.bulge_sed[the_ix]]
        else:
            return self.sed_names[self.disk_sed[the_ix]]


class Cmp(object):
    '''
    Handle writing of SED files and booking for either disk or bulge
    '''
    def __init__(self, cmp_name, obj_coll, output_dir, hp, n_seds, bins,
                 lookup_info):
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
        lookup_info  LookupInfo information
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

        (lmbda,f_lambda,magnorm) = _convert_tophat_sed(bins, sed_list)
        _write_sed_file(outpath, lmbda, f_lambda, wv_unit='nm')
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
