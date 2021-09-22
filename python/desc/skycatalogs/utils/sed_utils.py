from collections import namedtuple
import os
from astropy import units as u
import h5py
import pandas as pd

import numpy as np
import numpy.ma as ma
from numpy.random import default_rng
from lsst.sims.photUtils import Sed, Bandpass, CosmologyObject
from desc.skycatalogs.utils.config_utils import Tophat
from desc.skycatalogs.utils.common_utils import print_dated_msg
import GCRCatalogs

__all__ = ['SCALE_FACTOR', 'LookupInfo', 'Cmp']

# Can more or less use second equation under
# https://en.wikipedia.org/wiki/AB_magnitude#Expression_in_terms_of_f%CE%BB
# This is almost our case except for us f_nu and f_lam are not spectral flux
# densities; they're just "spectral flux".  So scale factor is something like
#    (1 / (units for tophat value)) * c    (c in cm/sec I guess)
#    (1/4.4659) * 10^(-13) * 2.99792 * 10^10  = 6.7129 * 10^(-4)
##### SCALE_FACTOR = 6.7129e-4
# or just let Sed class convert f_nu to f_lambda

# According to
# https://github.com/LSSTDESC/gcr-catalogs/blob/master/GCRCatalogs/SCHEMA.md
# Tophat value units are 4.4659e13 W/Hz
TOPHAT_SCALE = 4.4659e13

# Really these values should be be looked up from the galaxy catalog
# (e.g. cosmoDC2) If the catalog is loaded with GCRCatalogs (or, I assume,
# just opened with standard  hdf5 open), use
# cat.cosmology.H0.value and cat.cosmology.Om0
_H0 = 71.0
_Om0 = 0.2648

def _convert_tophat_sed(a_bins, f_nu, redshift=0, wavelen_step=5.0):
    '''
    Given a tophat SED and redshift, produce an equivalent SED as lists of
    wavelength and f_lambda.   Also return base magnorm and magnorm adjusted
    for redshift
    Parameters
    ----------
    a_bins: list of Tophat [tuples (start, width)] in Angstroms
    f_nu: list of values for the tophats
    redshift:   needed for computing distance modulus. Really should be
                redshift_hubble
    wavelen_step:   Used for resampling

    What to do about extrapolation?   Tophat lambdas (in nm) range from 100 to
    1740.2 + 259.6, so about 2000.
    A SED file has min lambda = 9.1, max = 160000 (but starting at 20000 bins
    are 20000 wide.     But LSST filters only cover a range comfortably within
    [300 nm, 1100 nm] so this shouldn't be an issue.   Can just start the
    SED file at 300 nm and run to 1100. In this case, wouldn't even use the
    first 5 tophats or the last  4.

    return
    ------
    arrays lambda, f_lambda where lambda is in nm and f_lambda is in
    erg / (cm**2 * s * nm)
    Also return base magnorm and final magnorm (base + distance modulus)
    '''

    lam_a = np.array([b.start for b in a_bins])
    lam_nm = 0.1 * lam_a
    lam_width_nm = 0.1 * np.array([b.width for b in a_bins])
    f_nu = TOPHAT_SCALE * np.array(f_nu)

    # Convert from f_nu to f_lambda:
    # Up to a constant - universal for all SEDs - all I need to do is divide
    # by lambda^2
    # But I'll let Sed class do ithis instead

    if (lam_nm[0] > lam_nm[1]):     # reverse
        lam_nm[:] = lam_nm[::-1]
        lam_width_nm = lam_width_nm[::1]
        f_nu[:] = f_nu[::-1]


    # Calculate magnorm.  Taken from SedFitter._create_library_one_sed in
    # sims_GCRCatSimInterface
    lam_min = lam_nm[0]
    lam_max = lam_nm[-1] + lam_width_nm[-1]
    base_spec = Sed(wavelen=lam_nm, fnu=f_nu)
    print(f'lam_min: {lam_min}  lam_max: {lam_max} wavelen_step: {wavelen_step}')
    base_spec.resampleSED(wavelen_min=lam_min, wavelen_max=lam_max,
                          wavelen_step=wavelen_step)

    print('resampled len: ', len(base_spec.wavelen))
    print(f'min and max: {base_spec.wavelen[0]}  {base_spec.wavelen[-1]}')

    imsim_bp = Bandpass()
    imsim_bp.imsimBandpass()
    magnorm_base = base_spec.calcMag(imsim_bp)

    # Now adjust for redshift,  a per-object adjustment.
    # Temporarily use redshift rather than redshift_true.  They're
    # not very different

    cosmo = CosmologyObject(H0=_H0, Om0=_Om0)
    distance_modulus = cosmo.distanceModulus(redshift=redshift)
    magnorm = magnorm_base + distance_modulus
    print(f'our base magnorm: {magnorm_base}   distmod:  {distance_modulus}  sum: {magnorm}')
    #print(f'redshift: {redshift}')

    return base_spec.wavelen, base_spec.flambda, magnorm_base, magnorm

def _write_sed_file(path, wv, f_lambda, wv_unit=None, f_lambda_unit=None):
    '''
    Write a two-column text file.  First column is wavelength,
    second is luminosity value
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
        self.cached = False

    def cache_info(self):
        if self.cached:   return

        with h5py.File(self.sed_lookup_file) as f:
            # Make a copy which will exist after file is closed
            self.sed_names = np.array(f['sed_names'])
            self.disk_sed = np.array(f['disk_sed'])
            self.bulge_sed = np.array(f['bulge_sed'])
            self.galaxy_id = np.array(f['galaxy_id'])
            self.redshift = np.array(f['galaxy_id'])
            self.cached = True

    def get_orig_sed_file(self, cmp, galaxy_id, min_ix=0):
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
            return (self.sed_names[self.bulge_sed[the_ix]]).decode("utf-8")
        else:
            return (self.sed_names[self.disk_sed[the_ix]]).decode("utf-8")


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
        obj_coll     object collection coming from sky catalog, typically all
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
        lookup_info.cache_info()
        self.lookup_info = lookup_info

    def _write_sed(self, outpath, sed_list, bins, redshift,
                   wavelen_step=5.0, summary_only=False):
        '''
        Convert cosmoDC2-style tophat SEDs to a file of the form expected by
        ImSim.
        Parameters
        ----------
        outpath   string              full path of output file
        sed_list  list of floats      list of values as they appear in
                                      cosmoDC2 catalog
        bins      list((start,width)) bin definitions
        redshift  -- for the object the sed file is associated with

        Return
        ------
        (magnorm, magnorm_adjust, norm_wv_value)  first is computed base magnorm
                   the second has been adjusted for redshift.
                   norm_wv_value is the sed value at or near 500 nm
        '''
        (lmbda,f_lambda,magnorm,
         magnorm_adjust) = _convert_tophat_sed(bins, sed_list,
                                               redshift=redshift,
                                               wavelen_step=wavelen_step)
        if not summary_only:
             _write_sed_file(outpath, lmbda, f_lambda, wv_unit='nm')
        start = (min([b.start for b in bins]))/10.0        # A to nm
        normwv_ix = int(np.floor((500 - start)/wavelen_step))
        return (magnorm, magnorm_adjust, f_lambda[normwv_ix])     # for now

    def _write_summary(self, ix, gal, sed, redshift, orig_magnorm, our_magnorm,
                       our_magnorm_adjust, orig_sed_file, norm_wv_value,
                       tp_sed_file):
        # Filepath.  Use same output dir.
        print_dated_msg(f'Entered _write_summary for component {self.cmp_name}')
        basename_csv = f'{self.cmp_name}_sed_hp{self.hp}_summary.csv'
        outpath_csv = os.path.join(self.output_dir, basename_csv)
        basename_pq = f'{self.cmp_name}_random_sed_hp{self.hp}_summary.parquet'
        outpath_pq = os.path.join(self.output_dir, basename_pq)

        # transpose SEDs
        ##sed_transposed = [[row[i] for row in sed] for i in range(len(sed[0]))]
        # form column names
        ##sed_col_names = [f'bin_{b.start}_{b.width}' for b in self.bins]

        out_dict = {}
        out_dict['chosen_ix'] = ix
        out_dict['gal_id'] = gal
        out_dict['redshift'] = redshift
        out_dict['orig_magnorm'] = orig_magnorm
        out_dict['our_magnorm'] = our_magnorm
        out_dict['our_magnorm_adjust'] = our_magnorm_adjust
        out_dict['orig_sed_file'] = orig_sed_file
        out_dict['norm_wv_value'] = norm_wv_value
        ##out_dict['tp_sed_file'] = tp_sed_file
        ##for ib in range(len(sed_col_names)):
        ##    out_dict[sed_col_names[ib]] = sed_transposed[ib]

        df = pd.DataFrame(data=out_dict)
        df.to_csv(path_or_buf=outpath_csv)
        df.to_parquet(outpath_pq)

    def create(self, count_start=0, summary_only=False):
        '''
        Create SED files as specified at init time and also table describing
        which tophat SEDs were used.
        count_start may be > 0 in case some of the required files have already
        been created and we just want to pick up where we left off.
        [But initial draft won't support this since there are complications]
        '''

        # For debugging predictability
        seed_dict = {}
        seed_dict['bulge'] = 135711 + 2 * self.hp
        seed_dict['disk'] = 890123 + 2 * self.hp

        print_dated_msg(f'Cmp.create called for component  {self.cmp_name}')
        sed_col = 'sed_val_' + self.cmp_name
        sed = np.array(self.coll.get_attribute(sed_col))
        magnorm_col = self.cmp_name + '_magnorm'
        magnorm = np.array(self.coll.get_attribute(magnorm_col))
        gal_id = np.array(self.coll.get_attribute('galaxy_id'))
        redshift = np.array(self.coll.get_attribute('redshift'))

        # For distance modulus calculation really should use redshift_true,
        # but the difference is small. It's not currently in Sky Catalogs but
        # plan to put it there, probably named redshift_hubble
        # Meanwhile could do the following
        # _GAL_CAT = 'cosmoDC2_v1.1.4_image'
        # cosmodc2 = GCRCatalogs.load_catalog(_GAL_CAT)
        # cosmo_df = cosmodc2.get_quantities(['galaxy_id', 'redshift_true'],
        #                        native_filters=f'healpix_pixel=={self.hp}')
        # for i in range(len(gal_id)):
        #     assert cosmo_df['galaxy_id'][i] == gal_id[i]
        # print(f'Up to len Sky Catalog for hp {self.hp} Sky Catalog and cosmoDC2 match')
        # redshift_true = cosmo_df['redshift_true'][:len(gal_id)]

        # mask off anything with magnorm infinite
        mask_inf = np.isinf(magnorm)
        good_sed = ma.array(sed, mask=mask_inf).compressed()
        good_gal_id = ma.array(gal_id, mask=mask_inf).compressed()
        good_magnorm = ma.array(magnorm, mask=mask_inf).compressed()
        good_redshift = ma.array(redshift, mask=mask_inf).compressed()

        # Choose entries at random
        rng = default_rng(seed_dict[self.cmp_name])
        ix_list = rng.integers(low=0, high=len(good_magnorm), size=self.n_seds)
        gal_chosen = [good_gal_id[i] for i in ix_list]
        sed_chosen = [good_sed[i] for i in ix_list]
        orig_magnorm_chosen = [good_magnorm[i] for i in ix_list]
        redshift_chosen = [good_redshift[i] for i in ix_list]
        our_magnorm = []
        our_magnorm_adjust = []
        norm_wv_value = []
        orig_sed_file = []
        tp_sed_file = []

        sed_rootdir = os.getenv('SIMS_SED_LIBRARY_DIR')
        for i in range(len(sed_chosen)):
            # Form output path
            filename = f'{self.cmp_name}_random_sed_{self.hp}_{i}.txt'
            outpath = os.path.join(self.output_dir, filename)

            (our_mag, our_mag_adjust,
             nm500) = self._write_sed(outpath, sed_chosen[i], self.bins,
                                      redshift_chosen[i], summary_only=summary_only)
            our_magnorm.append(our_mag)
            our_magnorm_adjust.append(our_mag_adjust)
            norm_wv_value.append(nm500)

            tp_sed_file.append(outpath)
            orig_sed = self.lookup_info.get_orig_sed_file(self.cmp_name,
                                                          gal_chosen[i],
                                                          min_ix=ix_list[i])
            orig_sed_file.append(os.path.join(sed_rootdir, orig_sed))

            print_dated_msg(f'Wrote file {i}')

        # Make summary table and write to a file
        self._write_summary(ix_list, gal_chosen, sed_chosen, redshift_chosen,
                            orig_magnorm_chosen, our_magnorm, our_magnorm_adjust,
                            norm_wv_value, orig_sed_file, tp_sed_file)
