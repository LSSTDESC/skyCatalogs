from collections import namedtuple
import os
import re
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
import h5py
import pandas as pd

import numpy as np
import numpy.ma as ma
from numpy.random import default_rng

from desc.skycatalogs.utils.common_utils import print_dated_msg

__all__ = ['LookupInfo', 'Cmp', 'MagNorm', 'convert_tophat_sed',
           'write_sed_file', 'NORMWV_IX', 'get_star_sed_path',
           'create_cosmology']

# Index for tophat bin containing 500 nm
NORMWV_IX = 13

def convert_tophat_sed(a_bins, f_nu_input, mag_norm_f, redshift=0,
                       wavelen_step=0.1):
    '''
    Given a tophat SED and redshift, produce an equivalent SED as lists of
    wavelength and f_lambda. Also compute magnorm

    Parameters
    ----------
    a_bins: list of Tophat [tuples (start, width)] in Angstroms
    f_nu: list of values for the tophats
    mag_norm_f: an instance of MagNorm
    redshift:   needed for computing distance modulus. Should be
                cosmoDC2 redshiftHubble, aka redshift_hubble in sky catalogs
    wavelen_step:    Re-cast tophat seds to use this bin width in nm (keeping
                     same step function in f_nu space).

    return
    ------
    arrays lambda, f_lambda where lambda is in nm and f_lambda is in
    erg / (cm**2 * s * nm)
    Also return final magnorm (including redshift adjustment) and f_nu value
    at 500 nm
    '''

    lam_nm = 0.1 * np.array([b.start + 0.5 * b.width for b in a_bins])
    lam_width_nm = 0.1 * np.array([b.width for b in a_bins])
    f_nu = 1.0 * np.array(f_nu_input)
    val_500nm = f_nu[NORMWV_IX]

    # Convert from f_nu to f_lambda:
    # In earlier versions tophats were in decreasing lambda order
    if (lam_nm[0] > lam_nm[1]):     # reverse
        lam_nm[:] = lam_nm[::-1]
        lam_width_nm[:] = lam_width_nm[::-1]
        f_nu[:] = f_nu[::-1]

    lam_min = lam_nm[0]
    lam_max = lam_nm[-1] + lam_width_nm[-1]

    # Keep the same step function but use fine bins instead of the
    # original tophat widths.
    n_bins = int((lam_max - lam_min) / wavelen_step)
    lam_fine = np.empty(n_bins)
    f_nu_fine = np.empty(n_bins)
    boundaries = list(lam_nm)
    boundaries.append(lam_max)
    b_ix = 0
    for i in range(n_bins):
        lam_fine[i] = lam_min + wavelen_step * i
        if (lam_fine[i] > boundaries[b_ix + 1]) :
            b_ix = b_ix + 1
        f_nu_fine[i] = f_nu[b_ix]

    # Convert fnu to flambda, ignoring constant factors.
    flambda = f_nu_fine/lam_fine**2

    # Normalize so flambda value at 500 nm is 1.0
    nm500_ix = int((500 - lam_min) / wavelen_step) + 1
    flambda_norm = flambda / flambda[nm500_ix]

    return lam_fine, flambda_norm, mag_norm_f(f_nu[NORMWV_IX],
                                             redshift), val_500nm

def write_sed_file(path, wv, f_lambda, wv_unit=None, f_lambda_unit=None):
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
            f.write(line)
    f.close()

_standard_dict = {'lte' : 'starSED/phoSimMLT',
                  'bergeron' : 'starSED/wDs',
                  'km|kp' : 'starSED/kurucz'}

def get_star_sed_path(filename, name_to_folder=_standard_dict):
    '''
    Return numpy array of full paths relative to SIMS_SED_LIBRARY_DIR,
    given filenames

    Parameters
    ----------
    filename       list of strings. Usually full filename but may be missing final ".gz"
    name_to_folder dict mapping regular expression (to be matched with
                   filename) to relative path for containing directory

    Returns
    -------
    Full path for file, relative to SIMS_SED_LIBRARY_DIR
    '''

    compiled = { re.compile(k) : v for (k, v) in name_to_folder.items()}

    path_list = []
    for f in filename:
        m = None
        matched = False
        for k,v in compiled.items():
            f = f.strip()
            m = k.match(f)
            if m:
                p = os.path.join(v, f)
                if not p.endswith('.gz'):
                    p = p + '.gz'
                path_list.append(p)
                matched = True
                break

        if not matched:
            raise ValueError(f'get_star_sed_path: Filename {f} does not match any known patterns')
    return np.array(path_list)


def create_cosmology(config):
    """
    Create a FlatLambdaCDM cosmology from a dictionary of input parameters.
    This code is based on/borrowed from
    https://github.com/LSSTDESC/gcr-catalogs/blob/master/GCRCatalogs/cosmodc2.py#L128
    """
    cosmo_astropy_allowed = FlatLambdaCDM.__init__.__code__.co_varnames[1:]
    cosmo_astropy = {k: v for k, v in config.items()
                     if k in cosmo_astropy_allowed}
    cosmology = FlatLambdaCDM(**cosmo_astropy)
    return cosmology


class MagNorm:
    def __init__(self, cosmology):
        """
        Parameters
        ----------
        cosmology : astropy.cosmology
            Cosmology object created from the gcr-catalogs galaxy catalog
            cosmology specification.
        """
        self.cosmology = cosmology

    def dl(self, z):
        """
        Return the luminosity distance in units of meters.
        """
        # Conversion factor from Mpc to meters (obtained from pyccl).
        MPC_TO_METER = 3.085677581491367e+22
        return self.cosmology.luminosity_distance(z).value*MPC_TO_METER

    def __call__(self, tophat_sed_value, redshift_hubble, one_maggy=4.3442e13):
        one_Jy = 1e-26  # W/Hz/m**2
        Lnu = tophat_sed_value*one_maggy    # convert from maggies to W/Hz
        Fnu = Lnu/4/np.pi/self.dl(redshift_hubble)**2
        with np.errstate(divide='ignore', invalid='ignore'):
            mgn = -2.5*np.log10(Fnu/one_Jy) + 8.90
            return mgn


class LookupInfo(object):
    '''
    Stash information from the lookup file for a particular hp which
    will be useful for Cmp class
    Also save tophat scale
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

# This class is no longer used. Consider deleting
class Cmp(object):
    '''
    Handle writing of SED files and booking for either disk or bulge
    '''
    def __init__(self, cmp_name, obj_coll, output_dir, hp, n_seds, bins,
                 lookup_info, mag_norm_f):
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
        lookup_info  LookupInfo information pertaining to a particular hp
        mag_norm_f   MagNorm   Used for computing mag norm

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
        self.mag_norm_f = mag_norm_f

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
        (magnorm, val_500nm)  magnorm is our computed magnorm value,
                   including adjustment for redshift.
                   val_500nm is the sed value at or near 500 nm
        '''
        (lmbda, f_lambda,
         magnorm, val_500nm) = convert_tophat_sed(bins, sed_list,
                                                  self.mag_norm_f,
                                                  redshift=redshift,
                                                  wavelen_step=wavelen_step)
        if not summary_only:
             write_sed_file(outpath, lmbda, f_lambda, wv_unit='nm')
        start = (min([b.start for b in bins]))/10.0        # A to nm
        return (magnorm, val_500nm)     # for now

    def _write_summary(self, ix, gal, sed, redshift, orig_magnorm, our_magnorm,
                       val_500nm, orig_sed_file, tp_sed_file):
        # Filepath.  Use same output dir.
        print_dated_msg(f'Entered _write_summary for component {self.cmp_name}')
        basename_csv = f'{self.cmp_name}_sed_hp{self.hp}_summary.csv'
        outpath_csv = os.path.join(self.output_dir, basename_csv)
        basename_csv_brief = f'{self.cmp_name}_sed_hp{self.hp}_brief.csv'
        outpath_csv_brief = os.path.join(self.output_dir, basename_csv_brief)
        basename_pq = f'{self.cmp_name}_sed_hp{self.hp}_summary.parquet'
        outpath_pq = os.path.join(self.output_dir, basename_pq)

        out_dict = {}
        out_dict['chosen_ix'] = ix
        out_dict['gal_id'] = gal
        out_dict['redshift'] = redshift
        out_dict['orig_magnorm'] = orig_magnorm
        out_dict['our_magnorm'] = our_magnorm
        out_dict['val_500nm'] = val_500nm
        df = pd.DataFrame(data=out_dict)

        # For convenience, output text file leaving off paths
        df.to_csv(path_or_buf=outpath_csv_brief)

        out_dict['orig_sed_file'] = orig_sed_file
        out_dict['tp_sed_file'] = tp_sed_file
        out_dict['tp_vals'] = sed

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
        seed_dict['bulge'] = 271423 + 2 * self.hp
        seed_dict['disk'] = 1780247 + 2 * self.hp

        print_dated_msg(f'Cmp.create called for component  {self.cmp_name}')
        ###  Really it should have _no_host_extinction suffix
        ###                       REALLY??
        ###   but for now schema is not using it
        ### sed_col = 'sed_val_' + self.cmp_name + '_no_host_extinction'
        sed_col = 'sed_val_' + self.cmp_name
        sed = np.array(self.coll.get_native_attribute(sed_col))
        magnorm_col = self.cmp_name + '_magnorm'
        magnorm = np.array(self.coll.get_native_attribute(magnorm_col))
        gal_id = np.array(self.coll.get_native_attribute('galaxy_id'))
        redshift = np.array(self.coll.get_native_attribute('redshift_hubble'))

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
        val_500nm = []
        orig_sed_file = []
        tp_sed_file = []

        sed_rootdir = os.getenv('SIMS_SED_LIBRARY_DIR')
        for i in range(len(sed_chosen)):
            # Form output path
            filename = f'{self.cmp_name}_random_sed_{self.hp}_{i}.txt'
            outpath = os.path.join(self.output_dir, filename)

            (our_mag, nm500) = self._write_sed(outpath, sed_chosen[i],
                                               self.bins, redshift_chosen[i],
                                               summary_only=summary_only)
            our_magnorm.append(our_mag)
            val_500nm.append(nm500)

            tp_sed_file.append(outpath)
            orig_sed = self.lookup_info.get_orig_sed_file(self.cmp_name,
                                                          gal_chosen[i],
                                                          min_ix=ix_list[i])
            orig_sed_file.append(os.path.join(sed_rootdir, orig_sed))

            if not summary_only:
                print_dated_msg(f'Wrote file {i}')

        # Make summary table and write to a file
        self._write_summary(ix_list, gal_chosen, sed_chosen, redshift_chosen,
                            orig_magnorm_chosen, our_magnorm,
                            val_500nm, orig_sed_file, tp_sed_file)
