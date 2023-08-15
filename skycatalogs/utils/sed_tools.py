import os
import re
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
import astropy.constants
import h5py
import pandas as pd

import numpy as np
import numpy.ma as ma
from dust_extinction.parameter_averages import F19
import galsim

__all__ = ['TophatSedFactory', 'MilkyWayExtinction', 'get_star_sed_path']
class TophatSedFactory:
    '''
    Used for modeling cosmoDC2 galaxy SEDs, which are represented with
    a small number of wide bins
    '''
    _clight = astropy.constants.c.to('m/s').value
    # Conversion factor below of cosmoDC2 tophat Lnu values to W/Hz comes from
    # https://github.com/LSSTDESC/gcr-catalogs/blob/master/GCRCatalogs/SCHEMA.md
    _to_W_per_Hz = 4.4659e13

    def __init__(self, th_definition, cosmology, delta_wl=0.001):
        # Get wavelength and frequency bin boundaries.
        if th_definition:
            bins = th_definition
            wl0 = [_[0] for _ in bins]

            # Find index of original bin which includes 500 nm == 5000 ang
            ix = -1
            for w in wl0:
                if w > 5000:
                    break
                ix += 1

            self._ix_500nm = ix

            wl0.append(bins[-1][0] + bins[-1][1])

            wl0 = 0.1*np.array(wl0)
            self.wl = np.array(wl0)
            self.nu = self._clight/(self.wl*1e-9)  # frequency in Hz

            # Also save version of wl where vertical rise is replaced by
            # steep slope
            wl_deltas = []
            for i in range(len(bins)):
                wl_deltas.extend((self.wl[i], self.wl[i+1] - delta_wl))

            # Prepend more bins which will have 0 value
            n_bins = int(wl0[0]) - 1
            pre_wl = [float(i) for i in range(n_bins)]

            wl_deltas = np.insert(wl_deltas, 0, pre_wl)

            # Also make a matching array of 0 values
            self.pre_val = [0.0 for i in range(n_bins)]

            self._wl_deltas = wl_deltas
            self._wl_deltas_u_nm = wl_deltas*u.nm

        # Create a FlatLambdaCDM cosmology from a dictionary of input
        # parameters.  This code is based on/borrowed from
        # https://github.com/LSSTDESC/gcr-catalogs/blob/master/GCRCatalogs/cosmodc2.py#L128
        cosmo_astropy_allowed = FlatLambdaCDM.__init__.__code__.co_varnames[1:]
        cosmo_astropy = {k: v for k, v in cosmology.items()
                         if k in cosmo_astropy_allowed}
        self.cosmology = FlatLambdaCDM(**cosmo_astropy)

    # Useful for getting magnorm from f_nu values
    @property
    def ix_500nm(self):
        return self._ix_500nm

    @property
    def wl_deltas(self):
        return self._wl_deltas

    @property
    def wl_deltas_u_nm(self):
        return self._wl_deltas_u_nm

    def dl(self, z):
        """
        Return the luminosity distance in units of meters.
        """
        # Conversion factor from Mpc to meters (obtained from pyccl).
        MPC_TO_METER = 3.085677581491367e+22
        return self.cosmology.luminosity_distance(z).value*MPC_TO_METER

    def create(self, Lnu, redshift_hubble, redshift, resolution=None):
        '''
        Given tophat values from cosmoDC2 produce redshifted sed.
        Does not apply extinction.
        '''
        # Compute Llambda in units of W/nm
        Llambda = (Lnu*self._to_W_per_Hz*(self.nu[:-1] - self.nu[1:])
                   /(self.wl[1:] - self.wl[:-1]))

        # Fill the arrays for the galsim.LookupTable.   Prepend
        # zero-valued bins down to mix extinction wl to handle redshifts z > 2.
        my_Llambda = []
        my_Llambda += self.pre_val
        for i in range(len(Llambda)):
            # Dealt with wl already in __init__
            my_Llambda.extend((Llambda[i], Llambda[i]))

        # Convert to (unredshifted) flux given redshift_hubble.
        flambda = np.array(my_Llambda)/(4.0*np.pi*self.dl(redshift_hubble)**2)

        # Convert to cgs units
        flambda *= (1e7/1e4)  # (erg/joule)*(m**2/cm**2)

        # Create the lookup table.
        lut = galsim.LookupTable(self.wl_deltas, flambda, interpolant='nearest')

        if resolution:
            wl_min = min(self.wl_deltas)
            wl_max = max(self.wl_deltas)
            wl_res = np.linspace(wl_min, wl_max, int((wl_max - wl_min)/resolution))
            flambda_res = [lut(wl) for wl in wl_res]
            lut = galsim.LookupTable(wl_res, flambda_res, interpolant='linear')

        # Create the SED object and apply redshift.
        sed = galsim.SED(lut, wave_type='nm', flux_type='flambda')\
                    .atRedshift(redshift)

        return sed

    def magnorm(self, tophat_values, z_H):
        one_Jy = 1e-26  # W/Hz/m**2
        Lnu = tophat_values[self.ix_500nm]*self._to_W_per_Hz  # convert to W/Hz
        Fnu = Lnu/4/np.pi/self.dl(z_H)**2
        with np.errstate(divide='ignore', invalid='ignore'):
            return -2.5*np.log10(Fnu/one_Jy) + 8.90

class MilkyWayExtinction:
    '''
    Applies extinction to a SED
    '''
    def __init__(self, delta_wl=0.1, mwRv=3.1, eps=1e-7):
        """
        Parameters
        ----------
        sed_factory: (remove)
        ext_bin_width: (rename to delta_wl?)
        delta_wl : float [0.1]
            Wavelength sampling of the extinction function in nm
        mwRv : float [3.1]
            Parameter describing the shape of the Milky Way extinction
            curve.
        eps : float [1e-7]
            Small numerical offset to avoid out-of-range errors in
            the wavelength array passed to the dust_extinction code.
        """
        # Wavelength sampling for the extinction function. F19.x_range
        # is in units of 1/micron so convert to nm.  The eps value
        # is needed to avoid numerical noise at the end points causing
        # out of range errors detected by the dust_extinction code.
        wl_min = 1e3/F19.x_range[1] + eps
        wl_max = 1e3/F19.x_range[0] - eps
        npts = int((wl_max - wl_min)/delta_wl)
        self.wls = np.linspace(wl_min, wl_max, npts)
        self.extinction = F19(Rv=mwRv)

    def extinguish(self, sed, mwAv):
        ext = self.extinction.extinguish(self.wls*u.nm, Av=mwAv)
        lut = galsim.LookupTable(self.wls, ext, interpolant='linear')
        mw_ext = galsim.SED(lut, wave_type='nm', flux_type='1')
        sed = sed*mw_ext
        return sed


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
