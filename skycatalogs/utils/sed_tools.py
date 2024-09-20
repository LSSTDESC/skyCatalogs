import os
import re
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
import astropy.constants
import h5py

import numpy as np
from pathlib import PurePath
from dust_extinction.parameter_averages import F19
import galsim

__all__ = ['TophatSedFactory', 'DiffskySedFactory', 'SsoSedFactory',
           'MilkyWayExtinction', 'get_star_sed_path', 'generate_sed_path',
           'normalize_sed']

_FILE_PATH = str(PurePath(__file__))
_SKYCATALOGS_DIR = _FILE_PATH[:_FILE_PATH.rindex('/skycatalogs')]


class TophatSedFactory:
    '''
    Used for modeling cosmoDC2 galaxy SEDs, which are represented with
    a small number of wide bins
    '''
    _clight = astropy.constants.c.to('m/s').value
    # Conversion factor below of cosmoDC2 tophat Lnu values to W/Hz comes from
    # https://github.com/LSSTDESC/gcr-catalogs/blob/master/GCRCatalogs/SCHEMA.md
    _to_W_per_Hz = 4.4659e13

    # def __init__(self, th_definition, cosmology, delta_wl=0.001):
    def __init__(self, th_definition, cosmology, delta_wl=0.001, knots=True):
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
                   / (self.wl[1:] - self.wl[:-1]))

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
        lut = galsim.LookupTable(self.wl_deltas, flambda, interpolant='linear')

        if resolution:
            wl_min = min(self.wl_deltas)
            wl_max = max(self.wl_deltas)
            wl_res = np.linspace(wl_min, wl_max,
                                 int((wl_max - wl_min)/resolution))
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


class DiffskySedFactory:
    '''
    Used for collecting diffsky galaxy SEDs
    '''

    def __init__(self, catalog_dir, file_template, cosmology):

        self._files = {}
        self._catalog_dir = catalog_dir
        self._file_template = file_template

        # Create a FlatLambdaCDM cosmology from a dictionary of input
        # parameters.  This code is based on/borrowed from
        # https://github.com/LSSTDESC/gcr-catalogs/blob/master/GCRCatalogs/cosmodc2.py#L128
        cosmo_astropy_allowed = FlatLambdaCDM.__init__.__code__.co_varnames[1:]
        cosmo_astropy = {k: v for k, v in cosmology.items()
                         if k in cosmo_astropy_allowed}
        self.cosmology = FlatLambdaCDM(**cosmo_astropy)

    def _load_file(self, pixel):

        if pixel not in self._files:
            # sed_filename = self._file_template.format(pixel)
            sed_filename = f'galaxy_sed_{pixel}.hdf5'
            sed_path = os.path.join(self._catalog_dir, sed_filename)
            self._files[pixel] = h5py.File(sed_path, 'r')

        if not hasattr(self, '_wave_list'):
            self._wave_list = self._files[pixel]['meta/wave_list'][:]

        return self._files[pixel]

    @property
    def wave_list(self):
        return self._wave_list

    def dl(self, z):
        """
        Return the luminosity distance in units of meters.
        """
        return self.cosmology.luminosity_distance(z).value

    def create(self, pixel, galaxy_id, redshift_hubble, redshift):
        '''
        Get stored diffsky SED from disk.
        Does not apply extinction.
        '''

        f = self._load_file(pixel)

        sed_array = f['galaxy/'+str(int(galaxy_id)//100000)+'/'+str(galaxy_id)][:].astype('float')
        sed_array /= (4.0*np.pi*(self.dl(redshift_hubble))**2)

        seds = {}
        for i, component in enumerate(['bulge', 'disk', 'knots']):
            lut = galsim.LookupTable(x=self._wave_list, f=sed_array[i, :],
                                     interpolant='linear')
            # Create the SED object and apply redshift.
            sed = galsim.SED(lut, wave_type='angstrom', flux_type='fnu')\
                        .atRedshift(redshift)
            seds[component] = sed

        return seds


class SsoSedFactory():
    '''
    Load the single SED used for SSO objects and make it available as galsim
    SED
    '''
    DEFAULT_SED_BNAME = 'solar_sed_thin.txt'
    def __init__(self, sed_path=None):
        '''
        Format of sed file is two-column text file, which galsim can
        read directly. Columns are
        "wavelength" (units angstroms) and "flux" (units flambda)
        '''
        if not sed_path:
            # Get directory for possible default sed files
            sed_path = os.path.join(_SKYCATALOGS_DIR, 'skycatalogs',
                                    'data', 'sso',
                                    SsoSedFactory.DEFAULT_SED_BNAME)
        wave_type = 'angstrom'
        flux_type = 'flambda'
        lut = galsim.LookupTable.from_file(sed_path, interpolant='linear')

        sed = galsim.SED(lut, wave_type=wave_type, flux_type=flux_type)
        self.sed = sed
        self._sed_path = sed_path    # In case we want to save it for posterity

    @property
    def sed_path(self):
        return self._sed_path

    def create(self):
        return self.sed


class MilkyWayExtinction:
    '''
    Applies extinction to a SED
    '''
    def __init__(self, delta_wl=1.0, mwRv=3.1, eps=1e-7):
        """
        Parameters
        ----------
        delta_wl : float [1.0]
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
        mw_ext = galsim.SED(lut, wave_type='nm', flux_type='1').thin()
        sed = sed*mw_ext
        return sed


_standard_dict = {'lte': 'starSED/phoSimMLT',
                  'bergeron': 'starSED/wDs',
                  'km|kp': 'starSED/kurucz'}


def get_star_sed_path(filename, name_to_folder=_standard_dict):
    '''
    Return numpy array of full paths relative to SIMS_SED_LIBRARY_DIR,
    given filenames

    Parameters
    ----------
    filename       list of strings. Usually full filename but may be
                   missing final ".gz"
    name_to_folder dict mapping regular expression (to be matched with
                   filename) to relative path for containing directory

    Returns
    -------
    Full path for file, relative to SIMS_SED_LIBRARY_DIR
    '''

    compiled = {re.compile(k): v for (k, v) in name_to_folder.items()}

    path_list = []
    for f in filename:
        m = None
        matched = False
        for k, v in compiled.items():
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


def generate_sed_path(ids, subdir, cmp):
    '''
    Generate paths (e.g. relative to SIMS_SED_LIBRARY_DIR) for galaxy component
    SED files
    Parameters
    ----------
    ids        list of galaxy ids
    subdir    user-supplied part of path
    cmp      component for which paths should be generated

    returns
    -------
    A list of strings.  The entries in the list have the form
    <subdir>/<cmp>_<id>.txt
    '''
    r = [f'{subdir}/{cmp}_{id}.txt' for id in ids]
    return r


def normalize_sed(sed, magnorm, wl=500*u.nm):
    """
    Set the normalization of a GalSim SED object given a monochromatic
    magnitude at a reference wavelength.

    Parameters
    ----------
    sed : galsim.SED
        The GalSim SED object.
    magnorm : float
        The monochromatic magnitude at the reference wavelength.
    wl : astropy.units.nm
        The reference wavelength.

    Returns
    -------
    galsim.SED : The renormalized SED object.
    """
    # Compute the flux density from magnorm in units of erg/cm^2/s/nm.
    fnu = (magnorm * u.ABmag).to_value(u.erg/u.s/u.cm**2/u.Hz)
    flambda = fnu * (astropy.constants.c/wl**2).to_value(u.Hz/u.nm)

    # GalSim expects the flux density in units of photons/cm^2/s/nm,
    # so divide flambda by the photon energy at the reference
    # wavelength.
    hnu = (astropy.constants.h * astropy.constants.c / wl).to_value(u.erg)

    flux_density = flambda/hnu

    return sed.withFluxDensity(flux_density, wl)
