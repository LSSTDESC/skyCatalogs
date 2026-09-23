import os
import re
import time
from collections import OrderedDict
from astropy import units as u
from astropy.coordinates import Distance
from astropy.cosmology import FlatLambdaCDM
import astropy.constants
import pandas as pd

import numpy as np
from pathlib import PurePath
from dust_extinction.parameter_averages import F19
import galsim

__all__ = ['TophatSedFactory', 'DiffskySedFactory', 'SsoSedFactory',
           'MilkyWayExtinction', 'TrilegalSedFactory', 'get_star_sed_path',
           'generate_sed_path', 'normalize_sed']

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
    """Compute current Diffsky component SEDs from packaged runtime state."""

    _flux_factor = 4.0204145742268754e-16
    DEFAULT_EMISSION_LINES_ANGSTROM = (
        3726.03, 3728.82, 4101.73, 4340.47, 4861.33,
        4958.91, 5006.84, 6548.05, 6562.80, 6583.45,
        6716.44, 6730.82,
    )

    def __init__(self, runtime_state_dir, cosmology, object_batch_size=256,
                 diffsky_batch_size=25, cache_size=8, rel_err=0.03,
                 wave_ang_min=500, wave_ang_max=100000,
                 pixel_cache_size=4, ssp_wave_min_micron=0.06,
                 ssp_wave_max_micron=2.34, sed_engine='fast',
                 sed_precision=None, thinning_mode='galsim',
                 thinning_selector=None, emission_lines_angstrom=None,
                 emission_line_half_width_angstrom=20.0):
        if object_batch_size < 1 or cache_size < 1 or pixel_cache_size < 1:
            raise ValueError(
                'Diffsky object and cache sizes must be positive integers')
        self._runtime_state_dir = os.path.abspath(runtime_state_dir)
        self._object_batch_size = object_batch_size
        self._diffsky_batch_size = diffsky_batch_size
        self._cache_size = cache_size
        self._sed_cache = OrderedDict()
        self._sed_engine = None
        self._sed_precision = None
        self.set_compute_options(sed_engine, sed_precision)
        self._sed_cache_capacity = cache_size * object_batch_size
        self._pixel_cache = OrderedDict()
        self._pixel_cache_size = pixel_cache_size
        self._aux_data = None
        self._profile_enabled = False
        self.reset_profile()
        self._thinning_options = (rel_err, wave_ang_min, wave_ang_max)
        self._thinning_mode = thinning_mode
        self._thinning_selector = thinning_selector
        self._emission_lines_angstrom = tuple(
            self.DEFAULT_EMISSION_LINES_ANGSTROM
            if emission_lines_angstrom is None else emission_lines_angstrom)
        self._emission_line_half_width_angstrom = \
            emission_line_half_width_angstrom
        self._validate_thinning_options()
        self._ssp_wave_bounds_micron = (
            ssp_wave_min_micron, ssp_wave_max_micron)
        self._validate_ssp_wave_bounds()

        # Create a FlatLambdaCDM cosmology from a dictionary of input
        # parameters.  This code is based on/borrowed from
        # https://github.com/LSSTDESC/gcr-catalogs/blob/master/GCRCatalogs/cosmodc2.py#L128
        cosmo_astropy_allowed = FlatLambdaCDM.__init__.__code__.co_varnames[1:]
        cosmo_astropy = {k: v for k, v in cosmology.items()
                         if k in cosmo_astropy_allowed}
        self.cosmology = FlatLambdaCDM(**cosmo_astropy)

    @property
    def sed_engine(self):
        """Name of the active component-SED calculation implementation."""
        return self._sed_engine

    @property
    def sed_precision(self):
        """Arithmetic precision used by the optimized SED implementation."""
        return self._sed_precision

    def set_compute_options(self, engine=None, precision=None):
        """Select the SED implementation before any spectra are cached."""
        engine = self._sed_engine if engine is None else engine
        if precision is None:
            if engine == 'reference':
                precision = 'float64'
            else:
                precision = self._sed_precision or 'float32'
        if engine not in ('fast', 'reference'):
            raise ValueError("SED engine must be 'fast' or 'reference'")
        if precision not in ('float32', 'float64'):
            raise ValueError(
                "SED precision must be 'float32' or 'float64'")
        if engine == 'reference' and precision != 'float64':
            raise ValueError(
                "The reference SED engine only supports float64 precision")
        if self._sed_cache and (
                engine != self._sed_engine
                or precision != self._sed_precision):
            raise RuntimeError(
                'SED compute options cannot change after spectra are cached')
        self._sed_engine = engine
        self._sed_precision = precision

    def reset_profile(self):
        """Reset optional fine-grained Diffsky runtime counters."""
        self._profile = {name: 0.0 for name in (
            'runtime_pixel', 'id_lookup', 'opencosmo_take_rows',
            'ssp_prepare', 'diffsky_compute', 'cache_insert',
            'galsim_sed_build', 'extinction')}
        self._profile.update({
            'requested': 0, 'cache_hits': 0, 'cache_misses': 0,
            'evictions': 0, 'diffsky_calls': 0})

    def profile(self, name, elapsed):
        if self._profile_enabled:
            self._profile[name] += elapsed

    @property
    def profile_stats(self):
        return dict(self._profile)

    def clear_sed_cache(self):
        """Release cached per-galaxy component arrays.

        This does not clear JAX executable caches or close packaged runtime
        pixels. It is intended for when you have completed a catalog
        output chunk and no longer benefit from retaining its spectra.
        """
        self._sed_cache.clear()

    def clear_runtime_caches(self, clear_compiled=False):
        """Release SED arrays, open runtime pixels, and optionally JAX code."""
        self.clear_sed_cache()
        self._pixel_cache.clear()
        if clear_compiled and self._sed_engine == 'fast':
            from skycatalogs.utils.diffsky_sed import clear_fast_sed_caches
            clear_fast_sed_caches()

    @property
    def thinning_mode(self):
        """Wavelength-selection policy used for every cached Diffsky SED."""
        return self._thinning_mode

    def _validate_thinning_options(self):
        valid_modes = ('galsim', 'none', 'emission_lines', 'custom')
        if self._thinning_mode not in valid_modes:
            raise ValueError(
                f'Diffsky thinning mode must be one of {valid_modes}')
        if self._thinning_mode == 'custom' and not callable(
                self._thinning_selector):
            raise ValueError(
                'Custom Diffsky thinning requires a callable selector')
        if (self._emission_line_half_width_angstrom is None
                or not np.isfinite(
                    self._emission_line_half_width_angstrom)
                or self._emission_line_half_width_angstrom < 0):
            raise ValueError(
                'Emission-line half-width must be non-negative')
        if not np.all(np.isfinite(self._emission_lines_angstrom)):
            raise ValueError('Emission-line wavelengths must be finite')

    def _validate_ssp_wave_bounds(self):
        lower, upper = self._ssp_wave_bounds_micron
        if (lower is None) != (upper is None):
            raise ValueError(
                'Both SSP wavelength bounds must be provided together')
        if lower is not None and (lower <= 0 or upper <= lower):
            raise ValueError(
                'SSP wavelength bounds must satisfy 0 < min < max')

    def set_ssp_wave_bounds(self, wave_min_micron, wave_max_micron):
        """Set the SSP wavelength slice before runtime data are loaded."""
        if self._aux_data is not None:
            raise RuntimeError(
                'SSP wavelength bounds cannot change after SED data load')
        self._ssp_wave_bounds_micron = (
            wave_min_micron, wave_max_micron)
        self._validate_ssp_wave_bounds()

    @staticmethod
    def _slice_ssp_data(ssp_data, lower_micron, upper_micron):
        wave = np.asarray(ssp_data.ssp_wave)
        lower_angstrom = lower_micron * 1.0e4
        upper_angstrom = upper_micron * 1.0e4
        if upper_angstrom < wave[0] or lower_angstrom > wave[-1]:
            raise ValueError(
                f'Wavelength bounds {lower_micron}--{upper_micron} micron '
                'do not overlap the available SSP grid')
        # Retain one grid sample outside either requested edge so GalSim's
        # linear interpolation covers the entire nonzero bandpass interval.
        start = max(0, np.searchsorted(wave, lower_angstrom) - 1)
        stop = min(len(wave), np.searchsorted(
            wave, upper_angstrom, side='right') + 1)
        if stop <= start:
            raise ValueError(
                f'Wavelength bounds {lower_micron}--{upper_micron} micron '
                'do not overlap the available SSP grid')
        return ssp_data._replace(
            ssp_wave=wave[start:stop],
            ssp_flux=np.asarray(ssp_data.ssp_flux)[..., start:stop])

    def _slice_ssp_wavelengths(self):
        lower, upper = self._ssp_wave_bounds_micron
        if lower is None:
            return
        self._aux_data['ssp_data'] = self._slice_ssp_data(
            self._aux_data['ssp_data'], lower, upper)

    def _ensure_aux_loaded(self, state_file):
        if self._aux_data is not None:
            return
        import h5py
        try:
            from diffsky.data_loaders.hacc_utils.lc_mock import (
                load_diffsky_param_collection_merging,
                load_diffsky_ssp_data,
                load_diffsky_t_table,
                load_diffsky_tcurves,
            )
        except ModuleNotFoundError as exc:
            raise ModuleNotFoundError(
                'Diffsky galaxy SED support requires the optional '
                'dependencies installed by `pip install skyCatalogs[diffsky]`'
            ) from exc

        with h5py.File(state_file) as handle:
            mock_name = handle['header']['catalog_info'].attrs[
                'mock_version_name']
        root = self._runtime_state_dir
        self._aux_data = {
            'ssp_data': load_diffsky_ssp_data(root, mock_name),
            'param_collection': load_diffsky_param_collection_merging(
                root, mock_name),
            'tcurves': load_diffsky_tcurves(root, mock_name),
            't_table': load_diffsky_t_table(root, mock_name),
        }
        self._slice_ssp_wavelengths()
        self._set_wavelengths()

    @staticmethod
    def _thinned_wavelengths(ssp_data, rel_err, wave_min, wave_max):
        wave = np.asarray(ssp_data.ssp_wave)
        use = (wave > wave_min) & (wave < wave_max)
        ssp_flux = np.asarray(ssp_data.ssp_flux)
        representative = np.sum(
            ssp_flux, axis=tuple(range(ssp_flux.ndim - 1)))
        wave_nm = wave[use] / 10.0
        sed = galsim.SED(
            galsim.LookupTable(wave_nm, representative[use]),
            wave_type='nm', flux_type='flambda')
        thinned = sed.thin(rel_err=rel_err, fast_search=False)
        candidates = np.flatnonzero(use)
        indices = candidates[np.isin(wave_nm, thinned.wave_list)]
        return indices, wave[indices]

    @staticmethod
    def _validate_wave_indices(indices, size):
        indices = np.asarray(indices)
        if indices.dtype == bool:
            if indices.ndim != 1 or len(indices) != size:
                raise ValueError(
                    'Custom wavelength mask must match the SSP grid')
            indices = np.flatnonzero(indices)
        elif not np.issubdtype(indices.dtype, np.integer):
            raise TypeError(
                'Custom wavelength selector must return indices or a mask')
        indices = np.unique(indices.astype(np.int64, copy=False))
        if len(indices) < 2:
            raise ValueError('At least two wavelength samples are required')
        if indices[0] < 0 or indices[-1] >= size:
            raise ValueError('Custom wavelength indices are out of range')
        return indices

    def _emission_line_indices(self, wave):
        """Return native samples protected around configured rest-frame lines."""
        protected = np.zeros(len(wave), dtype=bool)
        half_width = self._emission_line_half_width_angstrom
        for line in self._emission_lines_angstrom:
            if line < wave[0] - half_width or line > wave[-1] + half_width:
                continue
            inside = np.flatnonzero(np.abs(wave - line) <= half_width)
            if len(inside) == 0:
                inside = np.array([np.abs(wave - line).argmin()])
            lo = max(0, inside[0] - 1)
            hi = min(len(wave), inside[-1] + 2)
            protected[lo:hi] = True
        return np.flatnonzero(protected)

    def _set_wavelengths(self):
        """Select one immutable wavelength grid for this factory."""
        ssp_data = self._aux_data['ssp_data']
        wave = np.asarray(ssp_data.ssp_wave)
        if self._thinning_mode == 'none':
            indices = np.arange(len(wave), dtype=np.int64)
        elif self._thinning_mode == 'custom':
            indices = self._thinning_selector(wave, np.asarray(
                ssp_data.ssp_flux))
        else:
            indices, _ = self._thinned_wavelengths(
                ssp_data, *self._thinning_options)
            if self._thinning_mode == 'emission_lines':
                indices = np.union1d(
                    indices, self._emission_line_indices(wave))
        self._wave_indices = self._validate_wave_indices(indices, len(wave))
        self._wave_list = wave[self._wave_indices]

    @property
    def prefetch_batch_size(self):
        """Maximum requested-object batch recommended to callers."""
        return self._object_batch_size

    def _get_runtime_pixel(self, pixel):
        if pixel not in self._pixel_cache:
            from pathlib import Path
            try:
                import opencosmo as oc
            except ModuleNotFoundError as exc:
                raise ModuleNotFoundError(
                    'Diffsky galaxy SED support requires the optional '
                    'dependencies installed by '
                    '`pip install skyCatalogs[diffsky]`') from exc

            pixel_dir = Path(self._runtime_state_dir) / f'pixel_{pixel}'
            state_files = sorted(pixel_dir.glob('*.diffsky_gals.hdf5'))
            if not state_files:
                raise FileNotFoundError(
                    f'No packaged Diffsky state found for pixel {pixel} in '
                    f'{pixel_dir}')
            self._ensure_aux_loaded(state_files[0])
            # keep_top_host is required by OpenCosmo >=1.3 for the linked
            # host rows used in Diffsky SED calculations.
            catalog = oc.open(state_files, keep_top_host=True)
            galaxy_ids = np.atleast_1d(
                catalog.select('gal_id').get_data('numpy'))
            id_to_index = {int(galaxy_id): index
                           for index, galaxy_id in enumerate(galaxy_ids)}
            self._pixel_cache[pixel] = (catalog, id_to_index)
            while len(self._pixel_cache) > self._pixel_cache_size:
                self._pixel_cache.popitem(last=False)
        else:
            self._pixel_cache.move_to_end(pixel)
        return self._pixel_cache[pixel]

    def _cache_sed_info(self, batch_catalog, sed_info, wave_indices=None,
                        wave_list=None, requested_ids=None):
        """Thin and cache computed rest-frame component SED arrays."""
        batch_ids = np.atleast_1d(
            batch_catalog.select('gal_id').get_data('numpy'))
        if requested_ids is None:
            requested_ids = batch_ids
        else:
            requested_ids = np.atleast_1d(requested_ids)
        if wave_indices is None:
            wave_indices = self._wave_indices
        if wave_list is None:
            wave_list = self._wave_list
        storage_dtype = (np.float32 if self._sed_engine == 'fast'
                         and self._sed_precision == 'float32' else None)
        component_arrays = tuple(
            np.asarray(sed_info[name], dtype=storage_dtype)[:, wave_indices]
            for name in ('rest_sed_bulge', 'rest_sed_disk',
                         'rest_sed_knots')
        )
        batch_rows = {int(batch_id): row
                      for row, batch_id in enumerate(batch_ids)}
        for requested_id in requested_ids:
            galaxy_id = int(requested_id)
            row = batch_rows[galaxy_id]
            self._sed_cache[galaxy_id] = (
                np.asarray(wave_list),
                tuple(component[row] for component in component_arrays))
            self._sed_cache.move_to_end(galaxy_id)
        while len(self._sed_cache) > self._sed_cache_capacity:
            self._sed_cache.popitem(last=False)
            if self._profile_enabled:
                self._profile['evictions'] += 1

    def prefetch(self, galaxy_ids, partition_ids, redshifts=None):
        """Compute and cache SEDs for the explicitly requested galaxies.

        Requests are grouped by their SkyCatalog output pixel and selected
        with ``take_rows`` from the packaged native-state sidecar.
        """
        try:
            if self._sed_engine == 'fast':
                from skycatalogs.utils.diffsky_sed import (
                    compute_dbk_seds_fast as compute_seds,
                )
            else:
                from diffsky.data_loaders.opencosmo_utils import (
                    compute_dbk_seds_from_diffsky_mock as compute_seds,
                )
        except ModuleNotFoundError as exc:
            raise ModuleNotFoundError(
                'Diffsky galaxy SED support requires the optional '
                'dependencies installed by `pip install skyCatalogs[diffsky]`'
            ) from exc

        galaxy_ids = np.atleast_1d(galaxy_ids).astype(np.int64)
        partition_ids = np.atleast_1d(partition_ids).astype(np.int64)
        if len(galaxy_ids) != len(partition_ids):
            raise ValueError('galaxy_ids and partition_ids must align')
        if redshifts is not None:
            redshifts = np.atleast_1d(redshifts).astype(float)
            if len(galaxy_ids) != len(redshifts):
                raise ValueError('galaxy_ids and redshifts must align')

        missing = np.array(
            [int(galaxy_id) not in self._sed_cache
             for galaxy_id in galaxy_ids], dtype=bool)
        if self._profile_enabled:
            self._profile['requested'] += len(galaxy_ids)
            self._profile['cache_hits'] += int(np.count_nonzero(~missing))
            self._profile['cache_misses'] += int(np.count_nonzero(missing))
        if not np.any(missing):
            for galaxy_id in galaxy_ids:
                self._sed_cache.move_to_end(int(galaxy_id))
            return

        galaxy_ids = galaxy_ids[missing]
        partition_ids = partition_ids[missing]
        if redshifts is not None:
            redshifts = redshifts[missing]
            # Keep fixed-size Diffsky calls local in redshift so a
            # bandpass-limited SSP grid needs the narrowest interval.
            order = np.argsort(redshifts, kind='stable')
            galaxy_ids = galaxy_ids[order]
            partition_ids = partition_ids[order]
            redshifts = redshifts[order]
        for pixel in np.unique(partition_ids):
            started = time.perf_counter()
            pixel_catalog, id_to_index = self._get_runtime_pixel(int(pixel))
            self.profile('runtime_pixel', time.perf_counter() - started)
            pixel_ids = galaxy_ids[partition_ids == pixel]
            started = time.perf_counter()
            try:
                rows = np.array(
                    [id_to_index[int(galaxy_id)] for galaxy_id in pixel_ids],
                    dtype=np.int64)
            except KeyError as exc:
                raise KeyError(
                    f'Galaxy {exc.args[0]} not found in packaged Diffsky '
                    f'state for pixel {pixel}') from exc
            self.profile('id_lookup', time.perf_counter() - started)

            group_indices = np.arange(len(rows))
            for start in range(0, len(group_indices),
                                   self._object_batch_size):
                    batch_indices = group_indices[
                        start:start + self._object_batch_size]
                    requested_ids = pixel_ids[batch_indices]
                    # OpenCosmo requires monotonically increasing row indices.
                    batch_rows = np.sort(rows[batch_indices])
                    started = time.perf_counter()
                    batch_catalog = pixel_catalog.take_rows(batch_rows)
                    self.profile(
                        'opencosmo_take_rows', time.perf_counter() - started)
                    started = time.perf_counter()
                    if self._sed_engine == 'fast':
                        sed_info = compute_seds(
                            batch_catalog, self._aux_data,
                            batch_size=self._diffsky_batch_size,
                            precision=self._sed_precision)
                    else:
                        sed_info = compute_seds(
                            batch_catalog, self._aux_data, insert=False,
                            batch_size=self._diffsky_batch_size)
                    self.profile(
                        'diffsky_compute', time.perf_counter() - started)
                    if self._profile_enabled:
                        self._profile['diffsky_calls'] += 1
                    started = time.perf_counter()
                    # OpenCosmo may include linked top hosts needed by the SED
                    # calculation. Cache only the explicitly requested rows so
                    # helper hosts cannot evict objects the caller will render.
                    self._cache_sed_info(
                        batch_catalog, sed_info,
                        requested_ids=requested_ids)
                    self.profile('cache_insert', time.perf_counter() - started)

    @property
    def wave_list(self):
        if self._aux_data is None:
            raise RuntimeError(
                'Diffsky wavelengths are unavailable until a pixel is loaded')
        return self._wave_list

    def dl(self, z):
        """
        Return the luminosity distance in Mpc.
        """
        return self.cosmology.luminosity_distance(z).value

    def create(self, galaxy_id, partition_id, redshift_hubble, redshift):
        """Return lazily computed, unextincted component SEDs."""
        galaxy_id = int(galaxy_id)
        if galaxy_id not in self._sed_cache:
            self.prefetch(
                [galaxy_id], [partition_id], [redshift_hubble])
        try:
            wave_list, component_arrays = self._sed_cache[galaxy_id]
            sed_array = np.asarray(component_arrays)
        except KeyError as exc:
            raise KeyError(
                f'Galaxy {galaxy_id} missing after Diffsky SED computation') \
                    from exc
        self._sed_cache.move_to_end(galaxy_id)
        sed_array *= self._flux_factor
        sed_array /= (4.0*np.pi*(self.dl(redshift_hubble))**2)

        seds = {}
        started = time.perf_counter()
        for i, component in enumerate(['bulge', 'disk', 'knots']):
            lut = galsim.LookupTable(x=wave_list, f=sed_array[i, :],
                                     interpolant='linear')
            # Create the SED object and apply redshift.
            sed = galsim.SED(lut, wave_type='angstrom', flux_type='fnu')\
                        .atRedshift(redshift)
            seds[component] = sed
        self.profile('galsim_sed_build', time.perf_counter() - started)

        return seds


class TrilegalSedFactory():

    def __init__(self, object_type_config, logger):
        '''
        Parameters
        ----------
        object_type_config  dict containing section of the sky catalog
                            config pertaining to object type 'trilegal'
        '''
        self._pystellib = None
        self._logger = logger

        self._errors = 0

    @property
    def error_count(self):
        return self._errors

    def pystellib(self):
        if not self._pystellib:
            from pystellibs import BTSettl
            self._pystellib = BTSettl(medres=False)
        return self._pystellib

    def clear_errors(self):
        self._errors = 0

    def get_sed(self, tri):
        '''
        Parameters
        ----------
        tri  a trilegal object

        Returns
        -------
        galsim sed object.

        If SED file present, find row and wave length values in sed file
        Otherwise use quantities in main parquet file + pystellibs
        In either case return unextincted, unnormalized, unthinned SED,
        just as calculated by pystellib. The other transformations are
        handled elsewhere
        '''
        # if not self._pystellib:
        #     self._finish_init()

        # Get inputs from parquet file
        native = ['logT', 'logg', 'logL', 'Z']
        spec_inputs = [tri.get_native_attribute(x) for x in native]

        # May generate a runtime error if parameters are outside
        # interpolation range
        error_msg = None
        try:
            spectrum = self.pystellib().generate_stellar_spectrum(*spec_inputs)
        except RuntimeError:
            error_msg = 'Run-time error generating SED'

        if not error_msg:
            if spectrum is None:
                error_msg = 'No SED computed'
            elif any(np.isnan(spectrum)):
                error_msg = 'NaNs in SED'

        if error_msg:
            # Maybe keep track separately depending on value of evol_label?
            # evol_label = tri.get_native_attribute('evol_label')
            if not self._errors:
                self._logger.warning(error_msg)
            self._errors += 1
            return None

        # BTSettl.generate_stellar_spectrum produces spectra with
        # units erg/Angstrom/s, whereas galsim flambda units are
        # erg/Angstrom/cm^2/s, with wave_type='Angstrom'.  To convert
        # to flambda, we need to apply the inverse-square law.
        # We do this in log space to avoid underflows.

        # Compute distance to star in cm from the distance modulus.
        mu0 = tri.get_native_attribute("mu0")
        dist = Distance(distmod=mu0).to(u.cm).value
        # Apply 1/(4*pi*dist**2) dilution factor.
        log_dilution = np.log(4.0*np.pi) + 2.0*np.log(dist)
        index = np.where(spectrum > 0)  # avoid zeros
        spectrum[index] = np.exp(np.log(spectrum[index]) - log_dilution)
        sed_table = galsim.LookupTable(self.pystellib().wavelength, spectrum,
                                       interpolant='linear')
        sed = galsim.SED(sed_table, 'Angstrom', 'flambda')

        return sed

    def get_spectra_batch(self, pq_main, batch, l_bnd, u_bnd):
        '''
        Return spectra (still will need to be converted to observer SED)
        as computed by pystellibs for a subset (slice of a row group
        in the case of parquet input which for now is the only type
        supported).

        Parameters
        ----------
        pq_main     ParquetFile object for the "main" catalog for
                    healpixel of interest
        batch       row group for which spectra are to be returned
        l_bnd       Delimits slice
        u_bnd       Delimits slics

        Returns
        -------
        wavelength axis  numpy array of dimension n_wl
        array of galsim SED (not extincted)

        '''
        columns = ['id', 'logT', 'logg', 'logL', 'Z', 'mu0']
        a_dict = pq_main.read_row_group(batch, columns=columns).to_pydict()
        for k in a_dict.keys():
            a_dict[k] = a_dict[k][l_bnd: u_bnd]
        df = pd.DataFrame(a_dict)
        wl_axis, spectra = self.pystellib().generate_individual_spectra(df)
        spectra = np.array(spectra)

        dist = Distance(distmod=df['mu0'].to_numpy()).to(u.cm).value
        log_dilution = np.log(4.0*np.pi) + 2.0*np.log(dist)
        log_dilution = np.stack([log_dilution]*spectra.shape[1], axis=1)

        index = np.asarray(spectra > 0).nonzero()
        spectra[index] = np.exp(np.log(spectra[index]) - log_dilution[index])
        spectra_32 = spectra.astype(np.float32)
        seds = [galsim.SED(galsim.LookupTable(
            wl_axis, spectrum,
            interpolant='linear'), 'Angstrom', 'flambda')
                if not any(np.isnan(spectrum)) else None  for spectrum in spectra]

        del df
        del spectra
        del spectra_32

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
