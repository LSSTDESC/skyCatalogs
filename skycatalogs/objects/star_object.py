import os
import numpy as np
import galsim
from .base_object import BaseObject, load_lsst_bandpasses
from ..utils import normalize_sed


__all__ = ['StarObject']


class StarObject(BaseObject):

    _type_name = 'star'
    _mjd0 = 60400.0  # reference epoch for sinusoidal variability
    _lsst_bandpasses = load_lsst_bandpasses()

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not 'is_variable' in self.native_columns:
            self.is_variable = False
        else:
            self.is_variable = self.get_native_attribute('is_variable')

        if self.is_variable:
            # Get parameters governing sinusoidal variability.
            # The period has units of days
            self.period = self.get_native_attribute('period')
            # The phase has units of radians at self._mjd0
            self.phase = self.get_native_attribute('phase')
            # This is the amplitude of the sinusoidal modulation applied
            # to magnorm.
            self.mag_amplitude = self.get_native_attribute('mag_amplitude')

    def _mag_modulation(self, mjd):
        if not self.is_variable or mjd is None:
            return 0.0
        modulation = (self.mag_amplitude *
                      np.sin(2.*np.pi*(mjd - self._mjd0)/self.period + self.phase))
        return modulation

    def _get_sed(self, mjd=None, redshift=0):
        """Return the SED"""
        mag_norm = self.get_native_attribute('magnorm') + self._mag_modulation(mjd)
        fpath = os.path.join(os.getenv('SIMS_SED_LIBRARY_DIR'),
                             self.get_native_attribute('sed_filepath'))

        return normalize_sed(self._get_sed_from_file(fpath), mag_norm)

    def get_gsobject_components(self, gsparams=None, rng=None):
        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)
        return {'this_object': galsim.DeltaFunction(gsparams=gsparams)}

    def get_observer_sed_component(self, component, mjd=None):
        sed = self._get_sed(mjd=mjd)

        if sed is not None:
            sed = self._apply_component_extinction(sed)

        return sed

    def get_LSST_flux(self, band, sed=None, cache=True, mjd=None):
        if not self.is_variable:
            return super().get_LSST_flux(band, sed=sed, cache=cache, mjd=mjd)

        # Use the direct call to .get_flux to avoid unwanted flux caching
        # for variable stars.
        return self.get_flux(self._lsst_bandpasses[band], sed=sed, mjd=mjd)
