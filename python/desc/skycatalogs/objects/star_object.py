import os
import numpy as np
import galsim
from .base_object import BaseObject

__all__ = ['StarObject']

class StarObject(BaseObject):
    _type_name = 'star'
    def _get_sed(self, mjd=None, redshift=0):
        '''
        We'll need mjd when/if variable stars are supported. For now
        it's ignored.
        '''
        mag_norm = self.get_native_attribute('magnorm')
        rel_path = self.get_native_attribute('sed_filepath')

        fpath = os.path.join(os.getenv('SIMS_SED_LIBRARY_DIR'),
                             self.get_native_attribute('sed_filepath'))

        return self._get_sed_from_file(fpath), mag_norm

    def get_gsobject_components(self, gsparams=None, rng=None):
        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)
        return {'this_object': galsim.DeltaFunction(gsparams=gsparams)}

    def get_observer_sed_component(self, component, mjd=None):
        sed, magnorm = self._get_sed(mjd=mjd)

        # The SED is in units of photons/nm/cm^2/s
        # -0.9210340371976184 = -np.log(10)/2.5. Use to convert mag to flux
        flux_500 = np.exp(-0.9210340371976184 * magnorm)
        sed = sed.withMagnitude(0, self._bp500)
        sed = sed*flux_500

        if sed is not None:
            sed = self._apply_component_extinction(sed)
        return sed
