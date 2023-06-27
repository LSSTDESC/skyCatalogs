import os
import galsim
from .base_object import BaseObject

__all__ = ['StarObject']

class StarObject(BaseObject):
    _type_name = 'star'
    def _get_sed(self, mjd=None):
        '''
        We'll need mjd when/if variable stars are supported. For now
        it's ignored.
        '''
        mag_norm = self.get_native_attribute('magnorm')
        rel_path = self.get_native_attribute('sed_filepath')

        fpath = os.path.join(os.getenv('SIMS_SED_LIBRARY_DIR'),
                             self.get_native_attribute('sed_filepath'))
        sky_cat = self._belongs_to._sky_catalog

        return sky_cat.observed_sed_factory.create_pointsource(rel_path),mag_norm
    def get_gsobject_components(self, gsparams=None, rng=None):
        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)
        return {'this_object': galsim.DeltaFunction(gsparams=gsparams)}

    def get_observer_sed_component(self, component, mjd=None):
        sed, _ = self._get_sed(mjd=mjd)
        if sed is not None:
            sed = self._apply_component_extinction(sed)
        return sed
