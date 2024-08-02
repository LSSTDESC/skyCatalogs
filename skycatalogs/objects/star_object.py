import os
import galsim
from .base_object import BaseObject
from ..utils import normalize_sed
from .base_config_fragment import BaseConfigFragment

__all__ = ['StarObject', 'StarConfigFragment']


class StarObject(BaseObject):

    _type_name = 'star'

    def _get_sed(self, mjd=None, redshift=0):
        '''
        We'll need mjd when/if variable stars are supported. For now
        it's ignored.
        '''
        mag_norm = self.get_native_attribute('magnorm')
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


class StarConfigFragment(BaseConfigFragment):
    def __init__(self, prov, area_partition=None, data_file_type=None):
        super().__init__(prov, object_type_name='star',
                         area_partition=area_partition,
                         data_file_type=data_file_type)
