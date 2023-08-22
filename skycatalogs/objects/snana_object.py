import galsim
import h5py
from .base_object import BaseObject,ObjectCollection

__all__ = ['SnanaObject']

class SnanaObject(BaseObject):
    _type_name = 'snana'

    def _get_sed(self, mjd=None):
        if mjd is None:
            return None, 0.0
        model = self.get_native_attribute('model_name')
        param_names = self.get_native_attribute('model_param_names')
        param_values = self.get_native_attribute('model_param_values')
        mjd_start = self.get_native_attribute('mjd_start')
        mjd_end = self.get_native_attribute('mjd_start')
        if mjd < mjd_start or mjd > mjd_end:
            return None, 0.0
        # For now find the SED with mjd closest to ours and return that

        # What do we do about magnorm?


    def get_gsobject_components(self, gsparams=None, rng=None):
        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)
        return {'this_object': galsim.DeltaFunction(gsparams=gsparams)}

    def get_observer_sed_component(self, component, mjd=None):
        sed, _ = self._get_sed(mjd=mjd)
        if sed is not None:
            sed = self._apply_component_extinction(sed)
        return sed

    def get_LSST_flux(self, band, sed=None, mjd=None):
        if not band in LSST_BANDS:
            return None

        return self.get_flux(lsst_bandpasses[band], sed=sed, mjd=mjd)

    def _read_nearest_SED(self, mjd):
        pass   # until I figure it out

class SnanaCollection(ObjectCollection):
    '''
    This class (so far) differs from the vanilla ObjectCollection only
    in that it keeps track of where the file is which contains a library
    of SEDs for each sn
    '''
    def __init__(self, ra, dec, id, object_type, partition_id, sky_catalog,
                 region=None, mjd=None, mask=None, readers=None, row_group=0,
                 SED_file=None):
        self._SED_file = SED_file
        super().__init__(ra, dec, id, object_type, partition_id, sky_catalog,
                         region, mjd, mask, readers, row_group)
