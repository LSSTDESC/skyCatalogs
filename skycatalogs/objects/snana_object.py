import galsim
import h5py
import numpy as np
from .base_object import BaseObject,ObjectCollection

__all__ = ['SnanaObject', 'SnanaCollection']

class SnanaObject(BaseObject):
    _type_name = 'snana'

    def __init__(self, ra, dec, id, object_type, belongs_to, belongs_index,
                 redshift=None):
        super().__init__(ra, dec, id, object_type, belongs_to, belongs_index,
                         redshift)
        self._time_sampling = None
        self._lambda = None

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
        # What should we do about magnorm?
        return _read_nearest_SED(mjd)


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
        # Find row with closest mjd and return it along with
        # (for now) magnorm of 0.0

        f = h5py.File(self._belongs_to._SED_file, 'r')
        if self._time_sampling is None:
            self._time_sampling = np.array(f[self._id]['mjd'])
        if self._lambda is None:
            self._lambda = np.array(f[self._id]['lambda'])

        last_ix = len(self._time_sampling) - 1
        if mjd < self._time_sampling[0]:
            mjd_ix = 0
        elif mjd > self._time_sampling[last_ix]:
            mjd_ix = last_ix
        else:
            ixes = np.argmin(np.abs(self._time_sampling - mjd))
            if isinstance(ixes, list):
                mjd_ix = ixes[0]
            else:
                mjd_ix = ixes

        return f[self._id]['flambda'][mjd_ix]


class SnanaCollection(ObjectCollection):
    '''
    This class (so far) differs from the vanilla ObjectCollection only
    in that it keeps track of where the file is which contains a library
    of SEDs for each sn
    '''
    def set_SED_file(self, SED_file):
        self._SED_file = SED_file
