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
        self._mjds = None
        self._lambda = None

    def _get_sed(self, mjd=None):
        if mjd is None:
            return None, 0.0
        #model = self.get_native_attribute('model_name')
        #param_names = self.get_native_attribute('model_param_names')
        #param_values = self.get_native_attribute('model_param_values')
        mjd_start = self.get_native_attribute('start_mjd')
        mjd_end = self.get_native_attribute('end_mjd')
        if mjd < mjd_start or mjd > mjd_end:
            return None, 0.0
        # For now find the SED with mjd closest to ours and return that
        # What should we do about magnorm?
        return self._linear_interp_SED(mjd), 0.0


    def get_gsobject_components(self, gsparams=None, rng=None):
        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)
        return {'this_object': galsim.DeltaFunction(gsparams=gsparams)}

    def get_observer_sed_component(self, component, mjd=None):
        sed, _ = self._get_sed(mjd=mjd)
        # Can't do this now because don't have native attributes for
        # MW_av, MW_rv
        #if sed is not None:
        #    sed = self._apply_component_extinction(sed)
        return sed

    def _read_nearest_SED(self, mjd):
        # Find row with closest mjd and return galsim.SED generated
        # from it

        f = h5py.File(self._belongs_to._SED_file, 'r')
        if self._mjds is None:
            self._mjds = np.array(f[self._id]['mjd'])
        if self._lambda is None:
            self._lambda = np.array(f[self._id]['lambda'])

        last_ix = len(self._mjds) - 1
        if mjd < self._mjds[0]:
            mjd_ix = 0
        elif mjd > self._mjds[last_ix]:
            mjd_ix = last_ix
        else:
            ixes = np.argmin(np.abs(self._mjds - mjd))
            if isinstance(ixes, list):
                mjd_ix = ixes[0]
            else:
                mjd_ix = ixes

        lut = galsim.LookupTable(f[self._id]['lambda'],
                                 f[self._id]['flambda'][mjd_ix],
                                 interpolant='linear')
        return galsim.SED(lut, wave_type='A', flux_type='flambda')

    def _linear_interp_SED(self, mjd):
        f = h5py.File(self._belongs_to._SED_file, 'r')
        if self._mjds is None:
            self._mjds = np.array(f[self._id]['mjd'])
        if self._lambda is None:
            self._lambda = np.array(f[self._id]['lambda'])

        last_ix = len(self._mjds) - 1
        if mjd <= self._mjds[0]:
            flambda = f[self._id]['flambda'][0]
        elif mjd >= self._mjds[last_ix]:
            flambda = f[self._id]['flambda'][last_ix]
        else:
            # First index where condition is False
            mjd_ix = np.argmin((mjd - self._mjds ) > 0)
            mjds = self._mjds
            below = f[self._id]['flambda'][mjd_ix - 1]
            above = f[self._id]['flambda'][mjd_ix]
            ratio =  (mjd - mjds[mjd_ix - 1])/(mjds[mjd_ix] - mjds[mjd_ix - 1])
            flambda = below + ratio * (above - below)

        lut = galsim.LookupTable(f[self._id]['lambda'],
                                 flambda,
                                 interpolant='linear')
        return galsim.SED(lut, wave_type='A', flux_type='flambda')



class SnanaCollection(ObjectCollection):
    '''
    This class (so far) differs from the vanilla ObjectCollection only
    in that it keeps track of where the file is which contains a library
    of SEDs for each sn
    '''
    def set_SED_file(self, SED_file):
        self._SED_file = SED_file
