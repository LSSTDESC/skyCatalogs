import galsim
from skycatalogs.utils.sn_tools import SncosmoModel
from .base_object import BaseObject


__all__ = ['SncosmoObject']


class SncosmoObject(BaseObject):
    _type_name = 'sncosmo'

    def _get_sed(self, mjd=None):
        params = self.get_native_attribute('salt2_params')
        sn = SncosmoModel(params=params)

        if mjd < sn.mintime() or mjd > sn.maxtime():
            return None, 0.0
        # Already normalized so magnorm is zero
        return sn.get_sed(mjd), 0.0

    def get_gsobject_components(self, gsparams=None, rng=None):
        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)
        return {'this_object': galsim.DeltaFunction(gsparams=gsparams)}

    def get_observer_sed_component(self, component, mjd=None):
        if mjd is None:
            mjd = self._belongs_to._mjd
        if mjd is None:
            txt = 'SncosmoObject._get_sed: no mjd specified for this call\n'
            txt += 'nor when generating object list'
            raise ValueError(txt)
        sed, _ = self._get_sed(mjd=mjd)
        if sed is not None:
            sed = self._apply_component_extinction(sed)
        return sed

    def get_LSST_flux(self, band, sed=None, cache=False, mjd=None):
        # There is usually no reason to cache flux for SNe, in fact it could
        # cause problems
        return super().get_LSST_flux(band, sed=sed, cache=cache, mjd=mjd)
