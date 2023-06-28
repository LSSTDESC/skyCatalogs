import galsim
from desc.skycatalogs.utils.sn_tools import SNModel
from .base_object import BaseObject


__all__ = ['SNObject']

class SNObject(BaseObject):
    _type_name = 'sn'
    def _get_sed(self, mjd=None):
        params = self.get_native_attribute('salt2_params')
        sn = SNModel(params=params)
        if mjd < sn.mintime() or mjd > sn.maxtime():
            return None, 0.0
        # Already normalized so magnorm is zero
        return sn.get_sed(mjd), 0.0

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
