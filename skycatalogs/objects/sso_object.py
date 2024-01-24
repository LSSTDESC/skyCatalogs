import os
import numpy as np
import galsim
from .base_object import BaseObject

__all__ = ['SsoObject']

class SsoObject(BaseObject):
    _type_name = 'sso'

    def _get_sed(self, mjd=None):
        '''
        returns a SED and magnorm
        mjd is required
        '''

        # Always return the same SED.  Maybe need a class method
        # to load it?
        # For magnorm use the magnitude from Sorcha.  Can it be used
        # directly or are there other effects to be applied?
        # Have to find it by looking for entry for this id, this mjd
        # Do we look for specific entry or do we allow interpolation?
        pass

    def get_gsobject_components(self, gsparams=None, rng=None, exposure=15.0):
        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)
        # For a streak we use galsim.?? (galsim.Box?)
        # To get the dimensions of the box, use ra rate, dec rate and
        # exposure length.  The first two will be in the sky catalogs
        # parquet file; the last will be passed in.
        # For now start with the simpler thing: just a point source.
        return {'this_object': galsim.DeltaFunction(gsparams=gsparams)}

    def get_observer_sed_component(self, component, mjd=None, exposure=15.0):
        if mjd is None:
            mjd = self._belongs_to._mjd
        if mjd is None:
            txt = 'SsoObject.get_observer_sed_component: no mjd specified for this call\n'
            txt += 'nor when generating object list'
            raise ValueError(txt)
            
        sed, magnorm = self._get_sed(mjd=mjd)

        flux_500 = np.exp(-0.9210340371976184 * magnorm)
        sed = self.withMagnitude(0, self._bp500)
        sed = sed*flux_500

        # no extinction 
        return sed

    def get_flux(self, bandpass, sed=None, mjd=None, exposure=15.0):
        if not sed:
            sed = self.get_total_observer_sed(mjd=mjd, exposure=exposure)
        if sed is None:
            return 0.0

        flux = sed.calculateFlux(bandpass)

        return flux
            

    
