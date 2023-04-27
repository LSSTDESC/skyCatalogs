import numpy as np
from astropy import units as u
import sncosmo
import galsim


class SNObject(sncosmo.Model):
    def __init__(self, source='salt2-extended', params=None):
        '''
        params - dict of params suitable for the model

        See also https://sncosmo.readthedocs.io/en/stable/index.html

        '''
        # The following explicitly turns off host and Milky Way extinction.
        dust = sncosmo.F99Dust()
        super().__init__(source=source,
                         effects=[dust, dust],
                         effect_names=['host', 'mw'],
                         effect_frames=['rest', 'obs'])
        self.set(mwebv=0., hostebv=0.)
        if params:
            self.set(**params)

    def get_sed(self, mjd, npts=1000):
        """
        Return the SED in the observer frame at the requested time.
        """
        wl = np.linspace(self.minwave(), self.maxwave(), npts)
        flambda = self.flux(mjd, wl)
        lut = galsim.LookupTable(wl, flambda)
        return galsim.SED(lut, wave_type=u.Angstrom, flux_type='flambda')
