import numpy as np
from astropy import units as u
import sncosmo
import galsim

__all__ = ['SncosmoModel']


class SncosmoModel(sncosmo.Model):
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
        self.redshift = 0
        if params:
            self.set(**params)
            self.redshift = params['z']

    def get_sed(self, mjd, npts=1000):
        """
        Return the SED in the observer frame at the requested time.
        """
        wl = np.linspace(self.minwave(), self.maxwave(), npts)

        # prepend 0 bins
        n_bins = int(self.minwave()) - 1
        pre_wl = [float(i) for i in range(n_bins)]
        pre_val = [0.0 for i in range(n_bins)]
        flambda = self.flux(mjd, wl)
        wl = np.insert(wl, 0, pre_wl)
        flambda = np.insert(flambda, 0, pre_val)
        lut = galsim.LookupTable(wl, flambda, interpolant='linear')
        return galsim.SED(lut, wave_type=u.Angstrom, flux_type='flambda')
