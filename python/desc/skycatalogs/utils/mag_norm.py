'''
Compute magnorm from supplied tophat SED value

Taken from jchiang87/desc_simulation_tools repo, in particular the
notebook mag_norm_for_cosmodc2
'''

import numpy as np
import pyccl as ccl

__all__ = ['MagNorm']

class MagNorm:
    def __init__(self, Omega_c=0.2648, Omega_b=0.0448, h=0.71, sigma8=0.8,
                 n_s=0.963):
        self.cosmo = ccl.Cosmology(Omega_c=Omega_c, Omega_b=Omega_b, h=h,
                                   sigma8=sigma8, n_s=n_s)
    def dl(self, z):
        aa = 1/(1 + z)
        return ccl.luminosity_distance(self.cosmo, aa)*ccl.physical_constants.MPC_TO_METER
    def __call__(self, tophat_sed_value, redshift_hubble, one_maggy=4.3442e13):
        one_Jy = 1e-26  # W/Hz/m**2
        Lnu = tophat_sed_value*one_maggy    # convert from maggies to W/Hz
        Fnu = Lnu/4/np.pi/self.dl(redshift_hubble)**2
        return -2.5*np.log10(Fnu/one_Jy) + 8.90
