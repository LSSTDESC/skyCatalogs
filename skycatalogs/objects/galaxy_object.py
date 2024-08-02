import numpy as np
import galsim

from .base_object import BaseObject
from .base_config_fragment import BaseConfigFragment
from skycatalogs.utils.translate_utils import form_object_string

__all__ = ['GalaxyObject', 'GalaxyConfigFragment']


class GalaxyObject(BaseObject):
    _type_name = 'galaxy'

    def _get_sed(self, component=None, resolution=None):
        '''
        Return sed for a galaxy component
        Parameters
        ----------
        component    one of 'bulge', 'disk', 'knots' for now. Other components
                     may be supported.
        resolution   desired resolution of lambda in nanometers.

        Returns
        -------
        galsim.SED object
        '''
        if component not in ['disk', 'bulge', 'knots']:
            raise ValueError(f'Cannot fetch SED for component type {component}')

        th_val = self.get_native_attribute(f'sed_val_{component}')
        if th_val is None:   # values for this component are not in the file
            raise ValueError(f'{component} not part of this catalog')

        # if values are all zeros or nearly no point in trying to convert
        if max(th_val) < np.finfo('float').resolution:
            return None

        z_h = self.get_native_attribute('redshift_hubble')
        z = self.get_native_attribute('redshift')

        sky_cat = self._belongs_to._sky_catalog
        sed = sky_cat.observed_sed_factory.create(th_val, z_h, z,
                                                  resolution=resolution)
        return sed

    def get_knot_size(self, z):
        """
        Return the angular knot size. Knots are modelled as the same physical
        size
        """
        # Deceleration paramameter
        q = -0.55
        # Angular diameter scaling approximation in pc
        dA = (3e9/q**2)*(z*q+(q-1)*(np.sqrt(2*q*z+1)-1))/(1+z)**2*(1.4-0.53*z)
        # Using typical knot size 250pc, convert to sigma in arcmin
        if z < 0.6:
            return 206264.8*250/dA/2.355
        else:
            # Above z=0.6, fractional contribution to post-convolved size
            # is <20% for smallest Roman PSF size, so can treat as point source
            return None

    def get_wl_params(self):
        """Return the weak lensing parameters, g1, g2, mu."""
        gamma1 = self.get_native_attribute('shear_1')
        gamma2 = self.get_native_attribute('shear_2')
        kappa = self.get_native_attribute('convergence')
        # Compute reduced shears and magnification.
        g1 = gamma1/(1. - kappa)    # real part of reduced shear
        g2 = gamma2/(1. - kappa)    # imaginary part of reduced shear
        mu = 1./((1. - kappa)**2 - (gamma1**2 + gamma2**2))  # magnification
        return g1, g2, mu

    def get_total_observer_sed(self, mjd=None):
        """
        Return the SED summed over SEDs for the individual SkyCatalog
        components.
        """
        sed = super().get_total_observer_sed()

        if sed is None:
            return sed

        _, _, mu = self.get_wl_params()
        sed *= mu
        return sed

    def get_gsobject_components(self, gsparams=None, rng=None):

        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)

        if rng is None:
            rng = galsim.BaseDeviate(int(self.id))

        obj_dict = {}
        for component in self.subcomponents:
            # knots use the same major/minor axes as the disk component.
            my_component = 'disk' if component != 'bulge' else 'bulge'
            a = self.get_native_attribute(
                f'size_{my_component}_true')
            b = self.get_native_attribute(
                f'size_minor_{my_component}_true')
            assert a >= b
            hlr = (a*b)**0.5   # approximation for half-light radius

            e1 = self.get_native_attribute(
                f'ellipticity_1_{my_component}_true')
            e2 = self.get_native_attribute(
                f'ellipticity_2_{my_component}_true')
            # NOTE: The minus signs in the next executable line are
            # needed specifically for generating DC2-like from galaxy
            # catalogs in the cosmoDC2 format. They are included here in
            # order to reproduce the effect of adding 90 degrees to
            # position angle in the old code. Newer input galaxy catalogs
            # most likely will not need this adjustment
            shear = galsim.Shear(g1=-e1, g2=-e2)

            if component == 'knots':
                npoints = self.get_native_attribute('n_knots')
                assert npoints > 0
                obj = galsim.RandomKnots(npoints=npoints,
                                         profile=obj_dict['disk'], rng=rng,
                                         gsparams=gsparams)
                z = self.get_native_attribute('redshift')
                size = self.get_knot_size(z)  # get knot size
                if size is not None:
                    obj = galsim.Convolve(obj, galsim.Gaussian(sigma=size))
                obj_dict[component] = obj
            else:
                n = self.get_native_attribute(f'sersic_{component}')
                # Quantize the n values at 0.05 so that galsim can
                # possibly amortize sersic calculations from the previous
                # galaxy.
                n = round(n*20.)/20.
                obj = galsim.Sersic(n=n, half_light_radius=hlr,
                                    gsparams=gsparams)
                obj_dict[component] = obj._shear(shear)

        # Apply lensing
        g1, g2, mu = self.get_wl_params()
        for component in self.subcomponents:
            obj_dict[component] = obj_dict[component]._lens(g1, g2, mu)
        return obj_dict

    def get_observer_sed_component(self, component, mjd=None, resolution=None):
        sed = self._get_sed(component=component, resolution=resolution)
        if sed is not None:
            sed = self._apply_component_extinction(sed)

        return sed

    def get_instcat_entry(self, band='r', component=None):
        '''
        Return the string corresponding to instance catalog line
        Parameters:
            band       One of ['u', 'g', 'r', 'i', 'z', 'y']
            component  Required iff the object has subcomponents (i.e.,
                       object type is 'galaxy')
        Returns: A string formatted like a line in an instance catalog
        '''
        if component not in self.subcomponents:
            return ''
        return form_object_string(self, band, component)


class GalaxyConfigFragment(BaseConfigFragment):
    def __init__(self, prov, cosmology, tophat_bins,
                 area_partition=None, data_file_type=None):
        super().__init__(prov, object_type_name='galaxy',
                         area_partition=area_partition,
                         data_file_type=data_file_type)
        self._opt_dict['Cosmology'] = cosmology
        self._tophat_bins = tophat_bins

    def make_fragment(self):
        data = self.generic_create()
        data['tophat']['bins'] = self._tophat_bins
        return data
