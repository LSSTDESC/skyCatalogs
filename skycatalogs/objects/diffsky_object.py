import galsim
import numpy as np
from .base_object import BaseObject
from .base_config_fragment import BaseConfigFragment

__all__ = ['DiffskyObject', 'DiffskyConfigFragment']


class DiffskyObject(BaseObject):
    _type_name = 'diffsky_galaxy'

    # sersic values are constant for diffsky galaxies
    _sersic_disk = 1
    _sersic_bulge = 4

    def _get_sed(self, component=None, resolution=None):
        '''
        Return sed and mag_norm for a galaxy component or for a star
        Parameters
        ----------
        component    one of 'bulge', 'disk', 'knots' for now. Other components
                     may be supported.  Ignored for stars
        resolution   desired resolution of lambda in nanometers. Ignored
                     for stars.

        Returns
        -------
        galsim.SED object
        '''

        if component not in ['disk', 'bulge', 'knots']:
            raise ValueError(f'Cannot fetch SED for component type {component}')

        if not hasattr(self, '_seds'):
            z_h = self.get_native_attribute('redshiftHubble')
            z = self.get_native_attribute('redshift')
            pixel = self.partition_id

            sky_cat = self._belongs_to._sky_catalog
            self._seds = sky_cat.observed_sed_factory.create(pixel, self.id,
                                                             z_h, z)
        return self._seds[component]

    def get_knot_size(self, z):
        """
        Return the angular knot size. Knots are modelled as the same
        physical size
        """
        # Deceleration paramameter
        q = -0.5

        if z >= 0.6:
            # Above z=0.6, fractional contribution to post-convolved size
            # is <20% for smallest Roman PSF size, so can treat as point source
            # This also ensures sqrt in formula below has a
            # non-negative argument
            return None

        # Angular diameter scaling approximation in pc
        dA = (3e9/q**2)*(z*q+(q-1)*(np.sqrt(2*q*z+1)-1))/(1+z)**2*(1.4-0.53*z)
        # Using typical knot size 250pc, convert to sigma in arcmin
        return 206264.8*250/dA/2.355

    def get_knot_n(self, rng=None):
        """
        Return random value for number of knots based on galaxy sm.
        """
        if rng is not None:
            ud = galsim.UniformDeviate(rng)
        else:
            ud = galsim.UniformDeviate(int(self.id))
        sm = np.log10(self.get_native_attribute('um_source_galaxy_obs_sm'))
        m = (50-3)/(12-6)  # (knot_n range)/(logsm range)
        n_knot_max = m*(sm-6)+3
        n_knot = int(ud()*n_knot_max)  # random n up to n_knot_max
        if n_knot == 0:
            n_knot += 1  # need at least 1 knot
        return n_knot

    def get_wl_params(self):
        """Return the weak lensing parameters, g1, g2, mu."""
        gamma1 = self.get_native_attribute('shear1')
        gamma2 = self.get_native_attribute('shear2')
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
            my_cmp = 'disk' if component != 'bulge' else 'spheroid'
            hlr = self.get_native_attribute(f'{my_cmp}HalfLightRadiusArcsec')

            # Get ellipticities saved in catalog. Not sure they're what
            # we need
            e1 = self.get_native_attribute(f'{my_cmp}Ellipticity1')
            e2 = self.get_native_attribute(f'{my_cmp}Ellipticity2')
            shear = galsim.Shear(g1=e1, g2=e2)

            if component == 'knots':
                npoints = self.get_knot_n()
                assert npoints > 0
                knot_profile = galsim.Sersic(n=self._sersic_disk,
                                             half_light_radius=hlr/2.,
                                             gsparams=gsparams)
                knot_profile = knot_profile._shear(shear)
                obj = galsim.RandomKnots(npoints=npoints,
                                         profile=knot_profile, rng=rng,
                                         gsparams=gsparams)
                z = self.get_native_attribute('redshift')
                size = self.get_knot_size(z)  # get knot size
                if size is not None:
                    obj = galsim.Convolve(obj, galsim.Gaussian(sigma=size))
                obj_dict[component] = obj
            else:
                n = self._sersic_disk if component == 'disk' else self._sersic_bulge
                obj = galsim.Sersic(n=n, half_light_radius=hlr,
                                    gsparams=gsparams)
                obj_dict[component] = obj._shear(shear)

        # Apply lensing
        g1, g2, mu = self.get_wl_params()
        for component in self.subcomponents:
            obj_dict[component] = obj_dict[component]._lens(g1, g2, mu)
        return obj_dict

    def get_observer_sed_component(self, component, mjd=None, resolution=None):
        sed = self._get_sed(component)
        if sed is not None:
            sed = self._apply_component_extinction(sed)
        return sed


class DiffskyConfigFragment(BaseConfigFragment):
    def __init__(self, prov, cosmology,
                 area_partition=None, data_file_type=None):
        super().__init__(prov, object_type_name='diffsky_galaxy',
                         area_partition=area_partition,
                         data_file_type=data_file_type)
        self._opt_dict['Cosmology'] = cosmology
