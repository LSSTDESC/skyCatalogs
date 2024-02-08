from collections.abc import Iterable
import itertools
import numpy as np
import galsim
from .base_object import BaseObject, ObjectCollection


__all__ = ['SsoObject', 'SsoCollection']

class SsoObject(BaseObject):
    _type_name = 'sso'
    _solar_sed = None

    def __init__(self, ra, dec, id, object_type, belongs_to, belongs_index,
                 mjd):
        super().__init__(ra, dec, id, self._type_name, belongs_to,
                         belongs_index)
        self._mjd = mjd

    def _get_sed(self, mjd=None):
        '''
        returns a SED and magnorm
        mjd is required
        '''
        if SsoObject._solar_sed is None:
            SsoObject._solar_sed =\
                self._belongs_to._sky_catalog._sso_sed_factory.create()
            # For magnorm use the magnitude from Sorcha.  Can it be used
            # directly or are there other effects to be applied?
            # Have to find it by looking for entry for this id, this mjd
            # Do we look for specific entry or do we allow interpolation?
        return SsoObject._solar_sed, self.get_native_attribute('observedTrailedSourceMag')

    def get_gsobject_components(self, gsparams=None, rng=None, exposure=15.0):
        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)
        # For a streak we use galsim.?? (galsim.Box?)
        # To get the dimensions of the box, use ra rate, dec rate and
        # exposure length.  The first two will be in the sky catalogs
        # parquet file; the last will be passed in.
        # For now start with the simpler thing: just a point source.
        return {'this_object': galsim.DeltaFunction(gsparams=gsparams)}

    def get_observer_sed_component(self, component, mjd=None):
        if mjd is None:
            mjd = self._belongs_to._mjd
        if mjd is None:
            txt = 'SsoObject.get_observer_sed_component: no mjd specified for this call\n'
            txt += 'nor when generating object list'
            raise ValueError(txt)

        sed, magnorm = self._get_sed(mjd=mjd)

        flux_500 = np.exp(-0.9210340371976184 * magnorm)
        sed = sed.withMagnitude(0, self._bp500)
        sed = sed*flux_500

        # no extinction
        return sed

    def get_flux(self, bandpass, sed=None, mjd=None):
        if not sed:
            sed = self.get_total_observer_sed(mjd=mjd)
        if sed is None:
            return 0.0

        flux = sed.calculateFlux(bandpass)

        return flux


class SsoCollection(ObjectCollection):

    def __init__(self, ra, dec, id, hp, sky_catalog, mjd_individual=None,
                 region=None, mjd=None,
                 mask=None, readers=None, row_group=0):
        '''
        Parameters
        ----------
        ra, dec      array of float
        id           array of str
        sky_catalog  instance of SkyCatalog
        mjd_indiviual array of float or None
        region        Geometric region or string (representing file path)
        mjd_global    float or None
        mask          exclusion mask if cuts have been made due to
                      geometric region or mjd
        readers       parquet reader (in practice there is always only 1)
        row_group     int

        One of mjd_global, mjd_individual must not be None
        '''
        super().__init__(ra, dec, id, 'sso', hp, sky_catalog,
                         region=region, mjd=mjd, mask=mask,
                         readers=readers, row_group=row_group)
        self._mjds = np.array(mjd_individual)

    def __getitem__(self, key):
        '''
        Override implementation in base_object because sso objects
        must have an mjd

        Parameters
        ----------
        If key is an int (or int-like) return a single base object
        If key is a slice return a list of objects
        If key is a tuple where key[0] is iterable of int-like,
           return a list of objects
        '''

        if self._uniform_object_type:
            object_type = self._object_type_unique
        else:
            object_type = self._object_types[key]

        if isinstance(key, int) or isinstance(key, np.int64):
            return self._object_class(self._ra[key], self._dec[key],
                                      self._id[key], object_type, self, key,
                                      self._mjds[key])

        elif type(key) == slice:
            if key.start is None:
                key.start = 0
            ixdata = [i for i in range(min(key.stop, len(self._ra)))]
            ixes = itertools.islice(ixdata, key.start, key.stop, key.step)
            return [self._object_class(self._ra[i], self._dec[i], self._id[i],
                                       object_type, self, i, self._mjds[i])
                    for i in ixes]

        elif type(key) == tuple and isinstance(key[0], Iterable):
            #  check it's a list of int-like?
            return [self._object_class(self._ra[i], self._dec[i], self._id[i],
                                       object_type, self, i, self._mjds[i])
                    for i in key[0]]
