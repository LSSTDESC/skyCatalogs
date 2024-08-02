from collections.abc import Iterable
# import itertools
import numpy as np
import math
import galsim
from .base_object import BaseObject, ObjectCollection
from lsst.sphgeom import UnitVector3d, LonLat
from ..utils import normalize_sed
from .base_config_fragment import BaseConfigFragment

EXPOSURE_DEFAULT = 30.0         # seconds
__all__ = ['SsoObject', 'SsoCollection', 'SsoConfigFragment',
           'EXPOSURE_DEFAULT']

SECONDS_PER_DAY = 24.0*3600.0


class SsoObject(BaseObject):
    _type_name = 'sso'
    _solar_sed = None

    def __init__(self, ra, dec, id, object_type, belongs_to, belongs_index,
                 mjd):
        super().__init__(ra, dec, id, self._type_name, belongs_to,
                         belongs_index)
        self._mjd = mjd

    @property
    def mjd(self):
        return self._mjd

    def _get_sed(self, mjd=None):
        '''
        returns the SED
        mjd is required
        '''
        if SsoObject._solar_sed is None:
            SsoObject._solar_sed =\
                self._belongs_to._sky_catalog._sso_sed_factory.create()
            # For magnorm use the magnitude from Sorcha.  Can it be used
            # directly or are there other effects to be applied?
            # Have to find it by looking for entry for this id, this mjd
            # Do we look for specific entry or do we allow interpolation?
        magnorm = self.get_native_attribute('trailed_source_mag')
        return normalize_sed(SsoObject._solar_sed, magnorm)

    def get_gsobject_components(self, gsparams=None, rng=None,
                                streak=True):
        trail_width = 1.e-6   # in arcseconds
        exposure = self._belongs_to._exposure

        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)

        if streak:
            # For a streak make galsim.Box with length the angular distance
            # in arcseconds from initial ra, dec to final
            # Make the width very, very small
            # Then rotate appropriately
            ra = self.ra
            dec = self.dec
            # Convert from (arc?)degrees/day to degrees/sec
            ra_rate = self.get_native_attribute('ra_rate')/SECONDS_PER_DAY
            # Take out factor of cos(dec). np.cos expects radians
            ra_rate_raw = ra_rate/np.cos(math.radians(dec))
            dec_rate = self.get_native_attribute('dec_rate')/SECONDS_PER_DAY
            # ra_final is approximate since really ra_rate is a function
            # of dec, but average dec_rate is small so
            ra_final = ra + ra_rate_raw * exposure
            dec_final = dec + dec_rate * exposure

            init_v = UnitVector3d(LonLat.fromDegrees(ra, dec))
            final_v = UnitVector3d(LonLat.fromDegrees(ra_final, dec_final))
            cos_sep = init_v.dot(final_v)
            if cos_sep > 1.0:
                # Handle values like cos_sep = 1.0000000000000002
                length = 0.0
            else:
                length = np.degrees(np.arccos(cos_sep)) * 3600.0
            if np.isnan(length):
                # This will raise if cos_sep < -1.0. A value of cos_sep = -1.0
                # is almost certainly unphysical, so we want to catch these
                # cases.
                raise ValueError("SSO streak length is nan")

            if length * trail_width == 0:
                # Treat as a point source.
                return {'this_object': galsim.DeltaFunction(gsparams=gsparams)}

            gobj = galsim.Box(length, trail_width, gsparams=gsparams)

            # now rotate to direction of (ra_rate, dec_rate)
            # angle_rad = galsim.Angle(np.arctan2(dec_rate,
            #                                     (ra_rate * np.cos(dec))),
            #                          galsim.radians)
            # NOTE: cos(dec) has already been applied, so we want
            angle_rad = galsim.Angle(np.arctan2(dec_rate, ra_rate),
                                     galsim.radians)

            gobj = gobj.rotate(angle_rad)
            return {'this_object': gobj}
        else:
            # Treat as point source
            return {'this_object': galsim.DeltaFunction(gsparams=gsparams)}

    def get_observer_sed_component(self, component, mjd=None):
        if mjd is None:
            mjd = self._belongs_to._mjd
        if mjd is None:
            txt = 'SsoObject.get_observer_sed_component: no mjd specified for this call\n'
            txt += 'nor when generating object list'
            raise ValueError(txt)

        sed = self._get_sed(mjd=mjd)

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
                 mask=None, readers=None, row_group=0,
                 exposure=EXPOSURE_DEFAULT):
        '''
        Parameters
        ----------
        ra, dec        array of float
        id             array of str
        sky_catalog    instance of SkyCatalog
        mjd_individual array of float. Array of mjd values belonging
                       to the objects which will be in the new collection
        region         Geometric region
        mjd            float or None. The mjd value which was used (along with
                       region) to determine which objects should be in the
                       collection
        mask           exclusion mask if cuts have been made due to
                       geometric region or mjd
        readers        parquet reader (in practice there is always only 1)
        row_group      int

        '''
        super().__init__(ra, dec, id, 'sso', hp, sky_catalog,
                         region=region, mjd=mjd, mask=mask,
                         readers=readers, row_group=row_group)
        self._mjds = np.array(mjd_individual)
        self._exposure = exposure

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
            return [self.__getitem__(i) for i in range(self.__len__())[key]]

        elif type(key) == tuple and isinstance(key[0], Iterable):
            return [self.__getitem__(i) for i in key[0]]


class SsoConfigFragment(BaseConfigFragment):
    def __init__(self, prov, area_partition=None, data_file_type=None):
        super().__init__(prov, object_type_name='sso',
                         area_partition=area_partition,
                         data_file_type=data_file_type)
