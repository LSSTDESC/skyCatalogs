import bisect
import galsim
import h5py
import numpy as np
from .base_object import BaseObject, ObjectCollection
from .base_config_fragment import BaseConfigFragment
from skycatalogs.utils.exceptions import SkyCatalogsRuntimeError

__all__ = ['SnanaObject', 'SnanaCollection', 'SnanaConfigFragment']

_map_SNANA_bands = {'R062': 'R',
                    'Z087': 'Z',
                    'Y106': 'Y',
                    'J129': 'J',
                    'H158': 'H',
                    'F184': 'F',
                    'K213': 'K',
                    'W146': 'W',
                    'SNPrism': 'S',
                    'Grism_0thOrder': 'G0',
                    'Grism_1stOrder': 'G1',
                   }


class SnanaObject(BaseObject):
    _type_name = 'snana'

    def __init__(self, ra, dec, id, object_type, belongs_to, belongs_index):
        super().__init__(ra, dec, id, object_type, belongs_to, belongs_index)
        self._mjds = None
        self._lambda = None

        # indices of elements in mjd array bounding our mjd
        self._mjd_ix_l = None
        self._mjd_ix_u = None

    def _get_sed(self, mjd=None):
        mjd_start = self.get_native_attribute('start_mjd')
        mjd_end = self.get_native_attribute('end_mjd')
        if mjd < mjd_start or mjd > mjd_end:
            return None

        return self._linear_interp_SED(mjd)

    def _apply_flux_correction(self, flux, snana_band, mjd):
        def _flux_ratio(mag):
            # -0.9210340371976184 = -np.log(10)/2.5.
            return np.exp(-0.921034037196184 * mag)

        if flux < 0:
            raise SkyCatalogsRuntimeError('Negative flux')

        if flux == 0.0:
            return flux

        mjd_ix_l, mjd_ix_u, mjd_fraction = self._find_mjd_interval(mjd)

        with h5py.File(self._belongs_to._SED_file, 'r') as f:
            try:
                cors = f[self._id][snana_band]
            except KeyError:
                # nothing else to do
                return flux

            # interpolate corrections if we can.  Correction array
            # may include nans.
            if np.isnan(cors[mjd_ix_l]) or np.isnan(cors[mjd_ix_u]):
                txt = f'Cannot apply flux correction to SN {self._id} due to nan in correction array'
                self._logger.warn(txt)
                return flux

            if mjd_ix_l == mjd_ix_u:
                mag_cor = cors[mjd_ix_l]
            else:
                mag_cor = cors[mjd_ix_l] + mjd_fraction *\
                    (cors[mjd_ix_u] - cors[mjd_ix_l])

        # dbg = True
        dbg = False

        # Do everything in flux units
        flux_cor = _flux_ratio(mag_cor)
        corrected_flux = flux * flux_cor

        if dbg:
            print(f'Band {snana_band} uncorrected flux: {flux}')
            print(f'                  mag correction: {mag_cor}')
            print(f' multiplicative flux correction: {flux_cor}')

        return corrected_flux

    def get_gsobject_components(self, gsparams=None, rng=None):
        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)
        return {'this_object': galsim.DeltaFunction(gsparams=gsparams)}

    def get_observer_sed_component(self, component, mjd=None):
        if mjd is None:
            mjd = self._belongs_to._mjd
        if mjd is None:
            txt = 'SnananObject.get_observer_sed_component: no mjd specified for this call\n'
            txt += 'nor when generating object list'
            raise ValueError(txt)
        sed = self._get_sed(mjd=mjd)
        if sed is not None:
            sed = self._apply_component_extinction(sed)
        return sed

    def get_LSST_flux(self, band, sed=None, cache=False, mjd=None):
        # There is usually no reason to cache flux for SNe, in fact it could
        # cause problems. If flux has been cached and then this routine
        # is called again with a different value of mjd, it would
        # return the wrong answer.

        flux = super().get_LSST_flux(band, sed=sed, cache=cache, mjd=mjd)
        return self._apply_flux_correction(flux, f'magcor_{band}', mjd)

    def get_roman_flux(self, band, sed=None, cache=False, mjd=None):
        # There is usually no reason to cache flux for SNe, in fact it could
        # cause problems. If flux has been cached and then this routine
        # is called again with a different value of mjd, it would
        # return the wrong answer.

        flux = super().get_roman_flux(band, sed=sed, cache=cache, mjd=mjd)
        return self._apply_flux_correction(flux,
                                           f'magcor_{_map_SNANA_bands[band]}',
                                           mjd)

    def _find_mjd_interval(self, mjd=None):
        '''
        Find indices into mjd array of elements bounding our mjd
        Also compute and constant needed for interpolation.  If we're
        using "standard" mjd, also store these numbers.

        Parameters
        ----------
        mjd     float   If None use the one stored in our ObjectCollection

        Returns
        -------
        A tuple:     index below, index above, and fraction of distance
                     mjd is from entry below to entry above
        '''
        if not mjd:
            mjd = self._belongs_to._mjd
            store = True
        else:
            store = (mjd == self._belongs_to._mjd)

        if store:
            if self._mjd_ix_l is not None:
                # just return previously-computed values
                return self._mjd_ix_l, self._mjd_ix_u, self._mjd_fraction

        if self._mjds is None:
            with h5py.File(self._belongs_to._SED_file, 'r') as f:
                self._mjds = np.array(f[self._id]['mjd'])
        mjds = self._mjds

        mjd_fraction = None
        index = bisect.bisect(mjds, mjd)
        if index == 0:
            mjd_ix_l = mjd_ix_u = 0
        elif index == len(mjds):
            mjd_ix_l = mjd_ix_u = index - 1
        else:
            mjd_ix_l = index - 1
            mjd_ix_u = index
            mjd_fraction = (mjd - mjds[mjd_ix_l]) /\
                (mjds[mjd_ix_u] - mjds[mjd_ix_l])
        if store:
            self._mjd_ix_l = mjd_ix_l
            self._mjd_ix_u = mjd_ix_u
            self._mjd_fraction = mjd_fraction

        return mjd_ix_l, mjd_ix_u, mjd_fraction

    def _linear_interp_SED(self, mjd=None):
        '''
        Return galsim SED obtained by interpolating between SEDs
        for nearest mjds among the templates
        '''
        mjd_ix_l, mjd_ix_u, mjd_fraction = self._find_mjd_interval(mjd)

        with h5py.File(self._belongs_to._SED_file, 'r') as f:
            if self._mjds is None or self._lambda is None:
                self._mjds = np.array(f[self._id]['mjd'])
                self._lambda = np.array(f[self._id]['lambda'])

            if mjd_ix_l == mjd_ix_u:
                flambda = f[self._id]['flamba'][mjd_ix_l]
            else:
                mjd_ix = mjd_ix_u
                below = f[self._id]['flambda'][mjd_ix - 1]
                above = f[self._id]['flambda'][mjd_ix]
                flambda = below + mjd_fraction * (above - below)

            lut = galsim.LookupTable(f[self._id]['lambda'],
                                     flambda,
                                     interpolant='linear')
        return galsim.SED(lut, wave_type='A', flux_type='flambda')


class SnanaCollection(ObjectCollection):
    '''
    This class  differs from the vanilla ObjectCollection only
    in that
    * it keeps track of where the file is which contains a library
      of SEDs for each sn
    * it issues a warning if mjd is None
    '''
    def set_SED_file(self, SED_file):
        self._SED_file = SED_file

    def __init__(self, ra, dec, id, object_type, partition_id, sky_catalog,
                 region=None, mjd=None, mask=None, readers=None, row_group=0):
        # Normally mjd should be specified
        if mjd is None:
            sky_catalog._logger.warning('Creating SnanaCollection with no mjd value.')
            sky_catalog._logger.warning('Transient collections normally have non-None mjd')
        super().__init__(ra, dec, id, object_type, partition_id,
                         sky_catalog, region=region, mjd=mjd, mask=mask,
                         readers=readers, row_group=row_group)


class SnanaConfigFragment(BaseConfigFragment):
    def __init__(self, prov, area_partition=None, data_file_type=None):
        super().__init__(prov, object_type_name='snana',
                         area_partition=area_partition,
                         data_file_type=data_file_type)
