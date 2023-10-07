import galsim
import h5py
import numpy as np
from .base_object import BaseObject,ObjectCollection

__all__ = ['SnanaObject', 'SnanaCollection']

class SnanaObject(BaseObject):
    _type_name = 'snana'

    def __init__(self, ra, dec, id, object_type, belongs_to, belongs_index,
                 redshift=None):
        super().__init__(ra, dec, id, object_type, belongs_to, belongs_index,
                         redshift)
        self._mjds = None
        self._lambda = None

        # indices of elements in mjd array bounding our mjd
        self._mjd_ix_l = None
        self._mjd_ix_u = None

    def _get_sed(self, mjd=None):
        if mjd is None:
            mjd = self._belongs_to._mjd
        mjd_start = self.get_native_attribute('start_mjd')
        mjd_end = self.get_native_attribute('end_mjd')
        if mjd < mjd_start or mjd > mjd_end:
            return None, 0.0

        return self._linear_interp_SED(mjd), 0.0

    def get_gsobject_components(self, gsparams=None, rng=None):
        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)
        return {'this_object': galsim.DeltaFunction(gsparams=gsparams)}

    def get_observer_sed_component(self, component, mjd=None):
        sed, _ = self._get_sed(mjd=mjd)
        if sed is not None:
            sed = self._apply_component_extinction(sed)
        return sed

    def get_LSST_flux(self, bandpass, sed=None, cache=True, mjd=None):
        def _to_flux(m):
            return np.exp(-0.921034037196184 * m)
        def _to_mag(f):
            return np.log(flux)/(-0.921034037196184)

        flux = super().get_LSST_flux(bandpass, sed, mjd)

        mjd_ix_l, mjd_ix_u, mjd_fraction = self._find_mjd_interval(mjd)
        if flux > 0.0:
            with h5py.File(self._belongs_to._SED_file, 'r') as f:
                cors = f[self._id][f'magcor_{bandpass}']
                rough_mag = _to_mag(flux)

                # interpolate corrections
                if mjd_ix_l == mjd_ix_u:
                    mag_cor = cors[mjd_ix_l]
                else:
                    mag_cor = cors[mjd_ix_l] + mjd_fraction *\
                        (cors[mjd_ix_u] - cors[mjd_ix_l])
                corrected_flux = _to_flux(rough_mag + mag_cor)

                #dbg = True
                dbd = False
                if dbg:
                    print(f'Band {bandpass} uncorrected flux: {flux}')
                    print(f'                uncorrected mag: {rough_mag}')
                    print(f'                mag correction: {mag_cor}')

                if cache:
                    att = f'lsst_flux_{band}'
                    setattr(self, att, corrected_flux)
        return corrected_flux

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
        A triple:     index below, index above, and fraction of distance
                      mjd is from entry below to entry above
        '''
        if not mjd:
            mjd = self._belongs_to._mjd
            store = True
        else:
            store = (mjd == self._belongs_to._mjd)

        if store:
            if self._mjd_ix_l:
                # just return previously-computed values
                return self._mjd_ix_l, self._mjd_ix_u, self._mjd_fraction

        if not self._mjds:
            with h5py.File(self._belongs_to._SED_file, 'r') as f:
                self._mjds = np.array(f[self._id]['mjd'])

        last_ix = len(self._mjds) - 1
        if mjd <= self._mjds[0]:
            mjd_ix_l = 0
            mjd_ix_u = 0
            mjd_fraction = None
        elif mjd >= self._mjds[last_ix]:
            mjd_ix_l = last_ix
            mjd_ix_u = last_ix
            mjd_fraction = None
        else:
            mjd_ix_u = np.argmin((mjd - self._mjds ) > 0)
            mjd_ix_l = mjd_ix_u - 1
            mjds = self._mjds
            # Where is mjd relative to mjds on either side in the array?
            mjd_fraction = (mjd - mjds[mjd_ix_l]) /\
                (mjds[mjd_ix_u] - mjds[mjd_ix_l])
        if store:
            self._mjd_ix_l = mjd_ix_l
            self._mjd_ix_u = mjd_ix_u
            self._mjd_fraction = mjd_fraction

        return mjd_ix_l, mjd_ix_u, mjd_fraction

    def _read_nearest_SED(self, mjd=None):
        '''
        Find row with closest mjd and return galsim.SED generated
        from it
        '''

        if self._mjd_ix_l is None:
            mjd_ix_l, mjd_ix_u, mjd_fraction = self._find_mjd_interval()
        else:
            mjd_ix_l = self._mjd_ix_l
            mjd_ix_u = self._mjd_ix_u
            mjd_fraction = self._mjd_fraction

        f = h5py.File(self._belongs_to._SED_file, 'r')
        if self._lambda is None:
            self._lambda = np.array(f[self._id]['lambda'])

        if mjd_ix_l == mjd_ix_u:
            mjd_ix = mjd_ix_l

        else:
            if mjd_fraction < 0.5:
                mjd_ix = mjd_ix_l
            else:
                mjd_ix = mjd_ix_u

        lut = galsim.LookupTable(f[self._id]['lambda'],
                                 f[self._id]['flambda'][mjd_ix],
                                 interpolant='linear')
        return galsim.SED(lut, wave_type='A', flux_type='flambda')

    def _linear_interp_SED(self, mjd=None):
        '''
        Return galsim SED obtained by interpolating between SEDs
        for nearest mjds among the templates
        '''

        if not self._mjd_ix_l:
            mjd_ix_l,mjd_ix_u,mjd_fraction = self._find_mjd_interval(mjd)
        else:
            mjd_ix_l = self._mjd_ix_l
            mjd_ix_u = self._mjd_ix_u
            mjd_fraction = self._mjd_fraction

        with h5py.File(self._belongs_to._SED_file, 'r') as f:
            if self._mjds is None:
                self._mjds = np.array(f[self._id]['mjd'])
            if self._lambda is None:
                self._lambda = np.array(f[self._id]['lambda'])

            if self._mjd_ix_l == self._mjd_ix_u:
                flambda = f[slef._id]['flamba'][self._mjd_ix_l]
            else:
                mjd_ix = self._mjd_ix_u
                mjds = self._mjds
                below = f[self._id]['flambda'][mjd_ix - 1]
                above = f[self._id]['flambda'][mjd_ix]
                flambda = below + mjd_fraction * (above - below)

            lut = galsim.LookupTable(f[self._id]['lambda'],
                                     flambda,
                                     interpolant='linear')
        return galsim.SED(lut, wave_type='A', flux_type='flambda')


class SnanaCollection(ObjectCollection):
    '''
    This class (so far) differs from the vanilla ObjectCollection only
    in that it keeps track of where the file is which contains a library
    of SEDs for each sn
    '''
    def set_SED_file(self, SED_file):
        self._SED_file = SED_file
