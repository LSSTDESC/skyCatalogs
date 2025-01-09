# import os
import galsim
# import h5py
from .base_object import BaseObject, ObjectCollection, load_lsst_bandpasses
# from ..utils import normalize_sed
from .base_config_fragment import BaseConfigFragment

__all__ = ['TrilegalObject', 'TrilegalCollection', 'TrilegalConfigFragment']


class TrilegalObject(BaseObject):

    _type_name = 'trilegal'

    def __init__(self, ra, dec, id, object_type, belongs_to, belongs_index):
        super().__init__(ra, dec, id, self._type_name, belongs_to,
                         belongs_index)

    def _get_sed(self, mjd=None, redshift=0):
        '''
        '''
        factory = self._belongs_to._sky_catalog._trilegal_sed_factory
        return factory.get_sed(self)    # unextincted

    def get_gsobject_components(self, gsparams=None, rng=None):
        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)
        return {'this_object': galsim.DeltaFunction(gsparams=gsparams)}

    def get_observer_sed_component(self, component, mjd=None):
        '''
        Apply extinction, normalize
        '''
        sed = self._get_sed(mjd=mjd)

        if sed is not None:
            sed = self._apply_component_extinction(sed)
            sed = sed.thin()
            imag = self.get_native_attribute('imag')
            sed = sed.withMagnitude(imag,
                                    self._belongs_to._lsst_bandpasses['i'])
        return sed

    def _get_dust(self):
        "Return the Av, Rv parameters for internal and Milky Way extinction."
        internal_av = 0
        internal_rv = 1.

        galactic_av = self.get_native_attribute('av')
        # No native attribute for rv. Use standard value
        galactic_rv = 3.1
        return internal_av, internal_rv, galactic_av, galactic_rv


class TrilegalCollection(ObjectCollection):
    def __init__(self, ra, dec, id, object_type, hp, sky_catalog,
                 region=None, mjd=None,
                 mask=None, readers=None, row_group=0):
        '''
        Parameters
        ----------
        ra, dec        array of float
        id             array of str
        hp             int healpixel
        object_type    Should be 'trilegal'
        sky_catalog    instance of SkyCatalog
        region         Geometric region
        mjd            float or None. The mjd value which was used (along with
                       region) to determine which objects should be in the
                       collection
        mask           exclusion mask if cuts have been made due to
                       geometric region or mjd
        readers        parquet reader (in practice there is always only 1)
        row_group      int

        '''
        super().__init__(ra, dec, id, object_type, hp, sky_catalog,
                         region=region, mjd=mjd, mask=mask,
                         readers=readers, row_group=row_group)
        self._lsst_bandpasses = load_lsst_bandpasses()

        # See also classes TrilegalSedFactory, TrilegalSedFile, _SEDBatch in
        # sed_tools.py


        # Is this necessary? Probably not
        our_config = self.sky_catalog._config['object_types'][object_type]


class TrilegalConfigFragment(BaseConfigFragment):
    def __init__(self, prov, area_partition=None, data_file_type=None):
        super().__init__(prov, object_type_name='trilegal',
                         area_partition=area_partition,
                         data_file_type=data_file_type)
