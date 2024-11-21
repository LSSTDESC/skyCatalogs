import os
import galsim
import h5py
from .base_object import BaseObject, BaseCollection
from ..utils import normalize_sed
from .base_config_fragment import BaseConfigFragment

__all__ = ['TrilegalObject', 'TrilegalConfigFragment']


class TrilegalObject(BaseObject):

    _type_name = 'trilegal'

    def _get_sed(self, mjd=None, redshift=0):
        '''
        We'll need mjd when/if variable stars are supported. For now
        it's ignored.

        To be rewritten once form of SED files is known
        Will SEDs already be normalized?
        '''

        # mag_norm = self.get_native_attribute('magnorm')
        # fpath = os.path.join(os.getenv('SIMS_SED_LIBRARY_DIR'),
        #                      self.get_native_attribute('sed_filepath'))

        # return normalize_sed(self._get_sed_from_file(fpath), mag_norm)

        return None

    def get_gsobject_components(self, gsparams=None, rng=None):
        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)
        return {'this_object': galsim.DeltaFunction(gsparams=gsparams)}

    def get_observer_sed_component(self, component, mjd=None):
        sed = self._get_sed(mjd=mjd)

        if sed is not None:
            sed = self._apply_component_extinction(sed)
        return sed


class TrilegalCollection(BaseCollection):
    def __init__(self, ra, dec, id, hp, sky_catalog,
                 region=None, mjd=None,
                 mask=None, readers=None, row_group=0):
        '''
        Parameters
        ----------
        ra, dec        array of float
        id             array of str
        hp             int healpixel
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
        super().__init__(ra, dec, id, 'trilegal', hp, sky_catalog,
                         region=region, mjd=mjd, mask=mask,
                         readers=readers, row_group=row_group)

        # Assume for now SED files, if any, are in same directory as
        # parquet files. Make a SED reader which will manage SED data
        # for this collection
        our_config = self.sky_catalog._config['object_types'][self._object_type]
        self._sed_reader = _TrilegalSEDReader(hp, our_config,
                                              sky_catalog._cat_dir)


class _SEDBatch:
    def __init__(self, batch_name, hdf_file):
        self._batch_name = batch_name
        self._batch = hdf_file['batches'][batch_name]
        ids = self._batch['id']
        self._spectra = None
        self._hdf_file = hdf_file

        # Format for id is <truth_name>_hp<hp_id>_<count>
        # where <count> is the index of this source within the healpixel
        first_id = ids[0].decode()
        cmps = first_id.split('_')
        self._min_hp_ix = int(cmps[-1])
        # last_id = self._ids[-1].decode()
        # cmps = last_id.split('_')
        # self._max_hp_ix = int(cmps[-1])
        self._max_hp_ix = self._min_hp_ix + len(ids) - 1

    def get_sed_row(self, id):
        cmps = id.split('_')
        id_ix = int(cmps[-1])
        if (id_ix < self._min_hp_ix) or (id_ix > self._max_hp_ix):
            return None         # Not in this batch

        if not self._spectra:
            self._spectra = self._batch['spectra']

        return self._spectra[id_ix - self._min_hp_ix]


class _TrilegalSEDReader:
    def __init__(self, hp, object_type_config, cat_dir):
        tmpl = object_type_config.get('sed_file_template', None)
        self._hp = hp
        if not tmpl:
            raise Exception('Bad config file for Trilegal')
        self._fpath = os.path.join(cat_dir,
                                   tmpl.replace('(?P<healpix>\\d+)', str(hp)))
        self._file = h5py.File(self._fpath)
        self._batches = []
        self._wl = None

        # init the batchs
        for b in h5py.File['batches']:
            self._batches.append(_SEDBatch(b, self._file))

    def get_sed_row(self, id):
        cmps = id.split('_')
        if cmps[-2] != f'hp{self._hp}':
            raise Exception(f'This source is not in hp {self._hp}')
        for b in self._batches:
            sed_row = b.get_sed_row(id)
            if sed_row:
                return sed_row
        # Should have warning log message here
        print(f'No SED for id {id}')
        return None

    def get_wl_axis(self):
        if not self._wl:
            self._wl = self._File['wl_axis']


class TrilegalConfigFragment(BaseConfigFragment):
    def __init__(self, prov, area_partition=None, data_file_type=None):
        super().__init__(prov, object_type_name='trilegal',
                         area_partition=area_partition,
                         data_file_type=data_file_type)
