import os
import numpy as np
import astropy.modeling
from astropy import units as u
import galsim
import lsst.daf.butler as daf_butler
import lsst.geom
from lsst.meas.algorithms import LoadReferenceObjectsTask, \
    ReferenceObjectLoader
from desc.skycatalogs.skyCatalogs import SkyCatalog, open_catalog, Disk
from desc.skycatalogs.objects.base_object import BaseObject, ObjectCollection
#from desc.skycatalogs.objects.base_object import load_lsst_bandpasses


# Read in function for stellar temperature estimation given bp and rp
# fluxes
_TEMP_EST_DATA_FILE = os.path.join(os.environ['SKYCATALOGS_DIR'], 'data',
                                   'gaia_dr2',
                                   'gaia_dr2_temp_from_bp-rp_ratio.txt')
_TEMP_FUNC = galsim.LookupTable(*np.genfromtxt(_TEMP_EST_DATA_FILE))

# Read in the Gaia DR2 passband for the bp flux to use for setting the
# SED normalizations.
_GAIA_PB_FILE = os.path.join(os.environ['SKYCATALOGS_DIR'], 'data', 'gaia_dr2',
                             'GaiaDR2_RevisedPassbands.dat')
_GAIA_PBS = np.genfromtxt(_GAIA_PB_FILE, names=['wl', 'g', 'g_err', 'bp',
                                                'bp_err', 'rp', 'rp_err'])
_index = np.where(_GAIA_PBS['bp'] < 99)  # Omit sentinel values.
_GAIA_BP = galsim.Bandpass(
    galsim.LookupTable(_GAIA_PBS['wl'][_index], _GAIA_PBS['bp'][_index]),
    wave_type='nm').thin()


class GaiaObject(BaseObject):
    """
    Object class for Gaia DR2 data read from the LSST reference catalog files.
    """
    _stellar_temperature = _TEMP_FUNC
    _gaia_bp_bandpass = _GAIA_BP
    def __init__(self, obj_pars, parent_collection, index):
        """
        Parameters
        ----------
        obj_pars: dict-like
            Dictionary of object parameters as read directly from the
            reference catalog.
        """
        ra = np.degrees(obj_pars['coord_ra'])
        dec = np.degrees(obj_pars['coord_dec'])
        # Form the object id from the GAIA catalog id with the string
        # 'gaia_dr2_' prepended.
        obj_id = f"gaia_dr2_{obj_pars['id']}"
        super().__init__(ra, dec, obj_id, OBJECT_TYPES['gaia_star'],
                         belongs_to=parent_collecton, belongs_index=index)
        self.bp_flux = obj_pars['phot_bp_mean_flux']
        rp_flux = obj_pars['phot_rp_mean_flux']
        self.stellar_temp = self._stellar_temperature(self.bp_flux/rp_flux)

    def blambda(self, wl):
        clight = astropy.constants.c.value*1e2  # cm/s
        nu = clight/(wl*1e-7)  # Hz
        Bnu = astropy.modeling.physical_models.BlackBody(
            temperature=self.stellar_temp*u.K)
        return Bnu(nu*u.Hz).value*nu**2/clight/1e7  # erg/nm/cm^2/s

    def get_observer_sed_component(self, component):
        if component is not None:
            raise RuntimeError("Unknown SED component: %s", component)
        sed = galsim.SED(self.blambda, wave_type='nm', flux_type='flambda')
        return sed.withFlux(self.bp_flux, self._gaia_bp_bandpass)

class GaiaCollection(ObjectCollection):
    # Class methods
    def set_config(config=None):
        if config is None:
            config = dict()
            config['gaia_star'] = dict()
            config['gaia_star']['data_file_type'] = 'butler_refcat'
            config['gaia_star']['area_partition'] = 'None'
            butler_params = {'collections' : 'HSC/defaults',
                             'dstype' : 'gaia_dr2-20200414',
                             'repo' : '/repo/main'}
            config['gaia_star']['butler_parameters'] = butler_parameters
        GaiaCollection._gaia_config = config

    def load_collection(region, skycatalog, config=None):
        if not isinstance(region, Disk):
            raise TypeError('GaiaCollection.load_collection: region must be a Disk')
        if config is None:
            config = GaiaCollection._gaia_config
        butler_params = config['gaia_star']['butler_params']
        butler = daf_butler.Butler(butler_params['repo'],
                                   collections=butler_params['collections'])
        refs = butler.registry.queryDatasets(butler_params['dstype'])
        refCats = [daf_butler.DeferredDatasetHandle(butler, _, {})
                   for _ in refs]
        dataIds = [butler.registry.expandDataId(_.dataId) for _ in refs]
        config = LoadReferenceObjectsTask.ConfigClass()
        config.filterMap = {f'{_}': f'phot_{_}_mean' for _ in ('g', 'bp', 'rp')}
        ref_obj_loader = ReferenceObjectLoader(dataIds=dataIds,
                                               refCats=refCats,
                                               config=config)

        coord = lsst.geom.SpherePoint(lsst.geom.Angle(region.ra,
                                                      lsst.geom.degrees),
                                      lsst.geom.Angle(region.ra,
                                                      lsst.geom.degrees))
        rad = lsst.geom.Angle(region.radius_as, lsst.geom.arcseconds)
        band = 'bp'
        cat = ref_obj_loader.loadSkyCircle(coord, rad, band).refCat
        df =  cat.asAstropy().to_pandas()

        return GaiaCollection(df)

    def __init__(self, df, sky_catalog, region):
        self.df = df
        self._sky_catalog = sky_catalog
        self._partition_id = None
        self._id = np.array([f"gaia_dr2_{df.iloc[key]['id']}" for key in range(len(df))])
        self._mask = None
        self._object_type_unique = 'gaia_star'
        self._rdrs = []


    @property
    def native_columns(self):
        return set()

    def __getitem__(self, key):
        return GaiaObject(self.df.iloc[key], self, key)

    def __len__(self):
        return len(self.df)

if __name__ == '__main__':
    disk = Disk(60, -40, 0.17 * 360)
    GaiaCollection.set_config()
    skycatalog_root = os.getenv('SKYCATALOG_ROOT')
    config_path = os.path.join(skycatalog_root, 'reorg', 'skyCatalog.yaml')

    skycat = open_catalog(config_path, skycatalog_root=skycatalog_root)
    collection = GaiaCollection.load_collection(disk, skycat)
