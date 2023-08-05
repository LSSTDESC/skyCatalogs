import os
import sys
import warnings
import itertools
from pathlib import PurePath
import numpy as np
import erfa
import astropy.modeling
from astropy import units as u
import galsim
import lsst.daf.butler as daf_butler
import lsst.geom
from lsst.meas.algorithms import ReferenceObjectLoader
from skycatalogs.utils.shapes import Disk, PolygonalRegion
from skycatalogs.objects.base_object import BaseObject, ObjectCollection


__all__ = ['GaiaObject', 'GaiaCollection']


_FILE_PATH = str(PurePath(__file__))
_SKYCATALOGS_DIR = _FILE_PATH[:_FILE_PATH.rindex('/skycatalogs')]


# Read in function for stellar temperature estimation given bp and rp
# fluxes
_TEMP_EST_DATA_FILE = os.path.join(_SKYCATALOGS_DIR, 'skycatalogs',
                                   'data', 'gaia_dr2',
                                   'gaia_dr2_temp_from_bp-rp_ratio.txt')
_TEMP_FUNC = galsim.LookupTable(*np.genfromtxt(_TEMP_EST_DATA_FILE))

# Read in the Gaia DR2 passband for the bp flux to use for setting the
# SED normalizations.
_GAIA_PB_FILE = os.path.join(_SKYCATALOGS_DIR, 'skycatalogs',
                             'data', 'gaia_dr2',
                             'GaiaDR2_RevisedPassbands.dat')
_GAIA_PBS = np.genfromtxt(_GAIA_PB_FILE, names=['wl', 'g', 'g_err', 'bp',
                                                'bp_err', 'rp', 'rp_err'])
_index = np.where(_GAIA_PBS['bp'] < 99)  # Omit sentinel values.
# Create the Bandpass object for the bp passband, using the 'AB' zero point.
_GAIA_BP = galsim.Bandpass(
    galsim.LookupTable(_GAIA_PBS['wl'][_index], _GAIA_PBS['bp'][_index]),
    wave_type='nm').withZeropoint('AB').thin()


class GaiaObject(BaseObject):
    """
    Object class for Gaia DR2 data read from the LSST reference catalog files.
    """
    _stellar_temperature = _TEMP_FUNC
    _gaia_bp_bandpass = _GAIA_BP
    _wavelengths = np.arange(250, 1250, 5, dtype=float)
    def __init__(self, obj_pars, parent_collection, index):
        """
        Parameters
        ----------
        obj_pars: dict-like
            Dictionary of object parameters as read directly from the
            reference catalog.
        parent_collection: ObjectCollection
            Parent collection of this object.
        index: int
            Index of this object in the parent collection
        use_lut: bool [True]
            Flag to use a galsim LookupTable for representing the SED.
        """
        ra = np.degrees(obj_pars['coord_ra'])
        dec = np.degrees(obj_pars['coord_dec'])
        # Form the object id from the GAIA catalog id with the string
        # 'gaia_dr2_' prepended.
        obj_id = f"gaia_dr2_{obj_pars['id']}"
        super().__init__(ra, dec, obj_id, 'gaia_star',
                         belongs_to=parent_collection, belongs_index=index)
        self.use_lut = self._belongs_to._use_lut
        bp_flux = obj_pars['phot_bp_mean_flux']
        rp_flux = obj_pars['phot_rp_mean_flux']
        if rp_flux == 0.0 or bp_flux == 0.0:
            self.stellar_temp = None
        else:
            try:
                self.stellar_temp = self._stellar_temperature(bp_flux/rp_flux)
                # Convert from flux units of nJy to AB mag for the bp passband,
                # which we will use to normalize the SED.
                self.bp_mag = -2.5*np.log10(bp_flux*1e-9) + 8.90
            except galsim.errors.GalSimRangeError as ex:
                self.stellar_temp = None
            except RuntimeError as rex:
                self.stellar_temp = None

    def blambda(self, wl):
        clight = astropy.constants.c.value*1e2  # cm/s
        nu = clight/(wl*1e-7)  # Hz
        Bnu = astropy.modeling.physical_models.BlackBody(
            temperature=self.stellar_temp*u.K)
        return Bnu(nu*u.Hz).value*nu**2/clight/1e7  # erg/nm/cm^2/s

    def get_observer_sed_component(self, component, mjd=None):
        if component != 'this_object':
            raise RuntimeError("Unknown SED component: %s", component)
        if self.stellar_temp is None:
            return None
        if self.use_lut:
            flambda = self.blambda(self._wavelengths)
            lut = galsim.LookupTable(self._wavelengths, flambda)
            sed = galsim.SED(lut, wave_type='nm', flux_type='flambda')
        else:
            sed = galsim.SED(self.blambda, wave_type='nm', flux_type='flambda')
        return sed.withMagnitude(self.bp_mag, self._gaia_bp_bandpass)

    def get_gsobject_components(self, gsparams=None, rng=None):
        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)
        return {'this_object': galsim.DeltaFunction(gsparams=gsparams)}

    def set_use_lut(self, use_lut):
        self.use_lut = use_lut

class GaiaCollection(ObjectCollection):
    # Class methods
    _gaia_config = None
    def set_config(config):
        GaiaCollection._gaia_config = config

    def get_config():
        return GaiaCollection._gaia_config

    def load_collection(region, skycatalog, mjd=None):
        if isinstance(region, Disk):
            ra = lsst.geom.Angle(region.ra, lsst.geom.degrees)
            dec = lsst.geom.Angle(region.dec, lsst.geom.degrees)
            center = lsst.geom.SpherePoint(ra, dec)
            radius = lsst.geom.Angle(region.radius_as, lsst.geom.arcseconds)
            refcat_region = lsst.sphgeom.Circle(center.getVector(), radius)
        elif isinstance(region, PolygonalRegion):
            refcat_region = region._convex_polygon
        else:
            raise TypeError(f'GaiaCollection.load_collection: {region} not supported')

        source_type = 'gaia_star'

        if GaiaCollection.get_config() is None:
            gaia_section = skycatalog.raw_config['object_types']['gaia_star']
            GaiaCollection.set_config(gaia_section)

        butler_params = GaiaCollection.get_config()['butler_parameters']
        butler = daf_butler.Butler(butler_params['repo'],
                                   collections=butler_params['collections'])
        refs = set(butler.registry.queryDatasets(butler_params['dstype']))
        refCats = [daf_butler.DeferredDatasetHandle(butler, _, {})
                   for _ in refs]
        dataIds = [butler.registry.expandDataId(_.dataId) for _ in refs]
        config = ReferenceObjectLoader.ConfigClass()
        config.filterMap = {f'{_}': f'phot_{_}_mean' for _ in ('g', 'bp', 'rp')}
        ref_obj_loader = ReferenceObjectLoader(dataIds=dataIds,
                                               refCats=refCats,
                                               config=config)

        sed_method = GaiaCollection.get_config().get('sed_method', 'use_lut')
        use_lut = (sed_method.strip().lower() == 'use_lut')
        band = 'bp'
        cat = ref_obj_loader.loadRegion(refcat_region, band).refCat
        df =  cat.asAstropy().to_pandas()

        if mjd is None:
            raise RuntimeError("MJD needs to be provided for Gaia "
                               "proper motion calculations.")
        # The erfa.pmsafe code
        # (https://pyerfa.readthedocs.io/en/latest/api/erfa.pmsafe.html)
        # expects ra, dec in units of radians and pm_ra, pm_dec
        # (proper motion) in units of radians/year.  These are the
        # units used in the refcats files.  However, erfa.pmsafe
        # expects parallax in arcseconds so convert from the refcats
        # units of radians:
        arcsec_per_radian = (1.0*u.radian).to(u.arcsec).value
        df['parallax'] *= arcsec_per_radian
        # radial velocity is missing from the Gaia DR2 refcats, so pass
        # an array of zeros.
        rv1 = np.zeros(len(df))
        epNa = 2400000.5  # "part A" of starting and ending epochs for MJDs.
        ep2b = mjd
        pm_outputs = erfa.pmsafe(df['coord_ra'], df['coord_dec'],
                                 df['pm_ra'], df['pm_dec'], df['parallax'],
                                 rv1, epNa, df['epoch'], epNa, ep2b)
        # Update ra, dec values.
        df['coord_ra'] = pm_outputs[0]
        df['coord_dec'] = pm_outputs[1]

        return GaiaCollection(df, skycatalog, source_type, use_lut, mjd)

    def __init__(self, df, sky_catalog, source_type, use_lut, mjd):
        self.df = df
        self._sky_catalog = sky_catalog
        self._partition_id = None
        self._id = np.array([f"gaia_dr2_{df.iloc[key]['id']}" for key in range(len(df))])
        self._mask = None
        self._object_type_unique = source_type
        self._use_lut = use_lut
        self._rdrs = []
        self._object_class = GaiaObject
        self._use_lut = use_lut
        self._mjd = mjd

    @property
    def native_columns(self):
        return set()

    @property
    def use_lut(self):
        return self._use_lut

    @property
    def mjd(self):
        return self._mjd

    def __getitem__(self, key):
        if isinstance(key, int) or isinstance(key, np.int64):
            return GaiaObject(self.df.iloc[key], self, key)

        elif type(key) == slice:
            ixdata = [i for i in range(min(key.stop,len(self._id)))]
            ixes = itertools.islice(ixdata, key.start, key.stop, key.step)
            return [self._object_class(self.df.iloc[i], self, i) for i in ixes]

        elif type(key) == tuple and isinstance(key[0], Iterable):
            #  check it's a list of int-like?
            return [self._object_class(self.df.iloc[i], self, i) for i in key[0]]

    def __len__(self):
        return len(self.df)
