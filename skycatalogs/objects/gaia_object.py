import os
import warnings
from functools import wraps
from collections.abc import Iterable
from pathlib import PurePath
import numpy as np
import numpy.ma as ma
import pandas as pd
import erfa
import astropy.modeling
from astropy.io import fits
from astropy import units as u
import galsim

from lsst.sphgeom import HtmPixelization

from .base_object import BaseObject, ObjectCollection
from .base_config_fragment import BaseConfigFragment


__all__ = ['GaiaObject', 'GaiaCollection', 'GaiaConfigFragment']


def ignore_erfa_warnings(func):
    @wraps(func)
    def call_func(*args, **kwargs):
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', 'ERFA', erfa.ErfaWarning)
            return func(*args, **kwargs)
    return call_func


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
            complete object collection
        parent_collection: ObjectCollection
            Parent collection of this object.
        index: int
            Index of this object in the parent collection
        """
        # ra = np.degrees(obj_pars['coord_ra'])
        # dec = np.degrees(obj_pars['coord_dec'])
        ra = obj_pars['ra_deg']
        dec = obj_pars['dec_deg']
        # Form the object id from the GAIA catalog id with a string
        # like 'gaia_dr2_' prepended.
        id_prefix = GaiaCollection._id_prefix
        obj_id = f"{id_prefix}{obj_pars['id']}"
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
            except galsim.errors.GalSimRangeError:
                self.stellar_temp = None
            except RuntimeError:
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
            lut = galsim.LookupTable(self._wavelengths, flambda,
                                     interpolant='linear')
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


def _read_fits(htm_id, gaia_config, sky_root, out_dict, logger, region=None):
    '''
    Read data for columns in keys from fits file belonging to htm_id.
    Append to arrays in out_dict.  If region is not None, filter out
    entries not in region before appending.
    Note ra, dec must be in degrees for filtering, but
    coord_ra and coord_dec as stored in out_dict are in radians.

    Parameters
    ----------
    htm_id         int
    gaia_config    dict     gaia_star portion of skyCatalg config
    sky_root       string   absolute skycatalog_root path.  Data
                            is found relative to this
    out_dict       dict     out_dict.keys() are the columns to store
    region         lsst.sphgeom.Region or None
    '''
    f_dir = gaia_config['data_dir']
    f_name = gaia_config['basename_template'].replace('(?P<htm>\\d+)',
                                                      f'{htm_id}')
    if os.path.isabs(f_dir):
        f_path = os.path.join(f_dir, f_name)
    else:
        f_path = os.path.join(sky_root, f_dir, f_name)
    if not os.path.isfile(f_path):
        logger.info(f'No file for requested htm id {htm_id}: {f_path}')
        return

    tbl = fits.open(f_path)[1].data
    if region is None:
        for k in out_dict.keys():
            out_dict[k] += list(tbl[k])
        return

    # Otherwise start with ra, dec and create mask
    ra_full = tbl['coord_ra']
    dec_full = tbl['coord_dec']

    ra_full_deg = np.degrees(ra_full)
    dec_full_deg = np.degrees(dec_full)

    mask = region.compute_mask(ra_full_deg, dec_full_deg)

    if all(mask):
        return

    if any(mask):
        ra_compress = ma.array(ra_full, mask=mask).compressed()
        dec_compress = ma.array(dec_full, mask=mask).compressed()
    else:
        ra_compress = ra_full
        dec_compress = dec_full

    out_dict['coord_ra'] += list(ra_compress)
    out_dict['coord_dec'] += list(dec_compress)

    for k in out_dict.keys():
        if k in ('coord_ra', 'coord_dec'):
            continue
        full = tbl[k]
        if any(mask):
            out_dict[k] += list(ma.array(full, mask=mask).compressed())
        else:
            out_dict[k] += list(full)


class GaiaCollection(ObjectCollection):
    # Class methods
    _gaia_config = None

    @classmethod
    def set_config(cls, config):
        GaiaCollection._gaia_config = config
        if 'id_prefix' in config.keys():
            GaiaCollection._id_prefix = config['id_prefix']
        else:
            GaiaCollection._id_prefix = '_gaia'

    @classmethod
    def get_config(cls):
        return GaiaCollection._gaia_config

    @ignore_erfa_warnings
    @staticmethod
    def load_collection(region, skycatalog, mjd=None, exposure=None, object_type='gaia_star'):
        '''
        region      One of Disk, PolygonalRegion from skyCatalogs.utils.shapes.
                    Box is not currently supported
        skycatalog  An instance of the SkyCatalog class
        mjd         Time at which objects are to be assembled. Ignored for
                    Gaia stars
        exposure    exposure length.  Ignored for Gaia stars
        object_type The object type. This is 'gaia_star' by convention.
        '''
        if not skycatalog:
            raise ValueError('GaiaCollection.load_collection: skycatalog cannot be None')
        refcat_region = region.sphgeom_region()
        if GaiaCollection.get_config() is None:
            gaia_section = skycatalog.raw_config['object_types']['gaia_star']
            GaiaCollection.set_config(gaia_section)
        gaia_section = GaiaCollection.get_config()

        sed_method = GaiaCollection.get_config().get('sed_method', 'use_lut')
        use_lut = (sed_method.strip().lower() == 'use_lut')
        if gaia_section['data_file_type'] == 'butler_refcat':
            # ## Next two lines needed only for butler access
            import lsst.daf.butler as daf_butler
            from lsst.meas.algorithms import ReferenceObjectLoader

            butler_params = gaia_section['butler_parameters']
            butler = daf_butler.Butler(butler_params['repo'],
                                       collections=butler_params['collections'])
            refs = set(butler.registry.queryDatasets(butler_params['dstype']))
            refCats = [daf_butler.DeferredDatasetHandle(butler, _, {})
                       for _ in refs]
            dataIds = [butler.registry.expandDataId(_.dataId) for _ in refs]
            config = ReferenceObjectLoader.ConfigClass()
            config.filterMap = {f'{_}': f'phot_{_}_mean' for _ in ('g', 'bp',
                                                                   'rp')}
            ref_obj_loader = ReferenceObjectLoader(dataIds=dataIds,
                                                   refCats=refCats,
                                                   config=config)

            band = 'bp'
            cat = ref_obj_loader.loadRegion(refcat_region, band).refCat
            df = cat.asAstropy().to_pandas().sort_values('id')

        else:    # access fits files directly
            # find htms intersecting region
            level = gaia_section['area_partition']['level']
            htm_pix = HtmPixelization(level)

            overlap_ranges = htm_pix.envelope(refcat_region)
            interior_ranges = htm_pix.interior(refcat_region)
            partial_ranges = overlap_ranges - interior_ranges
            keys = ['id', 'coord_ra', 'coord_dec', 'parallax', 'pm_ra',
                    'pm_dec', 'epoch']
            keys += [f'phot_{_}_mean_flux' for _ in ('g', 'bp', 'rp')]
            out_dict = {k: [] for k in keys}

            config = GaiaCollection.get_config()
            for i_r in interior_ranges:
                for i_htm in i_r:
                    _read_fits(i_htm, config, skycatalog._sky_root, out_dict,
                               skycatalog._logger, region=None)
            for p_r in partial_ranges:
                for i_htm in p_r:
                    _read_fits(i_htm, config, skycatalog._sky_root,
                               out_dict, skycatalog._logger, region=region)
            df = pd.DataFrame(out_dict).sort_values('id')

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
        pm_outputs = erfa.pmsafe(np.array(df['coord_ra']),
                                 np.array(df['coord_dec']),
                                 np.array(df['pm_ra']), np.array(df['pm_dec']),
                                 np.array(df['parallax']),
                                 rv1, epNa, np.array(df['epoch']), epNa, ep2b)
        # Update ra, dec values.
        df['coord_ra'] = pm_outputs[0]
        df['coord_dec'] = pm_outputs[1]

        # and also store degrees version
        # is there any reason to keep coord_ra & coord_dec?
        df['ra_deg'] = np.degrees(pm_outputs[0])
        df['dec_deg'] = np.degrees(pm_outputs[1])

        return GaiaCollection(df, skycatalog, object_type, use_lut, mjd)

    def __init__(self, df, sky_catalog, object_type, use_lut, mjd):
        self.df = df
        self._sky_catalog = sky_catalog
        self._partition_id = None
        self._mask = None
        self._object_type_unique = object_type
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
        cols = ('id', 'ra_deg', 'dec_deg', 'phot_bp_mean_flux',
                'phot_rp_mean_flux')
        if isinstance(key, int) or isinstance(key, np.int64):
            row = {col: self.df[col][key] for col in cols}
            return GaiaObject(row, self, key)

        elif isinstance(key, slice):
            # ixdata = [i for i in range(min(key.stop, len(self.df['id'])))]
            # ixes = itertools.islice(ixdata, key.start, key.stop, key.step)
            # return [self.__getitem__(i) for i in ixes]
            return [self.__getitem__(i) for i in range(len(self.df))[key]]

        elif type(key) == tuple and isinstance(key[0], Iterable):
            return [self.__getitem__(i) for i in key[0]]

    def __len__(self):
        return len(self.df)


class GaiaConfigFragment(BaseConfigFragment):
    def __init__(self, prov, id_prefix=None, use_butler=False,
                 butler_parameters=None, area_partition=None,
                 data_file_type=None,
                 data_dir=None, basename_template=None):
        '''
        prov              dict       Provenance
        use_butler        boolean    Access data via Butler if True; else read
                                     FITS files directly
        id_prefix         string     Prepend to native (integer) id for
                                     each object
        butler_parameters dict       Used only if use_butler is true
        area_partition    dict       Used only if use_butler is False
        data_file_type    dict       Used only if use_butler is False
        data_dir          string     Used only if use_butler is False
        basename_template string     Used only if use_butler is False

        Parameters starting with id_prefix are optional.  They allow
        the caller to override the default information in the Gaia star
        fragment template
        '''
        if use_butler:
            template_name = 'gaia_star_butler_template.yaml'
        else:
            template_name = 'gaia_star_direct_template.yaml'

        super().__init__(prov, object_type_name='gaia_star',
                         template_name=template_name,
                         area_partition=area_partition,
                         data_file_type=data_file_type)

        if id_prefix:
            self._opt_dict['id_prefix'] = id_prefix
        if use_butler:
            if butler_parameters:
                self._opt_dict['butler_parameters'] = butler_parameters
        else:
            if basename_template:
                self._opt_dict['basename_template'] = basename_template
