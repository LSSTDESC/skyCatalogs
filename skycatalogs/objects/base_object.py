from collections.abc import Sequence, Iterable
from collections import namedtuple
from packaging import version
import os
import logging
import numpy as np
import galsim
from galsim.roman import longwave_bands as roman_longwave_bands
from galsim.roman import shortwave_bands as roman_shortwave_bands
from galsim.roman import getBandpasses as roman_getBandpasses

from skycatalogs.utils.translate_utils import form_object_string
from skycatalogs.utils.config_utils import Config

'''
Main object types.   There are also may be subtypes. For example,
there could be two subtypes for bulge components, differing in the
form of their associated SEDs
'''

__all__ = ['BaseObject', 'ObjectCollection', 'ObjectList',
           'LSST_BANDS', 'ROMAN_BANDS',
           'load_lsst_bandpasses', 'load_roman_bandpasses']

LSST_BANDS = ('ugrizy')
ROMAN_BANDS = roman_shortwave_bands+roman_longwave_bands

# global for easy access for code run within mp


def _load_lsst_bandpasses():
    '''
    Read in lsst bandpasses from standard place, trim, and store in global dict
    Returns: The dict and throughputs version
    '''
    global lsst_bandpasses
    lsst_bandpasses = dict()
    rubin_sim_dir = os.getenv('RUBIN_SIM_DATA_DIR', None)
    bp_dir = None
    if rubin_sim_dir:
        bp_dir = os.path.join(rubin_sim_dir, 'throughputs', 'baseline')

        if os.path.exists(bp_dir):
            BaseObject._bp_path = bp_dir
            # logger.info(f'Using rubin sim dir {rubin_sim_dir}')
        else:
            bp_dir = os.path.join(os.getenv('HOME'), 'rubin_sim_data',
                                  'throughputs', 'baseline')
            # logger.info(f'Using rubin sim dir rubin_sim_data under HOME')
            if os.path.exists(bp_dir):
                BaseObject._bp_path = bp_dir
            else:
                logger = logging.getLogger('skyCatalogs:load_lsst_bandpasses')
                logger.warning('Using galsim built-in bandpasses which may not be up to date')
                bp_dir = None
    if bp_dir:
        with open(os.path.join(bp_dir, 'version_info')) as f:
            version = f.read().strip()
    else:
        version = 'galsim_builtin'


    for band in LSST_BANDS:
        if bp_dir:
            bp_full_path = os.path.join(bp_dir, f'total_{band}.dat')
        else:
            bp_full_path = f'LSST_{band}.dat'
        bp = galsim.Bandpass(bp_full_path, 'nm')

        # Mirror behavior in imsim.RubinBandpass:
        # https://github.com/LSSTDESC/imSim/blob/main/imsim/bandpass.py#L9
        # Trim the edges to avoid 1.e-4 values out to very high and low
        # wavelengths.
        bp = bp.truncate(relative_throughput=1.e-3)
        # Remove wavelength values selectively for improved speed but
        # preserve flux integrals.
        bp = bp.thin()
        bp = bp.withZeropoint('AB')
        lsst_bandpasses[band] = bp

    return lsst_bandpasses, version

def load_lsst_bandpasses():
    '''
    Read in lsst bandpasses from standard place, trim, and store in global dict
    Returns
    -------
    The bandpasses
    '''
    return _load_lsst_bandpasses()[0]

def _load_roman_bandpasses(**kwargs):
    '''
    Read in Roman bandpasses from standard place, trim, and store in global dict
    Returns: The dict and version inforation
    '''
    global roman_bandpasses
    if version.parse(galsim.version) < version.parse('2.6.0'):
        kwargs.pop("include_all_bands", None)
    roman_bandpasses = roman_getBandpasses(**kwargs)
    return roman_bandpasses, 'galsim_builtin'

def load_roman_bandpasses(**kwargs):
    '''
    Read in Roman bandpasses from standard place, trim, and store in global dict
    Returns: The dict
    '''
    return _load_roman_bandpasses(**kwargs)[0]


class BaseObject(object):
    '''
    Abstract base class for static (in position coordinates) objects.
    Likely need a variant for SSO.
    '''

    def __init__(self, ra, dec, id, object_type, belongs_to, belongs_index):
        '''
        Save at least minimum info needed a fixed (not SSO) object to
        determine if it's in a region and discover all its other properties.
        Parameters
        ra, dec         float         in degrees
        id              string        int types are cast as str
        object_type     string        e.g. 'galaxy', 'star', ...
                                      must appear in catalog config file
        belongs_to      ObjectCollection  collection this object belongs to
        belongs_index   int           index of object within its collection
        '''
        self._ra = ra
        self._dec = dec
        self._id = str(id)
        self._object_type = object_type
        self._belongs_to = belongs_to
        self._belongs_index = belongs_index
        self._logger = belongs_to._sky_catalog._logger

    @property
    def ra(self):
        return self._ra

    @property
    def dec(self):
        return self._dec

    @property
    def id(self):
        return self._id

    @property
    def object_type(self):
        return self._object_type

    @property
    def partition_id(self):
        if self._belongs_to:
            return self._belongs_to.partition_id
        else:
            return None

    @property
    def belongs_to(self):
        return self._belongs_to

    @property
    def native_columns(self):
        '''
        Return set of all columns stored in parquet file for this object.
        May not include quantities which are constant for all objects
        of this type
        '''
        if self._belongs_to:
            return self._belongs_to.native_columns
        else:
            return None

    @property
    def subcomponents(self):
        '''
        Return list of all subcomponents (or None)
        '''
        if self._belongs_to:
            return self._belongs_to.subcomponents
        else:
            return None

    def get_native_attribute(self, attribute_name):
        if not self._belongs_to:
            raise ValueError('To fetch {attribute_name} object must be in an object collection')

        vals = self._belongs_to.get_native_attribute(attribute_name)
        return vals[self._belongs_index]

    def get_sed_metadata(self, **kwargs):
        '''
        E.g. list of wavelength or frequency intervals associated with sed values
        '''
        raise NotImplementedError

    def get_instcat_entry(self, band='r', component=None):
        '''
        Return the string corresponding to instance catalog line
        Parameters:
            band       One of ['u', 'g', 'r', 'i', 'z', 'y']
            component  Required iff the object has subcomponents (i.e.,
                       object type is 'galaxy')
        Returns: A string formatted like a line in an instance catalog
        '''
        return form_object_string(self, band, component)

    def _get_sed(self, component=None, resolution=None, mjd=None):
        '''
        Return sed and mag_norm for a galaxy component or for a star
        Parameters
        ----------
        component    one of 'bulge', 'disk', 'knots' for now. Other components
                     may be supported.  Ignored for stars
        resolution   desired resolution of lambda in nanometers. Ignored
                     for stars.
        mjd          ignored for static objects

        Returns
        -------
        galsim.SED object

        Must be implemented by subclass
        '''
        raise NotImplementedError('Must be implemented by BaseObject subclass if needed')

    def write_sed(self, sed_file_path, component=None, resolution=None,
                  mjd=None):
        sed = self._get_sed(component=component, resolution=None, mjd=None)

        wl = sed.wave_list
        flambda = [sed(w) for w in wl]

        with open(sed_file_path, 'w') as sed_file:
            for lam, flam in zip(wl, flambda):
                sed_file.write(f'{lam} {flam}\n')

    def _get_dust(self):
        """Return the Av, Rv parameters for internal and Milky Way extinction."""
        internal_av = 0
        internal_rv = 1.
        # Fails for object types without these native attributes
        galactic_av = self.get_native_attribute('MW_av')
        galactic_rv = self.get_native_attribute('MW_rv')
        return internal_av, internal_rv, galactic_av, galactic_rv

    def _apply_component_extinction(self, sed):
        # Apply Milky Way extinction.

        iAv, iRv, mwAv, mwRv = self._get_dust()

        sky_cat = self._belongs_to._sky_catalog
        sed = sky_cat.extinguisher.extinguish(sed, mwAv)

        return sed

    def _get_sed_from_file(self, fpath, redshift=0):
        sed = galsim.SED(fpath, wave_type='nm', flux_type='flambda')
        if redshift > 0:
            sed = sed.atRedshift(redshift)
        return sed

    def get_gsobject_components(self, gsparams=None, rng=None, exposure=None):
        """
        Return a dictionary of the GSObject components for the
        sky catalogs object, keyed by component name.
        Must be implemented by subclass
        """
        raise NotImplementedError('Subclass must implement get_gsobject_components')

    def get_observer_sed_component(self, component, mjd=None):
        """
        Return the SED for the specified subcomponent of the SkyCatalog
        object, applying internal extinction, redshift, and Milky Way
        extinction.

        For Milky Way extinction, the Fitzpatrick, et al. (2019) (F19)
        model, as implemented in the dust_extinction package, is used.

        The SEDs are computed assuming exposure times of 1 second.

        Must be implemented by subclass
        """

        raise NotImplementedError('get_observer_sed_component must be implemented by subclass')

    def get_observer_sed_components(self, mjd=None):
        """
        Return a dictionary of the SEDs, keyed by component name.
        """
        sed_components = {}
        subcomponents = [None] if not self.subcomponents \
            else self.subcomponents
        for component in subcomponents:
            sed = self.get_observer_sed_component(component, mjd=mjd)
            if sed is not None:
                sed_components[component] = sed
        return sed_components

    def get_total_observer_sed(self, mjd=None):
        """
        Return the SED summed over SEDs for the individual SkyCatalog
        components.
        """
        sed = None
        for sed_cmp in self.get_observer_sed_components(mjd=mjd).values():
            if sed is None:
                sed = sed_cmp
            else:
                sed += sed_cmp

        return sed

    def get_flux(self, bandpass, sed=None, mjd=None):
        """
        Return the total object flux integrated over the bandpass
        in photons/sec/cm^2
        Use supplied sed if there is one
        """

        if not sed:
            sed = self.get_total_observer_sed(mjd=mjd)

        if sed is None:
            return 0.0

        flux = sed.calculateFlux(bandpass)

        return flux

    def get_fluxes(self, bandpasses, mjd=None):
        # To avoid recomputing sed
        sed = self.get_total_observer_sed(mjd=mjd)
        if sed is None:
            return [0.0 for b in bandpasses]

        return [sed.calculateFlux(b) for b in bandpasses]

    def get_LSST_flux(self, band, sed=None, cache=True, mjd=None):
        if band not in LSST_BANDS:
            return None
        att = f'lsst_flux_{band}'

        # Check if it's already an attribute
        val = getattr(self, att, None)
        if val is not None:
            return val

        if att in self.native_columns:
            return self.get_native_attribute(att)

        val = self.get_flux(lsst_bandpasses[band], sed=sed, mjd=mjd)

        if cache:
            setattr(self, att, val)
        return val

    def get_LSST_fluxes(self, cache=True, as_dict=True, mjd=None):
        '''
        Return a dict of fluxes for LSST bandpasses or, if as_dict is False,
        just a list in the same order as LSST_BANDS
        '''
        fluxes = dict()
        sed = self.get_total_observer_sed(mjd=mjd)
        if sed is None:
            for band in LSST_BANDS:
                fluxes[band] = 0.0
        else:
            for band in LSST_BANDS:
                fluxes[band] = self.get_LSST_flux(band, sed=sed,
                                                  cache=cache, mjd=mjd)
        if as_dict:
            return fluxes
        else:
            return list(fluxes.values())

    def get_roman_flux(self, band, sed=None, cache=True, mjd=None):
        if band not in ROMAN_BANDS:
            return None
        att = f'roman_flux_{band}'

        # Check if it's already an attribute
        val = getattr(self, att, None)
        if val is not None:
            return val

        if att in self.native_columns:
            return self.get_native_attribute(att)

        val = self.get_flux(roman_bandpasses[band], sed=sed, mjd=mjd)

        if cache:
            setattr(self, att, val)
        return val

    def get_roman_fluxes(self, cache=True, as_dict=True, mjd=None):
        '''
        Return a dict of fluxes for Roman bandpasses or, if as_dict is False,
        just a list in the same order as ROMAN_BANDS
        '''
        fluxes = dict()
        sed = self.get_total_observer_sed(mjd=mjd)

        if sed is None:
            for band in ROMAN_BANDS:
                fluxes[band] = 0.0
        else:
            for band in ROMAN_BANDS:
                fluxes[band] = self.get_roman_flux(band, sed=sed,
                                                   cache=cache, mjd=mjd)
        if as_dict:
            return fluxes
        else:
            return list(fluxes.values())


class ObjectCollection(Sequence):
    '''
    Class for collection of fixed (not SSO) objects coming from the same
    source (e.g., file for particular healpix)
    Many of the methods look the same as for BaseObject but they return arrays
    rather than a single number.  There are some additional methods
    '''
    def __init__(self, ra, dec, id, object_type, partition_id, sky_catalog,
                 region=None, mjd=None, mask=None, readers=None, row_group=0):
        '''
        Parameters
        ra, dec      float, array-like of same length
        id           array-like, same length as ra and dec.   int or string
        object_type  single string or (if contained objects may be of
                     different types) array-like
        parition_id  int (e.g. healpix id)  if objects are partitioned by
                     location; else None
        sky_catalog  Instance of SkyCatalog class, typically obtained by
                     calling open_catalog
        region       maybe be used to determine which objects are in the
                     collection
        mjd          MJD of the observation epoch.  This is used by
                     time-varying objects, e.g., SNe, stars.
        mask         indices to be masked off, e.g. in case not all objects
                     in the partition (such as healpixel) are in the region
        readers      may be used to recover properties for objects in this
                     collection, one reader per relevant file
        row_group    used in case backing files are in parquet format

        '''
        self._ra = np.array(ra)
        self._dec = np.array(dec)
        self._id = np.array(id)
        self._rdrs = readers
        self._partition_id = partition_id
        self._sky_catalog = sky_catalog
        self._config = Config(sky_catalog.raw_config)
        self._mask = mask
        self._row_group = row_group

        self._object_class = sky_catalog.cat_cxt.lookup_source_type(object_type)

        # Maybe the following is silly and object_type should always be stored
        # as arrays
        if isinstance(object_type, list):
            self._object_types = np.array(object_type)
            self._object_type_unique = None
            self._uniform_object_type = False
        else:
            self._object_types = None
            self._object_type_unique = object_type
            self._uniform_object_type = True

        self._region = region
        self._mjd = mjd

    @property
    def mjd(self):
        return self._mjd

    @property
    def partition_id(self):
        return self._partition_id

    @property
    def config(self):
        return self._config

    @property
    def sky_catalog(self):
        return self._sky_catalog

    @property
    def native_columns(self):
        '''
        Return set of all columns stored for objects in this collection
        May not include quantities which are constant for all objects
        of this type
        '''
        columns = set()
        for rdr in self._rdrs:
            columns = columns.union(rdr.columns)
        return columns

    @property
    def subcomponents(self):
        '''
        Return list of all subcomponents for objects in this collection.
        Only galaxies have true subcomponents
        '''
        subs = []
        if self._object_type_unique == 'galaxy':
            for s in ['bulge', 'disk', 'knots']:
                if f'sed_val_{s}' in self.native_columns:
                    subs.append(s)
        elif self._object_type_unique == 'diffsky_galaxy':
            subs = ['bulge', 'disk', 'knots']
        else:
            return ['this_object']
        return subs

    def get_native_attribute(self, attribute_name, no_np=False):
        '''
        Retrieve a particular attribute for a collection
        If we already have it, just return it.  Otherwise attempt
        to fetch.   Reader should check whether the attribute actually
        exists.
        '''
        val = getattr(self, attribute_name, None)
        if val is not None:
            return val

        for r in self._rdrs:
            if attribute_name in r.columns:
                d = r.read_columns([attribute_name], self._mask,
                                   row_group=self._row_group, no_np=no_np)
                if not d:
                    return None
                val = d[attribute_name]
                setattr(self, attribute_name, val)
                return val

    def get_native_attributes(self, attribute_list):
        '''
        Return requested attributes as dict. Keys are column names.
        Use our mask if we have one
        '''
        remaining = set(attribute_list)
        final_d = dict()
        for r in self._rdrs:
            to_read = remaining.intersection(r.columns)
            d = r.read_columns(to_read, self._mask, row_group=self._row_group)
            remaining.difference_update(to_read)
            for k in d:
                final_d[k] = d[k]

        if len(remaining) > 0:
            # raise exception?   log error?
            print(f'Unknown column or columns {remaining}')
            return None

        return final_d

    def get_native_attributes_iterator(self, attribute_names):
        '''
        Return iterator for  list of attributes for a collection.  Most of
        the work probably happens in the Parquet reader

        Parameters
        ----------
        attribute_names  list of attribute names
        row_group

        Returns
        -------
        iterator which returns df for a chunk of values of the attributes
        '''
        pass        # for now

    # implement Sequence methods
    def __contains__(self, obj):
        '''
        Parameters
        ----------
        obj can be an (object id) or of type BaseObject
        '''
        if isinstance(obj, int):
            id = obj
        else:
            if isinstance(obj, BaseObject):
                id = obj.id
            else:
                raise TypeError
        return id in self._id

    def __len__(self):
        if self._id is not None:
            return len(self._id)
        else:
            return 0

    # def __iter__(self):   Standard impl based on __getitem__ should be ok
    # def __reversed__(self):   Default implementation ok

    def __getitem__(self, key):
        '''
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
                                      self._id[key], object_type, self, key)

        elif type(key) == slice:
            if key.start is None:
                key.start = 0
            return [self.__getitem__(i) for i in range(self.__len__())[key]]

        elif type(key) == tuple and isinstance(key[0], Iterable):
            return [self.__getitem__(i) for i in key[0]]

    def get_partition_id(self):
        return self._partition_id

    def count(self, obj):
        '''
        returns # of occurrences of obj.  It can only be 0 or 1
        '''
        if self.__contains__(obj):
            return 1
        return 0

    def index(self, obj):
        '''
        Use object id to find the index
        '''
        return self._id.index(obj.id)


LocatedCollection = namedtuple('LocatedCollection',
                               ['collection', 'first_index', 'upper_bound'])
'''
Describes collection included in a list of collections to be concatenated
   collection:     the data
   first_index: index in the concatenated list of first elt in this collection
   upper_bound: index + 1 in concatenated list of last elt in this collection
      so upper_bound for one collection = first_index in the next
'''


class ObjectList(Sequence):
    '''
    Keep track of a list of ObjectCollection objects, but from user
    perspective appears to be a list of objects.  Make use of ObjectCollection
    implementation of Sequence until it's time to move to the next
    ObjectCollection
    '''

    def __init__(self):
        # Elements of _located are named tuples:
        #   (collection, first_index, upper_bound)
        self._located = []
        self._total_len = 0

    def append_collection(self, coll):
        old = self._total_len
        if isinstance(coll, ObjectCollection):
            self._total_len += len(coll)
            self._located.append(LocatedCollection(coll, old, self._total_len))
        else:  # list of collections
            for c in coll:
                self._total_len += len(c)
                self._located.append(LocatedCollection(c, old, self._total_len))
                old = self._total_len

    def append_object_list(self, object_list):
        for e in object_list._located:
            self.append_collection(e.collection)

    def get_native_attribute(self, attribute_name):
        '''
        Retrieve a particular attribute for a source.
        The work is delegated to each of the constituent collections
        '''
        val = self._located[0].collection.get_native_attribute(attribute_name)
        for c in self._located[1:]:
            val = np.append(val,
                            c.collection.get_native_attribute(attribute_name))

        return val

    @property
    def collection_count(self):
        return len(self._located)

    def get_collections(self):
        '''
        return constituent ObjectCollection objects in a list
        '''
        collections = []
        for e in self._located:
            collections.append(e.collection)

        return collections

    # implement Sequence methods

    def __contains__(self, obj):
        '''
        Parameters
        ----------
        obj can be an (object id) or of type BaseObject
        '''

        for e in self._located:
            if obj in e.collection:
                return True

        return False

    def __len__(self):
        return self._total_len

    # def __iter__(self):   Standard impl based on __getitem__ should be ok??
    # def __reversed__(self):   Default implementation ok??

    def __getitem__(self, key):
        '''
        Parameters
        ----------
        If key is an int return a single base object
        If key is a slice return a list of object
        '''
        one_only = isinstance(key, int) or isinstance(key, np.int64)
        is_slice = isinstance(key, slice)
        is_list = isinstance(key, tuple) and isinstance(key[0], Iterable)

        if one_only:
            start = key
        elif is_slice:
            start = key.start
        elif is_list:
            start_ix = 0
            key_list = key[0]
            start = key_list[start_ix]
        else:
            raise ValueError(f'Unknown key type {type(key)}')

        to_return = []
        for e in self._located:
            if start >= e.first_index and start < e.upper_bound:
                rel_first_index = start - e.first_index
                if one_only:
                    return e.collection[rel_first_index]
                if is_slice:
                    if key.stop < e.upper_bound:
                        rel_stop_ix = key.stop - e.first_index
                        to_return += e.collection[slice(rel_first_index,
                                                        rel_stop_ix)]
                        return to_return
                    else:
                        to_return += e.collection[slice(rel_first_index, len(e.collection))]
                    start = e.upper_bound
                if is_list:
                    sub = [elem - e.first_index for elem in key_list if elem >= e.first_index and elem < e.upper_bound]
                    to_return += e.collection[(sub,)]

                    start_ix += len(sub)
                    if start_ix >= len(key_list):
                        break
                    start = key_list[start_ix]

        if len(to_return):
            return to_return
        else:
            raise IndexError
