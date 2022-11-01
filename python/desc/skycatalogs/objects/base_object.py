from collections.abc import Sequence, Iterable
from collections import namedtuple
import os
import gzip
import itertools
import numpy as np
import astropy.units as u
import warnings
from dust_extinction.parameter_averages import F19
import galsim
import logging
from galsim.errors import GalSimRangeError

from desc.skycatalogs.utils.translate_utils import form_object_string
from desc.skycatalogs.utils.config_utils import Config
from desc.skycatalogs.utils.sed_tools import ObservedSedFactory

'''
Main object types.   There are also may be subtypes. For example,
there could be two subtypes for bulge components, differing in the
form of their associated SEDs
'''

__all__ = ['BaseObject', 'ObjectCollection', 'ObjectList', 'OBJECT_TYPES',
           'LSST_BANDS', 'load_lsst_bandpasses', 'CatalogContext']
GALAXY=1
GALAXY_BULGE=2
GALAXY_DISK=3
GALAXY_KNOTS=4
STAR=5
AGN=6
SN=7

OBJECT_TYPES = {'galaxy' : GALAXY, 'bulge_basic' : GALAXY_BULGE,
                'disk_basic' : GALAXY_DISK, 'knots_basic' : GALAXY_KNOTS,
                'star' : STAR, 'agn' : AGN, 'sn' : SN}

LSST_BANDS = ('ugrizy')

# global for easy access for code run within mp

class CatalogContext:
    def __init__(self, the_sky_cat):
        global sky_cat
        self._sky_cat = the_sky_cat
        sky_cat = the_sky_cat

def load_lsst_bandpasses():
    '''
    Read in lsst bandpasses from standard place, trim, and store in global dict
    Returns: The dict

    '''
    global lsst_bandpasses
    lsst_bandpasses = dict()
    rubin_sim_dir = os.getenv('RUBIN_SIM_DATA_DIR', None)
    if rubin_sim_dir:
        bp_dir = os.path.join(rubin_sim_dir, 'throughputs', 'baseline')

        if os.path.exists(bp_dir):
            BaseObject._bp_path = bp_dir
            #logger.info(f'Using rubin sim dir {rubin_sim_dir}')
        else:
            bp_dir = os.path.join(os.getenv('HOME'), 'rubin_sim_data',
                                   'throughputs', 'baseline')
            #logger.info(f'Using rubin sim dir rubin_sim_data under HOME')
            if os.path.exists(bp_dir):
                BaseObject._bp_path = bp_dir
            else:
                #logger.info('Using galsim built-in bandpasses')
                bp_dir = None
                fname_fmt = 'LSST_{band}.dat'
        for band in LSST_BANDS:
            if bp_dir:
                bp_full_path = os.path.join(bp_dir, f'total_{band}.dat')
            else:
                bp_full_path = f'LSST_{band}.dat'
            lsst_bandpasses[band] = galsim.Bandpass(bp_full_path, 'nm').thin()

    return lsst_bandpasses

class BaseObject(object):
    '''
    Abstract base class for static (in position coordinates) objects.
    Likely need a variant for SSO.
    '''

    _bp500 = galsim.Bandpass(galsim.LookupTable([499, 500, 501],[0, 1, 0]),
                             wave_type='nm').withZeropoint('AB')
    def __init__(self, ra, dec, id, object_type, redshift=None,
                 belongs_to=None, belongs_index=None):
        '''
        Minimum information needed for static (not SSO) objects
        ra, dec needed to check for region containment
        belongs_to is object collection, if any, this object belongs to
        '''
        self._ra = ra
        self._dec = dec
        self._id = id
        self._object_type = object_type
        self._redshift = redshift
        self._belongs_to = belongs_to
        self._belongs_index = belongs_index

        # All objects also include redshift information. Also MW extinction,
        # but extinction is by subcomponent for galaxies

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
    def redshift(self):
        if self._redshift:        return self._redshift
        if self._belongs_to:
            self._redshift = self.get_native_attribute('redshift')
        return self._redshift

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

    def get_instcat_entry(self, band = 'r', component=None):
        '''
        Return the string corresponding to instance catalog line
        Parameters:
            band       One of ['u', 'g', 'r', 'i', 'z', 'y']
            component  Required iff the object has subcomponents (i.e.,
                       object type is 'galaxy')
        Returns: A string formatted like a line in an instance catalog
        '''
        if self.object_type == 'galaxy':
            if component not in self.subcomponents:
                return ''
        return form_object_string(self, band, component)

    def get_sed(self, component=None, resolution=None):
        '''
        Return sed and mag_norm for a galaxy component or for a star
        Parameters
        ----------
        component    one of 'bulge', 'disk', 'knots' for now. Other components
                     may be supported.  Ignored for stars
        resolution   desired resolution of lambda in nanometers. Ignored
                     for stars.

        Returns
        -------
        A pair (sed, mag_norm)
        '''

        if self._object_type == 'star':
            # read in the file; return data from it. Arguments are ignored
            if component:
                warnings.warn(f'get_sed: component argument ignored for object type {self._object_type}')
            if resolution:
                warnings.warn(f'get_sed: resolution argument ignored for object type {self._object_type}')
            mag_norm = self.get_native_attribute('magnorm')
            rel_path = self.get_native_attribute('sed_filepath')

            fpath = os.path.join(os.getenv('SIMS_SED_LIBRARY_DIR'),
                                 self.get_native_attribute('sed_filepath'))

            return sky_cat.observed_sed_factory.create_pointsource(rel_path), mag_norm
        if self._object_type != 'galaxy':
            raise ValueError('get_sed function only available for stars and galaxy components')
        if component not in ['disk', 'bulge', 'knots']:
            raise ValueError(f'Cannot fetch SED for component type {component}')

        th_val = self.get_native_attribute(f'sed_val_{component}')
        if  th_val is None:   #  values for this component are not in the file
            raise ValueError(f'{component} not part of this catalog')

        # if values are all zeros or nearly no point in trying to convert
        if max(th_val) < np.finfo('float').resolution:
            return None, 0.0

        z_h = self.get_native_attribute('redshift_hubble')
        z = self.get_native_attribute('redshift')

        sed = sky_cat.observed_sed_factory.create(th_val, z_h, z,
                                                     resolution=resolution)
        magnorm = sky_cat.observed_sed_factory.magnorm(th_val, z_h)

        return sed, magnorm

    def write_sed(self, sed_file_path, component=None, resolution=None):
        sed, _ = self.get_sed(component=component, resolution=None)

        wl = sed.wave_list
        flambda = [sed(w) for w in wl]

        with open(sed_file_path, 'w') as sed_file:
            for lam, flam in zip(wl, flambda):
                sed_file.write(f'{lam} {flam}\n')

    def get_dust(self):
        """Return the Av, Rv parameters for internal and Milky Way extinction."""
        internal_av = 0
        internal_rv = 1.
        galactic_av = self.get_native_attribute('MW_av')
        galactic_rv = self.get_native_attribute('MW_rv')
        return internal_av, internal_rv, galactic_av, galactic_rv

    def get_wl_params(self):
        """Return the weak lensing parameters, g1, g2, mu."""
        gamma1 = self.get_native_attribute('shear_1')
        gamma2 = self.get_native_attribute('shear_2')
        kappa =  self.get_native_attribute('convergence')
        # Compute reduced shears and magnification.
        g1 = gamma1/(1. - kappa)    # real part of reduced shear
        g2 = gamma2/(1. - kappa)    # imaginary part of reduced shear
        mu = 1./((1. - kappa)**2 - (gamma1**2 + gamma2**2)) # magnification
        return g1, g2, mu

    def get_gsobject_components(self, gsparams=None, rng=None):
        """
        Return a dictionary of the GSObject components for the
        sky catalogs object, keyed by component name.
        """
        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)
        if self.object_type == 'star':
            return {None: galsim.DeltaFunction(gsparams=gsparams)}
        if self.object_type != 'galaxy':
            raise RuntimeError("Do not know how to handle object type "
                               f"{self.object_type}")
        obj_dict = {}
        for component in self.subcomponents:
            # knots use the same major/minor axes as the disk component.
            my_component = 'disk' if component != 'bulge' else 'bulge'
            a = self.get_native_attribute(
                f'size_{my_component}_true')
            b = self.get_native_attribute(
                f'size_minor_{my_component}_true')
            assert a >= b
            pa = self.get_native_attribute('position_angle_unlensed')
            beta = float(90 + pa)*galsim.degrees
            hlr = (a*b)**0.5   # approximation for half-light radius
            if component == 'knots':
                npoints = self.get_native_attribute('n_knots')
                assert npoints > 0
                obj = galsim.RandomKnots(npoints=npoints,
                                         half_light_radius=hlr, rng=rng,
                                         gsparams=gsparams)
            else:
                n = self.get_native_attribute(f'sersic_{component}')
                # Quantize the n values at 0.05 so that galsim can
                # possibly amortize sersic calculations from the previous
                # galaxy.
                n = round(n*20.)/20.
                obj = galsim.Sersic(n=n, half_light_radius=hlr,
                                    gsparams=gsparams)
            shear = galsim.Shear(q=b/a, beta=beta)
            obj = obj._shear(shear)
            g1, g2, mu = self.get_wl_params()
            obj_dict[component] = obj._lens(g1, g2, mu)
        return obj_dict

    def get_observer_sed_component(self, component):
        """
        Return the SED for the specified subcomponent of the SkyCatalog
        object, applying internal extinction, redshift, and Milky Way
        extinction.

        For Milky Way extinction, the Fitzpatrick, et al. (2019) (F19)
        model, as implemented in the dust_extinction package, is used.

        The SEDs are computed assuming exposure times of 1 second.
        """
        # Create a galsim.SED for this subcomponent.
        if self._object_type == 'star':   # or other pointsource
            sed, magnorm = self.get_sed(component=component)
            # The SED is in units of photons/nm/cm^2/s
            # -0.9210340371976184 = -np.log(10)/2.5. Use to convert mag to flux
            flux_500 = np.exp(-0.9210340371976184 * magnorm)
            sed = sed.withMagnitude(0, self._bp500)
            sed = sed*flux_500
        else:
            # For galaxy components sed already has correct normalization
            sed, _ = self.get_sed(component=component)
            if sed is None:
                # This subcomponent has zero emission so return None.
                return None


        iAv, iRv, mwAv, mwRv = self.get_dust()
        if iAv > 0:
            # Apply internal extinction model, which is assumed
            # to be the same for all subcomponents.
            pass  #TODO add implementation for internal extinction.

        # redshift already applied by create or create_pointsource

        # Apply Milky Way extinction.
        sed = sky_cat.extinguisher.extinguish(sed, mwAv)
        return sed

    def get_observer_sed_components(self):
        """
        Return a dictionary of the SEDs, keyed by component name.
        """
        sed_components = {}
        subcomponents = [None] if not self.subcomponents \
            else self.subcomponents
        for component in subcomponents:
            sed = self.get_observer_sed_component(component)
            if sed is not None:
                sed_components[component] = sed
        return sed_components

    def get_total_observer_sed(self):
        """
        Return the SED summed over SEDs for the individual SkyCatalog
        components.
        """
        sed = None
        for sed_component in self.get_observer_sed_components().values():
            if sed is None:
                sed = sed_component
            else:
                sed += sed_component
        if 'shear_1' in self.native_columns:
            _, _, mu = self.get_wl_params()
            sed *= mu
        return sed

    def get_flux(self, bandpass, sed=None):
        """
        Return the total object flux integrated over the bandpass
        in photons/sec/cm^2
        Use supplied sed if there is one
        """

        if not sed:
            sed = self.get_total_observer_sed()

        flux = sed.calculateFlux(bandpass)

        return flux

    def get_fluxes(self, bandpasses):
        # To avoid recomputing sed
        sed = self.get_total_observer_sed()
        return [sed.calculateFlux(b) for b in bandpasses]

    def get_LSST_flux(self, band, sed=None, cache=True):
        if not band in LSST_BANDS:
            return None
        att = f'lsst_flux_{band}'
        # Check if it's already an attribute
        val = getattr(self, att, None)
        if val is not None:
            return val

        if att in self.native_columns:
            return self.get_native_attribute(att)

        val = self.get_flux(lsst_bandpasses[band], sed=sed)

        if cache:
            setattr(self, att, val)
        return val

    def get_LSST_fluxes(self, cache=True, as_dict=True):
        '''
        Return a dict of fluxes for LSST bandpasses or, if as_dict is False,
        just a list in the same order as LSST_BANDS
        '''
        fluxes = dict()
        sed = self.get_total_observer_sed()
        for band in LSST_BANDS:
            ##### temporary to debug error, now apparently fixed
            try:
                fluxes[band] = self.get_LSST_flux(band, sed=sed, cache=cache)
            except GalSimRangeError as e:
                logger = logging.getLogger('skyCatalogs.creator')
                galaxy_id = self.get_native_attribute('galaxy_id')
                logger.debug(f'galaxy_id: {galaxy_id}  band: {band}')
                redshift = self.get_native_attribute('redshift')
                logger.debug(f'redshift: {redshift}')
                logger.debug(f'GalSimRangeError: {e}')
                fluxes[band] = 0.0

        if as_dict:
            return fluxes
        else:
            return list(fluxes.values())

class ObjectCollection(Sequence):
    '''
    Class for collection of static objects coming from the same
    source (e.g., file for particular healpix)
    Many of the methods look the same as for BaseObject but they return arrays
    rather than a single number.  There are some additional methods
    '''
    def __init__(self, ra, dec, id, object_type, partition_id, sky_catalog,
                 region=None, mask=None, readers=None):
        '''
        Minimum information needed for static objects.
        specified, has already been used by the caller to generate mask.
        ra, dec, id must be array-like and the same length
        object_type  may be either single value or array-like (could be
        useful if more than one object type is stored in the same file)
        partition_id (e.g. healpix id)
        sky_catalog instance of SkyCatalog class
        (redshift should probably be ditched; no need for it)
        (similarly for region with current code structure. Information
        needed is encoded in mask)

        mask  indices to be masked off, e.g. because region does not
              include the entire healpix pixel
        readers one per file associated with the type(s). (Galaxy info
              is stored in two files)


        '''
        self._ra = np.array(ra)
        self._dec = np.array(dec)
        self._id = np.array(id)
        self._rdrs = readers
        self._partition_id = partition_id
        self._sky_catalog = sky_catalog   # might not need this
        self._config = Config(sky_catalog.raw_config)
        self._mask = mask

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
        Only galaxies have subcomponents
        '''
        subs = []
        if self._object_type_unique == 'galaxy':
            for s in ['bulge', 'disk', 'knots']:
                if f'sed_val_{s}' in self.native_columns:
                    subs.append(s)
        return subs

    def get_native_attribute(self, attribute_name):
        '''
        Retrieve a particular attribute for a collection
        If we already have it, just return it.  Otherwise attempt
        to fetch.   Reader should check whether the attribute actually
        exists.
        '''
        val = getattr(self, attribute_name, None)
        if val is not None: return val

        for r in self._rdrs:
            if attribute_name in r.columns:
                d = r.read_columns([attribute_name], self._mask)
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
            d = r.read_columns(to_read, self._mask)
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
        pass        #    for now

    # implement Sequence methods
    def __contains__(self, obj):
        '''
        Parameters
        ----------
        obj can be an (object id) or of type BaseObject
        '''
        if type(obj) == type(10):
            id = obj
        else:
            if isinstance(obj, BaseObject):
                id = obj.id
            else:
                raise TypeError
        return id in self._id

    def __len__(self):
        return len(self._id)

    #def __iter__(self):   Standard impl based on __getitem__ should be ok
    #def __reversed__(self):   Default implementation ok

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
            return BaseObject(self._ra[key], self._dec[key], self._id[key],
                              object_type, belongs_to=self,
                              belongs_index=key)

        elif type(key) == slice:
            ixdata = [i for i in range(min(key.stop,len(self._ra)))]
            ixes = itertools.islice(ixdata, key.start, key.stop, key.step)
            return [BaseObject(self._ra[i], self._dec[i], self._id[i],
                               object_type, belongs_to=self, belongs_index=i)
                    for i in ixes]

        elif type(key) == tuple and isinstance(key[0], Iterable):
            #  check it's a list of int-like?
            return [BaseObject(self._ra[i], self._dec[i], self._id[i],
                               object_type, belongs_to=self, belongs_index=i)
                    for i in key[0]]


    def get_partition_id(self):
        return self._partition_id

    def count(self, obj):
        '''
        returns # of occurrences of obj.  It can only be 0 or 1
        '''
        if self.__contains__(obj): return 1
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
        self._total_len += len(coll)
        self._located.append(LocatedCollection(coll, old, self._total_len))

    def append_object_list(self, object_list):
        for e in object_list._located:
            self.append_collection(e.collection)

    def redshifts(self):
        return self.get_native_attribute('redshift')

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
        for  e in self._located:
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

    #def __iter__(self):   Standard impl based on __getitem__ should be ok??
    #def __reversed__(self):   Default implementation ok??

    def __getitem__(self, key):
        '''
        Parameters
        ----------
        If key is an int return a single base object
        If key is a slice return a list of object
        '''
        one_only = isinstance(key, int) or isinstance(key, np.int64)
        is_slice = type(key) == slice
        is_list = type(key) == tuple and isinstance(key[0], Iterable)

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

                    start_ix +=len(sub)
                    if start_ix >= len(key_list):
                        break
                    start = key_list[start_ix]

        if len(to_return):
            return to_return
        else:
            raise IndexError
