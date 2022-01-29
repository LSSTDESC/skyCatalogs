from collections.abc import Sequence, Iterable
from collections import namedtuple, OrderedDict
import os
import gzip
import numpy as np
import itertools
import warnings
from desc.skycatalogs.utils.translate_utils import form_object_string
from desc.skycatalogs.utils.config_utils import Config
from desc.skycatalogs.utils.sed_utils import convert_tophat_sed

'''
Main object types.   There are also may be subtypes. For example,
there could be two subtypes for bulge components, differing in the
form of their associated SEDs
'''

__all__ = ['BaseObject', 'ObjectCollection', 'ObjectList', 'OBJECT_TYPES']
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

class BaseObject(object):
    '''
    Abstract base class for static (in position coordinates) objects.
    Likely need a variant for SSO.
    '''
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
        if self.object_type == 'galaxay':
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
                     for stars; defaults to 1.0 for galaxy components

        Returns
        -------
        A triple (lambda, f_lambda, mag_norm) where the first two are
                 parallel floating point arrays and the last is a scalar.
                 lambda is in nm.   f_lambda is in erg/(cm**2 * s * nm)
        '''

        if self._object_type == 'star':
            # read in the file; return data from it. Arguments are ignored
            if component:
                warnings.warn(f'get_sed: component argument ignored for object type {self._object_type}')
            if resolution:
                warnings.warn(f'get_sed: resolution argument ignored for object type {self._object_type}')
            mag_norm = self.get_native_attribute('magnorm')
            fpath = os.path.join(os.getenv('SIMS_SED_LIBRARY_DIR'),
                                 self.get_native_attribute('sed_filepath'))
            lmbda = []
            f_lambda = []
            with gzip.open(fpath, mode='rt') as f:
                ln = (f.readline()).strip()
                while ln:
                    if len(ln) == 0:
                        ln = (f.readline()).strip()
                        continue
                    if ln[0] == '#':
                        ln = (f.readline()).strip()
                        continue
                    atoms = ln.split()
                    lmbda.append(float(atoms[0]))
                    f_lambda.append(float(atoms[1]))
                    ln = (f.readline()).strip()
            return lmbda, f_lambda, mag_norm

        if self._object_type != 'galaxy':
            raise ValueError('get_sed function only available for stars and galaxy components')
        if not resolution:
            resolution = 1.0
        if component not in ['disk', 'bulge', 'knots']:
            raise ValueError(f'Cannot fetch SED for component type {component}')

        th_val = self.get_native_attribute(f'sed_val_{component}')
        if  th_val is None:   #  values for this component are not in the file
            raise ValueError(f'{component} not part of this catalog')
        r = self.get_native_attribute('redshift_hubble')
        mag_f = self._belongs_to.sky_catalog.mag_norm_f
        th_bins = self._belongs_to.config.get_tophat_parameters()

        lmbda,f_lambda,mag_norm,f_nu500 = convert_tophat_sed(th_bins, th_val,
                                          mag_f, redshift=r, wavelen_step=resolution)
        return lmbda, f_lambda, mag_norm


class ObjectCollection(Sequence):
    '''
    Class for collection of static objects coming from the same
    source (e.g., file for particular healpix)
    Many of the methods look the same as for BaseObject but they return arrays
    rather than a single number.  There are some additional methods
    '''
    def __init__(self, ra, dec, id, object_type, partition_id, sky_catalog,
                 region=None, mask=None, reader=None):
        '''
        Minimum information needed for static objects.
        specified, has already been used by the caller to generate mask.
        ra, dec, id must be array-like.
        object_type  may be either single value or array-like.
        partition_id (e.g. healpix id)
        sky_catalog instance of SkyCatalog class
        (redshift should probably be ditched; no need for it)
        (similarly for region with current code structure. Information
        needed is encoded in mask)
        mask  indices to be masked off, e.g. because region does not
              include the entire healpix pixel
        All arrays must be the same length


        '''
        self._ra = np.array(ra)
        self._dec = np.array(dec)
        self._id = np.array(id)
        self._rdr = reader
        self._partition_id = partition_id
        self._sky_catalog = sky_catalog   # might not need this
        self._config = Config(sky_catalog.raw_config)
        self._mag_norm_f = sky_catalog.mag_norm_f
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
        return self._rdr.columns

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

        d = self._rdr.read_columns([attribute_name], self._mask)
        if not d:
            return None
        val = d[attribute_name]
        setattr(self, attribute_name, val)
        return val

    def get_native_attributes(self, attribute_list):
        '''
        Return requested attributes as an OrderedDict. Keys are column names.
        Use our mask if we have one
        '''
        df = self._rdr.read_columns(attribute_list, self._mask)
        return df

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
