from collections.abc import Sequence, Iterable
# from collections import namedtuple
# import os
import numpy as np

'''
Minimal base classes to be subclassed when adding new, externally
maintained, source types.
'''

__all__ = ['MinObject', 'MinCollection']


class MinObject(object):
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
    def belongs_to(self):
        return self._belongs_to

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
        raise NotImplementedError('get_observer_sed_component must be implemented by subclass')

    def get_LSST_flux(self, band, sed=None, mjd=None):
        raise NotImplementedError('get_observer_sed_component must be implemented by subclass')


class MinCollection(Sequence):
    '''
    Class for collection of objects of the same type coming from the same
    place (e.g., region of the sky)
    Some methods look the same as for MinObject but they return arrays
    rather than a single number.  There are some additional methods
    '''
    def __init__(self, ra, dec, id, object_type, sky_catalog,
                 region=None, mjd=None):
        '''
        Parameters
        ra, dec      float, array-like of same length
        id           array-like, same length as ra and dec.   int or string
        object_type  single string or (if contained objects may be of
                     different types) array-like
        sky_catalog  Instance of SkyCatalog class, typically obtained by
                     calling open_catalog
        region       maybe be used to determine which objects are in the
                     collection
        mjd          MJD of the observation epoch.  This is used by
                     time-varying objects, e.g., SNe, stars.
        '''
        self._ra = np.array(ra)
        self._dec = np.array(dec)
        self._id = np.array(id)
        self._sky_catalog = sky_catalog
        self._object_class = sky_catalog.cat_cxt.lookup_source_type(object_type)
        self._object_type_unique = object_type
        self._region = region
        self._mjd = mjd

    @property
    def mjd(self):
        return self._mjd

    @property
    def sky_catalog(self):
        return self._sky_catalog

    @property
    def subcomponents(self):
        '''
        Return list of all subcomponents for objects in this collection.
        Typically only galaxy-like objects have true subcomponents.
        Others should return ['this_object']
        '''
        raise NotImplementedError('get_observer_sed_component must be implemented by subclass')

    # implement Sequence methods
    def __contains__(self, obj):
        '''
        Parameters
        ----------
        obj can be an (object id) or of type MinObject
        '''
        raise NotImplementedError('get_observer_sed_component must be implemented by subclass')

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

        object_type = self._object_type_unique

        if isinstance(key, int) or isinstance(key, np.int64):
            return self._object_class(self._ra[key], self._dec[key],
                                      self._id[key], object_type, self, key)

        elif type(key) is slice:
            if key.start is None:
                key.start = 0
            return [self.__getitem__(i) for i in range(self.__len__())[key]]

        elif type(key) is tuple and isinstance(key[0], Iterable):
            return [self.__getitem__(i) for i in key[0]]

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

    @staticmethod
    def load_collection(region, skycatalog, mjd=None, exposure=None):
        '''
        region      One of Disk, PolygonalRegion from skyCatalogs.utils.shapes.
                    Support of PolygonalRegion is required
        skycatalog  An instance of the SkyCatalog class
        mjd         Time at which objects are to be assembled.
        exposure    exposure length.

        Returns an instance of MinCollection
        '''
        raise NotImplementedError('get_observer_sed_component must be implemented by subclass')
