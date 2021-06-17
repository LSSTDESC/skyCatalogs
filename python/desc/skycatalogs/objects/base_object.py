from collections.abc import Sequence
import numpy as np
import itertools

'''
Main object types.   There are also may be subtypes. For example,
there could be two subtypes for bulge components, differing in the
form of their associated SEDs
'''

__all__ = ['BaseObject', 'BaseObjectCollection', 'OBJECT_TYPES']
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
                 hp_id=None, belongs_to=None, belongs_index=None):
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
        self._hp_id = hp_id
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
            self._redshift = self._belongs_to.redshifts()[self._belongs_index]
        return self._redshift

    @property
    def hp_id(self):
        if not self._hp_id:
            # Need to compute from ra, dec, but for now
            pass
        return self._hp_id

    def get_flux(self, date_time, band, noMW=False):
        '''
        Parameters
        ----------
        date_time   datetime object; time at which flux is requested
        band        specifies interval over which flux is to be integrated
                    (and filter characteristics?)
        noMW        If true, don't include Milky Way extinction

        Returns
        -------
        Flux of the object for specified time, band.  By default
        include Milky Way extinction.
        '''
        raise NotImplementedError

    def get_sed(self, **kwargs):
        '''
        For galaxies may want to specify subcomponent(s)
        '''
        raise NotImplementedError

    def get_sed_metadata(self, **kwargs):
        '''
        E.g. list of wavelength or frequency intervals associated with sed values
        '''
        raise NotImplementedError


class BaseObjectCollection(Sequence):
    '''
    Abstract base class for collection of static objects.
    Many of the methods look the same as for BaseObject but they return arrays
    rather than a single number.  There are some additional methods
    '''
    def __init__(self, ra, dec, id, object_type, include_mask = None,
                 redshift=None,
                 hp_id=None, indexes = None, region=None, reader=None):
        '''
        Minimum information needed for static objects.
        (Not sure redshift is necessary.  reader even less likely)

        ra, dec, id must be array-like.
        object_type  may be either single value or array-like.
        h_id  may be either single value or array-like or None
        All arrays must be the same length except for the mask, which may
        be longer. # zeros in the mask must equal length of other arrays.
        If mask array is None no filtering was done

        '''
        self._ra = np.array(ra)
        self._dec = np.array(dec)
        self._id = np.array(id)
        self._redshift = None
        self._rdr = reader

        #print("BaseObjectCollection constructor called with hp_id=")
        print(hp_id)
        #print("type is ", type(hp_id))

        # Save the mask in case we need to look up other columns later
        self._include_mask = include_mask

        # Maybe the following is silly and hp_id, object_type should always be stored
        # as arrays
        if isinstance(object_type, list):
            self._object_types = np.array(object_type)
            self._object_type_unique = None
            self._uniform_object_type = False
        else:
            self._object_types = None
            self._object_type_unique = object_type
            self._uniform_object_type = True

        if isinstance(hp_id, list):
            print("hp_id arg is a list")
            self._hp_ids = np.array(hp_id)
            self._hp_unique = None
            self._uniform_hp_id = False
        else:
            print("hp_id arg is a single item: ", hp_id)
            self._hp_ids = None
            self._hp_unique = hp_id
            self._uniform_hp_id = True

        if isinstance(redshift, list):
            self._redshifts = np.array(redshift)
        else:
            self._redshifts = redshift

        self._region = region

    def redshifts(self):
        if not self._redshift:
            # read from our file
            ##self._redshift = self._rdr.read_columns(['redshift'])['redshift']
            self._redshift = self.get_attribute('redshift')
        return self._redshift

    def get_attribute(self, attribute_name):
        '''
        Retrieve a particular attribute for a source.
        If we already have it, just return it.  Otherwise attempt
        to fetch.   Reader should check whether the attribute actually
        exists.
        '''
        val = getattr(self, attribute_name, None)
        if val is not None: return val

        val = self._rdr.read_columns([attribute_name])[attribute_name]
        if val is not None:
            setattr(self, attribute_name, val)
        return val



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
        If key is an int return a single base object
        If key is a slice return a list of object
        '''
        if self._uniform_hp_id:
            hp_id = self._hp_unique
        else:
            hp_id = self._hp_ids[key]

        if self._uniform_object_type:
            object_type = self._object_type_unique
        else:
            object_type = self._object_types[key]

        if type(key) == type(10):
            return BaseObject(self._ra[key], self._dec[key], self._id[key],
                              object_type, hp_id=hp_id, belongs_to=self,
                              belongs_index=key)

        else:
            ixdata = [i for i in range(min(key.stop,len(self._ra)))]
            ixes = itertools.islice(ixdata, key.start, key.stop, key.step)
            return [BaseObject(self._ra[i], self._dec[i], self._id[i],
                               object_type, hp_id, belongs_to=self,
                               belongs_index=i)
                    for i in ixes]

    def get_hpid(self):
        return self._hp_unique        # None if not uniform

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

    def add_columns(self, column_dict):
        '''
        Store more information about this collection
        '''
        pass
