from collections.abc import Sequence
import numpy as np

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
                 hp_id=None):
        '''
        Minimum information needed for static (not SSO) objects
        '''
        self._ra = ra
        self._dec = dec
        self._id = id
        self._object_type = object_type
        ##self._redshift = redshift    #  do we want this?
        self._hp_id = hp_id


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


class BaseObjectCollection(BaseObject, Sequence):
    '''
    Abstract base class for collection of static objects.
    Many of the methods look the same as for BaseObject but they return arrays
    rather than a single number.  There are some additional methods
    '''
    def __init__(self, ra, dec, id, object_type, include_mask = None,
                 redshift=None,
                 hp_id=None, region=None):
                 #reader=None):
        '''
        Minimum information needed for static objects.
        (Not sure redshift is necessary.  reader even less likely)

        ra, dec, id must be array-like.
        object_type  may be either single value or array-like.
        h_id  may be either single value or array-like or None
        All arrays must be the same length
        If collection is formed by filtering on region, save that; e.g.
        save mask array.
        '''
        self._ra = np.array(ra)
        self._dec = np.array(dec)
        self._id = np.array(id)

        # Maybe the following is silly and hp_id, object_type should always be stored
        # as arrays
        if isinstance(object_type, list):
            self._object_types = np.array(object_type)
            self._object_type_unique = None
        else:
            self._object_types = None
            self._object_type_unique = object_type

        if isinstance(hp_id, list):
            self._hp_ids = np.array(hp_id)
            self._hp_unique = None
        else:
            self._hp_ids = None
            self._hp_unique = hp_id

        if isinstance(redshift, list):
            self._redshift = np.array(redshift)
        else:
            self._redshift = redshift

        self._region = region
        ### self._reader = reader     # Not sure we need this

        self._uniform_object_type = (type(object_type) == type('galaxy'))
        self._uniform_hp_id = (type(hp_id) == type(10)) or (hp_id is None)

        #####if type(object_type) == type
        # do we need to do more here?

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
        return len(self.id)

    #def __iter__(self):   Standard impl based on __getitem__ should be ok
    #def __reversed__(self):   Default implementation ok

    def __getitem__(self, key):
        '''
        Parameters
        ----------
        If key is an int return a single base object
        If key is a slice return a BaseObjectCollection
        '''
        if self._uniform_hp_id:
            hp_id = self._hp_id
        else:
            hp_id = delf._hp_id[key]

        if self._uniform_object_type:
            object_type = self._object_type
        else:
            object_type = self._object_type[key]

        if type(key) == type(10):
            return BaseObject(self.ra[key], self.dec[key]. self.id[key],
                              object_type, hp_id, region, reader)

        else:
            return BaseObjectCollection(self.ra[key], self.dec[key].
                                        self.id[key], object_type, hp_id,
                                        region, reader)

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

###def read_objects(f, object_type=None, region=None):
