'''
Main object types.   There are also may be subtypes. For example,
there could be two subtypes for bulge components, differing in the
form of their associated SEDs
'''
GALAXY=1
GALAXY_BULGE=2
GALAXY_DISK=3
GALAXY_KNOTS=4
STAR=5
AGN=6
SN=7

def class BaseObject(object):
'''
Abstract base class for static objects. Might need a variant for SSO.
'''
    def __init__(self, ra, dec, id, object_type, hp_id=None, reader=None):
        '''
        Minimum information needed for static (not SSO) objects
        '''
        self._ra = ra
        self._dec = dec
        self._id = id
        self._object_type = object_type
        self._hp_id = hp_id

        # All objects also include redshift and extinction information


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
