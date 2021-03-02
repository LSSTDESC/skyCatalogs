from .base_object import BaseObject, GALAXY, GALAXY_BULGE, GALAXY_DISK, GALAXY_KNOTS
def class GalaxyObject(BaseObject):
'''
Galaxy object.
'''
    def __init__(self, ra, dec, id, hp_id):
        '''
        Minimum information needed for static (not SSO) objects
        Type of object_type should perhaps be an enumeration
        '''
        super().__init__(ra, dec, id, hp_id, BaseObject.GALAXY)
        self._cmps = {GALAXY_BULGE : BulgeObject(ra, dec, id, hp_id, self),
                      GALAXY_DISK : DiskObject(ra, dec, id, hp_id, self),
                      GALAXY_KNOTS : KnotsObject(ra,dec, id, hp_id, self)}

    def get_subcomponent(self, cmp_type):
        return self._cmps.get(cmp_type)

    def get_flux(self, date_time, band):
        '''
        Parameters
        ----------
        date_time   datetime object; time at which flux is requested
        band        specifies interval over which flux is to be integrated
                    (and filter characteristics?)
        '''
        raise NotImplementedError

    def get_sed(self, subcomponent_list=[GALAXY_BULGE, GALAXY_DISK, GALAXY_KNOTS], **kwargs):
        '''
        For galaxies may want to specify subcomponent(s)
        '''
        d = {}
        for c in subcomponent_list:
            d[c] = self._cmps[c].get_sed(**kwargs)

        return d

    def get_sed_metadata(self, **kwargs):
        '''
        Returns, e.g,. units (wavelength or frequency) and list intervals associated
        with sed values
        '''
        d = {}
        for c in subcomponent_list:
            d[c] = self._cmps[c].get_sed_metadata(**kwargs)

        return d

# If galaxy subcomponents all implement their methods the same way probably can
# get by with just one subcomponent class.  Does the bulge need to know it's a
# bulge, or is it sufficient that the parent knows?

def class GalaxyBulge(BaseObject):
    def __init__(self, ra, dec, id, parent):
        super().__init__(ra, dec, id, BaseObject.GALAXY_BULGE)
        self._parent = parent

def class GalaxyDisk(BaseObject):
    def __init__(self, ra, dec, id, parent):
        super().__init__(ra, dec, id, BaseObject.GALAXY_DISK)
        self._parent = parent

def class GalaxyKnots(BaseObject):
    def __init__(self, ra, dec, id, parent):
        super().__init__(ra, dec, id, BaseObject.GALAXY_KNOTS)
        self._parent = parent
