from .base_object import BaseObject, GALAXY, GALAXY_BULGE, GALAXY_DISK, GALAXY_KNOTS
def class GalaxyObject(BaseObject):
'''
Galaxy object.
'''
      def __init__(self, ra, dec, id, redshift=None, hp_id=None,
                   belongs_to=None, belongs_index=None)
        '''
        Minimum information needed for position static (not SSO) objects
        Type of object_type should perhaps be an enumeration
        '''
        super().__init__(ra, dec, id, BaseObject.GALAXY, redshift,
                         hp_id, belongs_to, belongs_index)

        self._cmps = {GALAXY_BULGE : GalaxySub(self, GALAXY_BULGE),
                      GALAXY_DISK : GalaxySub(self, GALAXY_DISK)}
        # Don't always have knots.  When we do, handling could be different
        #               GALAXY_KNOTS : GalaxySub(self, GALAXY_KNOTS)
        #    or         GALAXY_KNOTS : KnotsObject(self)

    def get_subcomponent(self, cmp_type):
        return self._cmps.get(cmp_type)

    def get_bulge_magnorm(self):
        if self._bulge_magnorm:
            return self._bulge_magnorm
        if self._belongs_to:
            self._bulge_magnorm = self._belongs_to.redshifts()[self._belongs_index]
            return self._bulge_magnorm

    def get_disk_magnorm(self):
        if self._disk_magnorm:
            return self._disk_magnorm
        if self._belongs_to:
            self._disk_magnorm = self._belongs_to.redshifts()[self._belongs_index]
        return self._disk_magnorm

    ##### -------------------------------------------
    # commented stuff below needs thought, likely redesign. Or, in the
    # case of get_flux, we probably just get rid of it

    # def get_flux(self, date_time, band):
    #     '''
    #     Parameters
    #     ----------
    #     date_time   datetime object; time at which flux is requested
    #     band        specifies interval over which flux is to be integrated
    #                 (and filter characteristics?)
    #     '''
    #     raise NotImplementedError

    # def get_sed(self, subcomponent_list=[GALAXY_BULGE, GALAXY_DISK, GALAXY_KNOTS], **kwargs):
    #     '''
    #     For galaxies may want to specify subcomponent(s)
    #     '''
    #     d = {}
    #     for c in subcomponent_list:
    #         d[c] = self._cmps[c].get_sed(**kwargs)

    #     return d

    # def get_sed_metadata(self, **kwargs):
    #     '''
    #     Returns, e.g,. units (wavelength or frequency) and list intervals associated
    #     with sed values
    #     '''
    #     d = {}
    #     for c in subcomponent_list:
    #         d[c] = self._cmps[c].get_sed_metadata(**kwargs)

    #     return d

# All galaxy information which is not indirected (e.g. if SEDs are stored
# in files there will be separate ones for different components) is read
# in at once. Provide a way to restrict information returned to a particular
# component if desired


def class GalaxySub(object):
    def __init__(self, parent, component_type):
        self._parent = parent
        self._cmp = component_type

    @property
    def ra(self):
        return self._parent.ra

    @property
    def dec(self):
        return self._parent.dec

    @property
    def id(self):
        return self._parent.id

    @property
    def object_type(self):
        return self._cmp

    @property
    def redshift(self):
        return self._parent.redshift

    def get_magnorm(self):
        if self._cmp == GALAXY_BULGE:
            return self._parent.get_bulge_magnorm()
        if self._cmp == GALAXY_DISK:
            return self._parent.get_disk_magnorm()

    def get_sed(self):
        pass

    def

def class GalaxyKnots(object):
    def __init__(self, parent):
        self._parent = parent
        self._cmp = GALAXY_KNOTS
