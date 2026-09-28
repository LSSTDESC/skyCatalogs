import numpy as np
import galsim
from skycatalogs.objects import BaseObject, ObjectCollection


__all__ = ["ExternalCollection", "register_objects"]


def register_objects(sky_catalog, object_type):
    '''
    When skyCatalogs encounters a third-party object type in its configuration
    it calls this routine to get the ball rolling.
    '''
    ExternalCollection.register(sky_catalog, object_type)


class ExternalObject(BaseObject):

    def __init__(self, ra, dec, id, object_type, belongs_to, belongs_index):
        '''
        Parameters
        ----------
        id     unique identifier
        object_type (string) will typicallly come from the configuration file
        belongs_to  is an instance of an ObjectCollection
        belongs_index is the index of this object in that collection
        '''
        super().__init__(self, ra, dec, id, object_type, belongs_to,
                         belongs_index)

    def get_observer_sed_component(self, component, mjd=None):
        '''
        Objects may have multiple components (e.g. "disk", "bulge" for
        a galaxy-like type) or only one. If only one, which is what we're
        assuming here, by default that component will be called "this_object".
        mjd would be needed to construct a SED for time-varying object types

        Returns
        -------
        galsim SED to which extinction has been applied
        '''
        # check that component is a known component type
        if component != "this_object":
            raise RuntimeError("Unknown SED component: %s", component)

        # Get the parametrization associated with this object type
        object_param = self._belongs_to._param

        # From the various inputs at our disposal construct a SED
        the_sed = ...    # some code needs to go here
        return the_sed


    def get_gsobject_components(self, gsparams=None, rng=None):
        '''
        Parameters
        ----------
        gsparams:     See galsim.GSParams
        rng:          Instance of galsim.BaseDeviate or subclass

        Returns
        -------
        dict of galsim objects (e.g. DeltaFunction, Sersic, etc.),
        lensed if appropriate, and indexed by component
        '''
        # For this example our object type is a point source
        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)
        return {"this_object": galsim.DeltaFunction(gsparams=gsparams)}


class ExternalCollection(ObjectCollection):

    def __init__(self, param, object_type,  sky_catalog, region):
        '''
        Parameters
        ----------
        param: value of "object_param" in the section of the configuration
               file for our object_type.   This is not a required feature
               for third-party types. It is used here to demonstrate
               how to pass information from the configuration file
               to the python module implementing the object type.
        object_type: string.   The object type name
        sky_catalog: instance of SkyCatalog class
        region: instance of Region class or None

        '''
        # generate ra, dec, id for this toy example.
        size = 10
        ra = np.random.uniform(0, 1, size=size)
        dec = np.random.uniform(0, 1, size=size)
        id = np.arange(size)

        # partition_id is used when input data for the catalog is partitioned,
        # e.g. by healpixel
        partition_id = None

        # Some optional arguments not appearing here may be useful when
        # inputs are in the form of parquet files (readers, row_group)
        # or for applying restrictions (mjd, mask)
        super().__init__(self, ra, dec, id, object_type, partition_id,
                         sky_catalog, region=region)

        # Note self._param will be available not only to the collection
        # but also to the object instances, so may be used to condition
        # behavior of, e.g., get_observer_sed_component
        self._param = param

    @staticmethod
    def register(sky_catalog, object_type):
        # This is how skyCatalogs knows which code to associate with
        # our object type
        sky_catalog.cat_cxt.register_source_type(
            object_type,
            object_class=ExternalObject,
            collection_class=ExternalCollection,
            custom_load=True,
        )

    @staticmethod
    def load_collection(region, sky_catalog, mjd=None, exposure=None,
                        object_type=None):
        # Getting access to that part of the yaml coonfiguration devoted
        # to our object type
        config = dict(sky_catalog.raw_config["object_types"][object_type])

        # Pass value of object_parm to the collection (and therefore also to
        # individual objects).
        return ExternalCollection(
            config["object_param"], object_type, sky_catalog, region
        )
