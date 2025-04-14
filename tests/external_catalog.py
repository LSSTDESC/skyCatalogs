import numpy as np
import galsim
from skycatalogs.objects import BaseObject, ObjectCollection


__all__ = ["ExternalCollection", "register_objects"]


def register_objects(sky_catalog, object_type):
    ExternalCollection.register(sky_catalog, object_type)


class ExternalObject(BaseObject):

    def __init__(self, *args, **kwds):
        super().__init__(self, *args[1:], **kwds)

    def get_observer_sed_component(self, component, mjd=None):
        if component != "this_object":
            raise RuntimeError("Unknown SED component: %s", component)
        return None

    def get_gsobject_components(self, gsparams=None, rng=None):
        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)
        return {"this_object": galsim.DeltaFunction(gsparams=gsparams)}


class ExternalCollection(ObjectCollection):

    def __init__(self, param, object_type, region, sky_catalog):
        self._param = param
        self._sky_catalog = sky_catalog
        self._object_type_unique = object_type
        self._object_class = ExternalObject
        self._uniform_object_type = True
        size = 10
        self._ra = np.random.uniform(0, 1, size=size)
        self._dec = np.random.uniform(0, 1, size=size)
        self._id = np.arange(size)

    @property
    def native_columns(self):
        return ()

    def __len__(self):
        return len(self._ra)

    @staticmethod
    def register(sky_catalog, object_type):
        sky_catalog.cat_cxt.register_source_type(
            object_type,
            object_class=ExternalObject,
            collection_class=ExternalCollection,
            custom_load=True,
        )

    @staticmethod
    def load_collection(region, sky_catalog, mjd=None, exposure=None, object_type=None):
        config = dict(sky_catalog.raw_config["object_types"][object_type])
        return ExternalCollection(
            config["object_param"], object_type, region, sky_catalog
        )
