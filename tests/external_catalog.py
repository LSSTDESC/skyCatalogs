import galsim
from skycatalogs.objects import BaseObject, ObjectCollection


__all__ = ["ExternalCollection", "register_objects"]


def register_objects(sky_catalog):
    ExternalCollection.register(sky_catalog)


class ExternalObject(BaseObject):

    def __init__(self):
        super().__init__(self)

    def get_observer_sed_component(self, component, mjd=None):
        if component != "this_object":
            raise RuntimeError("Unknown SED component: %s", component)
        return None

    def get_gsobject_components(self, gsparams=None, rng=None):
        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)
        return {'this_object': galsim.DeltaFunction(gsparams=gsparams)}


class ExternalCollection(ObjectCollection):
    _object_type = "external_objects"

    def __init__(self, region, sky_catalog):
        pass

    @property
    def native_columns(self):
        return ()

    def __len__(self):
        return 0

    @staticmethod
    def register(sky_catalog):
        sky_catalog.cat_cxt\
            .register_source_type(ExternalCollection._object_type,
                                  object_class=ExternalObject,
                                  collection_class=ExternalCollection,
                                  custom_load=True)

    @staticmethod
    def load_collection(region, sky_catalog, mjd=None, exposure=None):
        return ExternalCollection(region, sky_catalog)
