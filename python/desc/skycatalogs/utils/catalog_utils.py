#from desc.skycatalogs.objects.sncosmo_object import SncosmoObject
#from desc.skycatalogs.objects.galaxy_object import GalaxyObject
#from desc.skycatalogs.objects.star_object import StarObject

__all__ = ['CatalogContext']
class CatalogContext:
    def __init__(self, the_sky_cat):
        global sky_cat
        self._sky_cat = the_sky_cat
        sky_cat = the_sky_cat
        source_type_dict = {}

        # Initialize with known source types and corresponding object classes
        # Only include collection class if it's not "standard"
        # (that is, not ObjectCollection)
        # source_type_dict['star'] = {'object_class' : StarObject}
        # source_type_dict['sncosmo'] = {'object_class' : SncosmoObject}
        # source_type_dict['galaxy'] = {'object_class' : GalaxyObject}
        self._source_type_dict = source_type_dict

    def register_source_type(self, name, object_class, collection_class=None):
        self._source_type_dict[name] = {'object_class' : object_class,
                                        'collection_class' : collection_class}

    def lookup_source_type(self, name):
        if name in self._source_type_dict:
            return self._source_type_dict[name]['object_class']
        else:
            return None

    def lookup_collection_type(self, name):
        if name in self._source_type_dict:
            if 'collection_class' in self._source_type_dict[name]:
                return self._source_type_dict[name]['collection_class']
            else:
                return None
        else:
            return None
