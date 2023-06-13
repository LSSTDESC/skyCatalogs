import os
import sys
from desc.skycatalogs.utils.shapes import Disk
from desc.skycatalogs.skyCatalogs import open_catalog
from desc.skycatalogs.objects.gaia_object import GaiaCollection

if __name__ == '__main__':
    rad_degrees = 0.17
    disk = Disk(60, -40, rad_degrees * 360)

    config_path = os.path.join(os.environ['SKYCATALOGS_DIR'], 'local', 'gaia',
                               'gaia_config.yaml')

    skycat = open_catalog(config_path)
    collection = GaiaCollection.load_collection(disk, skycat)

    print('collection size: ', len(collection))
    obj = collection[0]

    print('For initial object:')
    print(f'id is {obj.id}   ra,dec are  {obj.ra}, {obj.dec}')

    print('id for objects in slice [0:2]:')
    for o in collection[0:2]:
        print(o.id)

    count = 0
    for o in collection:
        count = count + 1
        print(o._belongs_index)
        
    print(f'Counted {count} objects in the collection')
