import os
from pathlib import Path
from skycatalogs.skyCatalogs import open_catalog

PIXEL = 9556


def explore(cat, obj_type, ix_list=[0]):
    obj_list = cat.get_object_type_by_hp(PIXEL, obj_type)

    collects = obj_list.get_collections()
    print('Number of collections: ', len(collects))

    icoll = 0
    for c in collects:
        print(f'Object type: {c._object_type_unique}')
        print('Native columns:')
        print(c.native_columns)

        len_coll = len(c)
        obj0 = c[0]
        print(f"\nFor hpid {c.get_partition_id()} collection {icoll} found {len_coll} objects")
        print("First object: ")
        print('id=', obj0.id, ' ra=', obj0.ra, ' dec=', obj0.dec,
              ' belongs_index=', obj0._belongs_index)

        all_cols = obj0.native_columns

        for native in all_cols:
            print(native, '=', obj0.get_native_attribute(native))
        icoll += 1


skycatalog_root = os.path.join(Path(__file__).resolve().parents[1],
                               'skycatalogs', 'data')
config_path = os.path.join(skycatalog_root, 'ci_sample', 'skyCatalog.yaml')
cat = open_catalog(config_path, skycatalog_root=skycatalog_root)

print('Explore star collection')
explore(cat, 'star', ix_list=[10])
