"""
Unit tests for SkyCatalogs API
"""

import unittest
import os
from pathlib import Path

# Not currently used
#import pandas as pd
#import numpy as np

from desc.skycatalogs.skyCatalogs import SkyCatalog, open_catalog, Box
from desc.skycatalogs.skyCatalogs import _get_intersecting_hps

# This only works if my scratch area is accessible in the run environment
class APITester(unittest.TestCase):

    def setUp(self):
        '''
        Open the catalog
        '''
        skycatalog_root = os.path.join(Path(__file__).resolve().parents[1],
                                       'data')
        cfg_path = os.path.join(skycatalog_root, 'ci_sample', 'skyCatalog.yaml')
        self._cat = open_catalog(cfg_path, skycatalog_root=skycatalog_root)


    def tearDown(self):
        pass                  # nothing to do

    def testAPI_region(self):
        '''
        Should be broken down into separate small tests
        '''
        cat = self._cat


        hps = cat._find_all_hps()
        print('Found {} healpix pixels '.format(len(hps)))
        for h in hps: print(h)

        ra_min_tract = 55.736
        ra_max_tract = 57.564
        dec_min_tract = -37.190
        dec_max_tract = -35.702

        # Catalogs for testing are sparse, so make a large box
        ra_min_small = 55.8
        ra_max_small = 56.4
        dec_min_small = -36.5
        dec_max_small = -35.9

        rgn = Box(ra_min_small, ra_max_small, dec_min_small, dec_max_small)

        intersect_hps = _get_intersecting_hps('ring', 32, rgn)

        print("For region ", rgn)
        print("intersecting pixels are ", intersect_hps)

        print('Invoke get_objects_by_region with box region')
        object_list = cat.get_objects_by_region(rgn,
                                                obj_type_set=set(['galaxy']) )

        print('Number of collections returned:  ', object_list.collection_count)

        colls = object_list.get_collections()
        for c in colls:
            len_coll = len(c)
            print(f"For hpid {c.get_partition_id()} found {len_coll} objects")
            print("First object: ")
            print(c[0], '\nid=', c[0].id, ' ra=', c[0].ra, ' dec=', c[0].dec,
                  ' belongs_index=', c[0]._belongs_index)

            print("Slice [1:3]")
            slice13 = c[1:3]
            for o in slice13:
                print('id=',o.id, ' ra=',o.ra, ' dec=',o.dec, ' belongs_index=',
                      o._belongs_index)
                other = min(1000, len_coll - 1)
                print(f"Object {other}")
                print(c[other], '\nid=', c[other].id, ' ra=', c[other].ra, ' dec=',
                      c[other].dec,
                      ' belongs_index=', c[other]._belongs_index)
            slice_late = c[len_coll - 5:len_coll - 2]
            print(f'\nobjects indexed {len_coll - 5}  through {len_coll - 3}')

            for o in slice_late:
                print('id=',o.id, ' ra=',o.ra, ' dec=',o.dec, ' belongs_index=',
                      o._belongs_index)

        print('Total object count: ', len(object_list))

        obj = object_list[0]
        print("Type of element in object_list:", type(obj))

        redshift0 = object_list[0].get_native_attribute('redshift')
        print('First redshift: ', redshift0)

        sed_bulges = colls[0].get_native_attribute('sed_val_bulge')

        print("first bulge sed:")
        for v in sed_bulges[0]:
            print(v)


    def testAPI_hp(self):
        '''
        Try with and without specified object types
        '''
        def print_from_objects(label, obj_list):
            '''
        Print out information from supplied object list, accessed
        both from list and from collections it comes from. Values
        returned should be the same
        Parameters
        ----------
        label      string       Identify the list
        list       an object list returned by .get_objects_by...
        '''
            print(label)
            colls = obj_list.get_collections()

            print(f'Total objects: {len(obj_list)}  # collections {len(colls)}')

            obj0 = obj_list[0]
            obj0 = colls[0][0]
            native_columns = obj0.native_columns
            for att in ['object_type', 'redshift', 'id', 'galaxy_id']:

                if att in native_columns:
                    v0 = (obj_list[0]).get_native_attribute(att)
                    v0c0 = (colls[0][0]).get_native_attribute(att)
                    print(f'For attribute {att} obj_list 1st value: {v0}; coll 1st value: {v0c0}')
                    assert v0 == v0c0
                else:
                    print(f'No native attribute "{att}" for this object')
        ####   end of included function for printing

        hp = 9556

        cat = self._cat
        obj_types = {'galaxy'}
        specified = cat.get_objects_by_hp(hp, region=None, obj_type_set=obj_types)
        print_from_objects(f'From hp {hp} with object types {str(obj_types)}', specified)

        none_specified = cat.get_objects_by_hp(hp)
        print_from_objects(f'From hp {hp}, no specified object types', none_specified)

if __name__ == '__main__':
    unittest.main()
