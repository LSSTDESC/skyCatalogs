"""
Unit tests for SkyCatalogs API
"""

import unittest
import os
from pathlib import Path

# Not currently used
#import pandas as pd
#import numpy as np

from desc.skycatalogs.skyCatalogs import SkyCatalog, open_catalog, Box, Disk
from desc.skycatalogs.skyCatalogs import _get_intersecting_hps
from desc.skycatalogs.objects.base_object import BaseObject

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
        Exercise get_objects_by_region for box and disk and parts of the
        interfaces for BaseObject, ObjectCollection and ObjectList
        '''
        cat = self._cat


        hps = cat._find_all_hps()
        print('Found {} healpix pixels '.format(len(hps)))
        for h in hps: print(h)
        assert(set(hps) == {9556, 9683, 9684, 9812, 9813, 9940})

        # These hps are the ones covering tract 3828
        # For this tract ra is between 55.736 and 57.564
        # dec is between -37.190 and -35.702

        # Catalogs for testing are sparse, so make a large box
        ra_min_small = 55.8
        ra_max_small = 56.4
        dec_min_small = -36.5
        dec_max_small = -35.9

        rgn = Box(ra_min_small, ra_max_small, dec_min_small, dec_max_small)

        intersect_hps = _get_intersecting_hps('ring', 32, rgn)

        print("For region ", rgn)
        print("intersecting pixels are ", intersect_hps)
        assert(set(intersect_hps) == {9683, 9684, 9812})

        object_list = cat.get_objects_by_region(rgn,
                                                obj_type_set=set(['galaxy']) )

        print('Number of collections returned for box:  ',
              object_list.collection_count)
        assert(object_list.collection_count == 2)
        colls = object_list.get_collections()
        for c in colls:
            len_coll = len(c)
            print(f"For hpid {c.get_partition_id()} found {len_coll} objects")
            print("First object: ")
            print('id=', c[0].id, ' ra=', c[0].ra, ' dec=', c[0].dec,
                  ' belongs_index=', c[0]._belongs_index)
            fluxes = c[0].get_LSST_fluxes()
            assert(len(fluxes) == 6)
            print('Found fluxes')

            print("Slice [1:3]")
            slice13 = c[1:3]
            for o in slice13:
                print('id=',o.id, ' ra=',o.ra, ' dec=',o.dec, ' belongs_index=',
                      o._belongs_index)
                other = min(1000, len_coll - 1)
                print(f"Object {other}")
                print('id=', c[other].id, ' ra=', c[other].ra, ' dec=',
                      c[other].dec,
                      ' belongs_index=', c[other]._belongs_index)
            slice_late = c[len_coll - 5:len_coll - 2]
            print(f'\nobjects indexed {len_coll - 5}  through {len_coll - 3}')

            for o in slice_late:
                print('id=',o.id, ' ra=',o.ra, ' dec=',o.dec, ' belongs_index=',
                      o._belongs_index)

        print('Total object count: ', len(object_list))

        obj = object_list[0]
        assert(type(obj) == BaseObject)

        redshift0 = object_list[0].get_native_attribute('redshift')
        print('First redshift: ', redshift0)

        sed_bulges = colls[0].get_native_attribute('sed_val_bulge')

        #print("first bulge sed:")
        #for v in sed_bulges[0]:
        #    print(v)

        # Now disk.  Get all known object types
        rgn = Disk(56.6, -36.4, 1800)
        object_list = cat.get_objects_by_region(rgn)
        print('Number of collections returned for disk: ',
              object_list.collection_count)
        colls = object_list.get_collections()
        assert(len(object_list.get_collections()) == object_list.collection_count)
        assert(object_list.collection_count == 4)

        total_object_count = 0
        for c in colls:
            len_coll = len(c)
            total_object_count += len_coll
            print(f"For hpid {c.get_partition_id()} found {len_coll} objects")
            if c._object_type_unique:
                print(f"    of object type {c._object_type_unique}")
            print("First object: ")
            print('id=', c[0].id, ' ra=', c[0].ra, ' dec=', c[0].dec,
                  ' belongs_index=', c[0]._belongs_index, ' object_type=',
                  c[0].object_type)
            fluxes = c[0].get_LSST_fluxes()
            assert(len(fluxes) == 6)
            print('Found fluxes')

            assert(type(c[0]) == BaseObject)


        print('List len: ', len(object_list))
        print('Sum over collections: ', total_object_count)
        assert(total_object_count == len(object_list))



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
