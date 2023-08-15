"""
Unit tests for SkyCatalogs API
"""

import unittest
import os
from pathlib import Path
import numpy as np

# Not currently used
#import pandas as pd


from skycatalogs.skyCatalogs import SkyCatalog, open_catalog
from skycatalogs.skyCatalogs import Box, Disk, PolygonalRegion
from skycatalogs.skyCatalogs import _get_intersecting_hps
from skycatalogs.objects.base_object import BaseObject

class APITester(unittest.TestCase):

    def setUp(self):
        '''
        Open the catalog
        '''
        skycatalog_root = os.path.join(Path(__file__).resolve().parents[1],
                                       'skycatalogs', 'data')
        self._skycatalog_root = skycatalog_root
        cfg_path = os.path.join(skycatalog_root, 'ci_sample', 'skyCatalog.yaml')
        self._cat = open_catalog(cfg_path, skycatalog_root=skycatalog_root)


    def tearDown(self):
        pass                  # nothing to do

    def testAPI_region(self):
        '''
        Exercise get_objects_by_region for box and disk and parts of the
        interfaces for BaseObject, ObjectCollection and ObjectList
        '''
        import time
        cat = self._cat


        hps = cat._find_all_hps()
        print('Found {} healpix pixels '.format(len(hps)))
        for h in hps: print(h)
        assert(set(hps) == {9556, 9557, 9683, 9684, 9812, 9813, 9940})

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

        box_count = len(object_list)
        print('Total object count: ', box_count)

        obj = object_list[0]
        assert(isinstance(obj, BaseObject))

        redshift0 = object_list[0].get_native_attribute('redshift')
        print('First redshift: ', redshift0)

        sed_bulges = colls[0].get_native_attribute('sed_val_bulge')

        # Now make a polygon which describes same region as the box
        vertices = [(ra_min_small, dec_min_small),
                    (ra_min_small, dec_max_small),
                    (ra_max_small, dec_max_small),
                    (ra_max_small, dec_min_small)]
        rgn_poly = PolygonalRegion(vertices_radec=vertices)
        print("For polygonal region ", rgn_poly)
        print("intersecting pixels are ", intersect_hps)
        assert(set(intersect_hps) == {9683, 9684, 9812})

        t0 = time.time()
        object_list = cat.get_objects_by_region(rgn_poly,
                                                obj_type_set=set(['galaxy']) )
        t_end = time.time()
        print("Time to get box-like polygon objects: ", t_end - t0)


        print('Number of collections returned for polygon:  ',
              object_list.collection_count)
        assert(object_list.collection_count == 2)
        colls = object_list.get_collections()
        for c in colls:
            len_coll = len(c)
            print(f"For hpid {c.get_partition_id()} found {len_coll} objects")

        poly_count = len(object_list)
        print('Total object count in polygon: ', poly_count)

        assert(abs(poly_count - box_count) < 5)

        # Now try a convex polygon which is more of a diamond than a box
        top = ((ra_min_small + ra_max_small)/2, dec_max_small)
        bottom = ((ra_min_small + ra_max_small)/2, dec_min_small)
        left = (ra_min_small, (dec_min_small + dec_max_small)/2)
        right = (ra_max_small, (dec_min_small + dec_max_small)/2)
        vertices2 = [top, right, bottom, left]

        rgn_poly2 = PolygonalRegion(vertices_radec=vertices2)

        print("For polygonal region ", rgn_poly2)
        print("intersecting pixels are ", intersect_hps)

        t0 = time.time()
        object_list = cat.get_objects_by_region(rgn_poly2,
                                                obj_type_set=set(['galaxy']) )
        t_end = time.time()
        print("Time to get diamond polygon objects: ", t_end - t0)


        print('Number of collections returned for diamond polygon:  ',
              object_list.collection_count)
        #assert(object_list.collection_count == 2)
        colls = object_list.get_collections()
        for c in colls:
            len_coll = len(c)
            print(f"For hpid {c.get_partition_id()} found {len_coll} objects")

        poly2_count = len(object_list)
        print('Total object count in polygon: ', poly2_count)

        # Now disk.  Get all known object types
        rgn = Disk(56.6, -36.4, 1800)
        t0 = time.time()
        object_list = cat.get_objects_by_region(rgn)
        t_end = time.time()
        print("Time to get disk objects: ", t_end - t0)

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

            assert(isinstance(c[0], BaseObject))


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
        obj_type = 'galaxy'
        specified = cat.get_object_type_by_hp(hp, obj_type, region=None)
        print_from_objects(f'From hp {hp} with object types {obj_type}', specified)
    def testAPI_rowgroup(self):
        '''
        Compare results of get_objects_by_region when skyCatalog files
        each have a single row group or have several
        '''
        def compare_objects(no_rg, has_rg):
            no_rg_coll = no_rg.get_collections()
            has_rg_coll = has_rg.get_collections()

            no_rg_pieces = tuple([c.get_native_attribute('galaxy_id') for c in no_rg_coll])
            no_rg_concat = np.concatenate(no_rg_pieces)
            has_rg_pieces = tuple([c.get_native_attribute('galaxy_id') for c in has_rg_coll])
            has_rg_concat = np.concatenate(has_rg_pieces)

            assert(np.all(no_rg_concat == has_rg_concat))

        box = Box(56.0, 56.8, -36.8, -36.0)
        disk = Disk(56.4, -36.5, 3000)


        # Now try a convex polygon which is more of a diamond than a box
        ra_min_small = 55.8
        ra_max_small = 56.4
        dec_min_small = -36.5
        dec_max_small = -35.9

        top = ((ra_min_small + ra_max_small)/2, dec_max_small)
        bottom = ((ra_min_small + ra_max_small)/2, dec_min_small)
        left = (ra_min_small, (dec_min_small + dec_max_small)/2)
        right = (ra_max_small, (dec_min_small + dec_max_small)/2)
        vertices2 = [top, right, bottom, left]

        rhomb = PolygonalRegion(vertices_radec=vertices2)


        cat_one = self._cat
        rg_cfg_path = os.path.join(self._skycatalog_root, 'row_groups',
                                   'skyCatalog.yaml')
        cat_rg = open_catalog(rg_cfg_path, skycatalog_root=self._skycatalog_root)

        box_one = cat_one.get_objects_by_region(box,
                                                obj_type_set=set(['galaxy']))
        box_rg = cat_rg.get_objects_by_region(box,
                                              obj_type_set=set(['galaxy']))
        compare_objects(box_one, box_rg)

        disk_one = cat_one.get_objects_by_region(disk,
                                                 obj_type_set=set(['galaxy']))
        disk_rg = cat_rg.get_objects_by_region(disk,
                                               obj_type_set=set(['galaxy']))
        compare_objects(disk_one, disk_rg)

        rhomb_one = cat_one.get_objects_by_region(rhomb,
                                                  obj_type_set=set(['galaxy']))
        rhomb_rg = cat_rg.get_objects_by_region(rhomb,
                                                obj_type_set=set(['galaxy']))
        compare_objects(rhomb_one, rhomb_rg)


if __name__ == '__main__':
    unittest.main()
