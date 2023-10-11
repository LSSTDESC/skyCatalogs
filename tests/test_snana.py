import unittest
from pathlib import Path
import os
import numpy as np

from skycatalogs.skyCatalogs import SkyCatalog, open_catalog
from skycatalogs.skyCatalogs import Box, Disk, PolygonalRegion
from skycatalogs.skyCatalogs import _get_intersecting_hps
from skycatalogs.objects.base_object import BaseObject


PACKAGE_DIR = os.path.dirname(os.path.abspath(str(Path(__file__).parent)))
SKYCATALOG_ROOT = os.path.join(PACKAGE_DIR, "skycatalogs", "data")
CATALOG_DIR = os.path.join(SKYCATALOG_ROOT, "ci_snana_sample")
SKYCATALOG = os.path.join(CATALOG_DIR, "skyCatalog.yaml")

class SnanaTester(unittest.TestCase):
    def setUp(self):
        '''
        Open the catalog
        '''
        # To run as test on github, catalog must be accessible
        self._cat = None
        self._skip = ''
        try:
            self._cat = open_catalog(SKYCATALOG,
                                     skycatalog_root=SKYCATALOG_ROOT)
        except:
            pass

        if not self._cat:
            self._skycatalog_root = '/pscratch/sd/j/jrbogart/desc/skycatalogs'

            # Using parquet files with added MW extinction columns
            cfg_path = os.path.join(self._skycatalog_root, 'MW_columns_v3',
                                    'skyCatalog.yaml')
            try:
                self._cat = open_catalog(cfg_path)
            except:
                self._skip = 'Data source unavailable'

        self._hp = 9043
        ra_min = 52.4
        ra_max = 52.8
        dec_min = -28.4
        dec_max = -28.0
        self._rgn = Box(ra_min, ra_max, dec_min, dec_max)
        self._mjd = 62250.0

    def tearDown(self):
        pass                  # nothing to do


    def test_hp(self):
        '''
        Just ask for all objects in the hp
        '''
        if not self._cat:
            self.skipTest(self._skip)

        hps = self._cat._find_all_hps()
        self.assertEqual(len(hps), 5)

        obj_list = self._cat.get_object_type_by_hp(self._hp, 'snana')
        n_found = len(obj_list)
        self.assertEqual(len(obj_list), 20949)

    def test_mjd(self):
        if not self._cat:
            self.skipTest(self._skip)

        obj_list = self._cat.get_object_type_by_hp(self._hp, 'snana',
                                                   mjd=self._mjd)
        self.assertEqual(len(obj_list), 8795)

        coll = obj_list.get_collections()[0]
        first = True
        for o in coll[0:3]:
            sed_interp = o._linear_interp_SED(self._mjd)
            for band in 'ugrizy':
                flux = o.get_LSST_flux(band, cache=False, mjd=self._mjd)
                if first:
                    self.assertAlmostEqual(flux, 1.2835084851690354e-08)
                    first = False

    def test_region(self):
        if not self._cat:
            self.skipTest(self._skip)

        obj_list = self._cat.get_objects_by_region(self._rgn,
                                                   obj_type_set={'snana'})
        self.assertEqual(len(obj_list), 1345)

    def test_region_and_mjd(self):
        if not self._cat:
            self.skipTest(self._skip)

        obj_list = self._cat.get_objects_by_region(self._rgn,
                                                   obj_type_set={'snana'},
                                                   mjd=self._mjd)
        self.assertEqual(len(obj_list), 548)

if __name__ == '__main__':
    unittest.main()
