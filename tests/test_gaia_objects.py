import os
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
import lsst.daf.butler as daf_butler
from skycatalogs import skyCatalogs
from skycatalogs.objects.gaia_object import GaiaCollection
from skycatalogs.utils import Disk


PACKAGE_DIR = os.path.dirname(os.path.abspath(str(Path(__file__).parent)))
SKYCATALOG_ROOT = os.path.join(PACKAGE_DIR, "skycatalogs", "data")
CATALOG_DIR = os.path.join(PACKAGE_DIR, "skycatalogs", "data", "ci_sample")


CONFIG = {'area_partition': None,
          'butler_parameters':
          {'collections': 'refcats',
           'dstype': 'gaia_dr2_20200414',
           'repo': os.path.join(CATALOG_DIR, 'repo')},
          'id_prefix': 'gaia_dr2_',
          'data_file_type': 'butler_refcat'}


def get_gaia_data(butler_params):
    butler = daf_butler.Butler(butler_params['repo'],
                               collections=[butler_params['collections']])
    refs = set(butler.registry.queryDatasets(butler_params['dstype']))
    return pd.concat(butler.get(_).asAstropy().to_pandas() for _ in refs)


class GaiaObjectTestCase(unittest.TestCase):
    """TestCase class for GaiaObjects"""
    def setUp(self):
        '''
        Open the catalog; establish config
        '''
        skycatalog_root = os.path.join(Path(__file__).resolve().parents[1],
                                       'skycatalogs', 'data')
        self._skycatalog_root = skycatalog_root
        cfg_path = os.path.join(skycatalog_root, 'ci_sample', 'skyCatalog.yaml')
        self._cat = skyCatalogs.open_catalog(cfg_path,
                                             skycatalog_root=skycatalog_root)

        self.df = get_gaia_data(CONFIG['butler_parameters'])
        GaiaCollection.set_config(CONFIG)

    def tearDown(self):
        pass

    def test_proper_motion(self):
        ra, dec = 0, 0
        radius = 1
        region = Disk(ra, dec, radius*3600.0)

        # Use start of proper motion epoch, so ra, dec values should
        # be unchanged from refcat values.
        epoch_a_values = set(self.df['epoch'])
        assert len(epoch_a_values) == 1
        mjd0 = epoch_a_values.pop()
        object_list = GaiaCollection.load_collection(region, self._cat,
                                                     mjd=mjd0)
        for index in np.random.choice(range(len(object_list)), size=20):
            obj = object_list[index]
            gaia_id = obj.id.split('_')[-1]
            df = self.df.query(f"id=={gaia_id}")
            row = df.iloc[0]
            self.assertAlmostEqual(np.degrees(row.coord_ra), obj.ra, places=14)
            self.assertAlmostEqual(np.degrees(row.coord_dec), obj.dec,
                                   places=14)

        self.assertEqual(mjd0, object_list.mjd)

        # Use a plausible LSST observation date, and ensure that
        # calculated proper motion offsets are at least 10% larger
        # than the naive estimates.
        mjd = 60584.  # 2024-10-01 00:00:00
        object_list = GaiaCollection.load_collection(region, self._cat, mjd=mjd)
        for index in np.random.choice(range(len(object_list)), size=20):
            obj = object_list[index]
            gaia_id = obj.id.split('_')[-1]
            df = self.df.query(f"id=={gaia_id}")
            row = df.iloc[0]
            self.assertGreaterEqual(abs(np.radians(obj.ra) - row.coord_ra),
                                    abs(row['pm_ra']*(mjd - mjd0)/365.)/1.1)
            self.assertGreaterEqual(abs(np.radians(obj.dec) - row.coord_dec),
                                    abs(row['pm_dec']*(mjd - mjd0)/365.)/1.1)

        self.assertEqual(mjd, object_list.mjd)


if __name__ == '__main__':
    unittest.main()
