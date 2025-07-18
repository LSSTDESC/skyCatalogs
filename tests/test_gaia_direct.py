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

BUTLER_PARAMETERS = {'collections': 'refcats',
                     'dstype': 'gaia_dr2_20200414',
                     'repo': os.path.join(CATALOG_DIR, 'repo')}


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
        skycatalog_root = SKYCATALOG_ROOT
        cfg_path = os.path.join(skycatalog_root, 'ci_sample',
                                'gaia_skyCatalog.yaml')
        self._cat = skyCatalogs.open_catalog(cfg_path,
                                             skycatalog_root=skycatalog_root)

        self.df = get_gaia_data(BUTLER_PARAMETERS)

    def tearDown(self):
        pass

    def test_proper_motion(self):
        ra, dec = 0, 0
        radius = 1
        region = Disk(ra, dec, radius*3600.0)

        # make a quadrilateral in case we want to try out
        # masking for a polygonal region
        # ra_min = -0.7
        # ra_max =  0.7
        # dec_min = -0.8
        # dec_max = 0.6

        # vertices = [(ra_min, dec_min),
        #             (ra_min - 0.1, dec_max),
        #             (ra_max, dec_max),
        #             (ra_max, dec_min)]
        # rgn_poly = skyCatalogs.PolygonalRegion(vertices_radec=vertices)

        # Use start of proper motion epoch, so ra, dec values should
        # be unchanged from refcat values.
        epoch_a_values = set(self.df['epoch'])
        assert len(epoch_a_values) == 1
        mjd0 = epoch_a_values.pop()
        object_list = GaiaCollection.load_collection(region, self._cat,
                                                     mjd=mjd0)
        # First try a particular case
        # fixed = 2831
        # the_list.insert(0, fixed)
        # for index in np.random.choice(range(len(object_list)), size=20):
        the_list = list(np.random.choice(range(len(object_list)), size=20))
        for index in the_list:
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
        object_list = GaiaCollection.load_collection(region, self._cat,
                                                     mjd=mjd)
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
