import os
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
from skycatalogs import skyCatalogs
from galsim import Convolution,RandomKnots


PACKAGE_DIR = os.path.dirname(os.path.abspath(str(Path(__file__).parent)))
SKYCATALOG_ROOT = os.path.join(PACKAGE_DIR, "skycatalogs", "data")
CATALOG_DIR = os.path.join(PACKAGE_DIR, "skycatalogs", "data", "ci_sample")
SKYCATALOG = os.path.join(CATALOG_DIR, "skyCatalog.yaml")


class SkyCatalogsObjecInterfaceTestCase(unittest.TestCase):
    """
    TestCase class for skyCatalogs object interfaces.
    """

    def setUp(self):
        skycat = skyCatalogs.open_catalog(SKYCATALOG, skycatalog_root=SKYCATALOG_ROOT)
        hp = 9557
        object_type = "galaxy"
        self.objects = skycat.get_object_type_by_hp(hp, object_type)
        self.df = pd.read_parquet(os.path.join(CATALOG_DIR, f"galaxy_{hp}.parquet"))
        self.indexes = np.random.choice(range(len(self.objects)), size=100)

    def test_get_wl_params(self):
        """
        Test that weak lensing parameters match values computed from data
        retrieved directly from the galaxy parquet file.
        """
        for index in self.indexes:
            obj = self.objects[index]
            galaxy_id = obj.get_native_attribute("galaxy_id")
            row = self.df.query(f"galaxy_id == {galaxy_id}").iloc[0]
            g1, g2, mu = obj.get_wl_params()
            gamma1 = row["shear_1"]
            gamma2 = row["shear_2"]
            kappa = row["convergence"]
            self.assertAlmostEqual(g1, gamma1 / (1.0 - kappa))
            self.assertAlmostEqual(g2, gamma2 / (1.0 - kappa))
            self.assertAlmostEqual(
                mu, 1.0 / ((1.0 - kappa) ** 2 - (gamma1**2 + gamma2**2))
            )

    def test_get_dust(self):
        """
        Test that dust parameters match values retrieved directly
        from the galaxy parquet file.
        """
        for index in self.indexes:
            obj = self.objects[index]
            galaxy_id = obj.get_native_attribute("galaxy_id")
            row = self.df.query(f"galaxy_id == {galaxy_id}").iloc[0]
            iAv, iRv, gAv, gRv = obj._get_dust()
            # For galaxies, we use the SED values that have internal
            # extinction included, so should have iAv=0, iRv=1.
            self.assertEqual(iAv, 0)
            self.assertEqual(iRv, 1)
            self.assertEqual(gAv, row["MW_av"])
            self.assertEqual(gRv, row["MW_rv"])

    def test_get_gsobject_components(self):
        """Check some properties of the objects returned by the
        .get_gsobject_components function."""
        for index in self.indexes:
            obj = self.objects[index]
            galaxy_id = obj.get_native_attribute("galaxy_id")
            row = self.df.query(f"galaxy_id == {galaxy_id}").iloc[0]
            gs_objs = obj.get_gsobject_components()
            for component, gs_obj in gs_objs.items():
                if component in "disk bulge":
                    # Check sersic index
                    self.assertEqual(gs_obj.original.n, row[f"sersic_{component}"])
                elif component == "knots":
                    # Check number of knots
                    if type==Convolution:
                        if type(gs_obj.obj_list[0].original)==RandomKnots:
                            self.assertEqual(gs_obj.obj_list[0].original.npoints, row["n_knots"])
                        else:
                            self.assertEqual(gs_obj.obj_list[1].original.npoints, row["n_knots"])
                    else:
                        try:
                            self.assertEqual(gs_obj.original.npoints, row["n_knots"])

    def test_get_sed_components(self):
        """Check sed components."""
        for index in self.indexes:
            obj = self.objects[index]
            galaxy_id = obj.get_native_attribute("galaxy_id")
            row = self.df.query(f"galaxy_id == {galaxy_id}").iloc[0]
            seds = obj.get_observer_sed_components()
            for component, sed in seds.items():
                if sed is not None:
                    self.assertEqual(sed.redshift, row["redshift"])


if __name__ == "__main__":
    unittest.main()
