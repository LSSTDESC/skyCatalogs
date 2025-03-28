"""Unit tests for external catalog interface."""

import os
from pathlib import Path
import unittest
from skycatalogs import skyCatalogs
from skycatalogs.utils import Disk


PACKAGE_DIR = os.path.dirname(os.path.abspath(str(Path(__file__).parent)))


class ExternalCatalogsTestCase(unittest.TestCase):
    """
    TestCase class for third party catalog interface.
    """

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_skycatalog_yaml_loading(self):
        """Test loading of object types via skyCatalogs yaml config file."""
        yaml_file = os.path.join(PACKAGE_DIR, "tests", "external_skycat.yaml")
        skycat = skyCatalogs.open_catalog(yaml_file, skycatalog_root=".")
        region = Disk(0, 0, 3600.0)
        object_params = {"object_type_1": 111.0, "object_type_2": 2.1}
        for object_type in ("object_type_1", "object_type_2"):
            # Check that object types are registered by name.
            self.assertIn(object_type, skycat.raw_config["object_types"])
            self.assertIn(object_type, skycat.cat_cxt._source_type_dict)
            # Check that the object type parameters are set in the collection class.
            foo = skycat.get_object_type_by_region(region, object_type)
            self.assertEqual(foo[0].belongs_to._param, object_params[object_type])


if __name__ == "__main__":
    unittest.main()
