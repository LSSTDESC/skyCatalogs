"""Unit tests for external catalog interface."""
import os
from pathlib import Path
import unittest
from skycatalogs import skyCatalogs


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
        self.assertIn("external_objects", skycat.raw_config['object_types'])
        self.assertIn("external_objects", skycat.cat_cxt._source_type_dict)


if __name__ == "__main__":
    unittest.main()
