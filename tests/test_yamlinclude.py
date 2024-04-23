"""
Test yaml loader which handles !include tag
"""
import unittest
import yaml
import os
from pathlib import Path

from skycatalogs.utils.config_utils import YamlIncludeLoader

class skyCatalogsTestCase(unittest.TestCase):
    def setUp(self):
        # Get directory containing top-level file
        self._yaml_dir = os.path.join(Path(__file__).resolve().parents[1],
                                      'skycatalogs', 'data', 'ci_yamlinclude')
    def tearDown(self):
        pass

    def test_include(self):
        # Top level file references to
        #      file in the same directory
        #      file in subdirectory
        #      file which itself has an !include directive
        top_path = os.path.join(self._yaml_dir, 'top.yaml')

        # If the load succeeds, all required files have been found
        with open(top_path) as f:
            data = yaml.load(f, Loader=YamlIncludeLoader)

        # Verify that expected keys exist
        assert('object_types' in data.keys())
        objs = data['object_types']
        assert('gaia_star' in objs)
        assert('galaxy' in objs)
        assert('area_partition' in objs['galaxy'])
        assert('nside' in objs['galaxy']['area_partition'])

        print(data)

if __name__ == '__main__':
    unittest.main()
