"""
Example unit tests for skyCatalogs package
"""
import unittest
import desc.skycatalogs

class skyCatalogsTestCase(unittest.TestCase):
    def setUp(self):
        self.message = 'Hello, world'

    def tearDown(self):
        pass

    #def test_run(self):
    #    foo = desc.skycatalogs.skyCatalogs(self.message)
    #    self.assertEqual(foo.run(), self.message)

    def test_failure(self):
        #self.assertRaises(TypeError, desc.skycatalogs.skyCatalogs)
        #foo = desc.skycatalogs.skyCatalogs(self.message)
        #self.assertRaises(RuntimeError, foo.run, True)
        pass

if __name__ == '__main__':
    unittest.main()
