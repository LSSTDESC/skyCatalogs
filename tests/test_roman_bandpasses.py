import importlib.util
import unittest
from unittest import mock


@unittest.skipUnless(
    importlib.util.find_spec('roman_technical_information'),
    'roman-technical-information is an optional dependency')
class RomanBandpassTestCase(unittest.TestCase):

    def test_current_effective_area_curves_are_loaded(self):
        from skycatalogs.objects.base_object import ROMAN_BANDS
        from skycatalogs.objects.base_object import load_roman_bandpasses
        from skycatalogs.objects.base_object import \
            load_roman_bandpasses_version

        bandpasses = load_roman_bandpasses()

        self.assertEqual(tuple(bandpasses), ROMAN_BANDS)
        self.assertGreater(bandpasses['K213'].red_limit, 2300.0)
        self.assertLess(bandpasses['K213'].red_limit, 2340.0)
        self.assertIn(
            'roman-technical-information', load_roman_bandpasses_version())

    def test_legacy_names_map_to_all_current_curves(self):
        from skycatalogs.objects.base_object import load_roman_bandpasses

        bandpasses = load_roman_bandpasses(include_all_bands=True)

        self.assertEqual(len(bandpasses), 11)
        self.assertIn('R062', bandpasses)
        self.assertIn('W146', bandpasses)
        self.assertIn('K213', bandpasses)
        self.assertIn('SNPrism', bandpasses)

    def test_missing_optional_package_has_targeted_error(self):
        import builtins
        import skycatalogs.utils.roman_bandpasses as roman_bandpasses

        original_import = builtins.__import__
        cached_bandpasses = roman_bandpasses._bandpasses
        cached_version = roman_bandpasses._version

        def block_roman_import(name, *args, **kwargs):
            if name == 'roman_technical_information':
                raise ModuleNotFoundError(name)
            return original_import(name, *args, **kwargs)

        try:
            roman_bandpasses._bandpasses = None
            roman_bandpasses._version = None
            with mock.patch('builtins.__import__', block_roman_import):
                with self.assertRaisesRegex(
                        ModuleNotFoundError, r'skyCatalogs\[roman\]'):
                    roman_bandpasses.load_roman_bandpasses()
        finally:
            roman_bandpasses._bandpasses = cached_bandpasses
            roman_bandpasses._version = cached_version


if __name__ == '__main__':
    unittest.main()
