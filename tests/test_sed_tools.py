import unittest
from collections import namedtuple
import astropy.units as u
import astropy.constants
import numpy as np
import galsim
from skycatalogs.utils import normalize_sed
from skycatalogs.utils.diffsky_sed import create_diffsky_sed_factory
from skycatalogs.utils.sed_tools import DiffskySedFactory


class NormalizeSedTestCase(unittest.TestCase):
    """TestCase class for normalizing SEDs with magnorm."""

    def setUp(self):
        """No set up needed."""
        pass

    def tearDown(self):
        """Nothing to tear down."""
        pass

    def test_normalize_sed(self):
        """Test the normalize_sed function"""
        sed = galsim.SED(lambda wl: 1, "nm", "flambda", fast=False)
        wl = 400 * u.nm

        # Check that the ratio of differently normalized SEDs have the
        # expected magnitude difference at the reference wavelength:
        magnorm0 = 20.0
        magnorm1 = 25.0
        sed0 = normalize_sed(sed, magnorm0, wl=wl)
        sed1 = normalize_sed(sed, magnorm1, wl=wl)

        self.assertAlmostEqual(sed0(wl) / sed1(wl), 10 ** ((magnorm1 - magnorm0) / 2.5))

        # Check the implied magnitude of the renormalized SED at the
        # reference wavelength against the magnorm value.
        hnu = (astropy.constants.h * astropy.constants.c / wl).to_value(u.erg)
        fnu = sed0(wl) * hnu / (astropy.constants.c / wl**2).to_value(u.Hz / u.nm)
        mag = -2.5 * np.log10(fnu) - 48.60

        self.assertAlmostEqual(mag, magnorm0)


class DiffskySspSliceTestCase(unittest.TestCase):
    """Test optional restriction of the Diffsky SSP wavelength grid."""

    @staticmethod
    def make_factory(lower=0.755, upper=1.855):
        return DiffskySedFactory(
            '.', {'H0': 70, 'Om0': 0.3},
            ssp_wave_min_micron=lower,
            ssp_wave_max_micron=upper)

    def test_slice_updates_wave_and_flux_together(self):
        ssp_type = namedtuple('SSPData', ('ssp_wave', 'ssp_flux', 'other'))
        wave = np.array([7000., 7550., 10000., 18550., 19000.])
        flux = np.arange(2 * 3 * wave.size).reshape(2, 3, wave.size)
        factory = self.make_factory()
        factory._aux_data = {
            'ssp_data': ssp_type(wave, flux, 'preserved'),
        }

        factory._slice_ssp_wavelengths()

        sliced = factory._aux_data['ssp_data']
        np.testing.assert_array_equal(
            sliced.ssp_wave, [7000., 7550., 10000., 18550., 19000.])
        np.testing.assert_array_equal(sliced.ssp_flux, flux)
        self.assertEqual(sliced.other, 'preserved')

    def test_bounds_must_be_supplied_as_a_valid_pair(self):
        with self.assertRaises(ValueError):
            self.make_factory(lower=0.755, upper=None)
        with self.assertRaises(ValueError):
            self.make_factory(lower=1.855, upper=0.755)

    def test_default_bounds_cover_rubin_and_roman_through_z4(self):
        factory = DiffskySedFactory('.', {'H0': 70, 'Om0': 0.3})
        self.assertEqual(factory._ssp_wave_bounds_micron, (0.06, 2.34))

    @staticmethod
    def make_ssp_data(wave):
        ssp_type = namedtuple('SSPData', ('ssp_wave', 'ssp_flux'))
        flux = np.ones((2, 3, len(wave)), dtype=float)
        return ssp_type(np.asarray(wave, dtype=float), flux)

    def test_no_thinning_retains_full_sliced_ssp_grid(self):
        wave = np.array([3600., 3700., 3726., 3729., 3800., 5000.])
        factory = DiffskySedFactory(
            '.', {'H0': 70, 'Om0': 0.3}, thinning_mode='none')
        factory._aux_data = {'ssp_data': self.make_ssp_data(wave)}

        factory._set_wavelengths()

        np.testing.assert_array_equal(factory._wave_indices, np.arange(6))
        np.testing.assert_array_equal(factory.wave_list, wave)

    def test_emission_line_mode_unions_native_line_samples(self):
        wave = np.array([
            4800., 4900., 4950., 4958.91, 4970., 4990.,
            5000., 5006.84, 5020., 5100.,
        ])
        factory = DiffskySedFactory(
            '.', {'H0': 70, 'Om0': 0.3},
            thinning_mode='emission_lines',
            emission_lines_angstrom=[4958.91, 5006.84],
            emission_line_half_width_angstrom=5.0)
        factory._aux_data = {'ssp_data': self.make_ssp_data(wave)}
        factory._thinned_wavelengths = lambda *args: (
            np.array([0, len(wave) - 1]), wave[[0, -1]])

        factory._set_wavelengths()

        # Each line sample and its immediate native-grid neighbors survive.
        np.testing.assert_array_equal(
            factory._wave_indices, [0, 2, 3, 4, 6, 7, 8, 9])

    def test_custom_thinning_accepts_indices_or_boolean_mask(self):
        wave = np.array([1000., 2000., 3000., 4000.])
        for selector, expected in (
                (lambda w, f: [0, 2, 3], [0, 2, 3]),
                (lambda w, f: np.array([True, False, True, True]),
                 [0, 2, 3])):
            factory = DiffskySedFactory(
                '.', {'H0': 70, 'Om0': 0.3},
                thinning_mode='custom', thinning_selector=selector)
            factory._aux_data = {'ssp_data': self.make_ssp_data(wave)}
            factory._set_wavelengths()
            np.testing.assert_array_equal(factory._wave_indices, expected)

    def test_custom_thinning_requires_valid_callable_output(self):
        with self.assertRaises(ValueError):
            DiffskySedFactory(
                '.', {'H0': 70, 'Om0': 0.3}, thinning_mode='custom')
        factory = DiffskySedFactory(
            '.', {'H0': 70, 'Om0': 0.3}, thinning_mode='custom',
            thinning_selector=lambda wave, flux: [0])
        factory._aux_data = {
            'ssp_data': self.make_ssp_data([1000., 2000., 3000.])}
        with self.assertRaises(ValueError):
            factory._set_wavelengths()

    def test_catalog_config_selects_unthinned_factory(self):
        factory = create_diffsky_sed_factory(
            {'sed_thinning_mode': 'none'}, '.',
            {'H0': 70, 'Om0': 0.3})
        self.assertEqual(factory.thinning_mode, 'none')

    def test_sed_engine_and_precision_are_validated(self):
        with self.assertRaises(ValueError):
            DiffskySedFactory(
                '.', {'H0': 70, 'Om0': 0.3}, sed_engine='unknown')
        with self.assertRaises(ValueError):
            DiffskySedFactory(
                '.', {'H0': 70, 'Om0': 0.3}, sed_precision='float16')
        with self.assertRaises(ValueError):
            DiffskySedFactory(
                '.', {'H0': 70, 'Om0': 0.3}, sed_engine='reference',
                sed_precision='float32')
        factory = DiffskySedFactory('.', {'H0': 70, 'Om0': 0.3})
        self.assertEqual(factory.sed_engine, 'fast')
        self.assertEqual(factory.sed_precision, 'float32')
        factory.set_compute_options('reference')
        self.assertEqual(factory.sed_engine, 'reference')
        self.assertEqual(factory.sed_precision, 'float64')
        factory._sed_cache[1] = object()
        with self.assertRaises(RuntimeError):
            factory.set_compute_options('fast')

    def test_fast_default_caches_component_seds_as_float32(self):
        class FakeCatalog:
            def select(self, name):
                self.name = name
                return self

            def get_data(self, format):
                return np.array([17], dtype=np.int64)

        factory = DiffskySedFactory('.', {'H0': 70, 'Om0': 0.3})
        factory._wave_indices = np.array([0, 2])
        factory._wave_list = np.array([1000.0, 3000.0])
        sed_info = {
            name: np.ones((1, 3), dtype=np.float64)
            for name in ('rest_sed_bulge', 'rest_sed_disk',
                         'rest_sed_knots')
        }

        factory._cache_sed_info(FakeCatalog(), sed_info)

        wave, components = factory._sed_cache[17]
        self.assertEqual(wave.dtype, np.float64)
        self.assertTrue(all(item.dtype == np.float32 for item in components))

        factory.clear_sed_cache()
        self.assertEqual(len(factory._sed_cache), 0)

    def test_cache_sed_info_excludes_unrequested_linked_hosts(self):
        class FakeCatalog:
            def select(self, name):
                return self

            def get_data(self, format):
                return np.array([101, 202, 303], dtype=np.int64)

        factory = DiffskySedFactory('.', {'H0': 70, 'Om0': 0.3})
        factory._wave_indices = np.array([0, 1])
        factory._wave_list = np.array([1000.0, 2000.0])
        sed_info = {
            name: np.repeat(
                np.arange(3, dtype=float)[:, np.newaxis], 2, axis=1)
            for name in ('rest_sed_bulge', 'rest_sed_disk',
                         'rest_sed_knots')
        }

        factory._cache_sed_info(
            FakeCatalog(), sed_info, requested_ids=[303, 101])

        self.assertEqual(set(factory._sed_cache), {101, 303})
        np.testing.assert_array_equal(
            factory._sed_cache[303][1][0], [2.0, 2.0])
        np.testing.assert_array_equal(
            factory._sed_cache[101][1][0], [0.0, 0.0])

    def test_runtime_cache_cleanup_releases_seds_and_pixels(self):
        factory = DiffskySedFactory('.', {'H0': 70, 'Om0': 0.3})
        factory._sed_cache[1] = object()
        factory._pixel_cache[0] = object()

        factory.clear_runtime_caches()

        self.assertEqual(len(factory._sed_cache), 0)
        self.assertEqual(len(factory._pixel_cache), 0)


if __name__ == "__main__":
    unittest.main()
