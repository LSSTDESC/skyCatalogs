import unittest
import astropy.units as u
import astropy.constants
import numpy as np
import galsim
from skycatalogs.utils import normalize_sed


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


if __name__ == "__main__":
    unittest.main()
