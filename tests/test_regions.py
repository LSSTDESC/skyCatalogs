"""
Unit tests for Region classes.
"""

import unittest
import numpy as np
from skycatalogs.utils import Box, Disk, PolygonalRegion


class RegionTestCase(unittest.TestCase):
    """TestCase class for Region classes"""

    def setUp(self):
        # Bounds of a 0.5 x 0.5 degree region, centered on dec = 0.0
        self.ra0 = 40.0
        self.dec0 = 0.0
        self.half_size = 0.25
        self.ra_bounds = self.ra0 - self.half_size, self.ra0 + self.half_size
        self.dec_bounds = self.dec0 - self.half_size, self.dec0 + self.half_size

        # Create an array of ra, dec values inside that region.
        eps = 0.01
        ra_values = np.linspace(self.ra_bounds[0] + eps, self.ra_bounds[1] - eps, 10)
        dec_values = np.linspace(self.dec_bounds[0] + eps, self.dec_bounds[1] - eps, 10)
        ra, dec = np.meshgrid(ra_values, dec_values)
        self.ra_inside = ra.flatten()
        self.dec_inside = dec.flatten()

        # Create an array of ra, dec values just outside of the
        # perimeter of that region.
        nside = 10
        self.ra_outside = np.concatenate(
            [
                ra_values,
                np.ones(nside) * self.ra_bounds[1] + eps,
                ra_values,
                np.ones(nside) * self.ra_bounds[0] - eps,
            ]
        )
        self.dec_outside = np.concatenate(
            [
                np.ones(nside) * self.dec_bounds[0] - eps,
                dec_values,
                np.ones(nside) * self.dec_bounds[1] + eps,
                dec_values,
            ]
        )

    def tearDown(self):
        pass

    def _inclusion_test(self, region):
        """Compute the inclusion mask and apply it to the test grid points."""
        mask = region.compute_mask(self.ra_inside, self.dec_inside)
        ra = np.ma.array(self.ra_inside, mask=mask).compressed()
        dec = np.ma.array(self.dec_inside, mask=mask).compressed()

        # All points should remain.
        np.testing.assert_array_equal(ra, self.ra_inside)
        np.testing.assert_array_equal(dec, self.dec_inside)

    def _exclusion_test(self, region):
        """Compute the inclusion mask and apply it to the test grid points."""
        mask = region.compute_mask(self.ra_outside, self.dec_outside)
        ra = np.ma.array(self.ra_outside, mask=mask).compressed()
        dec = np.ma.array(self.dec_outside, mask=mask).compressed()

        # All points should be excluded.
        self.assertEqual(len(ra), 0)
        self.assertEqual(len(dec), 0)

    def testPolygonalRegion(self):
        """PolygonalRegion test for included and excluded points."""
        vertices = [
            (self.ra_bounds[0], self.dec_bounds[0]),
            (self.ra_bounds[1], self.dec_bounds[0]),
            (self.ra_bounds[1], self.dec_bounds[1]),
            (self.ra_bounds[0], self.dec_bounds[1]),
        ]
        region = PolygonalRegion(vertices_radec=vertices)
        self._inclusion_test(region)
        self._exclusion_test(region)

    def testBox(self):
        """Box region test for included and excluded points."""
        region = Box(*self.ra_bounds, *self.dec_bounds)
        self._inclusion_test(region)
        self._exclusion_test(region)

    def testDisk(self):
        """Disk region test for included and excluded points."""
        # Define a circular region that encloses the inside test points.
        radius = self.half_size * np.sqrt(2.0) * 3600.0  # radius in arcsec
        region = Disk(self.ra0, self.dec0, radius)
        self._inclusion_test(region)

        # Define a circular region that lies within the outside test points.
        radius = self.half_size * 3600.0  # radius in arcsec
        region = Disk(self.ra0, self.dec0, radius)
        self._exclusion_test(region)


if __name__ == "__main__":
    unittest.main()
