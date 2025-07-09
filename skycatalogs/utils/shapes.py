from collections import namedtuple
import numpy as np
import healpy
from astropy import units as u
import healpy
from lsst.sphgeom import Angle, ConvexPolygon, UnitVector3d, LonLat, Circle

__all__ = ['Box', 'Disk', 'PolygonalRegion']


class Region:
    """
    Base class for regions used to make object selections from
    catalogs.
    """
    def get_intersecting_hps(self, nside, hp_ordering):
        nest = (hp_ordering == "NESTED")
        return sorted(self._get_intersecting_hps(nside, nest))

    def _get_intersecting_hps(self, nside, nest):
        raise NotImplementedError

    def get_radec_bounds(self):
        raise NotImplementedError

    def compute_mask(self, ra, dec):
        raise NotImplementedError

    def sphgeom_region(self):
        raise NotImplementedError


class Box(Region):
    """
    Rectangular region in RA, Dec.
    """
    def __init__(self, ra_min, ra_max, dec_min, dec_max):
        self.ra_min = ra_min
        self.ra_max = ra_max
        self.dec_min = dec_min
        self.dec_max = dec_max

    def get_radec_bounds(self):
        """Return the bounds on ra, dec that enclose the region."""
        return self.ra_min, self.ra_max, self.dec_min, self.dec_max

    def _get_intersecting_hps(self, nside, nest):
        vec = healpy.pixelfunc.ang2vec([self.ra_min, self.ra_max,
                                        self.ra_max, self.ra_min],
                                       [self.dec_min, self.dec_min,
                                        self.dec_max, self.dec_max],
                                       lonlat=True)

        return healpy.query_polygon(nside, vec, inclusive=True, nest=nest)

    def compute_mask(self, ra, dec):
        '''
        Compute the mask for excluding entries in the ra, dec arrays.

        Parameters
        ----------
        ra, dec      parallel float arrays, units are degrees. Together
                     they define the list of points to be checked for
                     containment
        '''
        mask = np.logical_or((ra < self.ra_min),
                             (ra > self.ra_max))
        mask = np.logical_or(mask, (dec < self.dec_min))
        mask = np.logical_or(mask, (dec > self.dec_max))
        return mask


class Disk(Region):
    """
    Circular region on the sky.
    """
    def __init__(self, ra, dec, radius_as):
        self.ra = ra
        self.dec = dec
        self.radius = radius_as*u.arcsec

    def _get_intersecting_hps(self, nside, nest):
        center = healpy.pixelfunc.ang2vec(self.ra, self.dec, lonlat=True)
        return healpy.query_disc(nside, center, self.radius.to_value('radian'),
                                 inclusive=True, nest=nest)

    def get_radec_bounds(self):
        """Return the bounds on ra, dec that enclose the region."""
        radius = self.radius.to_value("degree")
        dec_min = self.dec - radius
        dec_max = self.dec + radius
        cos_dec_max = np.cos(np.radians(dec_max))
        ra_min = self.ra - radius/cos_dec_max
        ra_max = self.ra + radius/cos_dec_max
        return ra_min, ra_max, dec_min, dec_max

    def compute_mask(self, ra, dec):
        '''
        Compute the mask for excluding entries in the ra, dec arrays.

        Parameters
        ----------
        ra, dec      parallel float arrays, units are degrees. Together
                     they define the list of points to be checked for
                     containment
        '''
        # Use healpy rather than lsst.sphgeom because healpy takes
        # array inputs
        p_vec = healpy.pixelfunc.ang2vec(ra, dec, lonlat=True)

        c_vec = healpy.pixelfunc.ang2vec(self.ra,
                                         self.dec,
                                         lonlat=True)

        # Rather than comparing arcs, it is equivalent to compare chords
        # (or square of chord length)
        obj_chord_sq = np.sum(np.square(p_vec - c_vec), axis=1)

        # This is to be compared to square of chord for angle a corresponding
        # to disk radius.  That's 4(sin(a/2)^2)
        radius_rad = self.radius.to_value('radian')
        rad_chord_sq = 4 * np.square(np.sin(0.5 * radius_rad))
        return obj_chord_sq > rad_chord_sq

    def sphgeom_region(self):
        """Enclosing region expressed as lsst.sphgeom.Circle."""
        center = LonLat.fromDegrees(self.ra, self.dec)
        radius = Angle(self.radius.to_value("radian"))
        return Circle(UnitVector3d(center), radius)


class PolygonalRegion(Region):
    """
    Convex polygon region defined by a set of vertex positions on the sky.
    """
    def __init__(self, vertices_radec=None, convex_polygon=None):
        '''
        Supply either an object of type lsst.sphgeom.ConvexPolygon
        or a list of vertices, each a tuple (ra, dec) in degrees,
        which describe a convex polygon
        '''
        if convex_polygon:
            if isinstance(convex_polygon, ConvexPolygon):
                self._convex_polygon = convex_polygon
                return
        if vertices_radec:
            if not isinstance(vertices_radec, list):
                raise TypeError(f'PolygonalRegion: Argument {vertices_radec} '
                                'is not a list')
            vertices = [UnitVector3d(LonLat.fromDegrees(v_rd[0], v_rd[1]))
                        for v_rd in vertices_radec]
            self._convex_polygon = ConvexPolygon(vertices)
            return
        raise ValueError('PolygonalRegion: Either vertices_radec or '
                         'convex_polygon must have an acceptable value')

    def get_radec_bounds(self):
        """Return the bounds on ra, dec that enclose the region."""
        ra_vals, dec_vals = np.array(self.get_vertices_radec()).T
        return min(ra_vals), max(ra_vals), min(dec_vals), max(dec_vals)

    def get_vertices(self):
        '''
        Return vertices as list of 3d vectors
        '''
        return self._convex_polygon.getVertices()

    def get_vertices_radec(self):
        v3d = self.get_vertices()
        vertices_radec = []
        for v in v3d:
            vertices_radec.append((LonLat.longitudeOf(v).asDegrees(),
                                   LonLat.latitudeOf(v).asDegrees()))
        return vertices_radec

    def _get_bounding_disk(self):
        circle = self._convex_polygon.getBoundingCircle()
        center = circle.getCenter()
        ra_c = LonLat.longitudeOf(center).asDegrees()
        dec_c = LonLat.latitudeOf(center).asDegrees()
        # The opening angle seems a bit small, so include a 5% safety factor.
        rad_as = circle.getOpeningAngle().asDegrees() * 3600 * 1.05
        return Disk(ra_c, dec_c, rad_as)

    def compute_mask(self, ra, dec):
        '''
        Compute the mask for excluding entries in the ra, dec arrays.

        Parameters
        ----------
        ra, dec      parallel float arrays, units are degrees. Together
                     they define the list of points to be checked for
                     containment
        '''
        # Pre-filter using Disk region to do the bulk of the masking quickly.
        disk_region = self._get_bounding_disk()
        mask = disk_region.compute_mask(ra, dec)
        if all(mask):
            # Everything is masked, so no need to refine with convex polygon.
            return mask

        # Create masked arrays for the remaining positions.
        ra_compress = np.radians(np.ma.array(ra, mask=mask).compressed())
        dec_compress = np.radians(np.ma.array(dec, mask=mask).compressed())
        polygon_mask = self._convex_polygon.contains(ra_compress, dec_compress)

        # Indexes of unmasked by Disk.
        ixes = np.where(~mask)

        # Apply polygon region mask.
        mask[ixes] |= ~polygon_mask

        return mask

    def _get_intersecting_hps(self, nside, nest):
        return healpy.query_polygon(nside, self.get_vertices(),
                                    inclusive=True, nest=nest)

    def sphgeom_region(self):
        """Enclosing region expressed as lsst.sphgeom.ConvexPolygon."""
        return self._convex_polygon
