from collections import namedtuple
import numpy as np
from astropy import units as u
import healpy
from lsst.sphgeom import ConvexPolygon, UnitVector3d, LonLat

__all__ = ['Box', 'Disk', 'PolygonalRegion', 'compute_region_mask']

Box = namedtuple('Box', ['ra_min', 'ra_max', 'dec_min', 'dec_max'])

# radius is measured in arcseconds
Disk = namedtuple('Disk', ['ra', 'dec', 'radius_as'])


class PolygonalRegion:

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
                raise TypeError(f'PolygonalRegion: Argument {vertices_radec} is not a list')
            vertices = [UnitVector3d(LonLat.fromDegrees(v_rd[0], v_rd[1])) for v_rd in vertices_radec]
            self._convex_polygon = ConvexPolygon(vertices)
            return
        raise ValueError('PolygonalRegion: Either vertices_radec or convex_polygon must have an acceptable value')

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

    def get_containment_mask(self, ra, dec, included=True):
        '''
        Parameters
        ----------
        ra, dec      parallel float arrays, units are degrees. Together
                     they define the list of points to be checked for
                     containment
        included     boolean   If true, mask bit will be set to True for
                               contained points, else False.    Reverse
                               the settings if included is False
        '''
        # convert to radians
        ra = [(r * u.degree).to_value(u.radian) for r in ra]
        dec = [(d * u.degree).to_value(u.radian) for d in dec]

        mask = self._convex_polygon.contains(ra, dec)
        if included:
            return mask
        else:
            return np.logical_not(mask)


def compute_region_mask(region, ra, dec):
    '''
    Compute mask according to region for provided data
    Parameters
    ----------
    region         Supported shape (box, disk, PolygonalRegion)  or None
    ra,dec         Coordinates for data to be masked, in degrees
    Returns
    -------
    mask of elements in ra, dec arrays to be omitted

    '''
    mask = None
    if isinstance(region, Box):
        mask = np.logical_or((ra < region.ra_min),
                             (ra > region.ra_max))
        mask = np.logical_or(mask, (dec < region.dec_min))
        mask = np.logical_or(mask, (dec > region.dec_max))
    if isinstance(region, Disk):
        # Use healpy rather than lsst.sphgeom because healpy takes
        # array inputs
        p_vec = healpy.pixelfunc.ang2vec(ra, dec, lonlat=True)

        c_vec = healpy.pixelfunc.ang2vec(region.ra,
                                         region.dec,
                                         lonlat=True)
        radius_rad = (region.radius_as * u.arcsec).to_value('radian')

        # Rather than comparing arcs, it is equivalent to compare chords
        # (or square of chord length)
        obj_chord_sq = np.sum(np.square(p_vec - c_vec), axis=1)

        # This is to be compared to square of chord for angle a corresponding
        # to disk radius.  That's 4(sin(a/2)^2)
        rad_chord_sq = 4 * np.square(np.sin(0.5 * radius_rad))
        mask = obj_chord_sq > rad_chord_sq
    if isinstance(region, PolygonalRegion):
        mask = region.get_containment_mask(ra, dec, included=False)
    return mask
