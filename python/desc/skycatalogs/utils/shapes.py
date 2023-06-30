from collections import namedtuple
import numpy as np
from astropy import units as u


__all__ = ['Box', 'Disk', 'PolygonalRegion']

Box = namedtuple('Box', ['ra_min', 'ra_max', 'dec_min', 'dec_max'])

# radius is measured in arcseconds
Disk = namedtuple('Disk', ['ra', 'dec', 'radius_as'])

from lsst.sphgeom import ConvexPolygon, UnitVector3d, LonLat
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
