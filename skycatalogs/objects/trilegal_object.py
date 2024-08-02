import os
import itertools
from collections.abc import Iterable
import healpy
import astropy.units as u
import numpy as np
import numpy.ma as ma
import pandas as pd
import galsim

from skycatalogs.utils.shapes import Disk, PolygonalRegion, compute_region_mask
from skycatalogs.objects.base_object import BaseObject, ObjectCollection
from skycatalogs.objects.base_object import ObjectList

__all__ = ['TrilegalObject', 'TrilegalCollection']


def _get_disk_pixels(disk, nside):
    center = healpy.pixelfunc.ang2vec(disk.ra, disk.dec, lonlat=True)
    radius_rad = (disk.radius_as * u.arcsec).to_value('radian')

    pixels = healpy.query_disc(nside, center, radius_rad, inclusive=True,
                               nest=False)
    return pixels


def _execute_cone_query(disk, logger, hp_id=None):
    """
    Retrieve data for all stars in the disk (optionally only those in
    specified healpixel)
    disk     skycatalogs.shapes.Disk
    hp_id    if not None, restrict to objects in the specified healpixel
    """
    from dl import queryClient as qc

    nside = 32
    pixels = _get_disk_pixels(disk, nside)
    if hp_id:
        if hp_id not in pixels:
            logger.info(f'Healpixel {hp_id} does not intersect region')
            return None
        factor = 2**(-14)  # should be computed from config
        # also nest4096 below is column name available from config
        p_nest = healpy.ring2nest(nside, hp_id)
        hp_clause = 'FLOOR(nest4096 * {factor}) = {p_nest} and '
    else:
        hp_clause = ''

    radius_deg = disk.radius_as / 3600.0
    q = 'select ra, dec, random_id, nest4096 from lsst_sim.simdr2 where '
    q += hp_clause
    q += f'q3c_radial_query(ra, dec, {disk.ra}, {disk.dec}, {radius_deg})'
    results = qc.query(adql=q, fmt='pandas')

    return results


class TrilegalObject(BaseObject):
    """
    Object class for Trilegal simulated stars, read from Datalab db
    """

    def __init__(self, obj_pars, parent_collection, index):
        ra = obj_pars['ra']
        dec = obj_pars['dec']
        id = f"{TrilegalCollection._id_prefix}{obj_pars['id']}"

        super().__init__(ra, dec, id, parent_collection._object_type_unique,
                         belongs_to=parent_collection, belongs_index=index)

    def get_observer_sed_component(self, component, mjd=None):
        # No way to do this yet
        pass

    def get_gsobject_components(self, gsparams=None, rng=None):
        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)
        return {'this_object': galsim.DeltaFunction(gsparams=gsparams)}


class TrilegalCollection(ObjectCollection):
    # Class methods
    _trilegal_config = None

    @classmethod
    def set_config(cls, config):
        TrilegalCollection._trilegal_config = config
        TrilegalCollection._id_prefix = config['id_prefix']

    @classmethod
    def get_config(cls):
        return TrilegalCollection._trilegal_config

    @staticmethod
    def load_collection(region, skycatalog, mjd=None):
        '''
        region      One of Disk, PolygonalRegion from skyCatalogs.utils.shapes.
                    Box is not currently supported
        skycatalog  An instance of the SkyCatalog class
        mjd         Time at which objects are to be assembled. Ignored for
                    Gaia stars

        returns     ObjectList
        '''

        if not skycatalog:
            raise ValueError('GaiaCollection.load_collection: skycatalog cannot be None')
        if TrilegalCollection.get_config() is None:
            config = skycatalog.raw_config['object_types']['trilegal_star']
            TrilegalCollection.set_config(config)
        trilegal_section = TrilegalCollection.get_config()

        nside = 32
        need_mask = False
        if isinstance(region, Disk):
            disk = region
        elif isinstance(region, PolygonalRegion):
            need_mask = True
            # disk = circumscribed disk
            skycatalog._logger.error('Polygonal region not supported for Trilegal')
            return None
        else:
            skycatalog._logger.error('Unsupported region type')
            return None

        o_list = ObjectList()

        pixels = _get_disk_pixels(disk, nside)
        for p in pixels:
            df = _execute_cone_query(disk, skycatalog._logger, p)
            # if need_mask, compute here and prune df
            coll = None
            if df:
                coll = TrilegalCollection(df, skycatalog, 'trilegal_star', mjd,
                                          p)
            if coll:
                o_list.append_collection(coll)

        return o_list

    def __init__(self, df, skycatalog, source_type, mjd, pixel):
        self.df = df
        self._skycatalog = skycatalog
        self._partition_id = pixel
        id_prefix = self._config['id_prefix']
        self._id = np.array([f"{id_prefix}{df.iloc[key]['random_id']}" for key in range(len(df))])
        self._mask = None
        self._object_type_unique = source_type
        self._rdrs = []
        self._object_class = TrilegalObject
        self._mjd = mjd

    @property
    def native_columns(self):
        return set()

    @property
    def mjd(self):
        return self._mjd

    def __getitem__(self, key):
        if isinstance(key, int) or isinstance(key, np.int64):
            row = {col: self.df[col][key] for col in ('id', 'ra', 'dec')}
            return TrilegalObject(row, self, key)

        elif type(key) == slice:
            ixdata = [i for i in range(min(key.stop, len(self._id)))]
            ixes = itertools.islice(ixdata, key.start, key.stop, key.step)
            return [self.__getitem__(i) for i in ixes]

        elif type(key) == tuple and isinstance(key[0], Iterable):
            return [self.__getitem__(i) for i in key[0]]

    def __len__(self):
        return len(self.df)
