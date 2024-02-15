import os
import numpy as np
import sqlite3
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from .objects.base_object import LSST_BANDS

"""
Code for creating sky catalogs for sso objects
"""

__all__ = ['SsoCatalogCreator']

_DEFAULT_ROW_GROUP_SIZE = 100000      # Maybe could be larger


class SsoCatalogCreator:
    _sso_truth = '/sdf/home/j/jrb/rubin-user/sso/input/8feb2024'
    _sso_sed = '/sdf/home/j/jrb/rubin-user/sso/sed/solar_sed.db'

    def __init__(self, catalog_creator, output_dir,
                 sso_truth=None, sso_sed=None):
        '''
        Parameters
        ----------
        catalog_creator   instance of CatalogCreator
        output_dir        destination directory for generated catalogs
        sso_truth         path to input data directory
        sso_sed           path to solar sed
        '''
        self._catalog_creator = catalog_creator
        self._output_dir = catalog_creator._output_dir
        self._logger = catalog_creator._logger
        self._sso_truth = sso_truth
        if sso_truth is None:
            self._sso_truth = SsoCatalogCreator._sso_truth
        self._sso_sed = sso_sed
        if sso_sed is None:           # use default path
            self._sso_sed = SsoCatalogCreator._sso_sed

        self._row_group_size = _DEFAULT_ROW_GROUP_SIZE

        tbl = 'pp_results'
        mjd_c = 'FieldMJD_TAI'
        self._mjd_q = f'select min({mjd_c}), max({mjd_c}), count({mjd_c}) from {tbl}'
        self._dfhp_query = f'''select ObjID as id, {mjd_c} as mjd,
                 "AstRA(deg)" as ra, "AstDec(deg)" as dec,
                 "AstRARate(deg/day)" as ra_rate,
                 "AstDecRate(deg/day)" as dec_rate,
                 observedTrailedSourceMag from {tbl} where healpix = (?)
                 order by mjd, ObjID'''

    @property
    def sso_truth(self):
        return self._sso_truth

    @property
    def sso_sed(self):
        return self._sso_sed

    def _create_main_schema(self):

        fields = [
            pa.field('id', pa.string()),
            pa.field('mjd', pa.float64()),
            pa.field('ra', pa.float64()),
            pa.field('dec', pa.float64()),
            pa.field('observedTrailedSourceMag', pa.float64()),
            pa.field('ra_rate', pa.float64()),
            pa.field('dec_rate', pa.float64())]
        return pa.schema(fields)

    def _create_flux_schema(self):
        # id, mjd and 6 flux fields (for now.  Maybe later also Roman)
        fields = [
            pa.field('id', pa.string()),
            pa.field('mjd', pa.float64()),
            pa.field('lsst_flux_u', pa.float32(), True),
            pa.field('lsst_flux_g', pa.float32(), True),
            pa.field('lsst_flux_r', pa.float32(), True),
            pa.field('lsst_flux_i', pa.float32(), True),
            pa.field('lsst_flux_z', pa.float32(), True),
            pa.field('lsst_flux_y', pa.float32(), True)]
        return pa.schema(fields)

    def _get_hps(self, filepath):

        with sqlite3.connect(filepath) as conn:
            df = pd.read_sql_query('select distinct healpix from pp_results',
                                   conn)
        return set(df['healpix'])

    def _write_hp(self, hp, hps_by_file, arrow_schema):
        df_list = []
        for f in hps_by_file:
            if hp in hps_by_file[f]:
                conn = sqlite3.connect(f)
                one_df = pd.read_sql_query(self._dfhp_query, conn,
                                           params=(hp,))
                df_list.append(one_df)
        if df_list == []:
            return
        writer = pq.ParquetWriter(os.path.join(self._output_dir,
                                               f'sso_{hp}.parquet'),
                                  arrow_schema)

        df = pd.concat(df_list)
        df_sorted = df.sort_values('mjd')

        # Should be prepared to write multiple row groups here
        # depending on # of rows in the table.
        # For now only a single row group will be necessary
        tbl = pa.Table.from_pandas(df_sorted, schema=arrow_schema)
        writer.write_table(tbl)


    def create_sso_catalog(self):
        """
        Create the 'main' sso catalog, including everything except fluxes
        """
        #  Find all the db files from Sorcha.   They should all be in a single
        #  directory with no other files in that directory
        files = os.listdir(self._sso_truth)
        arrow_schema = self._create_main_schema()
        db_files = [os.path.join(self._sso_truth, f) for f in files if f.endswith('.db')]

        all_hps = set()
        hps_by_file = dict()
        for f in db_files:
            hps_by_file[f] = self._get_hps(f)
        all_hps = sorted(all_hps.union(hps_by_file[f]))

        todo = self._catalog_creator._parts
        if len(todo) == 0:
            todo = all_hps
        for h in todo:
            self._write_hp(h, hps_by_file, arrow_schema)

    def _create_sso_flux_pixel(self, pixel, arrow_schema):
        output_filename = f'sso_flux_{pixel}.parquet'
        output_path = os.path.join(self._catalog_creator._output_dir,
                                   output_filename)

        object_list = self._cat.get_object_type_by_hp(pixel, 'sso')
        colls = object_list.get_collections()
        writer = pq.ParquetWriter(output_path, arrow_schema)
        outs = dict()
        for c in colls:
            outs['id'] = c._id
            outs['mjd'] = c._mjds
            all_fluxes = [o.get_LSST_fluxes(as_dict=False, mjd=o._mjd) for o in c]
            all_fluxes_transpose = zip(*all_fluxes)
            for i, band in enumerate(LSST_BANDS):
                v = all_fluxes_transpose.__next__()
                outs[f'lsst_flux_{band}'] = v
            out_df = pd.DataFrame.from_dict(outs)
            out_table = pa.Table.from_pandas(out_df, schema=arrow_schema)
            writer.write_table(out_table)

        writer.close()

    def create_sso_flux_catalog(self):
        """
        Create sso flux catalog.  Includes id, mjd (since given sso normally
        has several observations) and fluxes.
        """
        from .skyCatalogs import open_catalog

        arrow_schema = self._create_flux_schema()
        config_file = self._catalog_creator.write_config(path_only=True)
        self._cat = open_catalog(config_file,
                                 skycatalog_root=self._catalog_creator._skycatalog_root)
        sso_config = self._cat.raw_config['object_types']['sso']
        self._flux_template = sso_config['flux_file_template']
        self._main_template = sso_config['file_template']

        self._logger.info('Creating sso flux files')

        # For each main file make a corresponding flux file
        # for f, info in self._cat._sso_files.items():
        #     if info['scope'] == 'main':
        #         self._create_sso_flux_file(info, arrow_schema)

        for p in self._catalog_creator._parts:
            self._create_sso_flux_pixel(p, arrow_schema)
