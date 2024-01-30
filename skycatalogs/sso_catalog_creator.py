import os
# import sys
# import re
# import logging
import numpy as np
# import numpy.ma as ma
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


def _partition_mjd(mins, maxes, counts, max_rowgroup_size):
    ''' Determine mjd intervals, each of which will be processed
            in one go'''
    total_rows = np.sum(np.array(counts))
    min_min = int(np.floor(np.min(np.array(mins))))
    max_max = int(np.ceil(np.max(np.array(maxes))))
    # n_rowgroup = (max_max + (max_rowgroup_size - 1) - min_min)/max_rowgroup_size
    n_rowgroup = int(np.ceil(total_rows/max_rowgroup_size))
    # for debugging
    if n_rowgroup < 2:
        n_rowgroup = 2
    mjd_interval = int(np.ceil((max_max - min_min)/n_rowgroup))
    mjds = [min_min]
    last = min_min
    for i in range(n_rowgroup):
        last += mjd_interval
        last = min(last, max_max)
        mjds.append(last)
        if last == max_max:
            break

    return mjds


class SsoCatalogCreator:
    _sso_truth = '/sdf/home/j/jrb/rubin-user/sso/input/19jan2024'
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
        self._df_query = f'''select ObjID as id, {mjd_c} as mjd,
                 "AstRA(deg)" as ra, "AstDec(deg)" as dec, optFilter as filter,
                 observedTrailedSourceMag from {tbl} where mjd >= (?)
                 and mjd < (?) order by mjd'''

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
            # pa.field('ra_rate', pa.float64()),
            # pa.field('dec_rate', pa.float64()),
            pa.field('filter', pa.string())]
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

    def _write_rg(self, writer, min_mjd, max_mjd, db_files, arrow_schema):
        df_list = []
        for f in db_files:
            conn = sqlite3.connect(f)
            df_list.append(pd.read_sql_query(self._df_query,
                                             conn,
                                             params=(min_mjd, max_mjd)))
        df = pd.concat(df_list)
        df_sorted = df.sort_values('mjd')
        tbl = pa.Table.from_pandas(df_sorted, schema=arrow_schema)
        writer.write_table(tbl)

    def create_sso_catalog(self):
        """
        Create the 'main' sso catalog, including everything except fluxes
        """
        #  Find all the db files from Sorcha.   They should all be in a single
        #  directory with no other files in that directory
        files = os.listdir(self._sso_truth)
        db_files = [os.path.join(self._sso_truth, f) for f in files if f.endswith('.db')]
        mins = []
        maxes = []
        counts = []
        for f in db_files:
            with sqlite3.connect(f) as conn:
                res = conn.execute(self._mjd_q)
                r = res.fetchone()
                mins.append(float(r[0]))
                maxes.append(float(r[1]))
                counts.append(int(r[2]))

        mjd_list = _partition_mjd(mins, maxes, counts, self._row_group_size)
        arrow_schema = self._create_main_schema()

        # ## Depending on length of mjd_list, may need more than one
        # ## output file in which case will need to add an outer loop
        # ## over output file. Decide on max row groups per file.
        mjd_min = min(mins)
        mjd_max = max(maxes)
        out_name = f'sso_{int(mjd_min)}_{int(np.ceil(mjd_max))}.parquet'
        writer = pq.ParquetWriter(os.path.join(self._output_dir, out_name),
                                  arrow_schema)
        for i in range(len(mjd_list) - 1):
            self._write_rg(writer, mjd_list[i], mjd_list[i+1], db_files,
                           arrow_schema)

        # close parquet file
        writer.close()

        #  In sso description in config come up with suitable re for filenames

    def _create_sso_flux_file(self, info, arrow_schema):
        '''
        Parameters
        ----------
        info          dict  information pertaining to an existing sso main file
        arrow_schema        to be used in creating the new flux file
        '''
        object_list = self._cat.get_object_type_by_region(None, 'sso',
                                                          mjd=None,
                                                          filepath=info['path'])
        colls = object_list.get_collections()
        outname = f"sso_flux_{info['mjd_min']}_{info['mjd_max']}.parquet"

        writer = pq.ParquetWriter(os.path.join(self._output_dir, outname),
                                  arrow_schema)
        outs = dict()
        for c in colls:
            outs['id'] = c._id
            outs['mjd'] = c._mjds
            all_fluxes = [o.get_LSST_fluxes(as_dict=False) for o in c]
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
        has several observations) and fluxes.   (Or just one flux and
        associated band?)
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
        for f, info in self._cat._sso_files.items():
            if info['scope'] == 'main':
                self._create_sso_flux_file(info, arrow_schema)
