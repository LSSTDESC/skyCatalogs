import os
import sys
from multiprocessing import Process, Pipe
import sqlite3
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import json
from .objects.base_object import LSST_BANDS
from .objects.sso_object import SsoConfigFragment
from .utils.config_utils import assemble_provenance
from .utils.config_utils import assemble_file_metadata


"""
Code for creating sky catalogs for sso objects
"""

__all__ = ['SsoMainCatalogCreator', 'SsoFluxCatalogCreator']

_DEFAULT_ROW_GROUP_SIZE = 100000      # Maybe could be larger


def _do_sso_flux_chunk(send_conn, sso_collection, instrument_needed,
                       l_bnd, u_bnd):
    '''
    end_conn         output connection
    star_collection  information from main file
    instrument_needed List of which calculations should be done
    l_bnd, u_bnd     demarcates slice to process

    returns
                    dict with keys id, lsst_flux_u, ... lsst_flux_y
    '''
    out_dict = {}

    o_list = sso_collection[l_bnd: u_bnd]
    out_dict['id'] = list(sso_collection._id[l_bnd: u_bnd])
    out_dict['mjd'] = list(sso_collection._mjds[l_bnd: u_bnd])
    if 'lsst' in instrument_needed:
        all_fluxes = [o.get_LSST_fluxes(as_dict=False, mjd=o._mjd) for o in o_list]
        all_fluxes_transpose = zip(*all_fluxes)
        for i, band in enumerate(LSST_BANDS):
            v = all_fluxes_transpose.__next__()
            out_dict[f'lsst_flux_{band}'] = v

    # if 'roman' in instrument_needed:
    #     all_fluxes = [o.get_roman_fluxes(as_dict=False) for o in o_list]
    #     all_fluxes_transpose = zip(*all_fluxes)
    #     for i, band in enumerate(ROMAN_BANDS):
    #         v = all_fluxes_transpose.__next__()
    #         out_dict[f'roman_flux_{band}'] = v

    if send_conn:
        send_conn.send(out_dict)
    else:
        return out_dict


class SsoMainCatalogCreator:
    # _sso_truth = '/sdf/home/j/jrb/rubin-user/sso/input/20feb2024'
    _sso_truth = '/sdf/data/rubin/shared/ops-rehearsals/ops-rehearsal-4/imSim_catalogs/inputs/sso'
    _sso_db_tbl = 'results'

    def __init__(self, catalog_creator, sso_truth=None, sso_db_tbl=None):
        '''
        Parameters
        ----------
        catalog_creator   instance of MainCatalogCreator
        sso_truth         path to input data directory
        sso_db_tbl        table to read from sso truth
        '''
        self._catalog_creator = catalog_creator
        self._output_dir = catalog_creator._output_dir
        self._logger = catalog_creator._logger
        self._sso_truth = sso_truth
        if sso_truth is None:
            self._sso_truth = SsoMainCatalogCreator._sso_truth

        self._sso_db_tbl = sso_db_tbl
        if sso_db_tbl is None:
            # correct for default input file
            self._sso_db_tbl = SsoMainCatalogCreator._sso_db_tbl

        self._row_group_size = _DEFAULT_ROW_GROUP_SIZE

        tbl = self._sso_db_tbl
        mjd_c = 'fieldMJD_TAI'
        self._mjd_q = f'select min({mjd_c}), max({mjd_c}), count({mjd_c}) from {tbl}'
        self._dfhp_query = f'''select ObjID as id, {mjd_c} as mjd,
                 "RA_deg" as ra, "Dec_deg" as dec,
                 "RARateCosDec_deg_day" as ra_rate,
                 "DecRate_deg_day" as dec_rate,
                 trailedSourceMag as trailed_source_mag from {tbl}
                 where healpix = (?)
                 order by mjd, ObjID'''

    @property
    def sso_truth(self):
        return self._sso_truth

    def _create_main_schema(self, metadata_input=None,
                            metadata_key='provenance'):

        fields = [
            pa.field('id', pa.string()),
            pa.field('mjd', pa.float64()),
            pa.field('ra', pa.float64()),
            pa.field('dec', pa.float64()),
            pa.field('trailed_source_mag', pa.float64()),
            pa.field('ra_rate', pa.float64()),
            pa.field('dec_rate', pa.float64())]
        if metadata_input:
            metadata_bytes = json.dumps(metadata_input).encode('utf8')
            final_metadata = {metadata_key: metadata_bytes}
        else:
            final_metadata = None

        return pa.schema(fields, metadata=final_metadata)

    def _get_hps(self, filepath):

        with sqlite3.connect(f'file:{filepath}?mode=ro', uri=True) as conn:
            df = pd.read_sql_query(f'select distinct healpix from {self._sso_db_tbl}',
                                   conn)
        return set(df['healpix'])

    def _write_hp(self, hp, hps_by_file, arrow_schema):
        df_list = []
        for f in hps_by_file:
            if hp in hps_by_file[f]:
                conn = sqlite3.connect(f'file:{f}?mode=ro', uri=True)
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
        inputs = {'sso_truth': self._sso_truth}
        run_options = self._catalog_creator._run_options
        file_metadata = assemble_file_metadata(self._catalog_creator._pkg_root,
                                               inputs=inputs,
                                               run_options=run_options)
        arrow_schema = self._create_main_schema(metadata_input=file_metadata)
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

        # Add config information for sso
        prov = assemble_provenance(
            self._catalog_creator._pkg_root,
            inputs={'sso_truth': self._sso_truth},
            run_options=self._catalog_creator._run_options)
        sso_fragment = SsoConfigFragment(prov)
        self._catalog_creator._config_writer.write_configs(sso_fragment)


class SsoFluxCatalogCreator:
    def __init__(self, catalog_creator):
        '''
        Parameters
        ----------
        catalog_creator   instance of FluxCatalogCreator
        sso_truth         path to input data directory
        '''
        self._catalog_creator = catalog_creator
        self._output_dir = catalog_creator._output_dir
        self._logger = catalog_creator._logger

    def _create_flux_schema(self, metadata_input=None,
                            metadata_key='provenance'):
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
        if metadata_input:
            metadata_bytes = json.dumps(metadata_input).encode('utf8')
            final_metadata = {metadata_key: metadata_bytes}
        else:
            final_metadata = None

        return pa.schema(fields, metadata=final_metadata)

    def _create_sso_flux_pixel(self, pixel, arrow_schema):
        '''
        Create parquet file for a single healpix pixel, containing
        id, mjd and fluxes
        Parameters
        pixel         int
        arrow_schema  schema for parquet file to be output
        Return
        ------
        None
        '''

        # global _sso_collection
        output_filename = f'sso_flux_{pixel}.parquet'
        output_path = os.path.join(self._catalog_creator._output_dir,
                                   output_filename)

        object_list = self._cat.get_object_type_by_hp(pixel, 'sso')
        n_parallel = self._catalog_creator._flux_parallel

        colls = object_list.get_collections()
        writer = pq.ParquetWriter(output_path, arrow_schema)
        fields_needed = arrow_schema.names
        instrument_needed = ['lsst']
        rg_written = 0
        writer = None

        # There is only one row group per main healpixel file currently
        # so the loop could be dispensed with
        for c in colls:
            _sso_collection = c
            l_bnd = 0
            u_bnd = len(c)
            if (u_bnd - l_bnd) < 5 * n_parallel:
                n_parallel = 1
            if n_parallel == 1:
                n_per = u_bnd - l_bnd
            else:
                n_per = int((u_bnd - l_bnd + n_parallel)/n_parallel)

            lb = l_bnd
            u = min(l_bnd + n_per, u_bnd)
            readers = []
            if n_parallel == 1:
                out_dict = _do_sso_flux_chunk(None, c, instrument_needed,
                                              l_bnd, u_bnd)
            else:
                out_dict = {}
                for field in fields_needed:
                    out_dict[field] = []

                tm = max(int((n_per*60)/500), 10)
                self._logger.info(f'Using timeout value {tm} for {n_per} sources')
                p_list = []
                for i in range(n_parallel):
                    conn_rd, conn_wrt = Pipe(duplex=False)
                    readers.append(conn_rd)

                    # For debugging call directly
                    proc = Process(target=_do_sso_flux_chunk,
                                   name=f'proc_{i}',
                                   args=(conn_wrt, _sso_collection,
                                         instrument_needed, lb, u))
                    proc.start()
                    p_list.append(proc)
                    lb = u
                    u = min(lb + n_per, u_bnd)
                self._logger.debug('Processes started')
                for i in range(n_parallel):
                    ready = readers[i].poll(tm)
                    if not ready:
                        self._logger.error(f'Process {i} timed out after {tm} sec')
                        sys.exit(1)
                    dat = readers[i].recv()
                    for field in fields_needed:
                        out_dict[field] += dat[field]
                for p in p_list:
                    p.join()

            out_df = pd.DataFrame.from_dict(out_dict)
            out_table = pa.Table.from_pandas(out_df,
                                             schema=arrow_schema)

            if not writer:
                writer = pq.ParquetWriter(output_path, arrow_schema)
            writer.write_table(out_table)

            rg_written += 1

        writer.close()
        self._logger.debug(f'# row groups written to flux file: {rg_written}')

    def create_sso_flux_catalog(self):
        """
        Create sso flux catalog.  Includes id, mjd (since given sso normally
        has several observations) and fluxes.
        """

        self._cat = self._catalog_creator._cat
        sso_config = self._cat.raw_config['object_types']['sso']
        self._flux_template = sso_config['flux_file_template']
        self._main_template = sso_config['file_template']

        thru_v = {'lsst_throughputs_version': self._cat._lsst_thru_v}
        file_metadata = assemble_file_metadata(
            self._catalog_creator._pkg_root,
            inputs={'sso_sed_path': self._cat._sso_sed_factory.sed_path},
            run_options=self._catalog_creator._run_options,
            flux_file=True,
            throughputs_versions=thru_v)

        arrow_schema = self._create_flux_schema(metadata_input=file_metadata)

        self._logger.info('Creating sso flux files')

        for p in self._catalog_creator._parts:
            self._create_sso_flux_pixel(p, arrow_schema)
