import os
import sys
from time import sleep
from multiprocessing import Process, Pipe
from datetime import datetime
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import healpy
import numpy as np
import json
import galsim
from skycatalogs.objects.base_object import LSST_BANDS, load_lsst_bandpasses
from skycatalogs.objects.trilegal_object import TrilegalConfigFragment
from skycatalogs.utils.config_utils import assemble_provenance
from skycatalogs.utils.config_utils import assemble_file_metadata
from skycatalogs.utils.creator_utils import get_trilegal_hp_nrows

"""
Code for creating sky catalogs for trilegal stars
"""

__all__ = ['TrilegalMainCatalogCreator', 'TrilegalFluxCatalogCreator']

_DEFAULT_ROW_GROUP_SIZE = 10_000_000  # since fewer columns than galaxies have
_DEFAULT_TRUTH_CATALOG = 'lsst_sim.simdr2'
_DEFAULT_START_EPOCH = 2000


class TrilegalMainCatalogCreator:
    def __init__(self, catalog_creator, truth_catalog,
                 start_epoch=_DEFAULT_START_EPOCH):
        '''
        Parameters
        ----------
        catalog_creator  instance of MainCatalogCreator
        input_catalog    name of Trilegal catalog to be queried
        '''
        self._catalog_creator = catalog_creator
        if not truth_catalog:
            self._truth_catalog = _DEFAULT_TRUTH_CATALOG
        else:
            self._truth_catalog = truth_catalog

        self._start_epoch = start_epoch
        self._output_dir = catalog_creator._output_dir
        self._logger = catalog_creator._logger
        # self._stride = _DEFAULT_ROW_GROUP_SIZE
        self._stride = self._catalog_creator._stride

    @property
    def trilegal_truth(self):
        return self._truth_catalog

    def _create_main_schema(self, metadata_input=None,
                            metadata_key='provenance'):

        fields = [
            pa.field('id', pa.string()),
            pa.field('ra', pa.float64()),
            pa.field('dec', pa.float64()),
            pa.field('av', pa.float64()),
            pa.field('pmdec', pa.float32()),  # proper motion in dec
            pa.field('pmracosd', pa.float32()),  # proper motion in cos(ra)
            pa.field('vrad', pa.float32()),  # radial velocity
            pa.field('mu0', pa.float32()),  # true distance modulus
            pa.field('evol_label', pa.int32()),  # evolutionary phase

            #  The following are used in SED generation
            pa.field('logT', pa.float32()),  # log effective temp. (K)
            pa.field('logg', pa.float32()),  # log surface gravity (cgs)
            pa.field('logL', pa.float32()),  # log luminosity  (L_sun)
            pa.field('Z', pa.float32()),  # heavy element abund.

            # Add in magnitudes for verification and normalization (imag)
            pa.field('umag', pa.float32()),
            pa.field('gmag', pa.float32()),
            pa.field('rmag', pa.float32()),
            pa.field('imag', pa.float32()),
            pa.field('zmag', pa.float32()),
            pa.field('ymag', pa.float32()),
            ]

        if metadata_input:
            metadata_bytes = json.dumps(metadata_input).encode('utf8')
            final_metadata = {metadata_key: metadata_bytes}
        else:
            final_metadata = None

        return pa.schema(fields, metadata=final_metadata)

    def _write_hp(self, hp, arrow_schema):
        '''
        Write out parquet file for specified healpixel

        Parameters
        ----------
        hp             int             the healpixel
        arrow_schema   pyarrow.schema  schema to use

        Returns
        -------
        Number of files written (0 or 1)

        '''

        ## MAX_QUERY_ROWS = 10_000_000
        MAX_QUERY_ROWS = 1_000_000
        #   MAX_QUERY_ROWS = 100_000            # temp for testing

        # Note dl is not part of desc-python or LSST Pipelines; it must be
        # separately installed
        from dl import queryClient as qc
        _NSIDE = 32
        _TRILEGAL_RING_NSIDE = 256
        _TRILEGAL_NEST_NSIDE = 4096

        nrows = get_trilegal_hp_nrows(hp, nside=_NSIDE)
        n_query = 1
        if nrows > MAX_QUERY_ROWS:
            # break up into several queries
            if nrows > 50 * MAX_QUERY_ROWS:
                n_query = 64
            elif nrows > 16 * MAX_QUERY_ROWS:
                n_query = 16
            else:
                n_query = 4

        if hp == 9246:   # 64 is not fine enough.  Query times out
            n_query = 256
        elif hp == 9119:
            n_query = 64

        self._logger.debug(f'MAX_QUERY_ROWS: {MAX_QUERY_ROWS}')
        self._logger.debug(f'Rows in healpix {hp}: {nrows}')
        self._logger.debug(f'Queries to db: {n_query}')
        def _next_level(pixel):
            return [4 * pixel, 4 * pixel + 1, 4 * pixel + 2, 4 * pixel + 3]

        # Find all subpixels contained in our pixel with nside 256,
        # which is what the Trilegal catalog uses
        current_nside = _NSIDE
        pixels = [healpy.ring2nest(_NSIDE, hp)]

        use_ring = False
        if n_query <= 64:
            use_column = 'ring256'
            while current_nside < _TRILEGAL_RING_NSIDE:
                pixels = [pix for p in pixels for pix in _next_level(p)]
                current_nside = current_nside * 2

            query_pixels = [healpy.nest2ring(_TRILEGAL_RING_NSIDE, p) for p in pixels]
        else:
            use_column = 'nest4096'
            while current_nside < _TRILEGAL_NEST_NSIDE:
                pixels = [pix for p in pixels for pix in _next_level(p)]
                current_nside = current_nside * 2

            query_pixels = pixels
        # Form the queries and issue them
        to_select = ['ra', 'dec', 'av', 'pmracosd', 'pmdec', 'vrad', 'mu0',
                     'label as evol_label', 'logte as logT', 'logg',
                     'logl as logL', 'z as Z',
                     'umag', 'gmag', 'rmag', 'imag', 'zmag', 'ymag']
        # all_results = []
        writer = None
        rg_written = 0
        per_query = int(len(query_pixels) / n_query)
        so_far = 0
        for iq in range(n_query):
            in_pixels = ','.join(str(p) for p in query_pixels[iq*per_query:(iq+1)*per_query])

            q = 'select ' + ','.join(to_select)

            q += f' from {self._truth_catalog} where {use_column} in ({in_pixels})'
            self._logger.debug(f'column {use_column} pixels {in_pixels}')
            # 600 seconds is max timeout allowed for synchronous query
            # 300 is generous.  Returning 6 million rows took 80 sec.
            # Hardly any of these queries return that many rows.
            try:
                results = qc.query(adql=q, fmt='pandas', timeout=300)
            except dl.queryClient.queryClientError as e:
                self._logger.debug(str(e))
                self._logger.debug('Sleep and retry')
                sleep(10)
                results = qc.query(adql=q, fmt='pandas', timeout=600)

            n_row = len(results['ra'])
            if not n_row:
                continue
            self._logger.info(f'rows returned: {n_row}')

            # generate id
            id_prefix = f'{self._truth_catalog}_hp{hp}_'
            results['id'] = [f'{id_prefix}{n}' for n in range(so_far, so_far + n_row)]
            so_far += n_row
            l_bnd = 0
            u_bnd = min(n_row, self._stride)
            outpath = ''

            if not writer:
                outpath = os.path.join(self._output_dir,
                                       f'trilegal_{hp}.parquet')
                writer = pq.ParquetWriter(outpath, arrow_schema)

            while u_bnd > l_bnd:
                out_dict = {k: results[k][l_bnd: u_bnd] for k in results}
                out_df = pd.DataFrame.from_dict(out_dict)
                out_table = pa.Table.from_pandas(out_df, schema=arrow_schema)

                # Parquet default max rows in a row group is 1M. Since
                # trilegal has a small number of columns, we can afford
                # to have more rows.
                writer.write_table(out_table, row_group_size=self._stride)
                l_bnd = u_bnd
                u_bnd = min(l_bnd + self._stride, n_row)
                rg_written += 1
                del out_dict
                del out_df
                del out_table

        writer.close()
        self._logger.debug(f'# row groups written to {outpath}: {rg_written}')
        return rg_written

    def create_catalog(self, hps):
        inputs = {'trilegal_truth': self._truth_catalog}
        file_metadata = assemble_file_metadata(
            self._catalog_creator._pkg_root,
            inputs=inputs,
            run_options=self._catalog_creator._run_options)
        schema = self._create_main_schema(metadata_input=file_metadata,
                                          metadata_key='provenance')
        written = 0
        for hp in hps:
            written += self._write_hp(hp, schema)

        if written == 0:
            return
        # Add config information for trilegal
        prov = assemble_provenance(
            self._catalog_creator._pkg_root,
            inputs=inputs,
            run_options=self._catalog_creator._run_options)
        trilegal_fragment = TrilegalConfigFragment(prov)
        self._catalog_creator._config_writer.write_configs(trilegal_fragment)


def _do_trilegal_flux_chunk(send_conn, collection, instrument_needed,
                            l_bnd, u_bnd,  main_path, row_group, debug=False):
    '''
    send_conn         output connection.  If none return output
    collection       object collection we're processing
    instrument_needed indicates which fluxes need to be computed
                     Ignored for now.  Just do LSST fluxes
    l_bnd, u_bnd     demarcates slice to process
    main_path        Path main skyCatalogs file for
                     current healpixel
    row_group        row group this chunk belongs to.
                     Row group # = (SED file) batch #

    returns
                    dict with keys id, lsst_flux_u, ... lsst_flux_y
    '''
    global tri_lsst_bandpasses
    out_dict = {}
    if debug:
        now = datetime.now().isoformat()[:19]
        print(f'{now}  Entering _do_trilegal_flux_chunk, l_bnd={l_bnd}, u_bnd={u_bnd}, row_group={row_group}', flush=True)

    if l_bnd >= u_bnd:
        if debug:
            print("_do_trilegal_flux_chunk: nothing to do", flush=True)
        if send_conn:
            send_conn.send(out_dict)
            return
        else:
            return out_dict

    skycat = collection._sky_catalog
    factory = skycat._trilegal_sed_factory
    extinguisher = skycat._extinguisher

    pq_main = pq.ParquetFile(main_path)

    wl, spectra = factory.get_spectra_batch(pq_main, row_group, l_bnd, u_bnd)
    if debug:
        now = datetime.now().isoformat()[:19]
        print(f'{now} Spectra computed', flush=True)
    av = collection.get_native_attribute('av')
    id = collection.get_native_attribute('id')
    imag = collection.get_native_attribute('imag')
    out_dict['id'] = id[l_bnd: u_bnd]
    fluxes = []
    for ix in range(l_bnd, u_bnd):
        obj_fluxes = []
        lut = galsim.LookupTable(wl, spectra[ix - l_bnd], interpolant='linear')
        sed = galsim.SED(lut, wave_type='nm', flux_type='flambda')
        sed = extinguisher.extinguish(sed, av[ix])

        if sed is None:
            obj_fluxes = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        else:
            # Normalize; calculate fluxes
            sed = sed.withMagnitude(imag[ix], tri_lsst_bandpasses['i'])
            for band in LSST_BANDS:
                obj_fluxes.append(collection[ix].get_LSST_flux(
                    band, sed=sed, cache=False))
        fluxes.append(obj_fluxes)

    colnames = [f'lsst_flux_{band}' for band in LSST_BANDS]
    fluxes_transpose = zip(*fluxes)
    flux_dict = dict(zip(colnames, fluxes_transpose))
    out_dict.update(flux_dict)
    del fluxes_transpose

    if debug:
        now = datetime.now().isoformat()[:19]
        print(f'{now}  Leaving _do_trilegal_flux_chunk, l_bnd={l_bnd}, row_group={row_group}', flush=True)

    if send_conn:
        send_conn.send(out_dict)
    else:
        return out_dict


class TrilegalFluxCatalogCreator:
    def __init__(self, catalog_creator):
        '''
        Parameters
        ----------
        catalog_creator   instance of FluxCatalogCreator
        '''
        self._catalog_creator = catalog_creator
        self._output_dir = catalog_creator._output_dir
        self._logger = catalog_creator._logger
        global tri_lsst_bandpasses
        tri_lsst_bandpasses = load_lsst_bandpasses()

    def _create_flux_schema(self, metadata_input=None,
                            metadata_key='provenance'):
        # id and 6 flux fields (for now.  Maybe later also Roman)
        fields = [
            pa.field('id', pa.string()),
            # pa.field('mjd', pa.float64()),
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

    def _create_trilegal_flux_pixel(self, pixel, arrow_schema):
        output_filename = f'trilegal_flux_{pixel}.parquet'
        output_path = os.path.join(self._catalog_creator._output_dir,
                                   output_filename)

        main_filename = f'trilegal_{pixel}.parquet'
        self._logger.info(f'Main file (input) is {main_filename}')
        main_path = os.path.join(self._catalog_creator._output_dir,
                                 main_filename)

        n_parallel = self._catalog_creator._flux_parallel

        writer = pq.ParquetWriter(output_path, arrow_schema)

        instrument_needed = ['lsst']
        rg_written = 0
        writer = None
        fields_needed = arrow_schema.names

        # Get all the objects in the pixel
        # For test pixel there is only one row group so only one collection
        # In general may have to iterate over row groups
        obj_list = self._cat.get_object_type_by_hp(pixel, 'trilegal')
        if len(obj_list) == 0:
            self._logger.warning(f'Cannot create flux file for pixel {pixel} because main file does not exist or is empty')
            return

        for rg, c in enumerate(obj_list.get_collections()):
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

            if n_parallel == 1:
                # For debugging call directly
                out_dict = _do_trilegal_flux_chunk(None, c, instrument_needed,
                                                   l_bnd, u_bnd, main_path,
                                                   rg, debug=True)
            else:
                out_dict = {}
                readers = []
                for field in fields_needed:
                    out_dict[field] = []

                tm = max(int((n_per*60)/500), 10)
                self._logger.info(f'Using timeout value {tm} for {n_per} sources')
                p_list = []
                for i in range(n_parallel):
                    conn_rd, conn_wrt = Pipe(duplex=False)
                    readers.append(conn_rd)

                    proc = Process(target=_do_trilegal_flux_chunk,
                                   name=f'proc_{i}',
                                   args=(conn_wrt, c,
                                         instrument_needed, lb, u,
                                         main_path, rg, False))
                    proc.start()
                    p_list.append(proc)
                    lb = u
                    u = min(lb + n_per, u_bnd)
                self._logger.debug('Processes started')
                for i in range(n_parallel):
                    ready = readers[i].poll(tm)
                    if not ready:
                        self._logger.error(
                            f'Process {i} timed out after {tm} sec')
                        sys.exit(1)
                    dat = readers[i].recv()
                    if len(dat.keys()) > 0:
                        for field in fields_needed:
                            if len(out_dict[field]) == 0:
                                out_dict[field] = dat[field]
                            else:
                                out_dict[field] = np.concatenate([out_dict[field],
                                                                  dat[field]])
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

    def create_trilegal_flux_catalog(self):
        """
        Create trilegal flux catalog.  Includes id, and fluxes.
        """

        self._cat = self._catalog_creator._cat
        trilegal_config = self._cat.raw_config['object_types']['trilegal']
        self._flux_template = trilegal_config['flux_file_template']
        self._main_template = trilegal_config['file_template']

        thru_v = {'lsst_throughputs_version': self._cat._lsst_thru_v}
        file_metadata = assemble_file_metadata(
            self._catalog_creator._pkg_root,
            run_options=self._catalog_creator._run_options,
            flux_file=True,
            throughputs_versions=thru_v)

        arrow_schema = self._create_flux_schema(metadata_input=file_metadata)

        self._logger.info('Creating trilegal flux files')

        for p in self._catalog_creator._parts:
            self._create_trilegal_flux_pixel(p, arrow_schema)
