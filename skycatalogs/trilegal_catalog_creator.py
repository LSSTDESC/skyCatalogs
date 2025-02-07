import os
# import sys                        May be needed later
# from multiprocessing import Process, Pipe    May be needed later
from datetime import datetime
import logging
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import healpy
#  ## import h5py
import numpy as np
import json
import galsim
from skycatalogs.objects.base_object import LSST_BANDS, load_lsst_bandpasses
from skycatalogs.objects.trilegal_object import TrilegalConfigFragment
from skycatalogs.utils.config_utils import assemble_provenance
from skycatalogs.utils.config_utils import assemble_file_metadata

"""
Code for creating sky catalogs for trilegal stars
"""

__all__ = ['TrilegalMainCatalogCreator',  # 'TrilegalSEDGenerator',
           'TrilegalFluxCatalogCreator']

_DEFAULT_ROW_GROUP_SIZE = 1000000
_DEFAULT_TRUTH_CATALOG = 'lsst_sim.simdr2'
_DEFAULT_START_EPOCH = 2000


class TrilegalMainCatalogCreator:
    # Note dl is not part of desc-python or LSST Pipelines; it must be
    # separately installed
    from dl import queryClient as qc

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
        from dl import queryClient as qc
        _NSIDE = 32
        _TRILEGAL_RING_NSIDE = 256

        def _next_level(pixel):
            return [4 * pixel, 4 * pixel + 1, 4 * pixel + 2, 4 * pixel + 3]

        # Find all subpixels contained in our pixel with nside 256,
        # which is what the Trilegal catalog uses
        current_nside = _NSIDE
        pixels = [healpy.ring2nest(_NSIDE, hp)]

        while current_nside < _TRILEGAL_RING_NSIDE:
            pixels = [pix for p in pixels for pix in _next_level(p)]
            current_nside = current_nside * 2

        ring_pixels = [healpy.nest2ring(_TRILEGAL_RING_NSIDE, p) for p in pixels]
        in_pixels = ','.join(str(p) for p in ring_pixels)

        # Form query
        to_select = ['ra', 'dec', 'av', 'pmracosd', 'pmdec', 'vrad', 'mu0',
                     'label as evol_label', 'logte as logT', 'logg',
                     'logl as logL', 'z as Z',
                     'umag', 'gmag', 'rmag', 'imag', 'zmag', 'ymag']
        q = 'select ' + ','.join(to_select)
        q += f' from {self._truth_catalog} where ring256 in ({in_pixels})'

        #   q += ' and label = 1'  # main sequence only

        results = qc.query(adql=q, fmt='pandas')
        n_row = len(results['ra'])
        if not n_row:
            return 0

        self._logger.info(f'rows returned: {n_row}')

        # generate id
        id_prefix = f'{self._truth_catalog}_hp{hp}_'
        results['id'] = [f'{id_prefix}{n}' for n in range(n_row)]

        l_bnd = 0
        u_bnd = min(n_row, self._stride)
        rg_written = 0
        outpath = ''

        writer = None
        if not writer:
            outpath = os.path.join(self._output_dir,
                                   f'trilegal_{hp}.parquet')
            writer = pq.ParquetWriter(outpath, arrow_schema)

        while u_bnd > l_bnd:
            out_dict = {k: results[k][l_bnd: u_bnd] for k in results}
            out_df = pd.DataFrame.from_dict(out_dict)
            out_table = pa.Table.from_pandas(out_df, schema=arrow_schema)
            writer.write_table(out_table)
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
        schema = self._create_main_schema()
        written = 0
        for hp in hps:
            written += self._write_hp(hp, schema)

        if written == 0:
            return
        # Add config information for trilegal
        prov = assemble_provenance(
            self._catalog_creator._pkg_root,
            inputs={'trilegal_truth': self._truth_catalog},
            run_options=self._catalog_creator._run_options)
        trilegal_fragment = TrilegalConfigFragment(prov)
        self._catalog_creator._config_writer.write_configs(trilegal_fragment)


# class TrilegalSEDGenerator:
#     '''
#     This class is responsible for calculating SEDs for Trilegal sources
#     and storing them in an h5df file (one file per healpixel). The
#     strucure of the file is as follows:
#     metadata (group)   one or more key/value pairs, stored as the attributes
#                        of the top-level group "metadata". The group has no
#                        other content
#     wl_axis (dataset)  All SEDs use the same binning. It's stored here.
#                        Array of 4-byte floats
#     batches (group)    Contains subgroups with names like "batch_X" where
#                        X is the batch number
#     batch_X (group)    Has metadata (batch number, count of sources in the
#                        batch, maybe lower and upper limits of count within the
#                        full healpixel).  It also contains two datasets, "id"
#                        and "SED"
#     id (dataset)       1-d array of (string) ids in the batch
#     spectra (dataset)  2-d array indexed by [source_number, lambda]. 4-byte
#                        floats
#     '''
#     def __init__(self, input_dir, output_dir=None, lib='BTSettl',
#                  log_level='INFO'):
#         '''
#         Parameters
#         ----------
#         input_dir      Where to find parquet main files
#         output_dir     Where to write output.  Defaults to input_dir
#         lib            Name of pystelllibs library used to generate SEDs
#         '''
#         from pystellibs import BTSettl

#         self._pystellib = BTSettl(medres=False)
#         self._wl = self._pystellib._wavelength / 10  # Convert to nm

#         self._input_dir = input_dir
#         self._output_dir = output_dir
#         if not output_dir:
#             self._output_dir = input_dir

#         self._logger = logging.getLogger('trilegal_SED')
#         if not self._logger.hasHandlers():
#             self._logger.setLevel(log_level)
#             ch = logging.StreamHandler()
#             ch.setLevel(log_level)
#             formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
#             ch.setFormatter(formatter)
#             self._logger.addHandler(ch)

#     def _write_SED_batch(self, outfile_path, batch, id, wl_axis, spectra):
#         '''
#         Write a collection of SEDs, typically corresponding to a single row
#         group in the parquet main file, to output.
#         Parameters
#         ----------
#         outfile_path  File to which SEDs will be appended
#         batch         E.g. row group number
#         id            Array of ids for the sources whose SEDs will be written
#         wl_axis       Array of wavelengths (nm) used for spectra
#         spectra       Flux values, units of flambda (erg/cm**2/s/nm)

#         '''
#         with h5py.File(outfile_path, 'a') as f:
#             if 'wl_axis' not in f.keys():
#                 axis_ds = f.create_dataset('wl_axis', shape=(len(wl_axis),),
#                                            dtype='f4', data=np.array(wl_axis))
#                 axis_ds.attrs.create('wl_units', 'nm')
#             if 'batches' not in f.keys():
#                 _ = f.create_group('batches')

#             batch_g = f.create_group(f'batches/batch_{batch}')
#             batch_g.attrs.create('source_count', len(id))
#             max_id_len = max([len(entry) for entry in id])
#             id_dat = np.array(id, dtype=f'S{max_id_len}')
#             id_chunk_size = min(50000, len(id_dat))
#             _ = batch_g.create_dataset('id', data=id_dat,
#                                        chunks=(id_chunk_size),
#                                        compression='gzip')

#             # Larger chunks speed up i/o but require more memory.
#             # Optimal value is TBD
#             data_chunk_nrow = min(1000, len(id_dat))
#             _ = batch_g.create_dataset('spectra',
#                                        data=np.array(spectra),
#                                        chunks=(data_chunk_nrow, len(wl_axis)),
#                                        compression='gzip',
#                                        dtype='f4')

#     def _generate_hp(self, hp):
#         # open parquet main file.
#         # For now use hardcoded templates.  Should read from config
#         # Modify to be suitable for creating rather than parsing a name
#         main_templ = 'trilegal_(?P<healpix>\d+).parquet'
#         sed_templ = 'trilegal_sed_(?P<healpix>\d+).hdf5'
#         in_fname = main_templ.replace('(?P<healpix>\\d+)', str(hp))
#         out_fname = sed_templ.replace('(?P<healpix>\\d+)', str(hp))

#         # Open input and output files
#         # For now read main file ourselves, not via skyCatalogs API
#         in_path = os.path.join(self._input_dir, in_fname)
#         pq_file = pq.ParquetFile(in_path)
#         n_gp = pq_file.metadata.num_row_groups
#         hp5_path = os.path.join(self._output_dir, out_fname)
#         with h5py.File(hp5_path, 'w') as f:
#             f.create_group('metadata')
#             f['metadata'].attrs.create('input_path', in_path)
#             f['metadata'].attrs.create('n_batch', n_gp)  # was in brackets

#         for batch in range(n_gp):
#             self._logger.info(f'Starting batch {batch}')
#             columns = ['id', 'logT', 'logg', 'logL', 'Z', 'mu0']
#             df = pq_file.read_row_group(batch, columns=columns).to_pandas()
#             wl_axis, spectra = self._pystellib.generate_individual_spectra(df)
#             self._logger.info('Computed spectra')
#             # Convert wl_axis, spectra from A to nm.
#             wl_axis = wl_axis / 10
#             spectra = spectra * 10
#             spectra_32 = spectra.astype(np.float32)
#             del spectra
#             # print(f'len(wavelength axis): {len(wl_axis)}')
#             self._debug.info(f'shape of spectra array: {spectra_32.shape}')
#             self._write_SED_batch(hp5_path, batch, df['id'],
#                                   wl_axis, spectra_32)
#             self._logger.info(f'Batch {batch} written')

#     def generate(self, hps):

#         for hp in hps:
#             self._generate_hp(hp)


def _do_trilegal_flux_chunk(send_conn, collection, instrument_needed,
                            l_bnd, u_bnd,  main_path, row_group):
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
    now = datetime.now().isoformat()[:19]
    print(f'{now}  Entering _do_trilegal_flux_chunk, l_bnd={l_bnd}, row_group={row_group}', flush=True)

    # hp = collection._partition_id
    skycat = collection._sky_catalog
    factory = skycat._trilegal_sed_factory
    #  ## sedfile = factory.get_hp_sedfile(hp, silent=False)  # should be there
    extinguisher = skycat._extinguisher

    pq_main = pq.ParquetFile(main_path)

    # SEDs we need should be in group named "batch_XX" where XX is row
    # group number

    # ## spectra = sedfile.get_sed_batch(f"batch_{row_group}")
    # ## wl = sedfile.wl_axis
    wl, spectra = factory.get_spectra_batch(pq_main, row_group)
    now = datetime.now().isoformat()[:19]
    print(f'{now} Spectra computed', flush=True)
    av = collection.get_native_attribute('av')
    id = collection.get_native_attribute('id')
    imag = collection.get_native_attribute('imag')
    # label = collection.get_native_attribute('evol_label')
    out_dict['id'] = id[l_bnd: u_bnd]
    fluxes = []
    for ix in range(l_bnd, u_bnd):
        obj_fluxes = []
        lut = galsim.LookupTable(wl, spectra[ix], interpolant='linear')
        sed = galsim.SED(lut, wave_type='nm', flux_type='flambda')
        sed = extinguisher.extinguish(sed, av[ix])

        if sed is None:
            obj_fluxes = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        else:
            # Normalize; calculate fluxes
            # sed = sed.thin()
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
        # self._sed_factory = catalog_creator._trilegal_sed_factory

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
        main_path = os.path.join(self._catalog_creator._output_dir,
                                 main_filename)

        n_parallel = self._catalog_creator._flux_parallel

        writer = pq.ParquetWriter(output_path, arrow_schema)
        #  fields_needed = arrow_schema.names
        instrument_needed = ['lsst']
        rg_written = 0
        writer = None

        # Get all the objects in the pixel
        # For test pixel there is only one row group so only one collection
        # In general may have to iterate over row groups
        obj_list = self._cat.get_object_type_by_hp(pixel, 'trilegal')

        for ix, c in enumerate(obj_list.get_collections()):
            # av = c.get_native_attribute('av')
            l_bnd = 0
            u_bnd = len(c)
            if (u_bnd - l_bnd) < 5 * n_parallel:
                n_parallel = 1
            if n_parallel == 1:
                n_per = u_bnd - l_bnd
            else:
                n_per = int((u_bnd - l_bnd + n_parallel)/n_parallel)

            lb = l_bnd
            #  u = min(l_bnd + n_per, u_bnd)
            #  readers = []
            if n_parallel == 1:
                out_dict = _do_trilegal_flux_chunk(None, c, instrument_needed,
                                                   l_bnd, u_bnd, main_path, ix)
            else:
                raise Exception('Parallel flux computation for trilegal NYI')
                # out_dict = {}
                # for field in fields_needed:
                #     out_dict[field] = []

                # tm = max(int((n_per*60)/500), 10)
                # self._logger.info(f'Using timeout value {tm} for {n_per} sources')
                # p_list = []
                # for i in range(n_parallel):
                #     conn_rd, conn_wrt = Pipe(duplex=False)
                #     readers.append(conn_rd)

                #     # For debugging call directly
                #     proc = Process(target=_do_trilegal_flux_chunk,
                #                    name=f'proc_{i}',
                #                    args=(conn_wrt, _trilegal_collection,
                #                          instrument_needed, lb, u,
                #                          main_path, ix))
                #     proc.start()
                #     p_list.append(proc)
                #     lb = u
                #     u = min(lb + n_per, u_bnd)
                # self._logger.debug('Processes started')
                # for i in range(n_parallel):
                #     ready = readers[i].poll(tm)
                #     if not ready:
                #         self._logger.error(f'Process {i} timed out after {tm} sec')
                #         sys.exit(1)
                #     dat = readers[i].recv()
                #     for field in fields_needed:
                #         out_dict[field] += dat[field]
                # for p in p_list:
                #     p.join()

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
            # inputs={'trilegal_sed_path': self._cat.??},
            run_options=self._catalog_creator._run_options,
            flux_file=True,
            throughputs_versions=thru_v)

        arrow_schema = self._create_flux_schema(metadata_input=file_metadata)

        self._logger.info('Creating trilegal flux files')

        for p in self._catalog_creator._parts:
            self._create_trilegal_flux_pixel(p, arrow_schema)
