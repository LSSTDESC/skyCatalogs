import os
# import sys                        May be needed later
# from multiprocessing import Process, Pipe    May be needed later
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import healpy
import h5py
import numpy as np
import json
from skycatalogs.objects.base_object import LSST_BANDS
from skycatalogs.objects.trilegal_object import TrilegalConfigFragment
from skycatalogs.utils.config_utils import assemble_provenance
from skycatalogs.utils.config_utils import assemble_file_metadata

"""
Code for creating sky catalogs for trilegal stars
"""

__all__ = ['TrilegalMainCatalogCreator', 'TrilegalSEDGenerator',
           'TrilegalFluxCatalogCreator']

_DEFAULT_ROW_GROUP_SIZE = 1000000
_DEFAULT_TRUTH_CATALOG = 'lsst_sim.simdr2'
_DEFAULT_START_EPOCH = 2000


def _do_trilegal_flux_chunk(send_conn, trilegal_collection, instrument_needed,
                            l_bnd, u_bnd):
    '''
    end_conn         output connection
    trilegal_collection  information from main file
    instrument_needed List of which calculations should be done
    l_bnd, u_bnd     demarcates slice to process

    returns
                    dict with keys id, lsst_flux_u, ... lsst_flux_y
    '''
    out_dict = {}

    # Not implemented yet
    return out_dict


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
        self._stride = _DEFAULT_ROW_GROUP_SIZE

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
            #  The following are used in SED generation
            pa.field('logT', pa.float32()),  # log effective temp. (K)
            pa.field('logg', pa.float32()),  # log surface gravity (cgs)
            pa.field('logL', pa.float32()),  # log luminosity  (L_sun)
            pa.field('Z', pa.float32()),  # heavy element abund.
            ]
        # The truth catalog supplies only mu0, not parallax
        # from https://en.wikipedia.org/wiki/Distance_modulus
        # luminous distance in parsecs =
        #  10**(1 + mu0/5)
        #
        # Formula for distance using parallax:
        #     d = E-sun/tan(0.5*theta) where E-sun is one AU
        # and 0.5*theta is parallax.   Or, approximately for small angles
        #    d (in parsecs) ~= 1/parallax

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
                     'logte as logT', 'logg', 'logl as logL', 'z as Z']
        q = 'select ' + ','.join(to_select)
        q += f' from {self._truth_catalog} where ring256 in ({in_pixels})'

        results = qc.query(adql=q, fmt='pandas')
        n_row = len(results['ra'])
        if not n_row:
            return 0

        print(f'rows returned: {n_row}')

        # generate id
        id_prefix = f'{self._truth_catalog}_hp{hp}_'
        results['id'] = [f'{id_prefix}{n}' for n in range(n_row)]

        l_bnd = 0
        u_bnd = min(n_row, self._stride)
        rg_written = 0
        outpath = ''

        while u_bnd > l_bnd:
            out_dict = {k: results[k][l_bnd: u_bnd] for k in results}
            out_df = pd.DataFrame.from_dict(out_dict)
            out_table = pa.Table.from_pandas(out_df, schema=arrow_schema)
            writer = None
            if not writer:
                outpath = os.path.join(self._output_dir,
                                       f'trilegal_{hp}.parquet')
                writer = pq.ParquetWriter(outpath, arrow_schema)
            writer.write_table(out_table)
            l_bnd = u_bnd
            u_bnd = min(l_bnd + self._stride, n_row)
            rg_written += 1

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


class TrilegalSEDGenerator:
    '''
    This class is responsible for calculating SEDs for Trilegal sources
    and storing them in an h5df file (one file per healpixel). The
    strucure of the file is as follows:
    metadata (group)   one or more key/value pairs, stored as the attributes
                       of the top-level group "metadata". The group has no
                       other content
    wl_axis (dataset)  All SEDs use the same binning. It's stored here.
                       Array of 4-byte floats
    batches (group)    Contains subgroups with names like "batch_X" where
                       X is the batch number
    batch_X (group)    Has metadata (batch number, count of sources in the
                       batch, maybe lower and upper limits of count within the
                       full healpixel).  It also contains two datasets, "id"
                       and "SED"
    id (dataset)       1-d array of (string) ids in the batch
    spectra (dataset)  2-d array indexed by [source_number, lambda]. 4-byte
                       floats
    '''
    def __init__(self, input_dir, output_dir=None, lib='BTSettl'):
        '''
        Parameters
        ----------
        input_dir      Where to find parquet main files
        output_dir     Where to write output.  Defaults to input_dir
        lib            Name of pystelllibs library used to generate SEDs
        '''
        from pystellibs import BTSettl

        self._pystellib = BTSettl(medres=False)
        self._wl = self._pystellib._wavelength / 10  # Convert to nm

        self._input_dir = input_dir
        self._output_dir = output_dir
        if not output_dir:
            self._output_dir = input_dir
        # self._logger = catalog_creator._logger

    def _write_SED_batch(self, outfile_path, batch, id, wl_axis, spectra):
        '''
        Write a collection of SEDs, typically corresponding to a single row
        group in the parquet main file, to output.
        Parameters
        ----------
        outfile_path  File to which SEDs will be appended
        batch         E.g. row group number
        id            Array of ids for the sources whose SEDs will be written
        wl_axis       Array of wavelengths (nm) used for spectra
        spectra       Flux values, units of flambda (erg/cm**2/s/nm)

        '''
        with h5py.File(outfile_path, 'a') as f:
            if 'wl_axis' not in f.keys():
                axis_ds = f.create_dataset('wl_axis', shape=(len(wl_axis),),
                                           dtype='f4', data=np.array(wl_axis))
                axis_ds.attrs.create('wl_units', ['nm'])
            if 'batches' not in f.keys():
                batches = f.create_group('batches')

            batch_g = f.create_group(f'batches/batch_{batch}')
            batch_g.attrs.create('source_count', [len(id)])
            max_id_len = max([len(entry) for entry in id])
            id_dat = np.array(id, dtype=f'S{max_id_len}')
            id_dataset = batch_g.create_dataset('id', data=id_dat,
                                                chunks=(50000),
                                                compression='gzip')
            spectra_dataset = batch_g.create_dataset('spectra',
                                                     data=np.array(spectra),
                                                     chunks=(200, len(wl_axis)),
                                                     compression='gzip',
                                                     dtype='f4')

    def _generate_hp(self, hp):
        PSEC_TO_CM = 3.085677581e16 * 100
        FOUR_PI = 4 * np.pi
        # open parquet main file.
        # For now use hardcoded templates.  Should read from config
        # Modify to be suitable for creating rather than parsing a name
        main_templ = 'trilegal_(?P<healpix>\d+).parquet'
        sed_templ = 'trilegal_sed_(?P<healpix>\d+).hdf5'
        in_fname = main_templ.replace('(?P<healpix>\\d+)', str(hp))
        out_fname = sed_templ.replace('(?P<healpix>\\d+)', str(hp))

        # Open input and output files
        # For now read main file ourselves, not via skyCatalogs API
        in_path = os.path.join(self._input_dir, in_fname)
        pq_file = pq.ParquetFile(in_path)
        hp5_path = os.path.join(self._output_dir, out_fname)
        with h5py.File(hp5_path, 'w') as f:
            f.create_group('metadata')
            f['metadata'].attrs.create('input_path', [in_path])

        n_gp = pq_file.metadata.num_row_groups

        for batch in range(n_gp):
            columns = ['id', 'logT', 'logg', 'logL', 'Z', 'mu0']
            # tbl = pq_file.read_row_group(i, columns=columns)
            py_dict = pq_file.read_row_group(batch, columns=columns).to_pydict()
            for c in ['logT', 'logg', 'logL', 'Z', 'mu0']:
                py_dict[c] = np.array(py_dict[c], dtype=np.float64)
            wl_axis, spectra = self._pystellib.generate_individual_spectra(py_dict)
            # Convert wl_axis from A to nm. Also need to multiply spectra by 10
            wl_axis = wl_axis / 10
            # Convert spectra from erg/s/A to erg/s/nm/cm**2
            dl = 10**(1 + py_dict['mu0']/5) * PSEC_TO_CM
            divisor = FOUR_PI * dl**2 * 0.1
            spectra = np.transpose(np.transpose(spectra.magnitude)/divisor)
            spectra_32 = spectra.astype(np.float32)
            del spectra
            print(f'len(wavelength axis): {len(wl_axis)}')
            print(f'shape of spectra array: {spectra_32.shape}')
            self._write_SED_batch(hp5_path, batch, py_dict['id'],
                                  wl_axis, spectra_32)

    def generate(self, hps):

        for hp in hps:
            self._generate_hp(hp)


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

    ####################
    def create_trilegal_flux_pixel(self, pixel, arrow_schema):
        output_filename = f'trilegal_flux_{pixel}.parquet'
        output_path = os.path.join(self._catalog_creator._output_dir,
                                   output_filename)

        object_list = self._cat.get_object_type_by_hp(pixel, 'trilegal')
        n_parallel = self._catalog_creator._flux_parallel

        colls = object_list.get_collections()
        writer = pq.ParquetWriter(output_path, arrow_schema)
        #  fields_needed = arrow_schema.names
        instrument_needed = ['lsst']
        rg_written = 0
        writer = None

        # Get all the objects in the pixel
        obj_collection = self._catalog_creator._cat.get_object_type_by_hp(
            pixel, 'trilegal').get_collections()[0]  # there only is one

        for c in colls:
            # trilegal_collection = c
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
                                                   l_bnd, u_bnd)
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
                #                          instrument_needed, lb, u))
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
