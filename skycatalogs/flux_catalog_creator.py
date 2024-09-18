import os
import sys
import logging
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from multiprocessing import Process, Pipe
from .utils.config_utils import assemble_file_metadata
from .utils.parquet_schema_utils import make_galaxy_flux_schema
from .utils.parquet_schema_utils import make_star_flux_schema
from .objects.base_object import LSST_BANDS
from .objects.base_object import ROMAN_BANDS
from .sso_catalog_creator import SsoFluxCatalogCreator

"""
Code to create flux sky catalogs for particular object types
"""

__all__ = ['FluxCatalogCreator']

_MW_rv_constant = 3.1
_nside_allowed = 2**np.arange(15)


# Collection of galaxy objects for current row group, current pixel
# Used while doing flux computation
def _do_flux_chunk(send_conn, object_collection, instrument_needed,
                   l_bnd, u_bnd, id_column):
    '''
    send_conn           output connection
    object_collection   objects whose fluxes are to be computed
    instrument_needed   determines which fluxes (LSST only or also Roman)
                        to compute
    l_bnd, u_bnd        demarcates slice to process
    id_column           column name

    returns
                    dict with keys id, lsst_flux_u, ... lsst_flux_y
    '''
    out_dict = {}

    o_list = object_collection[l_bnd: u_bnd]
    out_dict[id_column] = [o.get_native_attribute(id_column) for o in o_list]
    if 'lsst' in instrument_needed:
        all_fluxes = [o.get_LSST_fluxes(as_dict=False) for o in o_list]
        all_fluxes_transpose = zip(*all_fluxes)
        colnames = [f'lsst_flux_{band}' for band in LSST_BANDS]
        flux_dict = dict(zip(colnames, all_fluxes_transpose))
        out_dict.update(flux_dict)

    if 'roman' in instrument_needed:
        all_fluxes = [o.get_roman_fluxes(as_dict=False) for o in o_list]
        all_fluxes_transpose = zip(*all_fluxes)
        colnames = [f'roman_flux_{band}' for band in ROMAN_BANDS]
        flux_dict = dict(zip(colnames, all_fluxes_transpose))
        out_dict.update(flux_dict)

    if send_conn:
        send_conn.send(out_dict)
    else:
        return out_dict


class FluxCatalogCreator:
    def __init__(self, object_type, parts,
                 skycatalog_root=None,
                 catalog_dir='.',
                 config_path=None,
                 catalog_name='skyCatalog',
                 logname='skyCatalogs.creator',
                 pkg_root=None,
                 skip_done=False,
                 flux_parallel=16,
                 include_roman_flux=False,
                 sso_sed=None,
                 run_options=None):
        """
        Store context for catalog creation

        Parameters
        ----------
        object_type     'cosmodc2_galaxy', 'diffsky_galaxy', 'star', or 'sso'
        parts           Segments for which catalog is to be generated. If
                        partition type is HEALpix, parts will be a collection
                        of HEALpix pixels
        skycatalog_root Typically absolute directory containing one or
                        more subdirectories for sky catalogs. Defaults
                        to current directory
        catalog_dir     Directory relative to skycatalog_root where catalog
                        will be written.  Defaults to '.'
        config_path     Where to read config file. Default is data
                        directory.
        catalog_name    If a config file is written this value is saved
                        there and is also part of the filename of the
                        config file.
        logname         logname for Python logger
        pkg_root        defaults to one level up from __file__
        skip_done       If True, skip over files which already exist. Otherwise
                        (by default) overwrite with new version.
                        Output info message in either case if file exists.
        flux_parallel   Number of processes to divide work of computing fluxes
        # dc2             Whether to adjust values to provide input comparable
        #                to that for the DC2 run
        include_roman_flux Calculate and write Roman flux values
        sso_sed         Path to sed file to be used for all SSOs
        run_options     The options the outer script (create_sc.py) was
                        called with

        Might want to add a way to specify template for output file name
        and template for input sedLookup file name.
        """
        from .skyCatalogs import open_catalog

        self._object_type = object_type
        if object_type.endswith('_galaxy'):
            self._galaxy_type = object_type[:-7]  # len('_galaxy')

        if pkg_root:
            self._pkg_root = pkg_root
        else:
            self._pkg_root = os.path.join(os.path.dirname(__file__), '..')

        self._parts = parts
        if skycatalog_root:
            self._skycatalog_root = skycatalog_root
        else:
            self._skycatalog_root = './'
        if catalog_dir:
            self._catalog_dir = catalog_dir
        else:
            self._catalog_dir = '.'

        self._output_dir = os.path.join(self._skycatalog_root,
                                        self._catalog_dir)

        self._config_path = config_path
        self._catalog_name = catalog_name

        self._cat = open_catalog(self.get_config_file_path(),
                                 skycatalog_root=self._skycatalog_root)

        self._logname = logname
        self._logger = logging.getLogger(logname)
        self._skip_done = skip_done
        self._flux_parallel = flux_parallel
        self._include_roman_flux = include_roman_flux
        self._obs_sed_factory = None
        self._sso_creator = SsoFluxCatalogCreator(self)
        self._run_options = run_options
        self._tophat_sed_bins = None
        self._sed_gen = None

    def create(self):
        """
        Create catalog for our object_type, using stored context.

        Return
        ------
        None
        """
        object_type = self._object_type
        if object_type in {'cosmodc2_galaxy', 'diffsky_galaxy'}:
            self.create_galaxy_flux_catalog()
        elif object_type == ('star'):
            self.create_pointsource_flux_catalog()
        elif object_type == ('sso'):
            self._sso_creator.create_sso_flux_catalog()

        else:
            raise NotImplementedError(
                f'FluxCatalogCreator.create: unsupported object type {object_type}')

    def create_galaxy_flux_catalog(self, config_file=None):
        '''
        Create a second file per healpixel containing just galaxy id and
        LSST fluxes.  Use information in the main file to compute fluxes

        Parameters
        ----------
        Path to config created in first stage so we can find the main
        galaxy files.

        Return
        ------
        None
        '''

        # Throughput versions for fluxes included
        thru_v = {'lsst_throughputs_version': self._cat._lsst_thru_v}
        if self._include_roman_flux:
            thru_v['roman_throughputs_version'] = self._cat._roman_thru_v

        file_metadata = assemble_file_metadata(
            self._pkg_root,
            run_options=self._run_options,
            flux_file=True,
            throughputs_versions=thru_v)

        self._gal_flux_schema =\
            make_galaxy_flux_schema(self._logname, self._galaxy_type,
                                    include_roman_flux=self._include_roman_flux,
                                    metadata_input=file_metadata)
        self._gal_flux_needed = [field.name for field in self._gal_flux_schema]

        if self._galaxy_type == 'diffsky':
            type_field = 'diffsky_galaxy'
        else:
            type_field = 'galaxy'
        self._flux_template = self._cat.raw_config['object_types'][type_field]['flux_file_template']

        self._logger.info('Creating galaxy flux files')
        for p in self._parts:
            self._logger.info(f'Starting on pixel {p}')
            self._create_galaxy_flux_pixel(p)
            self._logger.info(f'Completed pixel {p}')

    def _get_needed_flux_attrs(self):
        if self._galaxy_type == 'diffsky':
            return ['galaxy_id', 'shear1', 'shear2', 'convergence',
                    'redshiftHubble', 'MW_av', 'MW_rv']
        else:
            return ['galaxy_id', 'shear_1', 'shear_2', 'convergence',
                    'redshift_hubble', 'MW_av', 'MW_rv', 'sed_val_bulge',
                    'sed_val_disk', 'sed_val_knots']

    def _create_galaxy_flux_pixel(self, pixel):
        '''
        Create a parquet file for a single healpix pixel containing only
        galaxy id and LSST fluxes

        Parameters
        ----------
        Pixel         int

        Return
        ------
        None
        '''

        # Would be better to obtain output filename from config or
        # at least from object_type
        output_filename = f'galaxy_flux_{pixel}.parquet'
        output_path = os.path.join(self._output_dir, output_filename)

        if os.path.exists(output_path):
            if not self._skip_done:
                self._logger.info(f'Overwriting file {output_path}')
            else:
                self._logger.info(f'Skipping regeneration of {output_path}')
                return

        # If there are multiple row groups, each is stored in a separate
        # object collection. Need to loop over them
        object_list = self._cat.get_object_type_by_hp(pixel, self._object_type)
        if len(object_list) == 0:
            self._logger.warning(f'Cannot create flux file for pixel {pixel} because main file does not exist or is empty')
            return

        if self._galaxy_type == 'diffsky':
            # Generate SEDs if necessary
            sed_output_path = os.path.join(self._output_dir,
                                           f'galaxy_sed_{pixel}.hdf5')
            if not os.path.exists(sed_output_path):
                if not self._sed_gen:
                    from .diffsky_sedgen import DiffskySedGenerator
                    # Default values are ok for all the diffsky-specific
                    # parameters: include_nonLSST_flux, sed_parallel, auto_loop,
                    # wave_ang_min, wave_ang_max, rel_err, n_per
                    self._sed_gen = DiffskySedGenerator(
                        logname=self._logname,
                        galaxy_truth=self._galaxy_truth,
                        output_dir=self._output_dir,
                        skip_done=True,
                        sky_cat=self._cat)

                self._sed_gen.generate_pixel(pixel)

        writer = None
        _instrument_needed = []
        rg_written = 0
        for field in self._gal_flux_needed:
            if 'lsst' in field and 'lsst' not in _instrument_needed:
                _instrument_needed.append('lsst')
            if 'roman' in field and 'roman' not in _instrument_needed:
                _instrument_needed.append('roman')

        for object_coll in object_list.get_collections():
            _galaxy_collection = object_coll
            # prefetch everything we need.
            for att in self._get_needed_flux_attrs():
                _ = object_coll.get_native_attribute(att)
            l_bnd = 0
            u_bnd = len(object_coll)

            self._logger.debug(f'Handling range {l_bnd} up to {u_bnd}')

            out_dict = {}
            for field in self._gal_flux_needed:
                out_dict[field] = []

            n_parallel = self._flux_parallel

            if n_parallel == 1:
                n_per = u_bnd - l_bnd
            else:
                n_per = int((u_bnd - l_bnd + n_parallel)/n_parallel)
            lb = l_bnd
            u = min(l_bnd + n_per, u_bnd)
            readers = []

            if n_parallel == 1:
                # For debugging call directly
                out_dict = _do_flux_chunk(None, _galaxy_collection,
                                          _instrument_needed, lb, u, 'galaxy_id')
            else:
                # Expect to be able to do about 1500/minute/process
                tm = max(int((n_per*60)/500), 5)  # Give ourselves a cushion
                self._logger.info(
                    f'Using timeout value {tm} for {n_per} sources')
                p_list = []
                for i in range(n_parallel):
                    conn_rd, conn_wrt = Pipe(duplex=False)
                    readers.append(conn_rd)
                    proc = Process(target=_do_flux_chunk,
                                   name=f'proc_{i}',
                                   args=(conn_wrt, _galaxy_collection,
                                         _instrument_needed, lb, u, 'galaxy_id'))
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
                    for field in self._gal_flux_needed:
                        out_dict[field] += dat[field]
                for p in p_list:
                    p.join()

            out_df = pd.DataFrame.from_dict(out_dict)
            out_table = pa.Table.from_pandas(out_df,
                                             schema=self._gal_flux_schema)

            if not writer:
                writer = pq.ParquetWriter(output_path, self._gal_flux_schema)
            writer.write_table(out_table)

            rg_written += 1

        writer.close()
        self._logger.debug(f'# row groups written to flux file: {rg_written}')

    def create_pointsource_flux_catalog(self, config_file=None):
        '''
        Create a second file per healpixel containing just id and
        LSST fluxes.  Use information in the main file to compute fluxes

        Parameters
        ----------
        Path to config created in first stage so we can find the main
        star files

        Return
        ------
        None
        '''
        thru_v = {'lsst_throughputs_version': self._cat._lsst_thru_v}
        file_metadata = assemble_file_metadata(
            self._pkg_root,
            run_options=self._run_options,
            flux_file=True,
            throughputs_versions=thru_v)

        self._ps_flux_schema = make_star_flux_schema(self._logname,
                                                     metadata_input=file_metadata)

        self._flux_template = self._cat.raw_config['object_types']['star']['flux_file_template']

        self._logger.info('Creating pointsource flux files')
        for p in self._parts:
            self._logger.info(f'Starting on pixel {p}')
            self._create_pointsource_flux_pixel(p)
            self._logger.info(f'Completed pixel {p}')

    def _create_pointsource_flux_pixel(self, pixel):
        '''
        Create a parquet file for a single healpix pixel containing only
        pointsource id and LSST fluxes

        Parameters
        ----------
        pixel         int

        Return
        ------
        None
        '''

        # For main catalog use self._cat
        # For schema use self._ps_flux_schema
        # output_template should be derived from value for flux_file_template
        #  in main catalog config.  Cheat for now
        output_filename = f'pointsource_flux_{pixel}.parquet'
        output_path = os.path.join(self._output_dir, output_filename)

        if os.path.exists(output_path):
            if not self._skip_done:
                self._logger.info(f'Overwriting file {output_path}')
            else:
                self._logger.info(f'Skipping regeneration of {output_path}')
                return
        n_parallel = self._flux_parallel

        object_list = self._cat.get_object_type_by_hp(pixel, 'star')
        writer = None
        instrument_needed = ['lsst']       # for now
        rg_written = 0
        fields_needed = self._ps_flux_schema.names

        for i in range(object_list.collection_count):
            _star_collection = object_list.get_collections()[i]

            l_bnd = 0
            u_bnd = len(_star_collection)
            out_dict = {}

            out_dict = {}
            for field in fields_needed:
                out_dict[field] = []

            if n_parallel == 1:
                n_per = u_bnd - l_bnd
            else:
                n_per = int((u_bnd - l_bnd + n_parallel)/n_parallel)

            lb = l_bnd
            u = min(l_bnd + n_per, u_bnd)
            readers = []

            if n_parallel == 1:
                # For debugging call directly
                out_dict = _do_flux_chunk(None, _star_collection,
                                          instrument_needed, lb, u, 'id')
            else:
                # Expect to be able to do about 1500/minute/process

                tm = max(int((n_per*60)/500), 5)  # Give ourselves a cushion
                self._logger.info(f'Using timeout value {tm} for {n_per} sources')
                p_list = []
                for i in range(n_parallel):
                    conn_rd, conn_wrt = Pipe(duplex=False)
                    readers.append(conn_rd)
                    proc = Process(target=_do_flux_chunk,
                                   name=f'proc_{i}',
                                   args=(conn_wrt, _star_collection,
                                         instrument_needed, lb, u, 'id'))
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
                                             schema=self._ps_flux_schema)

            if not writer:
                writer = pq.ParquetWriter(output_path, self._ps_flux_schema)
            writer.write_table(out_table)
            rg_written += 1

        writer.close()
        self._logger.debug(f'# row groups written to flux file: {rg_written}')

    def get_config_file_path(self):
        '''
        Return full path to config file.
        (self._config_path is path to *directory* containing config file)
        '''
        if not self._config_path:
            self._config_path = self._output_dir

        return os.path.join(self._config_path, self._catalog_name + '.yaml')
