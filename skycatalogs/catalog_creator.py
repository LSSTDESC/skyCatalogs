import os
import sys
import re
import math
import logging
import yaml
import numpy as np
import numpy.ma as ma
import healpy
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from multiprocessing import Process, Pipe
from astropy.coordinates import SkyCoord
import sqlite3
from .utils.common_utils import print_date
from .utils.sed_tools import TophatSedFactory, get_star_sed_path
from .utils.config_utils import create_config, assemble_SED_models
from .utils.config_utils import assemble_MW_extinction, assemble_cosmology, assemble_object_types, assemble_provenance, write_yaml
from .utils.parquet_schema_utils import make_galaxy_schema, make_galaxy_flux_schema, make_star_flux_schema, make_pointsource_schema
from .utils.creator_utils import make_MW_extinction_av, make_MW_extinction_rv
from .objects.base_object import LSST_BANDS
from .objects.base_object import ObjectCollection

# from dm stack
from dustmaps.sfd import SFDQuery

"""
Code to create a sky catalog for particular object types
"""

__all__ = ['CatalogCreator']

_MW_rv_constant = 3.1

def _generate_sed_path(ids, subdir, cmp):
    '''
    Generate paths (e.g. relative to SIMS_SED_LIBRARY_DIR) for galaxy component
    SED files
    Parameters
    ----------
    ids        list of galaxy ids
    subdir    user-supplied part of path
    cmp      component for which paths should be generated

    returns
    -------
    A list of strings.  The entries in the list have the form
    <subdir>/<cmp>_<id>.txt
    '''
    r = [f'{subdir}/{cmp}_{id}.txt' for id in ids]
    return r

def _get_tophat_info(columns):
    '''
    Parameters
    ----------
    columns    list of column names including the ones with per-tophat information

    Returns
    -------
    sed_bins        List of  the tophat bins, sorted by "start" (left edge) value
    sed_bulge_names To be fetched from input catalog
    sed_disk_names  To be fetched from input catalog
    '''
    tophat_bulge_re = r'sed_(?P<start>\d+)_(?P<width>\d+)_bulge'
    tophat_disk_re = r'sed_(?P<start>\d+)_(?P<width>\d+)_disk'

    # Save all the start, width values
    sed_bulge_names = [i for i in columns if (i.startswith('sed') and
                                              i.endswith('bulge'))]
    sed_disk_names = [i for i in columns if (i.startswith('sed') and
                                             i.endswith('disk'))]

    sed_bins = [[int(re.match(tophat_bulge_re, s)['start']),
                 int(re.match(tophat_bulge_re, s)['width'])] for s in sed_bulge_names]

    # Sort sed by value for start. Finishes off work on sed_bins
    def _bin_start_key(start_width):
        return start_width[0]
    sed_bins.sort(key=_bin_start_key)

    # Moving on to tophat_fetch
    def _sed_bulge_key(s):
        return int(re.match(tophat_bulge_re, s)['start'])
    def _sed_disk_key(s):
        return int(re.match(tophat_disk_re, s)['start'])

    # Sort into increaing order by start wavelength
    sed_bulge_names.sort(key=_sed_bulge_key)
    sed_disk_names.sort(key=_sed_disk_key)

    return sed_bins, sed_bulge_names, sed_disk_names

_nside_allowed = 2**np.arange(15)

def _find_subpixels(pixel, subpixel_nside, pixel_nside=32, nest=False):
    '''
    Return list of pixels of specified nside inside a given pixel
    Parameters
    ----------
    pixel           int       the id of the input pixel
    subpixel_nside  int       nside for subpixels
    pixel_nside     int       nside of original pixel (default=32)
    nest            boolean   True if pixel ordering for original pixel is nested
                              (default = False)
    Returns
    -------
    List of subpixel ids (nested ordering iff original was).  If subpixel resolution
    is no better than original, just return original pixel id
    '''
    if pixel_nside not in _nside_allowed:
        raise ValueError(f'Disallowed pixel nside value {pixel_nside}')
    if subpixel_nside not in _nside_allowed:
        raise ValueError(f'Disallowed subpixel nside value {subpixel_nside}')
    if pixel_nside >= subpixel_nside:
        return [pixel]

    if not nest:
        nest_pixel = healpy.ring2nest(pixel_nside, pixel)
    else:
        nest_pixel = pixel

    def _next_level(pixel):
        return [4 * pixel, 4 * pixel + 1, 4 * pixel + 2, 4 * pixel + 3]

    pixels = [nest_pixel]
    current_nside = pixel_nside
    while current_nside < subpixel_nside:
        pixels = [pix for p in pixels for pix in _next_level(p)]
        current_nside = current_nside * 2

    if nest:
        return pixels
    else:
        return [healpy.nest2ring(subpixel_nside, p) for p in pixels]

def _generate_subpixel_masks(ra, dec, subpixels, nside=32):
    '''
    Given ra, dec values for objects within a particular pixel and a list of its
    subpixels for some greater value of nside, return dict with subpixel ids as
    keys and values a mask which masks off all values except those belonging to subpixel

    Parameters
    ----------
    ra         float array   ra for all objects in a particular pixel
    dec        float array   dec for all objects in a particular pixel
    subpixels  int array     pixels for which masks should be generated
    nside      int           healpix ordering parameter

    Returns
    -------
    masks      dict          mask for each subpixel, keyed by subpixel id
    '''

    pix_id = np.array(healpy.pixelfunc.ang2pix(nside, ra, dec, lonlat=True))
    masks = dict()
    for p in subpixels:
        m = np.array(pix_id != p)
        masks[p] = m

    return masks


# Collection of galaxy objects for current row group, current pixel
# Used while doing flux computation

def _do_galaxy_flux_chunk(send_conn, galaxy_collection, l_bnd, u_bnd):
    '''
    output connection
    l_bnd, u_bnd     demarcates slice to process

    returns
                    dict with keys id, lsst_flux_u, ... lsst_flux_y
    '''
    out_dict = {}

    o_list = galaxy_collection[l_bnd : u_bnd]
    out_dict['galaxy_id'] = [o.get_native_attribute('galaxy_id') for o in o_list]
    all_fluxes =  [o.get_LSST_fluxes(as_dict=False) for o in o_list]
    all_fluxes_transpose = zip(*all_fluxes)
    for i, band in enumerate(LSST_BANDS):
        v = all_fluxes_transpose.__next__()
        out_dict[f'lsst_flux_{band}'] = v

    if send_conn:
        send_conn.send(out_dict)
    else:
        return out_dict

class CatalogCreator:
    def __init__(self, parts, area_partition=None, skycatalog_root=None,
                 catalog_dir='.', galaxy_truth=None,
                 star_truth=None, sn_truth=None,
                 config_path=None, catalog_name='skyCatalog',
                 output_type='parquet', mag_cut=None,
                 sed_subdir='galaxyTopHatSED', knots_mag_cut=27.0,
                 knots=True, logname='skyCatalogs.creator',
                 pkg_root=None, skip_done=False, flux_only=False,
                 main_only=False, flux_parallel=16, galaxy_nside=32,
                 galaxy_stride=1000000, provenance=None,
                 dc2=False, sn_object_type='sncosmo'):
        """
        Store context for catalog creation

        Parameters
        ----------
        parts           Segments for which catalog is to be generated. If
                        partition type is HEALpix, parts will be a collection
                        of HEALpix pixels
        area_partition  Dict characterizing partition globally (e.g. HEALpix,
                        nside=<something>) or None.  Defaults to None,
                        meaning partition is on per-source-type basis
        skycatalog_root Typically absolute directory containing one or
                        more subdirectories for sky catalogs. Defaults
                        to current directory
        catalog_dir     Directory relative to skycatalog_root where catalog
                        will be written.  Defaults to '.'
        galaxy_truth    GCRCatalogs name for galaxy truth (e.g. cosmoDC2)
        config_path     Where to write config file. Default is data
                        directory.
        catalog_name    If a config file is written this value is saved
                        there and is also part of the filename of the
                        config file.
        output_type     A format.  For now only parquet is supported
        mag_cut         If not None, exclude galaxies with mag_r > mag_cut
        sed_subdir      In instcat entry, prepend this value to SED filename
        knots_mag_cut   No knots for galaxies with i_mag > cut
        knots           If True include knots
        logname         logname for Python logger
        pkg_root        defaults to one level up from __file__
        skip_done       If True, skip over files which already exist. Otherwise
                        (by default) overwrite with new version.
                        Output info message in either case if file exists.
        flux_only       Only create flux files, not main files
        main_only       Only create main files, not flux files
        flux_parallel   Number of processes to divide work of computing fluxes
        galaxy_nside    Healpix configuration value "nside" for galaxy output
        galaxy_stride   Max number of rows per galaxy row group
        provenance      Whether to write per-output-file git repo provenance
        dc2             Whether to adjust values to provide input comparable
                        to that for the DC2 run
        sn_object_type  Which object type to use for SNe.

        Might want to add a way to specify template for output file name
        and template for input sedLookup file name.
        """

        _cosmo_cat = 'cosmodc2_v1.1.4_image_addon_knots'
        _star_db = '/global/cfs/cdirs/lsst/groups/SSim/DC2/dc2_stellar_healpixel.db'
        _sn_db = '/global/cfs/cdirs/lsst/groups/SSim/DC2/cosmoDC2_v1.1.4/sne_cosmoDC2_v1.1.4_MS_DDF_healpix.db'

        self._galaxy_stride = galaxy_stride
        if pkg_root:
            self._pkg_root = pkg_root
        else:
            self._pkg_root = os.path.join(os.path.dirname(__file__),
                                          '..')

        self._global_partition = area_partition
        if area_partition is not None:
            if area_partition['type'] != 'healpix':
                raise NotImplementedError(f'CatalogCreator: Unknown partition type {area_partition["type"]} ')

        if output_type != 'parquet':
            raise NotImplementedError(f'CatalogCreator: Output type {output_type} not supported')

        self._galaxy_truth = galaxy_truth
        if galaxy_truth is None:
            self._galaxy_truth = _cosmo_cat

        self._sn_truth = sn_truth
        if self._sn_truth is None:
            self._sn_truth = _sn_db

        self._sn_object_type = sn_object_type

        self._star_truth = star_truth
        if self._star_truth is None:
            self._star_truth = _star_db
        self._cat = None

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

        self._written_config = None
        self._config_path = config_path
        self._catalog_name = catalog_name

        self._output_type = output_type
        self._mag_cut = mag_cut
        self._sed_subdir = sed_subdir
        self._knots_mag_cut = knots_mag_cut
        self._knots = knots
        self._logname = logname
        self._logger = logging.getLogger(logname)
        self._skip_done = skip_done
        self._flux_only = flux_only
        self._main_only = main_only
        self._flux_parallel = flux_parallel
        self._galaxy_nside = galaxy_nside
        self._provenance = provenance
        self._dc2 = dc2

        self._obs_sed_factory = None

    def _make_tophat_columns(self, dat, names, cmp):
        '''
        Create columns sed_val_cmp, cmp_magnorm where cmp is one of "disk", "bulge", "knots"
        Parameters
        ----------
        dat          Data read from input galaxy catalog. Includes keys for everything in names
                     plus entry for redshiftHubble
        names        Names of SED columns for this component
        cmp          Component name

        Returns
        -------
        Add keys  sed_val_cmp, cmp_magnorm to dat.  Values are associated data. Return dat.
        '''
        sed_vals = (np.array([dat[k] for k in names]).T).tolist()
        dat['sed_val_' + cmp] = sed_vals
        dat[cmp + '_magnorm'] = [self._obs_sed_factory.magnorm(s, z) for (s, z)\
                                 in zip(sed_vals, dat['redshiftHubble'])]
        for k in names:
            del(dat[k])
        return dat

    def create(self, catalog_type):
        """
        Create catalog of specified type, using stored context.

        Parameters
        ----------
        catalog_type   string    Currently 'galaxy' and 'pointsource' are
                                 the only values allowed
        Return
        ------
        None
        """
        if catalog_type == ('galaxy'):
            if not self._flux_only:
                self.create_galaxy_catalog()
            if not self._main_only:
                self.create_galaxy_flux_catalog()
        elif catalog_type == ('pointsource'):
            if not self._flux_only:
                self.create_pointsource_catalog()
            if not self._main_only:
                self.create_pointsource_flux_catalog()
        else:
            raise NotImplemented(f'CatalogCreator.create: unsupported catalog type {catalog_type}')

    def create_galaxy_catalog(self):
        """
        Create the 'main' galaxy catalog, including everything except
        LSST fluxes

        Returns
        -------
        None

        """
        import GCRCatalogs

        self._cat = None

        gal_cat = GCRCatalogs.load_catalog(self._galaxy_truth)

        # Save cosmology in case we need to write parameters out later
        self._cosmology = gal_cat.cosmology

        arrow_schema = make_galaxy_schema(self._logname,
                                              self._sed_subdir,
                                              self._knots)

        for p in self._parts:
            self._logger.info(f'Starting on pixel {p}')
            self.create_galaxy_pixel(p, gal_cat, arrow_schema)
            self._logger.info(f'Completed pixel {p}')

        # Now make config.   We need it for computing LSST fluxes for
        # the second part of the galaxy catalog
        if self._skip_done:
            config_path = self.write_config(path_only=True)
            if os.path.exists(config_path):
                self._logger.info('Will not overwrite existing config file')
                return
        self.write_config()

    def _write_subpixel(self, dat=None, output_path=None, arrow_schema=None,
                        to_rename=dict(), stride=100000):
        '''
        Write out data for a single healpixel, single source type
        Parameters
        ----------
        dat           dict    The data to write out. Column names are the keys;
                              values are numpy array
        output_path   string  path to output file
        arrow_schema          Schema for output parquet file
        stride        int     number of rows to include in a row group
        to_rename     dict    Associate input column name with output name
                              if they differ
        '''
        dlen = 0
        for val in dat.values():
            dlen = len(val)
            break
        if dlen == 0: return

        last_row_ix = dlen - 1
        u_bnd = min(stride, dlen)
        l_bnd = 0
        rg_written = 0
        writer = None

        while u_bnd > l_bnd:
            out_dict = {k : dat[k][l_bnd : u_bnd] for k in dat if k not in to_rename}
            for k in to_rename:
                out_dict[to_rename[k]] = dat[k][l_bnd : u_bnd]
            out_df = pd.DataFrame.from_dict(out_dict)
            out_table = pa.Table.from_pandas(out_df, schema=arrow_schema)
            if not writer:
                writer = pq.ParquetWriter(output_path, arrow_schema)

            writer.write_table(out_table)
            rg_written += 1
            l_bnd = u_bnd
            u_bnd = min(l_bnd + stride, last_row_ix + 1)

        writer.close()
        self._logger.debug(f'# row groups written to {output_path}: {rg_written}')

    def create_galaxy_pixel(self, pixel, gal_cat, arrow_schema):
        """
        Parameters
        ----------
        pixel           Pixel for which catalog(s) is(are) to be generated.
                        Note: input pixels are assumed to use nside=32. Output
                        pixels may be finer
        gal_cat         GCRCatalogs-loaded galaxy truth (e.g. cosmoDC2)
        arrow_schema    schema to use for output file
        """

        # The typical galaxy input file to date (cosmoDC2 or diffsky)
        # is partitioned into nside=32 healpixels. This code will only
        # output pixels of of the same size or smaller.
        if self._galaxy_nside > 32:
            out_pixels = _find_subpixels(pixel, self._galaxy_nside)
            self._logger.debug(f'For nside={self._galaxy_nside} subpixels are')
            self._logger.debug(out_pixels)
        else:
            out_pixels = [pixel]
        self._out_pixels = out_pixels
        skip_count = 0
        for p in out_pixels:
            output_path = os.path.join(self._output_dir, f'galaxy_{p}.parquet')
            if os.path.exists(output_path):
                if self._skip_done:
                    self._logger.info(f'Skipping regeneration of {output_path}')
                    skip_count = skip_count + 1

        if skip_count == len(out_pixels):
            return

        # Number of rows to include in a row group
        stride = self._galaxy_stride

        hp_filter = [f'healpix_pixel=={pixel}']
        if self._mag_cut:
            r_mag_name = 'mag_r_lsst'
            mag_cut_filter = [f'{r_mag_name} < {self._mag_cut}']

        # to_fetch = all columns of interest in gal_cat
        non_sed = ['galaxy_id', 'ra', 'dec', 'redshift', 'redshiftHubble',
                   'peculiarVelocity', 'shear_1', 'shear_2',
                   'convergence',
                   'size_bulge_true', 'size_minor_bulge_true', 'sersic_bulge',
                   'size_disk_true', 'size_minor_disk_true', 'sersic_disk']
        if self._dc2:
            non_sed += ['ellipticity_1_disk_true_dc2',
                        'ellipticity_2_disk_true_dc2',
                        'ellipticity_1_bulge_true_dc2',
                        'ellipticity_2_bulge_true_dc2']
        else:
            non_sed += ['ellipticity_1_disk_true',
                        'ellipticity_2_disk_true',
                        'ellipticity_1_bulge_true',
                        'ellipticity_2_bulge_true']

        if self._knots:
            non_sed += ['knots_flux_ratio', 'n_knots', 'mag_i_lsst']

        # Find sed bin definition and all the tophat quantities needed
        q = gal_cat.list_all_quantities()

        sed_bins, sed_bulge_names, sed_disk_names = _get_tophat_info(q)
        self._sed_bins = sed_bins
        self._obs_sed_factory = TophatSedFactory(self._sed_bins,
                                                 assemble_cosmology(self._cosmology))

        #Fetch the data
        to_fetch = non_sed + sed_bulge_names + sed_disk_names

        # df is not a dataframe!  It's just a dict
        if not self._mag_cut:
            df = gal_cat.get_quantities(to_fetch, native_filters=hp_filter)
        else:
            df = gal_cat.get_quantities(to_fetch + [r_mag_name],
                                        native_filters=hp_filter,
                                        filters=mag_cut_filter)

        df['MW_rv'] = make_MW_extinction_rv(df['ra'], df['dec'])
        df['MW_av'] = make_MW_extinction_av(df['ra'], df['dec'])
        self._logger.debug('Made extinction')

        # Some columns need to be renamed
        to_rename = {'redshiftHubble' : 'redshift_hubble',
                     'peculiarVelocity' : 'peculiar_velocity'}
        if self._dc2:
            to_rename['ellipticity_1_disk_true_dc2'] = 'ellipticity_1_disk_true'
            to_rename['ellipticity_2_disk_true_dc2'] = 'ellipticity_2_disk_true'
            to_rename['ellipticity_1_bulge_true_dc2'] = 'ellipticity_1_bulge_true'
            to_rename['ellipticity_2_bulge_true_dc2'] = 'ellipticity_2_bulge_true'

        if self._sed_subdir:
            #  Generate full paths for disk and bulge SED files, even though
            #  we don't actually write the files here
            df['bulge_sed_file_path'] =\
                _generate_sed_path(df['galaxy_id'], self._sed_subdir, 'bulge')
            df['disk_sed_file_path'] =\
                _generate_sed_path(df['galaxy_id'], self._sed_subdir, 'disk')

        if self._knots:
            # adjust disk sed; create knots sed
            sed_knot_names = [i.replace('disk', 'knots') for i in sed_disk_names]
            eps = np.finfo(np.float32).eps
            mag_mask = np.where(np.array(df['mag_i_lsst']) > self._knots_mag_cut,0, 1)
            self._logger.debug(f'Count of mags <=  cut (so adjustment performed: {np.count_nonzero(mag_mask)}')

            for d_name, k_name in zip(sed_disk_names, sed_knot_names):
                df[k_name] = mag_mask * np.clip(df['knots_flux_ratio'], None, 1-eps) * df[d_name]
                df[d_name] = np.where(np.array(df['mag_i_lsst']) > self._knots_mag_cut, 1,
                                      np.clip(1 - df['knots_flux_ratio'], eps, None)) * df[d_name]

        if len(self._out_pixels) > 1:
            subpixel_masks = _generate_subpixel_masks(df['ra'], df['dec'],  self._out_pixels, nside=self._galaxy_nside)
        else:
            subpixel_masks = {pixel : None}


        for p, val in subpixel_masks.items():
            output_path = os.path.join(self._output_dir, f'galaxy_{p}.parquet')
            if os.path.exists(output_path):
                if not self._skip_done:
                    self._logger.info(f'Overwriting file {output_path}')
                else:
                    continue

            if val is not None:
                compressed = dict()
                for k in df:
                    compressed[k] = ma.array(df[k], mask=val).compressed()

                compressed = self._make_tophat_columns(compressed, sed_disk_names, 'disk')
                compressed = self._make_tophat_columns(compressed, sed_bulge_names, 'bulge')
                if self._knots:
                    compressed = self._make_tophat_columns(compressed, sed_knot_names, 'knots')

                self._write_subpixel(dat=compressed, output_path=output_path,
                                     arrow_schema=arrow_schema,
                                     stride=stride, to_rename=to_rename)
            else:
                df = self._make_tophat_columns(df, sed_disk_names, 'disk')
                df = self._make_tophat_columns(df, sed_bulge_names, 'bulge')
                if self._knots:
                    df = self._make_tophat_columns(df, sed_knot_names, 'knots')
                self._write_subpixel(dat=df, output_path=output_path,
                                     arrow_schema=arrow_schema,
                                     stride=stride, to_rename=to_rename)

        if self._provenance == 'yaml':
            self.write_provenance_file(output_path)

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

        from .skyCatalogs import open_catalog, SkyCatalog

        self._gal_flux_schema = make_galaxy_flux_schema(self._logname)

        if not config_file:
            config_file = self.write_config(path_only=True)
        if not self._cat:
            self._cat = open_catalog(config_file,
                                     skycatalog_root=self._skycatalog_root)

        # Might also open a main file and read its metadata to be sure
        # the value for stride is correct.
        # Especially need to do this if we want to support making flux
        # files in a separate job activation.

        self._flux_template = self._cat.raw_config['object_types']['galaxy']['flux_file_template']

        self._logger.info('Creating galaxy flux files')
        for p in self._parts:
            self._logger.info(f'Starting on pixel {p}')
            self._create_galaxy_flux_pixel(p)
            self._logger.info(f'Completed pixel {p}')


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
        object_list = self._cat.get_object_type_by_hp(pixel, 'galaxy')
        writer = None
        global _galaxy_collection

        for object_coll in object_list.get_collections():
            _galaxy_collection = object_coll
            # prefetch everything we need.
            for att in ['galaxy_id', 'shear_1', 'shear_2', 'convergence',
                        'redshift_hubble', 'MW_av', 'MW_rv', 'sed_val_bulge',
                        'sed_val_disk', 'sed_val_knots'] :
                v = object_coll.get_native_attribute(att)
            l_bnd = 0
            u_bnd = len(object_coll)
            rg_written = 0

            self._logger.debug(f'Handling range {l_bnd} up to {u_bnd}')

            out_dict = {'galaxy_id': [], 'lsst_flux_u' : [],
                        'lsst_flux_g' : [], 'lsst_flux_r' : [],
                        'lsst_flux_i' : [], 'lsst_flux_z' : [],
                        'lsst_flux_y' : []}

            n_parallel = self._flux_parallel

            if n_parallel == 1:
                n_per = u_bnd - l_bnd
            else:
                n_per = int((u_bnd - l_bnd + n_parallel)/n_parallel)
            l = l_bnd
            u = min(l_bnd + n_per, u_bnd)
            readers = []

            if n_parallel == 1:
                out_dict = _do_galaxy_flux_chunk(None, _galaxy_collection,
                                                 l, u)
            else:
                # Expect to be able to do about 1500/minute/process
                tm = max(int((n_per*60)/500), 5)  # Give ourselves a cushion
                self._logger.info(f'Using timeout value {tm} for {n_per} sources')
                p_list = []
                for i in range(n_parallel):
                    conn_rd, conn_wrt = Pipe(duplex=False)
                    readers.append(conn_rd)

                    # For debugging call directly
                    proc = Process(target=_do_galaxy_flux_chunk,
                                   name=f'proc_{i}',
                                   args=(conn_wrt, _galaxy_collection,l, u))
                    proc.start()
                    p_list.append(proc)
                    l = u
                    u = min(l + n_per, u_bnd)

                self._logger.debug('Processes started') # outside for loop
                for i in range(n_parallel):
                    ready = readers[i].poll(tm)
                    if not ready:
                        self._logger.error(f'Process {i} timed out after {tm} sec')
                        sys.exit(1)
                    dat = readers[i].recv()   # lines up with if
                    for k in ['galaxy_id', 'lsst_flux_u', 'lsst_flux_g',
                              'lsst_flux_r', 'lsst_flux_i', 'lsst_flux_z',
                              'lsst_flux_y']:
                        out_dict[k] += dat[k]
                for p in p_list:    # indent same as "for i in range(.."
                    p.join()

            out_df = pd.DataFrame.from_dict(out_dict) # outdent from for
            out_table = pa.Table.from_pandas(out_df,
                                             schema=self._gal_flux_schema)

            if not writer:
                writer = pq.ParquetWriter(output_path, self._gal_flux_schema)
            writer.write_table(out_table)

            rg_written +=1

        writer.close()
        self._logger.debug(f'# row groups written to flux file: {rg_written}')
        if self._provenance == 'yaml':
            self.write_provenance_file(output_path)

    def create_pointsource_catalog(self):

        """
        Parameters
        ----------
        star_truth      Where to find star parameters. If None use default
        sn_truth       Where to find SN parameters.  If None, use default

        Might want to add a way to specify template for output file name

        Returns
        -------
        None
        """
        arrow_schema = make_pointsource_schema()
        #  Need a way to indicate which object types to include; deal with that
        #  later.  For now, default is stars + sn
        for p in self._parts:
            self._logger.debug(f'Point sources. Starting on pixel {p}')
            self.create_pointsource_pixel(p, arrow_schema,
                                          star_cat=self._star_truth,
                                          sn_cat=self._sn_truth)
            self._logger.debug(f'Completed pixel {p}')

    def create_pointsource_pixel(self, pixel, arrow_schema, star_cat=None,
                             sn_cat=None):
        if not star_cat and not sn_cat:
            self._logger.info('No point source inputs specified')
            return

        output_filename = f'pointsource_{pixel}.parquet'
        output_path = os.path.join(self._output_dir, output_filename)

        if os.path.exists(output_path):
            if not self._skip_done:
                self._logger.info(f'Overwriting file {output_path}')
            else:
                self._logger.info(f'Skipping regeneration of {output_path}')
                return

        writer = pq.ParquetWriter(output_path, arrow_schema)

        if star_cat:
            # Get data for this pixel
            cols = ','.join(['format("%s",simobjid) as id', 'ra', 'decl as dec',
                             'magNorm as magnorm', 'mura', 'mudecl as mudec',
                             'radialVelocity as radial_velocity', 'parallax',
                             'sedFilename as sed_filepath', 'ebv'])
            q = f'select {cols} from stars where hpid={pixel} '
            with sqlite3.connect(star_cat) as conn:
                star_df = pd.read_sql_query(q, conn)

            star_df['sed_filepath'] = get_star_sed_path(star_df['sed_filepath'])

            nobj = len(star_df['id'])
            self._logger.debug(f'Found {nobj} stars')
            star_df['object_type'] = np.full((nobj,), 'star')
            star_df['host_galaxy_id'] = np.zeros((nobj,), np.int64())

            star_df['MW_rv'] = np.full((nobj,), _MW_rv_constant, np.float32())

            # NOTE MW_av calculation for stars does not use SFD dust map
            star_df['MW_av'] = star_df['ebv'] * _MW_rv_constant

            star_df['variability_model'] = np.full((nobj,), '')
            star_df['salt2_params'] = np.full((nobj,), None)
            out_table = pa.Table.from_pandas(star_df, schema=arrow_schema)
            self._logger.debug('Created arrow table from star dataframe')

            #writer = pq.ParquetWriter(output_path, arrow_schema)
            writer.write_table(out_table)

        if sn_cat:
            # Get data for this pixel
            cols = ','.join(['snid_in as id', 'snra_in as ra',
                             'sndec_in as dec', 'galaxy_id as host_galaxy_id'])

            params = ','.join(['z_in as z', 't0_in as t0, x0_in as x0',
                             'x1_in as x1', 'c_in as c'])

            q1 = f'select {cols} from sne_params where hpid={pixel} '
            q2 = f'select {params} from sne_params where hpid={pixel} '
            with sqlite3.connect(sn_cat) as conn:
                sn_df = pd.read_sql_query(q1, conn)
                params_df = pd.read_sql_query(q2, conn)

            nobj = len(sn_df['ra'])
            sn_df['object_type'] = np.full((nobj,), self._sn_object_type)

            sn_df['MW_rv'] = make_MW_extinction_rv(sn_df['ra'], sn_df['dec'])
            sn_df['MW_av'] = make_MW_extinction_ra(np.array(sn_df['ra']),
                                                   np.array(sn_df['dec']))

            # Add fillers for columns not relevant for sn
            sn_df['sed_filepath'] = np.full((nobj),'')
            sn_df['magnorm'] = np.full((nobj,), None)
            sn_df['mura'] = np.full((nobj,), None)
            sn_df['mudec'] = np.full((nobj,), None)
            sn_df['radial_velocity'] = np.full((nobj,), None)
            sn_df['parallax'] = np.full((nobj,), None)
            sn_df['variability_model'] = np.full((nobj,), 'salt2_extended')

            # Form array of struct from params_df
            sn_df['salt2_params'] = params_df.to_records(index=False)
            out_table = pa.Table.from_pandas(sn_df, schema=arrow_schema)
            self._logger.debug('Created arrow table from sn dataframe')

            writer.write_table(out_table)

        writer.close()
        if self._provenance == 'yaml':
            self.write_provenance_file(output_path)

        return

    def create_pointsource_flux_catalog(self, config_file=None):
        '''
        Create a second file per healpixel containing just id and
        LSST fluxes.  Use information in the main file to compute fluxes

        Parameters
        ----------
        Path to config created in first stage so we can find the main
        galaxy files.

        Return
        ------
        None
        '''

        from .skyCatalogs import open_catalog, SkyCatalog

        self._ps_flux_schema = make_star_flux_schema(self._logname)
        if not config_file:
            config_file = self.write_config(path_only=True)

        # Always open catalog. If it was opened for galaxies earlier
        # it won't know about star files.
        self._cat = open_catalog(config_file,
                                 skycatalog_root=self._skycatalog_root)

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

        object_list = self._cat.get_object_type_by_hp(pixel, 'star')
        last_row_ix = len(object_list) - 1
        writer = None

        # Write out as a single rowgroup as was done for main catalog
        l_bnd = 0
        u_bnd = last_row_ix + 1

        o_list = object_list[l_bnd : u_bnd]
        self._logger.debug(f'Handling range {l_bnd} up to {u_bnd}')
        out_dict = {}
        out_dict['id'] = [o.get_native_attribute('id') for o in o_list]
        all_fluxes =  [o.get_LSST_fluxes(as_dict=False) for o in o_list]
        all_fluxes_transpose = zip(*all_fluxes)
        for i, band in enumerate(LSST_BANDS):
            self._logger.debug(f'Band {band} is number {i}')
            v = all_fluxes_transpose.__next__()
            out_dict[f'lsst_flux_{band}'] = v
            if i == 1:
                self._logger.debug(f'Len of flux column: {len(v)}')
                self._logger.debug(f'Type of flux column: {type(v)}')
        out_df = pd.DataFrame.from_dict(out_dict)
        out_table = pa.Table.from_pandas(out_df,
                                         schema=self._ps_flux_schema)
        if not writer:
            writer = pq.ParquetWriter(output_path, self._ps_flux_schema)
        writer.write_table(out_table)

        writer.close()
        ##self._logger.debug(f'# row groups written to flux file: {rg_written}')
        if self._provenance == 'yaml':
            self.write_provenance_file(output_path)

    def write_config(self, overwrite=False, path_only=False):
        '''
        Parameters
        ----------
        overwrite   boolean default False.   If true, overwrite existing
                    config of the same name
        path_only   If true, just return the path; don't write anything

        Returns
        -------
        Path to would-be config file if path_only is True;
        else None

        Side-effects
        ------------
        Save path to config file written as instance variable

        '''
        if not self._config_path:
            self._config_path = self._output_dir

        if path_only:
            return os.path.join(self._config_path, self._catalog_name + '.yaml')


        config = create_config(self._catalog_name, self._logname)
        if self._global_partition is not None:
            config.add_key('area_partition', self._area_partition)
        config.add_key('skycatalog_root', self._skycatalog_root)
        config.add_key('catalog_dir' , self._catalog_dir)

        config.add_key('SED_models',
                       assemble_SED_models(self._sed_bins))
        config.add_key('MW_extinction_values', assemble_MW_extinction())
        config.add_key('Cosmology', assemble_cosmology(self._cosmology))
        config.add_key('object_types',
                       assemble_object_types(self._pkg_root,
                                             galaxy_nside=self._galaxy_nside))

        config.add_key('galaxy_magnitude_cut', self._mag_cut)
        config.add_key('knots_magnitude_cut', self._knots_mag_cut)

        inputs = {'galaxy_truth' : self._galaxy_truth}
        if self._sn_truth:
            inputs['sn_truth'] = self._sn_truth
        if self._star_truth:
            inputs['star_truth'] = self._star_truth
        config.add_key('provenance', assemble_provenance(self._pkg_root,
                                                         inputs=inputs))

        self._written_config = config.write_config(self._config_path,
                                                   overwrite=overwrite)

    def write_provenance_file(self, datafile_path):
        '''
        Write git provenance to a yaml file with name derived from a
        just-written datafile name
        '''
        outpath = datafile_path.rsplit('.', 1)[0] + '_provenance.yaml'

        prov = assemble_provenance(self._pkg_root, inputs=None)
        write_yaml(prov, outpath)
