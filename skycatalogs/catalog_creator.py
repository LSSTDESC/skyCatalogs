import os
import sys
import re
import logging
import numpy as np
import numpy.ma as ma
import healpy
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from multiprocessing import Process, Pipe
import sqlite3
from .utils.sed_tools import TophatSedFactory, get_star_sed_path
from .utils.sed_tools import generate_sed_path
from .utils.config_utils import create_config, assemble_SED_models
from .utils.config_utils import assemble_MW_extinction, assemble_cosmology
from .utils.config_utils import assemble_object_types, assemble_provenance
from .utils.config_utils import write_yaml
from .utils.star_parquet_input import _star_parquet_reader
from .utils.parquet_schema_utils import make_galaxy_schema
from .utils.parquet_schema_utils import make_galaxy_flux_schema
from .utils.parquet_schema_utils import make_star_flux_schema
from .utils.parquet_schema_utils import make_star_schema
from .utils.creator_utils import make_MW_extinction_av, make_MW_extinction_rv
from .objects.base_object import LSST_BANDS
from .objects.base_object import ROMAN_BANDS
from .sso_catalog_creator import SsoCatalogCreator

"""
Code to create a sky catalog for particular object types
"""

__all__ = ['CatalogCreator']

_MW_rv_constant = 3.1
_nside_allowed = 2**np.arange(15)


def _get_tophat_info(columns):
    '''
    Parameters
    ----------
    columns    list of column names including the ones with per-tophat info

    Returns
    -------
    sed_bins        List of  the tophat bins, sorted by "start" (left edge)
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


def _find_subpixels(pixel, subpixel_nside, pixel_nside=32, nest=False):
    '''
    Return list of pixels of specified nside inside a given pixel
    Parameters
    ----------
    pixel           int       the id of the input pixel
    subpixel_nside  int       nside for subpixels
    pixel_nside     int       nside of original pixel (default=32)
    nest            boolean   True if pixel ordering for original pixel is
                              nested (default = False)
    Returns
    -------
    List of subpixel ids (nested ordering iff original was).  If subpixel
    resolution is no better than original, just return original pixel id
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
    Given ra, dec values for objects within a particular pixel and a list of
    its subpixels for some greater value of nside, return dict with subpixel
    ids as keys and values a mask which masks off all values except those
    belonging to subpixel

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

def _do_galaxy_flux_chunk(send_conn, galaxy_collection, instrument_needed,
                          l_bnd, u_bnd):
    '''
    output connection
    l_bnd, u_bnd     demarcates slice to process
    instrument_needed List of which calculations should be done

    returns
                    dict with keys id, lsst_flux_u, ... lsst_flux_y
    '''
    out_dict = {}

    o_list = galaxy_collection[l_bnd: u_bnd]
    out_dict['galaxy_id'] = [o.get_native_attribute('galaxy_id') for o in o_list]
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


def _do_star_flux_chunk(send_conn, star_collection, instrument_needed,
                        l_bnd, u_bnd):
    '''
    send_conn         output connection, used to send results to
                      parent process
    star_collection   ObjectCollection. Information from main skyCatalogs
                      star file
    instrument_needed List of which calculations should be done. Currently
                      supported instrument names are 'lsst' and 'roman'
    l_bnd, u_bnd      demarcates slice to process

    returns
                    dict with keys id, lsst_flux_u, ... lsst_flux_y
    '''
    out_dict = {}

    o_list = star_collection[l_bnd: u_bnd]
    out_dict['id'] = [o.get_native_attribute('id') for o in o_list]
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


class CatalogCreator:
    def __init__(self, parts, area_partition=None, skycatalog_root=None,
                 catalog_dir='.', galaxy_truth=None,
                 star_truth=None,
                 config_path=None, catalog_name='skyCatalog',
                 output_type='parquet', mag_cut=None,
                 sed_subdir='galaxyTopHatSED', knots_mag_cut=27.0,
                 knots=True, logname='skyCatalogs.creator',
                 pkg_root=None, skip_done=False, no_main=False,
                 no_flux=False, flux_parallel=16, galaxy_nside=32,
                 galaxy_stride=1000000, provenance=None,
                 dc2=False, sn_object_type='sncosmo', galaxy_type='cosmodc2',
                 include_roman_flux=False, star_input_fmt='sqlite',
                 sso_truth=None, sso_sed=None, sso_partition='healpixel',
                 run_options=None):
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
        galaxy_truth    GCRCatalogs name or abs. path for galaxy truth
                        (e.g. cosmoDC2). Each known galaxy type
                        has a corresponding default value.
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
        no_main         Do not create main files
        no_flux         Do not create flux files
        flux_parallel   Number of processes to divide work of computing fluxes
        galaxy_nside    Healpix configuration value "nside" for galaxy output
        galaxy_stride   Max number of rows per galaxy row group
        provenance      Whether to write per-output-file git repo provenance
        dc2             Whether to adjust values to provide input comparable
                        to that for the DC2 run
        sn_object_type  Which object type to use for SNe.
        galaxy_type     Currently allowed values are cosmodc2 and diffsky
        include_roman_flux Calculate and write Roman flux values
        star_input_fmt  May be either 'sqlite' or 'parquet'
        sso_truth       Directory containing Sorcha output
        sso_sed         Path to sed file to be used for all SSOs
        sso_partition   Whether to partition by time or by healpixels

        Might want to add a way to specify template for output file name
        and template for input sedLookup file name.
        """

        _cosmo_cat = 'cosmodc2_v1.1.4_image_addon_knots'
        _diffsky_cat = 'roman_rubin_2023_v1.1.2_elais'
        _star_db = '/global/cfs/cdirs/lsst/groups/SSim/DC2/dc2_stellar_healpixel.db'
#        _star_parquet = '/global/cfs/cdirs/descssim/postDC2/UW_star_catalog'
        _star_parquet = '/sdf/data/rubin/shared/ops-rehearsal-3/imSim_catalogs/UW_stars'

        self._galaxy_stride = galaxy_stride

        # Temporary. Should add a separate star_stride argument or change name
        # e.g. galaxy_stride --> catalog_stride
        self._star_stride = galaxy_stride
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
        self._galaxy_type = galaxy_type
        self._galaxy_truth = galaxy_truth
        if galaxy_truth is None:
            if galaxy_type == 'cosmodc2':
                self._galaxy_truth = _cosmo_cat
            else:
                self._galaxy_truth = _diffsky_cat

        self._sn_object_type = sn_object_type

        self._star_truth = star_truth
        self._star_input_fmt = star_input_fmt
        if self._star_truth is None:
            if self._star_input_fmt == 'sqlite':
                self._star_truth = _star_db
            elif self._star_input_fmt == 'parquet':
                self._star_truth = _star_parquet

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
        self._no_main = no_main
        self._no_flux = no_flux
        self._flux_parallel = flux_parallel
        self._galaxy_nside = galaxy_nside
        self._provenance = provenance
        self._dc2 = dc2
        self._include_roman_flux = include_roman_flux
        self._obs_sed_factory = None
        self._sso_sed_factory = None               # do we need this?
        self._sso_creator = SsoCatalogCreator(self, sso_truth, sso_sed)
        self._sso_truth = self._sso_creator.sso_truth
        self._sso_sed = self._sso_creator.sso_sed
        self._sso_partition = sso_partition
        self._run_options = run_options

        # Do we need this?
        self._sed_bins = None

    def _make_tophat_columns(self, dat, names, cmp):
        '''
        Create columns sed_val_cmp, cmp_magnorm where cmp is one of "disk",
        "bulge", "knots"

        Parameters
        ----------
        dat          Data read from input galaxy catalog. Includes keys for
                     everything in names plus entry for redshiftHubble
        names        Names of SED columns for this component
        cmp          Component name

        Returns
        -------
        Add keys  sed_val_cmp, cmp_magnorm to input dat. Then return dat.
        '''
        sed_vals = (np.array([dat[k] for k in names]).T).tolist()
        dat['sed_val_' + cmp] = sed_vals
        dat[cmp + '_magnorm'] = [self._obs_sed_factory.magnorm(s, z) for (s, z)
                                 in zip(sed_vals, dat['redshiftHubble'])]
        for k in names:
            del dat[k]
        return dat

    def create(self, catalog_type):
        """
        Create catalog of specified type, using stored context.

        Parameters
        ----------
        catalog_type   string    Currently 'galaxy', 'pointsource'
                                 and 'sso' are the only values allowed
        Return
        ------
        None
        """
        if catalog_type == ('galaxy'):
            if not self._no_main:
                self.create_galaxy_catalog()
            if not self._no_flux:
                self.create_galaxy_flux_catalog()
        elif catalog_type == ('pointsource'):
            if not self._no_main:
                self.create_pointsource_catalog()
            if not self._no_flux:
                self.create_pointsource_flux_catalog()
        elif catalog_type == ('sso'):
            if not self._no_main:
                self._sso_creator.create_sso_catalog()
            if not self._no_flux:
                self._sso_creator.create_sso_flux_catalog()

        else:
            raise NotImplementedError(f'CatalogCreator.create: unsupported catalog type {catalog_type}')

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
        self._gal_cat = gal_cat

        # Save cosmology in case we need to write parameters out later
        self._cosmology = gal_cat.cosmology

        arrow_schema = make_galaxy_schema(self._logname,
                                          sed_subdir=self._sed_subdir,
                                          knots=self._knots,
                                          galaxy_type=self._galaxy_type)

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
        if dlen == 0:
            return
        last_row_ix = dlen - 1
        u_bnd = min(stride, dlen)
        l_bnd = 0
        rg_written = 0
        writer = None

        while u_bnd > l_bnd:
            out_dict = {k: dat[k][l_bnd: u_bnd] for k in dat if k not in to_rename}
            for k in to_rename:
                out_dict[to_rename[k]] = dat[k][l_bnd: u_bnd]
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

        if self._galaxy_type == 'cosmodc2':

            # to_fetch = all columns of interest in gal_cat
            non_sed = ['galaxy_id', 'ra', 'dec', 'redshift', 'redshiftHubble',
                       'peculiarVelocity', 'shear_1', 'shear_2',
                       'convergence', 'size_bulge_true',
                       'size_minor_bulge_true', 'sersic_bulge',
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
            all_q = gal_cat.list_all_quantities()
            sed_bins, sed_bulge_names, sed_disk_names = _get_tophat_info(all_q)
            self._sed_bins = sed_bins

            th_fact = TophatSedFactory(sed_bins,
                                       assemble_cosmology(self._cosmology))
            self._obs_sed_factory = th_fact

            to_fetch = non_sed + sed_bulge_names + sed_disk_names

        elif self._galaxy_type == 'diffsky':
            to_fetch = ['galaxy_id', 'ra', 'dec', 'redshift', 'redshiftHubble',
                        'peculiarVelocity', 'shear1', 'shear2',
                        'convergence', 'diskEllipticity1', 'diskEllipticity2',
                        'spheroidEllipticity1', 'spheroidEllipticity2',
                        'spheroidHalfLightRadiusArcsec',
                        'diskHalfLightRadiusArcsec', 'um_source_galaxy_obs_sm']

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

        # For cosmodc2 input some columns need to be renamed and there is
        # special handling for knots
        to_rename = dict()
        if self._galaxy_type == 'cosmodc2':
            to_rename = {'redshiftHubble': 'redshift_hubble',
                         'peculiarVelocity': 'peculiar_velocity'}
            if self._dc2:
                to_rename['ellipticity_1_disk_true_dc2'] = 'ellipticity_1_disk_true'
                to_rename['ellipticity_2_disk_true_dc2'] = 'ellipticity_2_disk_true'
                to_rename['ellipticity_1_bulge_true_dc2'] = 'ellipticity_1_bulge_true'
                to_rename['ellipticity_2_bulge_true_dc2'] = 'ellipticity_2_bulge_true'

            if self._sed_subdir:
                #  Generate full paths for disk and bulge SED files even though
                #  we don't actually write the files here
                df['bulge_sed_file_path'] =\
                    generate_sed_path(df['galaxy_id'], self._sed_subdir,
                                      'bulge')
                df['disk_sed_file_path'] =\
                    generate_sed_path(df['galaxy_id'], self._sed_subdir,
                                      'disk')

            if self._knots:
                # adjust disk sed; create knots sed
                sed_knot_names = [i.replace('disk', 'knots') for i in sed_disk_names]
                eps = np.finfo(np.float32).eps
                mag_mask = np.where(np.array(df['mag_i_lsst']) > self._knots_mag_cut, 0, 1)
                self._logger.debug(f'Count of mags <=  cut (so adjustment performed: {np.count_nonzero(mag_mask)}')

                for d_name, k_name in zip(sed_disk_names, sed_knot_names):
                    df[k_name] = mag_mask * np.clip(df['knots_flux_ratio'],
                                                    None, 1-eps) * df[d_name]
                    df[d_name] = np.where(np.array(df['mag_i_lsst']) > self._knots_mag_cut, 1,
                                          np.clip(1 - df['knots_flux_ratio'],
                                                  eps, None)) * df[d_name]

        if len(self._out_pixels) > 1:
            subpixel_masks = _generate_subpixel_masks(df['ra'], df['dec'],
                                                      self._out_pixels,
                                                      nside=self._galaxy_nside)
        else:
            subpixel_masks = {pixel: None}

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

                if self._galaxy_type == 'cosmodc2':
                    compressed = self._make_tophat_columns(compressed,
                                                           sed_disk_names,
                                                           'disk')
                    compressed = self._make_tophat_columns(compressed,
                                                           sed_bulge_names,
                                                           'bulge')
                    if self._knots:
                        compressed = self._make_tophat_columns(compressed,
                                                               sed_knot_names,
                                                               'knots')

                self._write_subpixel(dat=compressed, output_path=output_path,
                                     arrow_schema=arrow_schema,
                                     stride=stride, to_rename=to_rename)
            else:
                if self._galaxy_type == 'cosmodc2':
                    df = self._make_tophat_columns(df, sed_disk_names, 'disk')
                    df = self._make_tophat_columns(df, sed_bulge_names, 'bulge')
                    if self._knots:
                        df = self._make_tophat_columns(df, sed_knot_names,
                                                       'knots')
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

        from .skyCatalogs import open_catalog
        self._sed_gen = None

        self._gal_flux_schema =\
            make_galaxy_flux_schema(self._logname, self._galaxy_type,
                                    include_roman_flux=self._include_roman_flux)
        self._gal_flux_needed = [field.name for field in self._gal_flux_schema]

        if not config_file:
            config_file = self.write_config(path_only=True)
        if not self._cat:
            self._cat = open_catalog(config_file,
                                     skycatalog_root=self._skycatalog_root)

        # Might also open a main file and read its metadata to be sure
        # the value for stride is correct.
        # Especially need to do this if we want to support making flux
        # files in a separate job activation.
        self.object_type = 'galaxy'
        if self._galaxy_type == 'diffsky':
            self.object_type = 'diffsky_galaxy'
            from .diffsky_sedgen import DiffskySedGenerator
            # Default values are ok for all the diffsky-specific
            # parameters: include_nonLSST_flux, sed_parallel, auto_loop,
            # wave_ang_min, wave_ang_max, rel_err, n_per
            self._sed_gen = DiffskySedGenerator(logname=self._logname,
                                                galaxy_truth=self._galaxy_truth,
                                                output_dir=self._output_dir,
                                                skip_done=True,
                                                sky_cat=self._cat)

        self._flux_template = self._cat.raw_config['object_types'][self.object_type]['flux_file_template']

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
        object_list = self._cat.get_object_type_by_hp(pixel, self.object_type)
        if len(object_list) == 0:
            self._logger.warning(f'Cannot create flux file for pixel {pixel} because main file does not exist or is empty')
            return

        if self._galaxy_type == 'diffsky':
            # Generate SEDs if necessary
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
                out_dict = _do_galaxy_flux_chunk(None, _galaxy_collection,
                                                 _instrument_needed, lb, u)
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
                                   args=(conn_wrt, _galaxy_collection,
                                         _instrument_needed, lb, u))
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
        arrow_schema = make_star_schema()
        #  Need a way to indicate which object types to include; deal with that
        #  later.  For now, default is stars + sn
        for p in self._parts:
            self._logger.debug(f'Point sources. Starting on pixel {p}')
            self.create_pointsource_pixel(p, arrow_schema,
                                          star_cat=self._star_truth)
            self._logger.debug(f'Completed pixel {p}')

    def create_pointsource_pixel(self, pixel, arrow_schema, star_cat=None):
        if not star_cat:
            self._logger.info('No star input specified')
            return

        output_filename = f'pointsource_{pixel}.parquet'
        output_path = os.path.join(self._output_dir, output_filename)
        stride = self._star_stride

        if os.path.exists(output_path):
            if not self._skip_done:
                self._logger.info(f'Overwriting file {output_path}')
            else:
                self._logger.info(f'Skipping regeneration of {output_path}')
                return

        writer = pq.ParquetWriter(output_path, arrow_schema)

        # Get data for this pixel
        if self._star_input_fmt == 'sqlite':
            cols = ','.join(['format("%s",simobjid) as id', 'ra',
                             'decl as dec', 'magNorm as magnorm', 'mura',
                             'mudecl as mudec',
                             'radialVelocity as radial_velocity',
                             'parallax',
                             'sedFilename as sed_filepath', 'ebv'])
            q = f'select {cols} from stars where hpid={pixel} '
            with sqlite3.connect(star_cat) as conn:
                star_df = pd.read_sql_query(q, conn)
        elif self._star_input_fmt == 'parquet':
            star_df = _star_parquet_reader(self._star_truth, pixel,
                                           arrow_schema)
        nobj = len(star_df['id'])
        self._logger.debug(f'Found {nobj} stars')
        if nobj == 0:
            return
        star_df['sed_filepath'] = get_star_sed_path(star_df['sed_filepath'])
        star_df['object_type'] = np.full((nobj,), 'star')
        star_df['host_galaxy_id'] = np.zeros((nobj,), np.int64())

        star_df['MW_rv'] = np.full((nobj,), _MW_rv_constant, np.float32())

        # NOTE MW_av calculation for stars does not use SFD dust map
        star_df['MW_av'] = star_df['ebv'] * _MW_rv_constant

        star_df['variability_model'] = np.full((nobj,), '')
        star_df['salt2_params'] = np.full((nobj,), None)

        last_row_ix = nobj - 1
        u_bnd = min(stride, nobj)
        l_bnd = 0
        rg_written = 0

        while u_bnd > l_bnd:
            out_dict = {k: star_df[k][l_bnd: u_bnd] for k in star_df.columns}
            out_df = pd.DataFrame.from_dict(out_dict)

            out_table = pa.Table.from_pandas(out_df, schema=arrow_schema)
            self._logger.debug('Created arrow table from star dataframe')

            # write a row broup
            writer.write_table(out_table)
            rg_written += 1
            l_bnd = u_bnd
            u_bnd = min(l_bnd + stride, last_row_ix + 1)

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

        from .skyCatalogs import open_catalog

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
                out_dict = _do_star_flux_chunk(None, _star_collection,
                                               instrument_needed, lb, u)
            else:
                # Expect to be able to do about 1500/minute/process

                tm = max(int((n_per*60)/500), 5)  # Give ourselves a cushion
                self._logger.info(f'Using timeout value {tm} for {n_per} sources')
                p_list = []
                for i in range(n_parallel):
                    conn_rd, conn_wrt = Pipe(duplex=False)
                    readers.append(conn_rd)

                    # For debugging call directly
                    proc = Process(target=_do_star_flux_chunk,
                                   name=f'proc_{i}',
                                   args=(conn_wrt, _star_collection,
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
                                             schema=self._ps_flux_schema)

            if not writer:
                writer = pq.ParquetWriter(output_path, self._ps_flux_schema)
            writer.write_table(out_table)
            rg_written += 1

        writer.close()
        self._logger.debug(f'# row groups written to flux file: {rg_written}')
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
            return os.path.join(self._config_path,
                                self._catalog_name + '.yaml')

        config = create_config(self._catalog_name, self._logname)
        if self._global_partition is not None:
            config.add_key('area_partition', self._area_partition)
        # These will be covered in run_options so no need to add
        # separately
        # config.add_key('skycatalog_root', self._skycatalog_root)
        # config.add_key('catalog_dir', self._catalog_dir)

        #   Do we need this?
        config.add_key('active_skycatalog_root', self._skycatalog_root)

        if self._galaxy_type == 'cosmodc2':
            config.add_key('SED_models',
                           assemble_SED_models(self._sed_bins))
        config.add_key('MW_extinction_values', assemble_MW_extinction())
        config.add_key('Cosmology', assemble_cosmology(self._cosmology))
        config.add_key('object_types',
                       assemble_object_types(self._pkg_root,
                                             galaxy_nside=self._galaxy_nside))

        # Following will be in run_options
        # config.add_key('galaxy_magnitude_cut', self._mag_cut)
        # config.add_key('knots_magnitude_cut', self._knots_mag_cut)

        inputs = {'galaxy_truth': self._galaxy_truth}
        if self._star_truth:
            inputs['star_truth'] = self._star_truth
        if self._sso_truth:
            inputs['sso_truth'] = self._sso_truth
            inputs['sso_sed'] = self._sso_sed
        config.add_key('provenance',
                       assemble_provenance(self._pkg_root, inputs=inputs,
                                           run_options=self._run_options))

        self._written_config = config.write_config(self._config_path,
                                                   overwrite=overwrite)

    def write_provenance_file(self, datafile_path):
        '''
        Write git provenance to a yaml file with name derived from a
        just-written datafile name
        '''
        outpath = datafile_path.rsplit('.', 1)[0] + '_provenance.yaml'

        prov = assemble_provenance(self._pkg_root, inputs=None,
                                   run_options=self._run_options)
        write_yaml(prov, outpath)
