import os
import re
import logging
import numpy as np
import numpy.ma as ma
import healpy
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import sqlite3
from .utils.sed_tools import TophatSedFactory, get_star_sed_path
from .utils.config_utils import assemble_cosmology
from .utils.config_utils import assemble_provenance
from .utils.config_utils import assemble_file_metadata
from .utils.config_utils import ConfigWriter
from .utils.star_parquet_input import _star_parquet_reader
from .utils.parquet_schema_utils import make_galaxy_schema
from .utils.parquet_schema_utils import make_star_schema
from .utils.creator_utils import make_MW_extinction_av, make_MW_extinction_rv
from .objects.star_object import StarConfigFragment
from .objects.galaxy_object import GalaxyConfigFragment
from .objects.diffsky_object import DiffskyConfigFragment
from .sso_catalog_creator import SsoMainCatalogCreator

"""
Code to create a sky catalog for particular object types
"""

__all__ = ['MainCatalogCreator']

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


class MainCatalogCreator:
    def __init__(self, object_type, parts, skycatalog_root=None,
                 catalog_dir='.', truth=None,
                 config_path=None, catalog_name='skyCatalog',
                 mag_cut=None,
                 knots_mag_cut=27.0,
                 knots=True, logname='skyCatalogs.creator',
                 pkg_root=None, skip_done=False,
                 nside=32, stride=1000000, dc2=False,
                 star_input_fmt='sqlite', sso_sed=None,
                 run_options=None):
        """
        Store context for catalog creation

        Parameters
        ----------
        object_type     One of {"star", "cosmodc_galaxy", "diffsky_galaxy",
                        "sso"}
        parts           Segments for which catalog is to be generated. If
                        partition type is HEALpix, parts will be a collection
                        of HEALpix pixels
        skycatalog_root Typically absolute directory containing one or
                        more subdirectories for sky catalogs. Defaults
                        to current directory
        catalog_dir     Directory relative to skycatalog_root where catalog
                        will be written.  Defaults to '.'
        truth           GCRCatalogs name or abs. path for input truth catalog
                        Default depends on object type
        config_path     Where to write config file. Default is data
                        directory.
        catalog_name    If a config file is written this value is saved
                        there and is also part of the filename of the
                        config file.
        mag_cut         If not None, exclude galaxies with mag_r > mag_cut
        knots_mag_cut   No knots for galaxies with i_mag > cut
        knots           If True include knots
        logname         logname for Python logger
        pkg_root        defaults to one level up from __file__
        skip_done       If True, skip over files which already exist. Otherwise
                        (by default) overwrite with new version.
                        Output info message in either case if file exists.
        nside           Healpix configuration value "nside" for output
        stride          Max number of rows per row group
        dc2             Whether to adjust values to provide input comparable
                        to that for the DC2 run
        star_input_fmt  May be either 'sqlite' or 'parquet'
        run_options     The options the outer script (create_sc.py) was
                        called with

        Might want to add a way to specify template for output file name
        and template for input sedLookup file name.
        """

        self._object_type = object_type
        self._stride = stride
        if pkg_root:
            self._pkg_root = pkg_root
        else:
            self._pkg_root = os.path.join(os.path.dirname(__file__), '..')

        self._truth = truth
        self._star_input_fmt = star_input_fmt

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

        self._mag_cut = mag_cut
        self._knots_mag_cut = knots_mag_cut
        self._knots = knots
        self._logname = logname
        self._logger = logging.getLogger(logname)
        self._skip_done = skip_done
        self._nside = nside
        self._dc2 = dc2
        self._obs_sed_factory = None
        if object_type == 'sso':
            self._sso_creator = SsoMainCatalogCreator(self, self._truth)
        self._run_options = run_options
        self._tophat_sed_bins = None

        self._config_writer = ConfigWriter(self._skycatalog_root,
                                           self._catalog_dir,
                                           self._catalog_name,
                                           not self._skip_done,
                                           self._logname)

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

    def create(self):
        """
        Create catalog of specified type, using stored context.


        Return
        ------
        None
        """
        object_type = self._object_type
        if object_type in {'cosmodc2_galaxy', 'diffsky_galaxy'}:
            self.create_galaxy_catalog()
        elif object_type == ('star'):
            self.create_pointsource_catalog()
        elif object_type == ('sso'):
            self._sso_creator.create_sso_catalog()
        else:
            raise NotImplementedError(
                f'MainCatalogCreator.create: unsupported object type {object_type}')

    def create_galaxy_catalog(self):
        """
        Create the 'main' galaxy catalog, including everything except
        fluxes

        Returns
        -------
        None

        """
        _cosmo_cat = 'cosmodc2_v1.1.4_image_addon_knots'
        _diffsky_cat = 'roman_rubin_2023_v1.1.2_elais'

        import GCRCatalogs

        if self._object_type == 'cosmodc2_galaxy':
            self._galaxy_type = 'cosmodc2'
            if self._truth is None:
                self._truth = _cosmo_cat
        else:    # only other possibility is diffsky
            self._galaxy_type = 'diffsky'
            if self._truth is None:
                self._truth = _diffsky_cat

        gal_cat = GCRCatalogs.load_catalog(self._truth)
        self._gal_cat = gal_cat

        # Save cosmology in case we need to write parameters out later
        self._cosmology = gal_cat.cosmology

        inputs = {'galaxy_truth': self._truth}
        file_metadata = assemble_file_metadata(self._pkg_root,
                                               inputs=inputs,
                                               run_options=self._run_options)

        arrow_schema = make_galaxy_schema(self._logname,
                                          knots=self._knots,
                                          galaxy_type=self._galaxy_type,
                                          metadata_input=file_metadata)

        for p in self._parts:
            self._logger.info(f'Starting on pixel {p}')
            self.create_galaxy_pixel(p, gal_cat, arrow_schema)
            self._logger.info(f'Completed pixel {p}')

        # Now make config.   We need it for computing LSST fluxes for
        # the second part of the galaxy catalog
        prov = assemble_provenance(self._pkg_root,
                                   inputs={'galaxy_truth': self._truth},
                                   run_options=self._run_options)
        cosmo = assemble_cosmology(self._cosmology)
        if self._galaxy_type == 'diffsky':
            fragment = DiffskyConfigFragment(prov, cosmo)
            self._config_writer.write_configs(fragment)
        else:
            fragment = GalaxyConfigFragment(prov, cosmo, self._tophat_sed_bins)
            self._config_writer.write_configs(fragment)

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
        if self._nside > 32:
            out_pixels = _find_subpixels(pixel, self._nside)
            self._logger.debug(f'For nside={self._nside} subpixels are')
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
        stride = self._stride

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
            self._tophat_sed_bins = sed_bins

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
                                                      nside=self._nside)
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
        _star_db = '/global/cfs/cdirs/lsst/groups/SSim/DC2/dc2_stellar_healpixel.db'
#        _star_parquet = '/global/cfs/cdirs/descssim/postDC2/UW_star_catalog'
        _star_parquet = '/sdf/data/rubin/shared/ops-rehearsal-3/imSim_catalogs/UW_stars'

        if self._truth is None:
            if self._star_input_fmt == 'sqlite':
                self._truth = _star_db
            else:              # must be parquet
                self._truth = _star_parquet

        inputs = {'star_truth': self._truth}
        file_metadata = assemble_file_metadata(self._pkg_root,
                                               inputs=inputs,
                                               run_options=self._run_options)

        arrow_schema = make_star_schema(metadata_input=file_metadata)

        for p in self._parts:
            self._logger.debug(f'Point sources. Starting on pixel {p}')
            self.create_pointsource_pixel(p, arrow_schema,
                                          star_cat=self._truth)
            self._logger.debug(f'Completed pixel {p}')

        prov = assemble_provenance(self._pkg_root,
                                   inputs={'star_truth': self._truth},
                                   run_options=self._run_options)
        fragment = StarConfigFragment(prov)
        self._config_writer.write_configs(fragment)

    def create_pointsource_pixel(self, pixel, arrow_schema, star_cat=None):
        if not star_cat:
            self._logger.info('No star input specified')
            return

        output_filename = f'pointsource_{pixel}.parquet'
        output_path = os.path.join(self._output_dir, output_filename)
        stride = self._stride

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
            star_df = _star_parquet_reader(self._truth, pixel,
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
        return
