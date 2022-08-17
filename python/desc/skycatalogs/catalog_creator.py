import os
import re
import math
import logging
import yaml
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from astropy.coordinates import SkyCoord
import sqlite3
from desc.skycatalogs.utils.common_utils import print_date
from desc.skycatalogs.utils.sed_utils import MagNorm, NORMWV_IX, get_star_sed_path
from desc.skycatalogs.utils.config_utils import create_config, assemble_SED_models
from desc.skycatalogs.utils.config_utils import assemble_MW_extinction, assemble_cosmology, assemble_object_types, assemble_provenance
from desc.skycatalogs.utils.parquet_schema_utils import make_galaxy_schema, make_galaxy_flux_schema, make_star_schema, make_star_flux_schema
from desc.skycatalogs.objects.base_object import LSST_BANDS, BaseObject

# from dm stack
from dustmaps.sfd import SFDQuery

"""
Code to create a sky catalog for particular object types
"""

#####pixels = [9556, 9683, 9684, 9812, 9813, 9940]

__all__ = ['CatalogCreator']

_MW_rv_constant = 3.1


_Av_adjustment = 2.742
def _make_MW_extinction(ra, dec):
    '''
    Given arrays of ra & dec, create a MW Av column corresponding to V-band
    correction.
    See "Plotting Dust Maps" example in
    https://dustmaps.readthedocs.io/en/latest/examples.html

    The coefficient _Av_adjustment comes Table 6 in Schlafly & Finkbeiner (2011)
    See http://iopscience.iop.org/0004-637X/737/2/103/article#apj398709t6

    Parameters
    ----------
    ra, dec - arrays specifying positions where Av is to be computed
    Return:
    Array of Av values
    '''

    sfd = SFDQuery()
    ebv_raw = np.array(sfd.query_equ(ra, dec))

    return _Av_adjustment * ebv_raw

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

class CatalogCreator:
    def __init__(self, parts, area_partition, skycatalog_root=None,
                 catalog_dir='.', galaxy_truth=None,
                 config_path=None, catalog_name='skyCatalog',
                 output_type='parquet', mag_cut=None,
                 sed_subdir='galaxyTopHatSED', knots_mag_cut=27.0,
                 knots=True, logname='skyCatalogs.creator',
                 pkg_root=None):
        """
        Store context for catalog creation

        Parameters
        ----------
        parts           Segments for which catalog is to be generated. If
                        partition type is HEALpix, parts will be a collection
                        of HEALpix pixels
        area_partition  Dict characterizing partition; e.g. HEALpix,
                        nside=<something>
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
        pkg_root        defaults to three levels up from __file__

        Might want to add a way to specify template for output file name
        and template for input sedLookup file name.
        """

        _cosmo_cat = 'cosmodc2_v1.1.4_image_addon_knots'
        self._galaxy_stride = 1000000
        if pkg_root:
            self._pkg_root = pkg_root
        else:
            self._pkg_root = os.path.join(os.path.dirname(__file__),
                                          '../../..')

        if area_partition['type'] != 'healpix':
            raise NotImplementedError(f'CatalogCreator: Unknown partition type {area_partition["type"]} ')

        if output_type != 'parquet':
            raise NotImplementedError(f'CatalogCreator: Output type {output_type} not supported')

        self._galaxy_truth = galaxy_truth
        if galaxy_truth is None:
            self._galaxy_truth = _cosmo_cat

        self._sn_truth = None
        self._star_truth = None
        self._cat = None

        self._parts = parts
        self._area_partition = area_partition
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
            self.create_galaxy_catalog()
            self.create_galaxy_flux_catalog()
        elif catalog_type == ('pointsource'):
            self.create_pointsource_catalog()
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

        self._mag_norm_f = MagNorm(self._cosmology)

        arrow_schema = make_galaxy_schema(self._logname, self._sed_subdir,
                                           self._knots)
        self._gal_flux_schema = make_galaxy_flux_schema(self._logname)


        for p in self._parts:
            self._logger.info(f'Starting on pixel {p}')
            self.create_galaxy_pixel(p, gal_cat, arrow_schema)
            self._logger.info(f'Completed pixel {p}')

        # Now make config.   We need it for computing LSST fluxes for
        # the second part of the galaxy catalog
        self.write_config()

    def create_galaxy_pixel(self, pixel, gal_cat, arrow_schema):
        """
        Parameters
        ----------
        pixel           Pixel for which catalog is to be generated.
        gal_cat         GCRCatalogs-loaded galaxy truth (e.g. cosmoDC2)
        arrow_schema    schema to use for output file
        """

        # Output filename template
        output_template = 'galaxy_{}.parquet'

        # Used to finding tophat parameters from cosmoDC2 column names
        tophat_bulge_re = r'sed_(?P<start>\d+)_(?P<width>\d+)_bulge'
        tophat_disk_re = r'sed_(?P<start>\d+)_(?P<width>\d+)_disk'

        # Number of rows to include in a row group
        #stride = 1000000
        stride = self._galaxy_stride

        hp_filter = [f'healpix_pixel=={pixel}']
        if self._mag_cut:
            r_mag_name = 'mag_r_lsst'
            mag_cut_filter = [f'{r_mag_name} < {self._mag_cut}']

        # to_fetch = all columns of interest in gal_cat
        non_sed = ['galaxy_id', 'ra', 'dec', 'redshift', 'redshiftHubble',
                   'peculiarVelocity', 'shear_1', 'shear_2',
                   'convergence', 'position_angle_true',
                   'size_bulge_true', 'size_minor_bulge_true', 'sersic_bulge',
                   'size_disk_true', 'size_minor_disk_true', 'sersic_disk']
        if self._knots:
                   non_sed += ['knots_flux_ratio', 'n_knots', 'mag_i_lsst']
        # Find all the tophat sed numbers
        q = gal_cat.list_all_quantities()
        sed_bulge_names = [i for i in q if (i.startswith('sed') and
                                            i.endswith('bulge'))]
        sed_disk_names = [i for i in q if (i.startswith('sed') and
                                           i.endswith('disk'))]
        # Save all the start, width values
        self._sed_bins = [[int(re.match(tophat_bulge_re, s)['start']), int(re.match(tophat_bulge_re, s)['width'])] for s in sed_bulge_names]

        # Sort by value for start
        def _bin_start_key(start_width):
            return start_width[0]
        self._sed_bins.sort(key=_bin_start_key)

        def _sed_bulge_key(s):
            return int(re.match(tophat_bulge_re, s)['start'])
        def _sed_disk_key(s):
            return int(re.match(tophat_disk_re, s)['start'])

        # Sort into increaing order by start wavelength
        sed_bulge_names.sort(key=_sed_bulge_key)
        sed_disk_names.sort(key=_sed_disk_key)
        if self._knots:
            sed_knot_names = [i.replace('disk', 'knots') for i in sed_disk_names]

        #Fetch the data
        to_fetch = non_sed + sed_bulge_names + sed_disk_names

        if not self._mag_cut:
            df = gal_cat.get_quantities(to_fetch, native_filters=hp_filter)
        else:
            df = gal_cat.get_quantities(to_fetch + [r_mag_name],
                                        native_filters=hp_filter,
                                        filters=mag_cut_filter)

        if self._sed_subdir:
            #  Generate full paths for disk and bulge SED files, even though
            #  we don't actually write the files here
            bulge_path = _generate_sed_path(df['galaxy_id'], self._sed_subdir,
                                            'bulge')
            disk_path =  _generate_sed_path(df['galaxy_id'], self._sed_subdir,
                                            'disk')

        if self._knots:
            # adjust disk sed; create knots sed
            eps = np.finfo(np.float32).eps
            mag_mask = np.where(np.array(df['mag_i_lsst']) > self._knots_mag_cut,0, 1)
            self._logger.debug(f'Count of mags <=  cut (so adjustment performed: {np.count_nonzero(mag_mask)}')

            for d_name, k_name in zip(sed_disk_names, sed_knot_names):
                df[k_name] = mag_mask * np.clip(df['knots_flux_ratio'], None, 1-eps) * df[d_name]
                df[d_name] = np.where(np.array(df['mag_i_lsst']) > self._knots_mag_cut, 1,
                                      np.clip(1 - df['knots_flux_ratio'], eps, None)) * df[d_name]
        # Re-form sed columns into two arrays
        bulge_seds = (np.array([df[sbn] for sbn in sed_bulge_names]).T).tolist()
        disk_seds = (np.array([df[sdn] for sdn in sed_disk_names]).T).tolist()
        if self._knots:
            knots_seds = (np.array([df[kdn] for kdn in sed_knot_names]).T).tolist()

        #  Compute mag_norm from TH sed and redshift
        bulge_magnorm = [self._mag_norm_f(s[NORMWV_IX], r) for (s, r) in zip(bulge_seds, df['redshiftHubble']) ]
        disk_magnorm = [self._mag_norm_f(s[NORMWV_IX], r) for (s, r) in zip(disk_seds, df['redshiftHubble']) ]
        if self._knots:
            knots_magnorm = [self._mag_norm_f(s[NORMWV_IX], r) for (s, r) in zip(knots_seds, df['redshiftHubble']) ]

        MW_rv = np.full_like(df['sersic_bulge'], _MW_rv_constant)
        MW_av = _make_MW_extinction(df['ra'], df['dec'])
        self._logger.debug('Made extinction')

        #  Write row groups of size stride (or less) until input is exhausted
        last_row_ix = len(df['galaxy_id']) - 1
        u_bnd = min(stride, last_row_ix + 1)
        l_bnd = 0
        rg_written = 0
        writer = None

        # Some columns need to be renamed
        to_modify = ['position_angle_true', 'redshiftHubble', 'peculiarVelocity']
        while u_bnd > l_bnd:
            out_dict = {k : df[k][l_bnd : u_bnd] for k in non_sed if k not in to_modify}
            out_dict['redshift_hubble'] = df['redshiftHubble'][l_bnd : u_bnd]
            out_dict['peculiar_velocity'] = df['peculiarVelocity'][l_bnd : u_bnd]
            out_dict['position_angle_unlensed'] = df['position_angle_true'][l_bnd : u_bnd]
            out_dict['sed_val_bulge'] = bulge_seds[l_bnd : u_bnd]
            out_dict['sed_val_disk'] = disk_seds[l_bnd : u_bnd]
            out_dict['bulge_magnorm'] = bulge_magnorm[l_bnd : u_bnd]
            out_dict['disk_magnorm'] = disk_magnorm[l_bnd : u_bnd]
            out_dict['MW_rv'] = MW_rv[l_bnd : u_bnd]
            out_dict['MW_av'] = MW_av[l_bnd : u_bnd]

            if self._knots:
                out_dict['sed_val_knots'] = knots_seds[l_bnd : u_bnd]
                out_dict['n_knots'] = df['n_knots'][l_bnd : u_bnd]
                out_dict['knots_magnorm'] = knots_magnorm[l_bnd : u_bnd]

            if self._sed_subdir:
                out_dict['bulge_sed_file_path'] = bulge_path[l_bnd : u_bnd]
                out_dict['disk_sed_file_path'] = disk_path[l_bnd : u_bnd]

            if self._logger.getEffectiveLevel() == logging.DEBUG:
                for kd,i in out_dict.items():
                    self._logger.debug(f'Key={kd}, type={type(i)}, len={len(i)}')
            print_col = False

            out_df = pd.DataFrame.from_dict(out_dict)
            out_table = pa.Table.from_pandas(out_df, schema=arrow_schema)
            if not writer:
                writer = pq.ParquetWriter(os.path.join(
                    self._output_dir, output_template.format(pixel)),
                                          arrow_schema)

            writer.write_table(out_table)
            rg_written += 1
            l_bnd = u_bnd
            u_bnd = min(l_bnd + stride, last_row_ix + 1)

        writer.close()
        self._logger.debug(f'# row groups written: {rg_written}')

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

        from desc.skycatalogs import open_catalog, SkyCatalog

        if not config_file:
            config_file = self._written_config
        if not self._cat:
            self._cat = open_catalog(config_file)

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

        # For main catalog use self._cat
        # For schema use self._gal_flux_schema
        # output_template should be derived from value for flux_file_template
        #  in main catalog config.  Cheat for now
        output_template = 'galaxy_flux_{}.parquet'
        stride = self._galaxy_stride

        object_list = self._cat.get_objects_by_hp(pixel,
                                                  obj_type_set={'galaxy'})
        last_row_ix = len(object_list) - 1
        writer = None
        l_bnd = 0
        u_bnd = min(last_row_ix + 1, stride)
        ##o_list = object_list[l_bnd : u_bnd]
        rg_written = 0
        while u_bnd > l_bnd:
            o_list = object_list[l_bnd : u_bnd]
            self._logger.debug(f'Handling range {l_bnd} up to {u_bnd}')
            out_dict = {}
            out_dict['galaxy_id'] = [o.get_native_attribute('galaxy_id') for o in o_list]
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
                                             schema=self._gal_flux_schema)
            if not writer:
                writer = pq.ParquetWriter(os.path.join(self._output_dir,
                                                       output_template.format(pixel)),
                                                       self._gal_flux_schema)
            writer.write_table(out_table)
            l_bnd = u_bnd
            u_bnd = min(l_bnd + stride, last_row_ix + 1)

            rg_written +=1

        writer.close()
        self._logger.debug(f'# row groups written to flux file: {rg_written}')


    def create_pointsource_catalog(self, star_truth=None, sn_truth=None):

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

        # For now fixed location for star, SNe parameter files.
        _star_db = '/global/cfs/cdirs/lsst/groups/SSim/DC2/dc2_stellar_healpixel.db'
        _sn_db = '/global/cfs/cdirs/lsst/groups/SSim/DC2/cosmoDC2_v1.1.4/sne_cosmoDC2_v1.1.4_MS_DDF_healpix.db'
        self._star_truth = _star_db
        self._sn_truth = _sn_db

        arrow_schema = make_star_schema()
        self._ps_flux_schema = make_star_flux_schema(self._logname)
        #  Need a way to indicate which object types to include; deal with that
        #  later.  For now, default is stars only.  Use default star parameter file.
        for p in self._parts:
            self._logger.debug(f'Point sources. Starting on pixel {p}')
            self.create_pointsource_pixel(p, arrow_schema, star_cat=_star_db)
            self._logger.debug(f'Completed pixel {p}')

        # If we already did galaxies don't need to write config.
        # And if we didn't, we don't have all the information needed
        # for config
        #if not self._written_config:
        #    self.write_config()

    def create_pointsource_pixel(self, pixel, arrow_schema, star_cat=None,
                             sn_cat=None):

        if not star_cat and not sn_cat:
            self._logger.info('No point source inputs specified')
            return

        output_template = 'pointsource_{}.parquet'

        if star_cat:
            # Get data for this pixel
            cols = ','.join(['simobjid as id', 'ra', 'decl as dec',
                             'magNorm as magnorm', 'sedFilename as sed_filepath'])
            q = f'select {cols} from stars where hpid={pixel} '
            with sqlite3.connect(star_cat) as conn:
                star_df = pd.read_sql_query(q, conn)

            star_df['sed_filepath'] = get_star_sed_path(star_df['sed_filepath'])

            nobj = len(star_df['id'])
            self._logger.debug(f'Found {nobj} stars')
            star_df['object_type'] = np.full((nobj,), 'star')
            star_df['host_galaxy_id'] = np.zeros((nobj,), np.int64())

            star_df['MW_rv'] = np.full((nobj,), _MW_rv_constant, np.float32())

            star_df['MW_av'] = _make_MW_extinction(np.array(star_df['ra']),
                                                   np.array(star_df['dec']))
            out_table = pa.Table.from_pandas(star_df, schema=arrow_schema)
            self._logger.debug('Created arrow table from dataframe')

            writer = pq.ParquetWriter(os.path.join(
                self._output_dir, output_template.format(pixel)), arrow_schema)
            writer.write_table(out_table)

            writer.close()

        if sn_cat:
            raise NotImplementedError('SNe not yet supported. Have a nice day.')

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

        from desc.skycatalogs import open_catalog, SkyCatalog

        if not config_file:
            config_file = self.write_config(path_only=True)

        # Always open catalog. If it was opened for galaxies earlier
        # it won't know about star files.
        self._cat = open_catalog(config_file)

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
        Pixel         int

        Return
        ------
        None
        '''

        # For main catalog use self._cat
        # For schema use self._ps_flux_schema
        # output_template should be derived from value for flux_file_template
        #  in main catalog config.  Cheat for now
        output_template = 'pointsource_flux_{}.parquet'

        object_list = self._cat.get_objects_by_hp(pixel,
                                                  obj_type_set={'star'})
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
            writer = pq.ParquetWriter(os.path.join(self._output_dir,
                                                   output_template.format(pixel)),
                                      self._ps_flux_schema)
        writer.write_table(out_table)

        writer.close()
        ##self._logger.debug(f'# row groups written to flux file: {rg_written}')

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
        config.add_key('area_partition', self._area_partition)
        config.add_key('skycatalog_root', self._skycatalog_root)
        config.add_key('catalog_dir' , self._catalog_dir)

        config.add_key('SED_models',
                       assemble_SED_models(self._sed_bins))
        config.add_key('MW_extinction_values', assemble_MW_extinction())
        config.add_key('Cosmology', assemble_cosmology(self._cosmology))
        config.add_key('object_types', assemble_object_types(self._pkg_root))

        inputs = {'galaxy_truth' : self._galaxy_truth}
        if self._sn_truth:
            inputs['sn_truth'] = self._sn_truth
        if self._star_truth:
            inputs['star_truth'] = self._star_truth
        config.add_key('provenance', assemble_provenance(self._pkg_root,
                                                         inputs=inputs))

        self._written_config = config.write_config(self._config_path,
                                                   overwrite=overwrite)
