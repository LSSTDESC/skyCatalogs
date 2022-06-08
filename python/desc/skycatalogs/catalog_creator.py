import os
import re
import math
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from astropy.coordinates import SkyCoord
import sqlite3
from desc.skycatalogs.utils.common_utils import print_date
from desc.skycatalogs.utils.sed_utils import MagNorm, NORMWV_IX, get_star_sed_path

# from dm stack
from dustmaps.sfd import SFDQuery

"""
Code to create a sky catalog for a particular object type
"""

#####pixels = [9556, 9683, 9684, 9812, 9813, 9940]

__all__ = ['CatalogCreator']

_Av_adjustment = 2.742
_MW_rv_constant = 3.1

# This schema is not the same as the one taken from the data,
# probably because of the indexing in the schema derived from a pandas df.
def _make_galaxy_schema(sed_subdir=False, knots=True):
    fields = [pa.field('galaxy_id', pa.int64()),
              pa.field('ra', pa.float64() , True),
##                       metadata={"units" : "radians"}),
              pa.field('dec', pa.float64() , True),
##                       metadata={"units" : "radians"}),
              pa.field('redshift', pa.float64(), True),
              pa.field('redshift_hubble', pa.float64(), True),
              pa.field('peculiar_velocity', pa.float64(), True),
              pa.field('shear_1', pa.float64(), True),
              pa.field('shear_2', pa.float64(), True),
              pa.field('convergence', pa.float64(), True),
              pa.field('size_bulge_true', pa.float32(), True),
              pa.field('size_minor_bulge_true', pa.float32(), True),
              pa.field('sersic_bulge', pa.float32(), True),
              pa.field('size_disk_true', pa.float32(), True),
              pa.field('size_minor_disk_true', pa.float32(), True),
              pa.field('sersic_disk', pa.float32(), True),
              pa.field('position_angle_unlensed', pa.float64(), True),
              pa.field('sed_val_bulge',
                       pa.list_(pa.float64()), True),
              pa.field('sed_val_disk',
                       pa.list_(pa.float64()), True),
              pa.field('bulge_magnorm', pa.float64(), True),
              pa.field('disk_magnorm', pa.float64(), True),
              pa.field('MW_rv', pa.float32(), True),
              pa.field('MW_av', pa.float32(), True)]
    if knots:
        print("knots requested")
        fields.append(pa.field('sed_val_knots',
                               pa.list_(pa.float64()), True))
        ### For sizes API can alias to disk sizes
        ###  position angle, shears and convergence are all
        ###  galaxy-wide quantities.
        fields.append(pa.field('n_knots', pa.float32(), True))
        fields.append(pa.field('knots_magnorm', pa.float64(), True))

    if sed_subdir:
        fields.append(pa.field('bulge_sed_file_path', pa.string(), True))
        fields.append(pa.field('disk_sed_file_path', pa.string(), True))

    for f in fields:
        print(f.name)
    return pa.schema(fields)


def _make_star_schema():
    '''
    Minimal schema for non-variable stars.  For variables will need to add fields
    to express variability.   Will also likely have to make changes to accomodate SNe.
    If AGN also go in this file will need to include gamma1, gamma2, kappa.
    Could add field for galactic extinction model, but currently it's always 'CCM'
    so will put it in config.
    '''
    fields = [pa.field('object_type', pa.string(), False),
              pa.field('id', pa.int64(), False),
              pa.field('ra', pa.float64(), False),
              pa.field('dec', pa.float64(), False),
              pa.field('host_galaxy_id', pa.int64(), True),
              pa.field('magnorm', pa.float64(), True),
              pa.field('sed_filepath', pa.string(), True),
              pa.field('MW_rv', pa.float32(), True),
              pa.field('MW_av', pa.float32(), True)]
    return pa.schema(fields)

def _make_MW_extinction(ra, dec):
    '''
    Given arrays of ra & dec, create a MW Av column corresponding to V-band
    correction.
    See "Plotting Dust Maps" example in
    https://dustmaps.readthedocs.io/en/latest/examples.html

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
                 catalog_dir='.', galaxy_truth=None, write_config=False,
                 config_path=None,
                 output_type='parquet', verbose=False,
                 mag_cut=None,  sed_subdir='galaxyTopHatSED', knots_mag_cut=27.0,
                 knots=True):
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
        write_config    If True, write config file. Default is False
        config_path     Where to write config file. Default is data
                        directory. Ignored if write_config
                        is False.
        output_type     A format.  For now only parquet is supported
        mag_cut         If not None, exclude galaxies with mag_r > mag_cut
        sed_subdir      In instcat entry, prepend this value to SED filename
        knots_mag_cut   No knots for galaxies with i_mag > cut
        knots           If True include knots

        Might want to add a way to specify template for output file name
        and template for input sedLookup file name.
        """

        _cosmo_cat = 'cosmodc2_v1.1.4_image_addon_knots'

        if area_partition['type'] != 'healpix':
            raise NotImplementedError(f'CatalogCreator: Unknown partition type {area_partition["type"]} ')

        if output_type != 'parquet':
            raise NotImplementedError(f'CatalogCreator: Output type {output_type} not supported')

        self._galaxy_truth = galaxy_truth
        if galaxy_truth is None:
            self._galaxy_truth = _cosmo_cat

        self._parts = parts
        self._area_partition = area_partition
        if skycatalog_root:
            self._skycatalog_root = skycatalog_root
        else:
            self._skycatalog_root = './'
        if catalog_dir:
            self._catalog_dir = catalog_dir
        else:
            catalog_dir = '.'

        self._output_dir = os.path.join(self._skycatalog_root,
                                        self._catalog_dir)

        self._output_type = output_type
        self._verbose = verbose
        self._mag_cut = mag_cut
        self._sed_subdir = sed_subdir
        self._knots_mag_cut = knots_mag_cut
        self._knots = knots

    def create(self, catalog_type):
        """
        Create catalog of specified type, using stored context.

        Parameters
        ----------
        catalog_type   string    Currently 'galaxy' and 'pointsource' are
                                 the only values allowed
        """
        if catalog_type == ('galaxy'):
            return self.create_galaxy_catalog()
        elif catalog_type == ('pointsource'):
            return self.create_pointsource_catalog()
        else:
            raise NotImplemented(f'CatalogCreator.create: unsupported catalog type {catalog_type}')

    def create_galaxy_catalog(self):
        """
        Returns
        -------
        Dict describing catalog produced
        """
        import GCRCatalogs

        gal_cat = GCRCatalogs.load_catalog(self._galaxy_truth)

        self._mag_norm_f = MagNorm(gal_cat.cosmology)

        arrow_schema = _make_galaxy_schema(self._sed_subdir, self._knots)

        for p in self._parts:
            print("Starting on pixel ", p)
            print_date()
            self.create_galaxy_pixel(p, gal_cat, arrow_schema)
            print("completed pixel ", p)
            print_date()

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
        stride = 1000000

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

        #Sort sed columns by start value, descending
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

        # Save the 'start' and 'width' values; they'll be needed for our output
        # config.  Though we only need to do this once, not once per pixel
        tophat_parms = []
        for s in sed_bulge_names:
            m = re.match(tophat_bulge_re, s)
            if m:
                tophat_parms.append((m['start'], m['width']))

        if not self._mag_cut:
            df = gal_cat.get_quantities(to_fetch, native_filters=hp_filter)
        else:
            df = gal_cat.get_quantities(to_fetch + [r_mag_name],
                                        native_filters=hp_filter,
                                        filters=mag_cut_filter)

        if self._sed_subdir:
            #  Generate full paths for disk and bulge SED files, even though
            #  we don't actually write the files here
            # Name of file can be <galaxy_id>_<cmp>.txt
            # prepend subdirectory supplied as command-line arg.
            bulge_path = _generate_sed_path(df['galaxy_id'], self._sed_subdir,
                                            'bulge')
            disk_path =  _generate_sed_path(df['galaxy_id'], self._sed_subdir,
                                            'disk')

        if self._knots:
            # adjust disk sed; create knots sed
            eps = np.finfo(np.float32).eps
            mag_mask = np.where(np.array(df['mag_i_lsst']) > self._knots_mag_cut,0, 1)
            if self._verbose:
                print("Count of mags <=  cut (so adjustment performed: ",
                      np.count_nonzero(mag_mask))

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
        if self._verbose: print("Made extinction")

        #  Write row groups of size stride (or less) until input is exhausted
        last_row = len(df['galaxy_id']) - 1
        u_bnd = min(stride, last_row)
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

            if self._verbose:
                for kd,i in out_dict.items():
                    print(f'Key={kd}, type={type(i)}, len={len(i)}')
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
            u_bnd = min(l_bnd + stride, last_row)

        writer.close()
        if self._verbose: print("# row groups written: ", rg_written)

    def create_pointsource_catalog(self, star_truth=None, sne_truth=None):

        """
        Parameters
        ----------
        star_truth      Where to find star parameters. If None, omit stars
        sne_truth       Where to find SNe parameters.  If None, omit SNe

        Might want to add a way to specify template for output file name

        Returns
        -------
        Dict describing catalog produced
        """

        # For now fixed location for star, SNe parameter files.
        _star_db = '/global/cfs/cdirs/lsst/groups/SSim/DC2/dc2_stellar_healpixel.db'
        _sn_db = '/global/cfs/cdirs/lsst/groups/SSim/DC2/cosmoDC2_v1.1.4/sne_cosmoDC2_v1.1.4_MS_DDF_healpix.db'

        arrow_schema = _make_star_schema()
        #  Need a way to indicate which object types to include; deal with that
        #  later.  For now, default is stars only.  Use default star parameter file.
        for p in self._parts:
            if self._verbose:
                print("Point sources. Starting on pixel ", p)
                print_date()
            self.create_pointsource_pixel(p, arrow_schema, star_cat=_star_db)
            if self._verbose:
                print("completed pixel ", p)
                print_date()

    def create_pointsource_pixel(self, pixel, arrow_schema, star_cat=None,
                             sn_cat=None):

        if not star_cat and not sn_cat:
            print("no point source inputs specified")
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
            print(f"Found {nobj} stars")
            star_df['object_type'] = np.full((nobj,), 'star')
            star_df['host_galaxy_id'] = np.zeros((nobj,), np.int64())

            star_df['MW_rv'] = np.full((nobj,), _MW_rv_constant, np.float32())

            star_df['MW_av'] = _make_MW_extinction(np.array(star_df['ra']),
                                                   np.array(star_df['dec']))
            out_table = pa.Table.from_pandas(star_df, schema=arrow_schema)
            if self._verbose: print("created arrow table from dataframe")

            writer = pq.ParquetWriter(os.path.join(
                self._output_dir, output_template.format(pixel)), arrow_schema)
            writer.write_table(out_table)

            writer.close()

        if sn_cat:
            raise NotImplementedError('SNe not yet supported. Have a nice day.')

        return
