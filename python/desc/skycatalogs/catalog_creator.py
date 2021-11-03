import os
import re
import math
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from astropy.coordinates import SkyCoord
import h5py
import sqlite3
import GCRCatalogs
import pyccl as ccl
from desc.skycatalogs.utils.common_utils import print_date
from desc.skycatalogs.utils.sed_utils import MagNorm, get_random_sed, NORMWV_IX

# from dm stack
from dustmaps.sfd import SFDQuery

"""
Code to create a sky catalog for a particular object type
"""

#####pixels = [9556, 9683, 9684, 9812, 9813, 9940]

__all__ = ['CatalogCreator']

'''
Dict of MW av column names and multipliers needed to create from ebv, MW_rv
Multipliers come from
https://iopscience.iop.org/article/10.1088/0004-637X/737/2/103#apj398709t6
appendix, table 6
'''
_MW_extinction_bands = {'MW_av_lsst_u' : 4.145, 'MW_av_lsst_g' : 3.237,
                       'MW_av_lsst_r' : 2.273, 'MW_av_lsst_i' : 1.684,
                       'MW_av_lsst_z' : 1.323, 'MW_av_lsst_y' : 1.088}

# This schema is not the same as the one taken from the data,
# probably because of the indexing in the schema derived from a pandas df.
def _make_galaxy_schema(random_seds=False):
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
              pa.field('sed_val_bulge_no_host_extinction',
                       pa.list_(pa.float64()), True),
              pa.field('sed_val_disk_no_host_extinction',
                       pa.list_(pa.float64()), True),
              pa.field('internalAv_bulge', pa.float64(), True),
              pa.field('internalRv_bulge', pa.float64(), True),
              pa.field('internalAv_disk', pa.float64(), True),
              pa.field('internalRv_disk', pa.float64(), True),
              pa.field('bulge_magnorm', pa.float64(), True),
              pa.field('disk_magnorm', pa.float64(), True),
              pa.field('MW_rv', pa.float32(), True),
              pa.field('MW_av_lsst_u', pa.float32(), True),
              pa.field('MW_av_lsst_g', pa.float32(), True),
              pa.field('MW_av_lsst_r', pa.float32(), True),
              pa.field('MW_av_lsst_i', pa.float32(), True),
              pa.field('MW_av_lsst_z', pa.float32(), True),
              pa.field('MW_av_lsst_y', pa.float32(), True)]
    if random_seds:
        fields.append(pa.field('bulge_sed_file_path', pa.string(), True))
        fields.append(pa.field('disk_sed_file_path', pa.string(), True))

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
              pa.field('MW_av_lsst_u', pa.float32(), True),
              pa.field('MW_av_lsst_g', pa.float32(), True),
              pa.field('MW_av_lsst_r', pa.float32(), True),
              pa.field('MW_av_lsst_i', pa.float32(), True),
              pa.field('MW_av_lsst_z', pa.float32(), True),
              pa.field('MW_av_lsst_y', pa.float32(), True)]
    return pa.schema(fields)

def _make_MW_extinction(ra, dec, MW_rv_constant, band_dict):
    '''
    Given arrays of ra & dec,  fixed Rv and per-band column names
    and multipliers, create a MW Av column for each band
    Parameters:
    ra, dec - arrays specifying positions where Av is to be computed
    MW_rv - single constant value for Rv
    band_dict - keys are column names; values are multipliers
    Return:
    dict with keys = column names. Value for each column is array of
    Av values for a particular band at the ra,dec positions
    '''

    sfd = SFDQuery()
    ebv_raw = np.array(sfd.query_equ(ra, dec))
    av_dict = {}
    for k,v in band_dict.items():
        av_dict[k] = MW_rv_constant * v * ebv_raw

    return av_dict


class CatalogCreator:
    def __init__(self, parts, area_partition, output_dir=None, galaxy_truth=None,
                 sedLookup_dir=None, output_type='parquet', verbose=False,
                 random_sed_dir=None):
        """
        Store context for catalog creation

        Parameters
        ----------
        parts           Segments for which catalog is to be generated. If
                        partition type is HEALpix, parts will be a collection
                        of HEALpix pixels
        area_partition  Dict characterizing partition; e.g. HEALpix,
                        nside=<something>
        output_dir      Where to put created sky catalog. Defaults to
                        current directory.
        galaxy_truth    GCRCatalogs name for galaxy truth (e.g. cosmoDC2)
        sedLookup_dir   Where to find files with some per-galaxy information
                        relevant to finding and using appropriate SED file
        output_type     A format.  For now only parquet is supported
        random_sed_dir  Indicates special mode where tophat SEDs from cosmoDC2
                        are replaced with ones from a pre-determined sample.
                        Also implies addition of output columns to sky catalog
                        which contain paths to equivalent sed file for each
                        galaxy component with substituted SED.

        Might want to add a way to specify template for output file name
        and template for input sedLookup file name.
        """

        # Following directory contains per-healpix pixel files, each of which
        # has some per-galaxy information, including internal Av, Rv for
        # disk and bulge components, fluxes, mags for lsst bands, appropriate
        # index into sed_names array (which is typically around a 1000 entries,
        # whereas #galaxies is of order 10 million) and so forth.
        # If multiprocessing probably should defer this to create_pixel
        _sedLookup_dir = '/global/cfs/cdirs/lsst/groups/SSim/DC2/cosmoDC2_v1.1.4/sedLookup'
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
        if output_dir:
            self._output_dir = output_dir
        else:
            self._output_dir = os.path.abspath(os.curdir)
        if not sedLookup_dir:
            self._sedLookup_dir = _sedLookup_dir
        else:
            self._sedLookup_dir = sedLookup_dir

        self._output_type = output_type
        self._verbose = verbose
        self._random_sed_dir = random_sed_dir

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

        gal_cat = GCRCatalogs.load_catalog(self._galaxy_truth)

        self._mag_norm_f = MagNorm(Omega_c=gal_cat.cosmology.Om0,
                         Omega_b=gal_cat.cosmology.Ob0, h=gal_cat.cosmology.h,
                         sigma8=gal_cat.cosmology.sigma8,
                         n_s=gal_cat.cosmology.n_s)

        arrow_schema = _make_galaxy_schema(self._random_sed_dir)
        #print(f'Number of columns in schema: {len(arrow_schema.names)}')
        #print(arrow_schema.names)

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

        # Filename templates: input (sedLookup) and our output.
        # Hardcode for now.
        sedLookup_template = 'sed_fit_{}.h5'
        output_template = 'galaxy_{}.parquet'
        tophat_bulge_re = r'sed_(?P<start>\d+)_(?P<width>\d+)_bulge_no_host_extinction'
        tophat_disk_re = r'sed_(?P<start>\d+)_(?P<width>\d+)_disk_no_host_extinction'

        # Number of rows to include in a row group
        stride = 1000000

        # to_fetch = all columns of interest in gal_cat
        non_sed = ['galaxy_id', 'ra', 'dec', 'redshift', 'redshiftHubble',
                   'peculiarVelocity', 'shear_1', 'shear_2',
                   'convergence', 'position_angle_true',
                   'size_bulge_true', 'size_minor_bulge_true', 'sersic_bulge',
                   'size_disk_true', 'size_minor_disk_true', 'sersic_disk']
        # Find all the tophat sed numbers
        q = gal_cat.list_all_quantities()
        sed_bulge_names = [i for i in q if (i.startswith('sed') and
                                            i.endswith('bulge_no_host_extinction'))]
        sed_disk_names = [i for i in q if (i.startswith('sed') and
                                           i.endswith('disk_no_host_extinction'))]
        #Sort sed columns by start value, descending
        def _sed_bulge_key(s):
            return int(re.match(tophat_bulge_re, s)['start'])
        def _sed_disk_key(s):
            return int(re.match(tophat_disk_re, s)['start'])

        # Sort into increaing order by start wavelength
        sed_bulge_names.sort(key=_sed_bulge_key)
        sed_disk_names.sort(key=_sed_disk_key)

        #Fetch the data
        if not self._random_sed_dir:
            to_fetch = non_sed + sed_bulge_names + sed_disk_names
        else:
            to_fetch = non_sed

        # Save the 'start' and 'width' values; they'll be needed for our output
        # config.  Though we only need to do this once, not once per pixel
        tophat_parms = []
        for s in sed_bulge_names:
            m = re.match(tophat_bulge_re, s)
            if m:
                tophat_parms.append((m['start'], m['width']))

        df = gal_cat.get_quantities(to_fetch,
                                    native_filters=f'healpix_pixel=={pixel}')

        if self._random_sed_dir:
            (bulge_seds, bulge_path) = get_random_sed('bulge',
                                                      self._random_sed_dir,
                                                      len(df['galaxy_id']),
                                                      pixel=pixel)
            (disk_seds, disk_path) = get_random_sed('disk', self._random_sed_dir,
                                                    len(df['galaxy_id']),
                                                    pixel=pixel)
        else:
            # Re-form sed columns into two arrays
            bulge_seds = (np.array([df[sbn] for sbn in sed_bulge_names]).T).tolist()
            disk_seds = (np.array([df[sdn] for sdn in sed_disk_names]).T).tolist()
            #bulge_seds = np.array([df[sbn] for sbn in sed_bulge_names])
            #disk_seds = np.array([df[sdn] for sdn in sed_disk_names])

        #  Compute mag_norm from TH sed and redshift
        bulge_magnorm = [self._mag_norm_f(s[NORMWV_IX], r) for (s, r) in zip(bulge_seds, df['redshiftHubble']) ]
        disk_magnorm = [self._mag_norm_f(s[NORMWV_IX], r) for (s, r) in zip(disk_seds, df['redshiftHubble']) ]
        # Look up internal A_v, R_v
        with  h5py.File(os.path.join(self._sedLookup_dir,
                                     sedLookup_template.format(pixel))) as lookup:
            lookup_gid = np.array(lookup['galaxy_id'])
            bulge_av = np.array(lookup['bulge_av'])
            bulge_rv = np.array(lookup['bulge_rv'])
            disk_av = np.array(lookup['disk_av'])
            disk_rv = np.array(lookup['disk_rv'])
            # Old stuff - reading magnorm from sed fit file
            # Note shape of bulge_magnorm, disk_magnorm is (6, #objects)
            # Pick a middle column to use
            # magnorm_col = 3
            # bulge_magnorm = np.array(lookup['bulge_magnorm'][magnorm_col])
            # disk_magnorm =  np.array(lookup['disk_magnorm'][magnorm_col])

            ##  New: Compute magnorm ourselves from tophat values instead

        # Check that galaxies match and are ordered the same way
        cosmo_gid = np.array(df['galaxy_id'])
        if cosmo_gid.shape != lookup_gid.shape:
            print('#lookup galaxies != #cosmodc2 galaxies')
            exit(1)
        if not (cosmo_gid == lookup_gid).all():
            print('lookup galaxies differ from cosmodc2 galaxies in content or ordering')
            exit(1)

        # Assume R(V) = 3.1.  Calculate A(V) from R(V), E(B-V). See "Plotting
        # Dust Maps" example in
        # https://dustmaps.readthedocs.io/en/latest/examples.html
        MW_rv_constant = 3.1
        MW_rv = np.full_like(df['sersic_bulge'], MW_rv_constant)
        MW_av_columns = _make_MW_extinction(df['ra'], df['dec'],
                                            MW_rv_constant,
                                            _MW_extinction_bands)
        #MW_av = 2.742 * ebv_raw
        if self._verbose: print("Made extinction")

        #  Write row groups of size stride (or less) until input is exhausted
        total_row = lookup_gid.shape[0] - 1
        u_bnd = min(stride, total_row)
        l_bnd = 0
        rg_written = 0
        writer = None

        # Some columns need to be renamed and, in one case, units changed
        to_modify = ['position_angle_true', 'redshiftHubble', 'peculiarVelocity']
        print_col = True
        while u_bnd > l_bnd:
            out_dict = {k : df[k][l_bnd : u_bnd] for k in non_sed if k not in to_modify}
            out_dict['redshift_hubble'] = df['redshiftHubble'][l_bnd : u_bnd]
            out_dict['peculiar_velocity'] = df['peculiarVelocity'][l_bnd : u_bnd]
            out_dict['position_angle_unlensed'] = np.radians(df['position_angle_true']
                                                             [l_bnd : u_bnd])
            out_dict['sed_val_bulge_no_host_extinction'] = bulge_seds[l_bnd : u_bnd]
            out_dict['sed_val_disk_no_host_extinction'] = disk_seds[l_bnd : u_bnd]
            out_dict['internalAv_bulge'] = bulge_av[l_bnd : u_bnd]
            out_dict['internalRv_bulge'] = bulge_rv[l_bnd : u_bnd]
            out_dict['internalAv_disk'] = disk_av[l_bnd : u_bnd]
            out_dict['internalRv_disk'] = disk_rv[l_bnd : u_bnd]
            out_dict['bulge_magnorm'] = bulge_magnorm[l_bnd : u_bnd]
            out_dict['disk_magnorm'] = disk_magnorm[l_bnd : u_bnd]
            out_dict['MW_rv'] = MW_rv[l_bnd : u_bnd]
            for k,v in  MW_av_columns.items():
                out_dict[k] = v[l_bnd : u_bnd]

            if self._random_sed_dir:
                out_dict['bulge_sed_file_path'] = bulge_path[l_bnd : u_bnd]
                out_dict['disk_sed_file_path'] = disk_path[l_bnd : u_bnd]

            if print_col:
                for kd,i in out_dict.items():
                    print(f'Key={kd}, type={type(i)}, len={len(i)}')
            print_col = False
            ##exit(0)                                      ######### DEBUG ONLY #######

            out_df = pd.DataFrame.from_dict(out_dict)
            out_table = pa.Table.from_pandas(out_df, schema=arrow_schema)
            if not writer:
                writer = pq.ParquetWriter(os.path.join(
                    self._output_dir, output_template.format(pixel)),
                                          arrow_schema)

            writer.write_table(out_table)
            rg_written += 1
            l_bnd = u_bnd
            u_bnd = min(l_bnd + stride, total_row)

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
            create_pointsource_pixel(self, p, arrow_schema, star_cat=_star_db)
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

            nobj = len(star_df['id'])
            print(f"Found {nobj} stars")
            star_df['object_type'] = np.full((nobj,), 'star')
            star_df['host_galaxy_id'] = np.zeros((nobj,), np.int64())

            MW_rv_constant = 3.1
            star_df['MW_rv'] = np.full((nobj,), MW_rv_constant, np.float32())

            MW_av_columns = _make_MW_extinction(np.array(star_df['ra']),
                                                np.array(star_df['dec']),
                                                MW_rv_constant,
                                                _MW_extinction_bands)
            for k,v in  MW_av_columns.items():
                # No need for multiple row groups. Data size is small.
                #star_df[k] = v[l_bnd : u_bnd]
                star_df[k] = v
            out_table = pa.Table.from_pandas(star_df, schema=arrow_schema)
            if self._verbose: print("created arrow table from dataframe")

            writer = pq.ParquetWriter(os.path.join(
                output_dir, output_template.format(pixel)), arrow_schema)
            writer.write_table(out_table)

            writer.close()

        if sn_cat:
            raise NotImplementedError('SNe not yet supported. Have a nice day.')

        return
