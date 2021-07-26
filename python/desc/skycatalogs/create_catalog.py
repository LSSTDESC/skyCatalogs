import os
import re
import argparse
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from astropy.coordinates import SkyCoord
import h5py
import GCRCatalogs
from desc.skycatalogs.utils.common_utils import print_date, print_callinfo

# from dm stack
from dustmaps.sfd import SFDQuery

"""
Code to create a sky catalog for a particular object type
"""

pixels = [9556, 9683, 9684, 9812, 9813, 9940]

'''
Dict of MW av column names and multipliers needed to create from ebv, MW_rv
Multipliers come from
https://iopscience.iop.org/article/10.1088/0004-637X/737/2/103#apj398709t6
appendix, table 6
'''
MW_extinction_bands = {'MW_av_lsst_u' : 4.145, 'MW_av_lsst_g' : 3.237,
                       'MW_av_lsst_r' : 2.273, 'MW_av_lsst_i' : 1.684,
                       'MW_av_lsst_z' : 1.323, 'MW_av_lsst_y' : 1.088}

def create_galaxy_catalog(parts, area_partition, galaxy_truth=None,
                          sedLookup_dir=None, output_type='parquet',
                          output_dir=None):
    """
    Parameters
    ----------
    parts           Segments for which catalog is to be generated. If partition
                    type is HEALpix, parts would be a collection of HEALpix pixels
    area_partition  Dict characterizing partition; e.g. HEALpix, nside=<something>
    galaxy_truth    GCRCatalogs name for galaxy truth (e.g. cosmoDC2)
    sedLookup_dir   Where to find files with some per-galaxy information relevant
                    to finding and using appropriate SED file
    output_type     A format.  For now only parquet is supported
    output_dir      Where to put created sky catalog. Defaults to
                    current directory.

    Might want to add a way to specify template for output file name
    and template for input sedLookup file name.

    Returns
    -------
    Dict describing catalog produced
    """

    # Following directory contains per-healpix pixel files, each of which
    # has some per-galaxy information, including internal Av, Rv for
    # disk and bulge components, fluxes, mags for lsst bands, appropriate
    # index into sed_names array (which is typically around a 1000 entries,
    # whereas #galaxies is of order 10 million) and so forth.
    _sedLookup_dir = '/global/cfs/cdirs/lsst/groups/SSim/DC2/cosmoDC2_v1.1.4/sedLookup'
    _cosmo_cat = 'cosmodc2_v1.1.4_image_addon_knots'

    if area_partition['type'] != 'healpix':
        raise NotImplementedError('Unknown partition type ')
    if output_type != 'parquet':
        raise NotImplementedError('Unknown partition type ')

    gal_truth = galaxy_truth
    if gal_truth is None:
        gal_truth = _cosmo_cat

    print('gal_truth is ', gal_truth)

    # If multiprocessing probably should defer this to create_pixel
    gal_cat = GCRCatalogs.load_catalog(gal_truth)

    lookup = sedLookup_dir
    if lookup is None:
        lookup = _sedLookup_dir

    for p in parts:
        print("Starting on pixel ", p)
        print_date()
        create_pixel(p, area_partition, gal_cat, lookup, output_type, output_dir)
        print("completed pixel ", p)
        print_date()

def create_pixel(pixel, area_partition, gal_cat, lookup_dir, output_type, output_dir):

    # Filename templates: input (sedLookup) and our output.  Hardcode for now.
    sedLookup_template = 'sed_fit_{}.h5'
    output_template = 'galaxy_{}.parquet'
    tophat_bulge_re = r'sed_(?P<start>\d+)_(?P<width>\d+)_bulge'
    tophat_disk_re = r'sed_(?P<start>\d+)_(?P<width>\d+)_disk'

    # Number of rows to include in a row group
    stride = 1000000

    # to_fetch = all columns of interest in gal_cat
    non_sed = ['galaxy_id', 'ra', 'dec', 'redshift', 'shear_1',
               'shear_2',
               'convergence', 'position_angle_true',
               'size_bulge_true', 'size_minor_bulge_true', 'sersic_bulge',
               'size_disk_true', 'size_minor_disk_true', 'sersic_disk']
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

    sed_bulge_names.sort(key=_sed_bulge_key, reverse=True)
    sed_disk_names.sort(key=_sed_disk_key, reverse=True)

    #Fetch the data
    to_fetch = non_sed + sed_bulge_names + sed_disk_names

    # Save the 'start' and 'width' values; they'll be needed for our output
    # config.  Though we only need to do this once, not once per pixel
    tophat_parms = []
    for s in sed_bulge_names:
        m = re.match(tophat_bulge_re, s)
        if m:
            tophat_parms.append((m['start'], m['width']))

    #for prm in tophat_parms:
    #    print('start={}, width={}'.format(prm[0], prm[1]))

    #exit(0)

    df = gal_cat.get_quantities(to_fetch,
                                native_filters=f'healpix_pixel=={pixel}')
    # Re-form sed columns into two arrays

    bulge_seds = (np.array([df[sbn] for sbn in sed_bulge_names]).T).tolist()
    disk_seds = (np.array([df[sdn] for sdn in sed_disk_names]).T).tolist()

    # Look up internal A_v, R_v

    with  h5py.File(os.path.join(lookup_dir,
                                 sedLookup_template.format(pixel))) as lookup:
        lookup_gid = np.array(lookup['galaxy_id'])
        bulge_av = np.array(lookup['bulge_av'])
        bulge_rv = np.array(lookup['bulge_rv'])
        disk_av = np.array(lookup['disk_av'])
        disk_rv = np.array(lookup['disk_rv'])
        # The following will be of interest when using file sed from our
        # lookup file
        # Note shape of bulge_magnorm, disk_magnorm is (6, #objects)
        # Pick a middle column to use
        magnorm_col = 3
        bulge_magnorm = np.array(lookup['bulge_magnorm'][magnorm_col])
        disk_magnorm =  np.array(lookup['disk_magnorm'][magnorm_col])
        print('bulge_magnorm shape: ', bulge_magnorm.shape)
        print('disk_magnorm shape: ', disk_magnorm.shape)

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
    MW_rv = np.full_like(disk_rv, MW_rv_constant)
    MW_av_columns = make_MW_extinction(df['ra'], df['dec'],
                                       MW_rv_constant,MW_extinction_bands)
    #MW_av = 2.742 * ebv_raw
    print("Made extinction")
    # Form arrow schema.  To get the type write, easiest to make a tiny
    # table and extract schema
    dummy_dict = { k : df[k][:10] for k in non_sed if k != 'position_angle_true'}
    dummy_dict['position_angle'] = np.radians(df['position_angle_true'][:10])
    dummy_dict['sed_val_bulge'] = bulge_seds[:10]
    dummy_dict['sed_val_disk'] = disk_seds[:10]
    dummy_dict['internalAv_bulge'] = bulge_av[:10]
    dummy_dict['internalRv_bulge'] = bulge_rv[:10]
    dummy_dict['internalAv_disk'] = disk_av[:10]
    dummy_dict['internalRv_disk'] = disk_rv[:10]

    dummy_dict['bulge_magnorm'] = bulge_magnorm[:10]
    dummy_dict['disk_magnorm'] = disk_magnorm[:10]
    dummy_dict['MW_rv'] = MW_rv[:10]

    for k,v in MW_av_columns.items():
        dummy_dict[k] = v[:10]


    dummy_df = pd.DataFrame.from_dict(dummy_dict)
    print("Created data frame from dummy_dict")
    dummy_table = pa.Table.from_pandas(dummy_df)
    print("Created pyarrow table from data frame")
    arrow_schema = dummy_table.schema

    writer = pq.ParquetWriter(os.path.join(output_dir, output_template.format(pixel)), arrow_schema)

    #  Write row groups of size stride (or less) until input is exhausted
    total_row = lookup_gid.shape[0] - 1
    u_bnd = min(stride, total_row)
    l_bnd = 0
    rg_written = 0

    while u_bnd > l_bnd:
        out_dict = { k : df[k][l_bnd : u_bnd] for k in non_sed if k != 'position_angle_true'}
        out_dict['position_angle'] = np.radians(df['position_angle_true'][l_bnd : u_bnd])
        out_dict['sed_val_bulge'] = bulge_seds[l_bnd : u_bnd]
        out_dict['sed_val_disk'] = disk_seds[l_bnd : u_bnd]
        out_dict['internalAv_bulge'] = bulge_av[l_bnd : u_bnd]
        out_dict['internalRv_bulge'] = bulge_rv[l_bnd : u_bnd]
        out_dict['internalAv_disk'] = disk_av[l_bnd : u_bnd]
        out_dict['internalRv_disk'] = disk_rv[l_bnd : u_bnd]
        out_dict['bulge_magnorm'] = bulge_magnorm[l_bnd : u_bnd]
        out_dict['disk_magnorm'] = disk_magnorm[l_bnd : u_bnd]
        out_dict['MW_rv'] = MW_rv[l_bnd : u_bnd]
        for k,v in  MW_av_columns.items():
            out_dict[k] = v[l_bnd : u_bnd]

        out_df = pd.DataFrame.from_dict(out_dict)
        out_table = pa.Table.from_pandas(out_df)
        writer.write_table(out_table)
        rg_written += 1
        l_bnd = u_bnd
        u_bnd = min(l_bnd + stride, total_row)

    writer.close()
    print("# row groups written: ", rg_written)

def make_MW_extinction(ra, dec, MW_rv_constant, band_dict):
    '''
    Given array of MW dust map values, fixed rv and per-band column names
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

# May want a base truth class for this

# class GalaxyTruth():
#         """
#         Responsible for reading from a source like CosmoDC2 catalog
#         and making available to GalaxySkyCatalog object
#         """

# Try it out
if __name__ == "__main__":
    '''
    Create sky catalogs for one or more healpixels. Invoke with --help
    for details
    '''
    # For now partitioning is fixed
    area_partition = {'type' : 'healpix', 'ordering' : 'ring', 'nside' : 32}

    parser = argparse.ArgumentParser(description='''
      Create Sky Catalogs. By default create a galaxy catalog for a
      single healpix pixel''')
    parser.add_argument('--pointsource', action='store_true',
                        help='if used, create point source catalog(s)')
    parser.add_argument('--no-galaxies', action='store_true',
                        help='if used galaxy catalogs will NOT be created')
    parser.add_argument('--pixels', type=int, nargs='*', default=[9556],
                        help='healpix pixels for which catalogs will be created')
    out_dir = os.path.join(os.getenv('SCRATCH'), 'desc', 'skycatalogs', 'test')
    parser.add_argument('--output-dir', help='directory for output files',
                        default=out_dir)
    args = parser.parse_args()

    print_callinfo('create_catalog', args)

    ##parts = pixels[0:1]
    ##output_dir='/global/cscratch1/sd/jrbogart/desc/skycatalogs/toy6'
    output_dir = args.output_dir
    parts = args.pixels
    print('Starting with healpix pixel ', parts[0])
    if not args.no_galaxies:
        print("Creating galaxy catalogs")
        create_galaxy_catalog(parts, area_partition, output_dir=output_dir)

    if args.pointsource:
        print("Creating point source catalogs")
        create_pointsource_Catalog(parts, area_partition, output_dir=output_dir)

    print('All done')
