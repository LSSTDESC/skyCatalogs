import os
import re
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import h5py
import GCRCatalogs

"""
Code to create a sky catalog for a particular object type
"""

pixels = [9556, 9683, 9684, 9812, 9813, 9940]
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
        print("Start on pixel ", p)
        create_pixel(p, area_partition, gal_cat, lookup, output_type, output_dir)
        print("completed pixel ", p)

def create_pixel(pixel, area_partition, gal_cat, lookup_dir, output_type, output_dir):

    # Filename templates: input (sedLookup) and our output.  Hardcode for now.
    sedLookup_template = 'sed_fit_{}.h5'
    output_template = 'galaxy_{}.parquet'
    tophat_bulge_re = r'sed_(?P<start>\d+)_(?P<width>\d+)_bulge'
    tophat_disk_re = r'sed_(?P<start>\d+)_(?P<width>\d+)_disk'

    # to_fetch = all columns of interest in gal_cat
    non_sed = ['galaxy_id', 'ra', 'dec', 'redshift', 'shear_1',
               'shear_2_phosim',
               'convergence', 'position_angle_true',
               'size_bulge_true', 'size_minor_bulge_true', 'sersic_bulge',
               'A_v_bulge', 'R_v_bulge',
               'size_disk_true', 'size_minor_disk_true', 'sersic_disk',
               'A_v_disk', 'R_v_disk']
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

    # Check that galaxies match and are ordered the same way
    cosmo_gid = np.array(df['galaxy_id'])
    if cosmo_gid.shape != lookup_gid.shape:
        print('#lookup galaxies != #cosmodc2 galaxies')
        exit(1)
    if not (cosmo_gid == lookup_gid).all():
        print('lookup galaxies differ from cosmodc2 galaxies in content or ordering')
        exit(1)

    # For first attempt make a parquet file in the simplest possible way
    out_dict = { k : df[k] for k in non_sed if k != 'position_angle_true'}
    out_dict['position_angle'] = np.radians(df['position_angle_true'])
    out_dict['sed_val_bulge'] = bulge_seds
    out_dict['sed_val_disk'] = disk_seds
    out_dict['internalAv_bulge'] = bulge_av
    out_dict['internalRv_bulge'] = bulge_rv
    out_dict['internalAv_disk'] = disk_av
    out_dict['internalRv_disk'] = disk_rv

    out_df = pd.DataFrame.from_dict(out_dict)
    out_table = pa.Table.from_pandas(out_df)

    pq.write_table(out_table, os.path.join(output_dir,
                                           output_template.format(pixel)))


# May want a base truth class for this

# class GalaxyTruth():
#         """
#         Responsible for reading from a source like CosmoDC2 catalog
#         and making available to GalaxySkyCatalog object
#         """

# Try it out
if __name__ == "__main__":
    area_partition = {'type' : 'healpix', 'ordering' : 'ring', 'nside' : 32}
    parts = [9556]
    output_dir='/global/cscratch1/sd/jrbogart/desc/skycatalogs/toy2'
    print('Starting with healpix pixel ', parts[0])
    create_galaxy_catalog(parts, area_partition, output_dir=output_dir)
    print('All done')
