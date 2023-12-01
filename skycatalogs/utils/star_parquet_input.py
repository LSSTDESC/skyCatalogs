# import pyarrow as pa
import os
import pyarrow.parquet as pq
import pandas as pd
import numpy as np
import numpy.ma as ma
import healpy


def _find_star_files(dirpath, pixel):
    '''
    Parameters
    ----------
    dirpath    string       Path to directory containing input parquet files
                            HTM partitioning
    pixel      int          Healpix pixel

    Returns
    -------
    List of parquet file paths for files which might have sources belonging
    to specified pixel
    '''
    # For now, hard-code to 3 files which have everything needed
    # for joint simulation region

    fnames = ['stars_chunk_8800000000000_8900000000000.parquet',
              'stars_chunk_9200000000000_9300000000000.parquet',
              'stars_chunk_9700000000000_9800000000000.parquet']

    return [os.path.join(dirpath, fname) for fname in fnames]


def _calculate_pixel_mask(ra, dec, pixel, nside=32):
    '''
    Parameters
    ----------
    ra          Array of float
    dec         Array of float
    pixel       int

    Returns
    -------
    Mask, set to 1 for all elements not in the pixel
    '''
    ra = np.array(ra)
    dec = np.array(dec)
    in_pix = healpy.pixelfunc.ang2pix(nside, ra, dec, nest=False, lonlat=True)
    mask = ma.logical_not(in_pix == pixel)
    return mask


def _star_parquet_reader(dirpath, pixel, output_arrow_schema):
    # Get requisite info from parquet files for sources in pixel.
    # Next do renames and calculation for magnorm
    to_read = ['simobjid', 'ra', 'decl', 'mura', 'mudecl', 'vrad',
               'parallax', 'sedfilename', 'flux_scale', 'ebv']
    rename = {'decl': 'dec', 'sedfilename': 'sed_filepath', 'mudecl': 'mudec',
              'vrad': 'radial_velocity'}
    out_fields = ['id', 'ra', 'dec', 'mura', 'mudec', 'radial_velocity',
                  'parallax', 'sed_filepath', 'magnorm', 'ebv']

    paths = _find_star_files(dirpath, pixel)
    out_dict = {}
    for k in out_fields:
        out_dict[k] = []
    for f in paths:
        f_dict = {}
        for k in out_fields:
            f_dict[k] = []
        meta = pq.read_metadata(f)
        pq_file = pq.ParquetFile(f)
        for rg in range(meta.num_row_groups):
            tbl = pq_file.read_row_group(rg, columns=to_read)
            msk = _calculate_pixel_mask(tbl['ra'], tbl['decl'], pixel)
            for col in to_read:
                if col in rename:
                    f_dict[rename[col]] += list(ma.array(np.array(tbl[col]),
                                                         mask=msk).compressed())
                elif col == 'flux_scale':
                    # compute magnorm from flux_scale
                    tmp = ma.array(np.array(tbl['flux_scale']),
                                   mask=msk).compressed()
                    f_dict['magnorm'] += list(-2.5*np.log(tmp)/np.log(10.0) - 18.402732642)
                elif col == 'simobjid':
                    # convert simobjid to string, change name to id
                    tmp = [str(sim_id) for sim_id in tbl['simobjid']]
                    f_dict['id'] += list(ma.array(np.array(tmp),
                                                  mask=msk).compressed())
                else:
                    f_dict[col] += list(ma.array(np.array(tbl[col]),
                                                 mask=msk).compressed())
        for k in out_fields:
            out_dict[k] += f_dict[k]

    return pd.DataFrame(out_dict)
