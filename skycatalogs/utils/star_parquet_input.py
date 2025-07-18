import os
import glob
import pyarrow.parquet as pq
import pandas as pd
import numpy as np
import numpy.ma as ma
import healpy
from esutil.htm import HTM


class UWStarFiles:
    _files = {}  # shared state

    def __init__(self, input_dir):
        self._index_files(input_dir)
        self.htm_indexer = HTM(depth=20)

    def _index_files(self, input_dir):
        if self._files:
            # The files have already been indexed.
            return
        files = sorted(glob.glob(os.path.join(input_dir,
                                              'stars_chunk_*.parquet')))
        for item in files:
            tokens = os.path.basename(item).split('_')
            imin = int(tokens[2])
            imax = int(tokens[3].split('.')[0])
            self._files[(imin, imax)] = item

    def find_files(self, pixel, nside=32, res_factor=16):
        """
        Find the UW input files that cover a given healpixel by sampling the
        interior of the healpixel and computing the htm indexes at those
        locations.

        Parameters
        ----------
        pixel : int
            Healpixel id.
        nside : int
            Resolution parameter of the healpixel map.
        res_factor : int
            Sampling factor of the healpix to ensure coverage of the healpixel
            by the UW files.

        Returns
        -------
        list : List of filenames.
        """
        corners = np.transpose(healpy.boundaries(nside, pixel))
        subpixels = healpy.query_polygon(nside*res_factor, corners)
        ra, dec = [], []
        for subpix in subpixels:
            coord = healpy.pix2ang(nside*res_factor, subpix, lonlat=True)
            ra.append(coord[0])
            dec.append(coord[1])
        files = set()
        indices = set(self.htm_indexer.lookup_id(ra, dec))
        for index in indices:
            for (imin, imax), item in self._files.items():
                if imin <= index <= imax:
                    files.add(item)
        return files


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


def _star_parquet_reader(dirpath, pixel, output_arrow_schema, nside=32):
    # Get requisite info from parquet files for sources in pixel.
    # Next do renames and calculation for magnorm
    to_read = ['simobjid', 'ra', 'decl', 'mura', 'mudecl', 'vrad',
               'parallax', 'sedfilename', 'flux_scale', 'ebv']
    rename = {'decl': 'dec', 'sedfilename': 'sed_filepath', 'mudecl': 'mudec',
              'vrad': 'radial_velocity'}
    out_fields = ['id', 'ra', 'dec', 'mura', 'mudec', 'radial_velocity',
                  'parallax', 'sed_filepath', 'magnorm', 'ebv']

    uw_files = UWStarFiles(dirpath)
    paths = uw_files.find_files(pixel, nside)
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
