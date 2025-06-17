import os
from pathlib import Path
import numpy as np
import pandas as pd
import healpy
from dustmaps.sfd import SFDQuery

_Av_adjustment = 2.742
_MW_rv_constant = 3.1


def make_MW_extinction_av(ra, dec):
    '''
    Given arrays of ra & dec, create a MW Av column corresponding to V-band
    correction.
    See "Plotting Dust Maps" example in
    https://dustmaps.readthedocs.io/en/latest/examples.html

    The coefficient _Av_adjustment comes from Table 6 in
    Schlafly & Finkbeiner (2011)
    See http://iopscience.iop.org/0004-637X/737/2/103/article#apj398709t6

    Parameters
    ----------
    ra, dec - arrays specifying positions where Av is to be computed
    Return:
    Array of Av values
    '''

    sfd = SFDQuery()
    ebv_raw = np.array(sfd.query_equ(np.array(ra), np.array(dec)))

    return _Av_adjustment * ebv_raw


def make_MW_extinction_rv(ra, dec):
    return np.full_like(np.array(ra), _MW_rv_constant)


PACKAGE_DIR = os.path.dirname(os.path.abspath(str(Path(__file__).parent.parent)))
DATA_DIR = os.path.join(PACKAGE_DIR, "skycatalogs", "data")


def get_trilegal_hp_nrows(hp, nside=32):
    counts_path = os.path.join(DATA_DIR, "trilegal", "star_counts.parquet")
    tbl = pd.read_parquet(counts_path)
    nrows = (tbl.query("hp_ring_id == @hp")["hp_star_count"]).values[0]
    return nrows


_NSIDE = 32
_TRILEGAL_RING_NSIDE = 256
_TRILEGAL_NEST_NSIDE = 4096


def find_trilegal_subpixels(hp, n_rows, in_nside=32, n_gps=None,
                            max_query_rows=1_000_000):
    '''
    Return subhealpixels belonging to specified healpixel for all row groups

    Parameters
    ----------
    hp             int      healpixel (ring naming) id
    in_nside       int      NSIDE value for original tiling (normally 32)
    max_query_rows int      Used to determine how many row groups will be
                            written

    Returns
    -------
    out_nside  int      32,  256 or 4096
    out_ring   boolean  True if out list ids use ring naming; false for nest
    out_list  hp ids for each row group.

    So, for example, if there are 4 row groups, out_list will have 4 elememts,
    each itself a list.   If there is only to be one row group, out_list will
    be   [ [hp] ]

    '''
    if n_rows <= max_query_rows:
        return in_nside, True, [[hp]]

    out_nside = _TRILEGAL_RING_NSIDE
    out_ring = True

    if n_rows > max_query_rows:
        # break up into several queries
        if n_rows > 60 * max_query_rows:
            n_group = 64
        elif n_rows > 30 * max_query_rows:
            n_group = 32
        elif n_rows > 16 * max_query_rows:
            n_group = 16
        else:
            n_group = 4

    if hp == 9246:
        # This healpixel, even though it is not the healpixel with the most
        # stars, apparently includes a very dense section so that 64
        # subqueries is not fine enough.  One times out.
        n_group = 256
        out_nside = _TRILEGAL_NEST_NSIDE
        out_ring = False
    elif hp == 9119:
        n_group = 64

    if n_gps:            # Use it if it was passed in
        n_group = n_gps


    def _next_level(pixel):
        return [4 * pixel, 4 * pixel + 1, 4 * pixel + 2, 4 * pixel + 3]

    pixels = [healpy.ring2nest(_NSIDE, hp)]
    current_nside = in_nside

    while current_nside < out_nside:
        pixels = [pix for p in pixels for pix in _next_level(p)]
        current_nside = current_nside * 2

    if n_group <= 64:
        subpixels = [healpy.nest2ring(_TRILEGAL_RING_NSIDE, p) for p in pixels]
    else:
        subpixels = pixels

    per_group = len(subpixels) // n_group
    out_list = []
    for i in range(n_group):
        out_list.append(subpixels[i * per_group: (i + 1) * per_group])

    return out_nside, out_ring, out_list
