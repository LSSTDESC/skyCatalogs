import os
from pathlib import Path
import numpy as np
import pandas as pd
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
