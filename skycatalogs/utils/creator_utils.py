import numpy as np
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
