from collections import namedtuple

from astropy import units as u
import numpy as np

__all__ = ['Tophat', 'convert_tophat_sed', 'SCALE_FACTOR', 'write_sed_file']

# Can more or less use second equation under
# https://en.wikipedia.org/wiki/AB_magnitude#Expression_in_terms_of_f%CE%BB
# This is almost our case except for us f_nu and f_lam are not spectral flux densities:
# they're just "spectral flux".  So the scale factor is something like
#    (1 / (units for tophat value)) * c    (c in cm/sec I guess)
#    (1/4.4659) * 10^(-13) * 2.99792 * 10^10  = 6.7129 * 10^(-4)

SCALE_FACTOR = 6.7129e-4

Tophat = namedtuple('Tophat', ['start', 'width'])


def convert_tophat_sed(a_bins, f_nu, granularity=None):
    '''
    Given a tophat SED, produce an equivalent SED as list of (wavelength, f_lambda)
    Parameters
    ----------
    a_bins: list of tuples (start, width) in Angstroms
    values: list of  f_nu
    granularity: spacing between SED values.  If None, use input start values/end values
                 plus midpoints of each bin, translated from Angstrom to nm

    What to do about extrapolation?   Tophat lambdas (in nm) range from 100 to 1740.2 + 259.6,
    so about 2000.
    A SED file has min lambda = 9.1, max = 160000 (but starting at 20000 bins are 20000 wide.
    But LSST filters only cover a range comfortably within [300 nm, 1100 nm] so this shouldn't be
    an issue.   Can just start the SED file at 300 nm and run to 1100. In this case,
    don't even use the first 5 tophats or the last  4.

    return: arrays lambda, f_lambda where lambda is in nm and f_lambda
            is in erg / (cm**2 * s * ang)
            Or return lambda, f_lambda needing to be scaled, and also the scale factor
    '''


    # For each tophat value
    # Choose lambda to use.  For first draft, just use midpoint
    lam_a = np.array([ b.start + 0.5*b.width for b in a_bins])     # Angstroms
    f_nu = np.array(f_nu)

    #        convert from f_nu to f_lambda:
    # Up to a constant - universal for all SEDs - all I need to do is divide by lambda^2
    f_lam = f_nu / (lam_a * lam_a)

    if (lam_a[0] > lam_a[1]):     # reverse
        lam_a[:] = lam_a[::-1]
        f_lam[:] = f_lam[::-1]

    magnorm = 1.0         # for now
    return 0.1 * lam_a, f_lam, magnorm


def write_sed_file(path, wv, f_lambda, wv_unit=None, f_lambda_unit=None):
    '''
    Write a two-column text file.  First column is wavelength, second is luminosity value
    If units are supplied, write a comment line at the top
    Parameters
    ----------
    path           Where to write the file and what to call it
    wv             List or array of wavelength values
    f_lambda       List or array of luminosities.  Must be the same length as wv
    wv_unit        String describing units for first column
    f_lambda_unit  String describing units for second column
    '''

    header = '#  '
    if wv_unit:
        header += wv_unit + ' '
    else:
        header += ' lambda unit unknown '
    if f_lambda_unit:
        header += f_lambda_unit
    else:
        header += ' f_lambda unit unknown'
    header += '\n'
    with open(path, mode="w") as f:
        f.write(header)
        for i in range(len(wv)):
            line = '{:8.2f}  {:g}\n'.format(wv[i], f_lambda[i])
            #line = f'{wv[i]}   {f_lambda[i]}\n'
            f.write(line)
    f.close()
