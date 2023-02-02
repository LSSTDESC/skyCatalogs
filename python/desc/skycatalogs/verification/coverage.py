import healpy as hp
import numpy as np
import logging

__all__ = ['locate_hp', 'accumulate', 'merge_dicts']

def locate_hp(ra, dec, nside=1024):
    '''
    Given lonlat coordinates find hp ids.
    nside of 1024 has resolution of about 0.057 degrees

    returns array of healpix ids
    '''
    return hp.pixelfunc.ang2pix(nside, ra, dec, lonlat=True)

def accumulate(input_dict, ids, logname='skyCatalogs.verifier'):
    '''
    Given array of (integer) ids and input dict (keyed by id;
    value = occurrence count), update dict taking into account ids array
    dict
    Return dict formed just from new ids
    '''
    unique, unique_counts = np.unique(ids, return_counts=True)
    new_dict = {unique[i] : unique_counts[i] for i in range(len(unique))}
    overlap = set(new_dict.keys()).intersection(set(input_dict.keys()))
    logger = logging.getLogger(logname)
    logger.debug(f'# overlap entries: {len(overlap)}')
    new_entries = {k : new_dict[k] for k in new_dict.keys() - input_dict.keys()}
    input_dict.update(new_entries)

    for k in overlap:
        input_dict[k] = input_dict[k] + new_dict[k]

    return new_dict

def merge_dicts(summary_dict, new_dict):
    overlap = set(new_dict.keys()).intersection(set(summary_dict.keys()))
    new_entries = {k : new_dict[k] for k in new_dict.keys() - summary_dict.keys()}
    summary_dict.update(new_entries)

    for k in overlap:
        summary_dict[k] = summary_dict[k] + new_dict[k]
