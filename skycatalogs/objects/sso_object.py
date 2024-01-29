import os
import re
import numpy as np
import numpy.ma as ma
# import pyarrow.parquet as pq
import galsim
from .base_object import BaseObject, ObjectCollection, ObjectList
from ..readers import ParquetReader


__all__ = ['SsoObject', 'find_sso_files']
_MJD_EPS = 0.00002    # about 1.7 seconds


def find_sso_files(sso_path, sso_config):
    '''
    Parameters
    ----------
    sso_path    string       where sso catalog files may be found
    sso_config  dict         configuration information, including filename
                             templates

    Returns
    -------
    dict.  Keys are files (fullpath?) Each value is itself a dict including
           keys for at least mjd_min, mjd_max, path and reader (of 
           type ParquetReader, defined in parquet_reader.py). But 
           reader key is only added when needed, not in this routine.
    '''
    files_dict = dict()
    files = os.listdir(sso_path)
    for f in files:
        match = re.fullmatch(sso_config['file_template'], f)
        if match:
            new_entry = {'mjd_min': int(match['mjd_min']),
                         'mjd_max': int(match['mjd_min']),
                         'scope': 'main',
                         'path' : os.path.join(sso_path, f)}
            files_dict[f] = new_entry
            continue
        match = re.fullmatch(sso_config['flux_file_template'], f_base)
        if match:
            new_entry = {'mjd_min': int(match['mjd_min']),
                         'mjd_max': int(match['mjd_min']),
                         'scope': 'flux',
                         'path' : os.path.join(sso_path, f)}
            files_dict[f] = new_entry
            continue
    return files_dict
            

class SsoObject(BaseObject):
    _type_name = 'sso'
    _solar_sed = None

    def __init__(self, ra, dec, id, object_type, belongs_to, belongs_index,
                 mjd):
        super.__init__(ra, dec, id, self._type_name, belongs_to, belongs_index)
        self._mjd = mjd

    def _get_sed(self, mjd=None):
        '''
        returns a SED and magnorm
        mjd is required
        '''
        if SsoObject._solar_sed is None:
            SsoObject._solar_sed =\
                self._belongs_to._sky_catalog._sso_sed_factory.create()
            # For magnorm use the magnitude from Sorcha.  Can it be used
            # directly or are there other effects to be applied?
            # Have to find it by looking for entry for this id, this mjd
            # Do we look for specific entry or do we allow interpolation?
        return SsoObject._solar_sed, self.get_native_attribute('observedTrailedSourceMag')

    def get_gsobject_components(self, gsparams=None, rng=None, exposure=15.0):
        if gsparams is not None:
            gsparams = galsim.GSParams(**gsparams)
        # For a streak we use galsim.?? (galsim.Box?)
        # To get the dimensions of the box, use ra rate, dec rate and
        # exposure length.  The first two will be in the sky catalogs
        # parquet file; the last will be passed in.
        # For now start with the simpler thing: just a point source.
        return {'this_object': galsim.DeltaFunction(gsparams=gsparams)}

    def get_observer_sed_component(self, component, mjd=None):
        if mjd is None:
            mjd = self._belongs_to._mjd
        if mjd is None:
            txt = 'SsoObject.get_observer_sed_component: no mjd specified for this call\n'
            txt += 'nor when generating object list'
            raise ValueError(txt)

        sed, magnorm = self._get_sed(mjd=mjd)

        flux_500 = np.exp(-0.9210340371976184 * magnorm)
        sed = self.withMagnitude(0, self._bp500)
        sed = sed*flux_500

        # no extinction
        return sed

    def get_flux(self, bandpass, sed=None, mjd=None):
        if not sed:
            sed = self.get_total_observer_sed(mjd=mjd)
        if sed is None:
            return 0.0

        flux = sed.calculateFlux(bandpass)

        return flux


class SsoCollection(ObjectCollection):
    @staticmethod
    def find_row_groups(reader, mjd):
        num_row_groups = reader._meta.num_row_groups
        if mjd is None:
            return [i for i in range(num_row_groups)]
        rg_list = []
        for i in range(num_row_groups):
            # mjds = pq_file.read_row_group(i, columns=['mjd'])['mjd']
            mjds = reader.read_columns(['mjd'], None, i)['mjd']
            if mjd < mjds[0] - _MJD_EPS:
                continue
            if mjd >= mjds[-1] + _MJD_EPS:
                break
            rg_list.append(i)
        return rg_list

    @staticmethod
    def form_row_group_collection(skycatalog, filepath, row_group,
                                  mjd=None, region=None):
        '''
        Parameters
        ----------
        filepath
        row_group  row group to use as input
        mjd        if not None, exclude rows with mjd not within _MJD_EPS
                   of this value
        region     if not None, exclude rows with (ra, dec) not in the region

        Returns
        -------
        An SsoCollection
        '''
        if region is not None:
            raise ValueError('SsoCollection.form_row_group_collection: region parameter not yet implemented')
        reader = ParquetReader(filepath)
        values = reader.read_row_group(row_group,
                                       columns=['id', 'mjd', 'ra', 'dec'])
        if mjd:     # make mask
            mask = np.logical_or((values['mjd'] >= mjd + _MJD_EPS),
                                 (values['mjd'] < mjd - _MJD_EPS))
            id_compress = ma.array(values['id'], mask=mask).compressed()
            mjd_compress = ma.array(values['mjd'], mask=mask).compressed()
            ra_compress = ma.array(values['ra'], mask=mask).compressed()
            dec_compress = ma.array(values['dec'], mask=mask).compressed()

            if all(mask):
                return None
            new_collection = SsoCollection(ra_compress, dec_compress,
                                           id_compress, skycatalog,
                                           mjd_individual=mjd_compress,
                                           region=None, mjd_global=mjd,
                                           mask=mask,
                                           readers=[reader],
                                           row_group=row_group)
        else:
            new_collection = SsoCollection(values['ra'], values['dec'],
                                           values['id'], skycatalog,
                                           mjd_individual=values['mjd'],
                                           region=None, mask=None,
                                           readers=[reader],
                                           row_group=row_group)
        return new_collection

    @staticmethod
    def load_collection(region, skycatalog, mjd=None, filepath=None):
        '''
        region      A standard region (Disk, PolygonalRegion, Box)
                    or None
        skycatalog  An instance of the SkyCatalog class
        mjd         mjd (single float) of objects to be loaded. May be
                    None if filepath is used, in which case all objects from the
                    file will be loaded
        filepath    Restrict search to a particular file.

        Returns:    A list of SsoCollection objects
        '''
        object_list = ObjectList()
        files_dict = skycatalog._sso_files
        if filepath:
            if filepath not in files_dict:
                return []
            f_entry = files_dict[filepath]
            if mjd is not None:   # parse file name to see if mjd is included.
                if mjd < f_entry['mjd_min'] - _MJD_EPS or mjd >= f_entry['mjd_max'] + _MJD_EPS:
                    return []
            if 'reader' not in f_entry:
                f_entry['reader'] = ParquetReader(filepath)
            row_groups = SsoObject.find_row_groups(f_entry['reader'], mjd)
            for r in row_groups:
                coll = SsoObject.form_row_group_collection(skycatalog,
                                                           region, r, mjd)
                if coll:
                    object_list.append_collection(coll)
            return object_list
        else:
            # There must be an mjd.  Look up all our files (or maybe this
            # is done in advance) to determine which one or two might have
            # objects with this mjd.
            # Find the one or two row groups that are relevant. For each
            # *  fetch id, mjd, ra, dec
            # *  form mask to exclude based on
            #       mjd (must be within epsilon of specified)
            #       and on (ra, dec) (must be within region)
            # *  compress out unwanted objects
            # *  return one or two collections
            pass
        # return object_list

    def __init__(self, ra, dec, id, sky_catalog, mjd_individual=None,
                 region=None, mjd_global=None,
                 mask=None, readers=None, row_group=0):
        '''
        Parameters
        ----------
        ra, dec      array of float
        id           array of str
        sky_catalog  instance of SkyCatalog
        mjd_indiviual array of float or None
        region        Geometric region or string (representing file path)
        mjd_global    float or None
        mask          exclusion mask if cuts have been made due to
                      geometric region or mjd
        readers       parquet reader (in practice there is always only 1)
        row_group     int

        One of mjd_global, mjd_individual must not be None
        '''
        super.__init__(ra, dic, id, 'sso', None, sky_catalog,
                       region=region, mjd=mjd_global, mask=mask,
                       readers=readers, row_group=row_group)
        self._mjds = np.array(mjd_individuals)
