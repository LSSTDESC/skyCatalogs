import os
from collections import namedtuple, OrderedDict
import numpy as np
import pandas as pd
from desc.skycatalogs.skyCatalogs import SkyCatalog, open_catalog

'''
Guts of the code to translate sky catalog to instance catalog(s)
'''

__all__ = ['Translator']

_DATA_SOURCE = 0
_CONFIG_SOURCE = 1
_FIXED = 2

_STAR_FMT = '{:s} {:d} {:.14f} {:.14f} {:.8f} {:s} {:d} {:d} {:d} {:d} {:d} {:d} {:s} {:s} {:s} {:.8f} {:f}\n'

column_finder = namedtuple('ColumnFinder', ['instance_name', 'source_type',
                                            'source_parm'])

def _check_file(path):
    '''Look for a file that should not exist'''
    try:
        f = open(path, mode='r')
    except FileNotFoundError as e:
        return

    raise ValueError(f'File for {path} already exists')

def _form_star_instance_columns(band):
    star_instance = [column_finder('prefix', _FIXED, ('object', np.dtype('U6'))),
                     column_finder('uniqueId', _DATA_SOURCE, 'id'),
                     column_finder('raPhoSim', _DATA_SOURCE, 'ra'),
                     column_finder('decPhoSim', _DATA_SOURCE, 'dec'),
                     column_finder('maskedMagNorm', _DATA_SOURCE, 'magnorm'),
                     column_finder('sedFilepath',_DATA_SOURCE, 'sed_filepath'),
                     column_finder('redshift', _FIXED, (0, int)),
                     column_finder('gamma1', _FIXED, (0, int)),
                     column_finder('gamma2', _FIXED, (0, int)),
                     column_finder('kappa', _FIXED, (0, int)),
                     column_finder('raOffset', _FIXED, (0, int)),
                     column_finder('decOffset', _FIXED, (0, int)),
                     column_finder('spatialmodel', _FIXED, ('point', np.dtype('U5'))),
                     column_finder('internalExtinctionModel', _FIXED, ('none', np.dtype('U4'))),
                     column_finder('galacticExtinctionModel', _CONFIG_SOURCE,
                                   'object_types/star/MW_extinction'),
                     column_finder('galactivAv', _DATA_SOURCE,
                                   f'MW_av_lsst_{band}'),
                     column_finder('galacticRv', _CONFIG_SOURCE,
                                   'MW_extinction_values/r_v/value')]
    return star_instance

def _write_to_instance(handle, ordered_data_dict, fmt):
    '''
    Parameters
    ----------
    handle             file handle to which text will be written
    ordered_data_dict  columns in the order they should be written
    '''
    col_list = [ordered_data_dict[k] for k in ordered_data_dict]
    for row in zip(*col_list):
        handle.write(fmt.format(*row))

class Translator:
    def __init__(self, visit, config_path, output_dir, object_types, band='r',
                 verbose=False, clear=False):
        self._visit = visit
        self._config = config_path
        self._output_dir = output_dir
        self._sky_cat = open_catalog(config_path)
        types_set = set(object_types)
        if 'galaxy' in types_set:
            types_set.remove('galaxy')
            types_set.update(set(['disk', 'bulge']))
        self._object_types = types_set
        self._band = band
        self._verbose = verbose
        self._clear = clear

    def translate_visit(self, pixels=[9556], region=None):
        '''
        For now the area of a sky associated with a "visit" is determined
        by a set of healpix pixels. Should also allow for a disk defined
        by ra,dec and radius (in which case value for region maybe supersedes
        value for pixels?)
        '''
        # Make a visit subdirectory. It's ok if it already exists
        d = os.path.join(self._output_dir, str(self._visit))
        try:
            os.mkdir(d)
        except FileExistsError as e:
            if self._clear:
                fnames = os.listdir(d)
                for f in fnames:
                    if f.endswith(f'_{visit}.txt'):
                        os.rm(f)
            pass

        # Open summary file and a file for each object type.
        visit = self._visit
        summ_path = os.path.join(self._output_dir, f'{visit}/instcat_{visit}.txt')
        if not self._clear: _check_file(summ_path)
        f_summ =  open(summ_path, mode='w')
        handle_dict = {'summary' : (summ_path, f_summ)}
        for cmp in self._object_types:
            if cmp in ['disk', 'bulge']:
                fpath = os.path.join(self._output_dir,
                                     f'{cmp}_gal_cat_{visit}.txt')
            else:
                fpath = os.path.join(self._output_dir, f'{visit}/{cmp}_cat_{visit}.txt')
            if not self._clear: _check_file(fpath)
            handle_dict[cmp] = (fpath, open(fpath, mode='w'))
        self._handle_dict = handle_dict

        if not region:
            for p in pixels:
                self.translate_pixel(p)
        else:
            self.translate_region(region)

        for v in handle_dict.values():
            v[1].close()

    def translate_pixel(self, pixel=9556):
        '''
        Make all instance catalog entries for the healpix pixel
        '''
        if 'star' in self._object_types:
            star_instance = _form_star_instance_columns(self._band)

            star_data_columns = [q.source_parm for q in star_instance if q.source_type == _DATA_SOURCE]

            star_config_columns = {q.instance_name : q.source_parm for q in star_instance if q.source_type == _CONFIG_SOURCE}

            #  Get columns from SkyCatalog
            collections = self._sky_cat.get_objects_by_hp(0, pixel,
                                                          obj_type_set=set(['star'])).get_collections()
            #print('Len of collections: ', len(collections))
            star_collection = collections[0]
            #print('Type of star_collection', type(star_collection))
            skydata_star = star_collection.get_attributes(star_data_columns)
            data_len = len(skydata_star['ra'])
            #print(f'skydata_star #columns: {len(skydata_star)} column length: {data_len}')
            extinction_model = self._sky_cat.get_config_value(star_config_columns['galacticExtinctionModel'])
            galacticRv = self._sky_cat.get_config_value(star_config_columns['galacticRv'])
            # Make a new ordered dict including everything
            star_write = OrderedDict()
            data_ix = 0
            key_list = [k for k in skydata_star.keys()]
            for c in star_instance:
                if c.source_type == _DATA_SOURCE:
                    star_write[c.instance_name] = skydata_star[key_list[data_ix]]
                    data_ix = data_ix + 1
                elif c.source_type == _FIXED:
                    star_write[c.instance_name] = np.full(data_len,
                                                          c.source_parm[0], dtype=c.source_parm[1])
                else:   # _CONFIG_SOURCE.  We only have two
                    val = self._sky_cat.get_config_value(c.source_parm)
                    if c.instance_name == 'galacticExtinctionModel':
                        star_write[c.instance_name] = np.full(data_len, val)
                    else:     # it's r_v
                        star_write[c.instance_name] = np.full(data_len, int(val))

            _write_to_instance(self._handle_dict['star'][1], star_write, _STAR_FMT)

            # Also write a line to the summary file
            self._handle_dict['summary'][1].write(f'includeobj star_cat_{self._visit}.txt.gz\n')

        if 'disk' in self._object_types or 'bulge' in self._object_types:
            self._translate_galaxy(pixel)

        # Not sure where I was going with this
        #self._generate
