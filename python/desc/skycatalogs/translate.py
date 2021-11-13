import os
import numpy as np
import pandas as pd
from desc.skycatalogs.skyCatalogs import SkyCatalog, open_catalog

'''
Guts of the code to translate sky catalog to instance catalog(s)
'''

__all__ = ['Translator']

class Translator:
    def __init__(self, visit, config_path, output_dir, object_types, band='r',
                 verbose=False):
        self._visit = visit
        self._config = config_path
        self._output_dir = output_dir
        self._sky_cat = open_catalog(config_path)
        types_set = set(object_types)
        if 'galaxy' in types_set:
            types_set.remove('galaxy')
            types_set.update(set(['disk', 'bulge']))
        self._object_types = types_set

    def translate_visit(self, pixels=[9556], region=None):
        '''
        For now the area of a sky associated with a "visit" is determined
        by a set of healpix pixels. Should also allow for a disk defined
        by ra,dec and radius (in which case value for region maybe supersedes
        value for pixels?)
        '''
        # Make a visit subdirectory. It's ok if it already exists
        try:
            os.mkdir(os.path.join(self._output_dir, str(self._visit)))
        except FileExistsError as e:
            pass

        exit(0)            ### DEBUG: See if we can get this far

        # Open summary file and a file for each object type.
        # Save handles in a dict
        f_summ =  open(os.path.join(self._output_dir,
                                    f'instcat{self._visit}.txt'),
                       mode='w')
        handle_dict = {'summary' : f_summ}
        visit = self._visit
        for cmp in types_set:
            if cmp in ['disk', 'bulge']:
                fpath = os.path.join(self._output_dir, f'{cmp}_gal_cat_{visit}.txt')
            else:
                fpath = os.path.join(self._output_dir, f'{cmp}_cat_{visit}.txt')
            handle_dict[cmp] = open(fpath, mode='w')



    def translate_pixel(self, pixel=9556):
        '''
        Make all instance catalog entries for the healpix pixel
        '''
        if 'star' in self._object_types:
            self._translate_pointsource(pixel)
        if 'disk' in self._object_types or 'bulge' in self._object_types:
            self._translate_galaxy(pixel)

        self._generate
