import os
import shutil
import gzip
from collections import namedtuple, OrderedDict
import numpy as np
import pandas as pd
from desc.skycatalogs.skyCatalogs import SkyCatalog, open_catalog
from desc.skycatalogs.utils.config_utils import Config, open_config_file
from desc.skycatalogs.utils.sed_utils import MagNorm, convert_tophat_sed, write_sed_file
from desc.skycatalogs.utils.translate_utils import column_finder, check_file, write_to_instance, SourceType
from desc.skycatalogs.utils.translate_utils import STAR_FMT, CMP_FMT, form_star_instance_columns, form_cmp_instance_columns

'''
Guts of the code to translate sky catalog to instance catalog(s)
'''

__all__ = ['Translator']


class Translator:
    def __init__(self, visit, config_path, output_dir, object_types, band='r',
                 verbose=False, clear=False):
        self._visit = visit
        self._config = open_config_file(config_path)         # store Config object
        self._output_dir = output_dir
        self._sky_cat = open_catalog(config_path)
        types_set = set(object_types)
        if 'galaxy' in types_set:
            #types_set.remove('galaxy')
            types_set.update(set(['disk', 'bulge']))
        if 'disk' in types_set or 'bulge' in types_set:
            types_set.update(set(['galaxy']))
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
        if not self._clear: check_file(summ_path)
        f_summ =  open(summ_path, mode='w')
        handle_dict = {'summary' : (summ_path, f_summ)}
        for cmp in self._object_types:
            if cmp == 'galaxy':
                continue
            elif cmp in ['disk', 'bulge']:
                fpath = os.path.join(self._output_dir,
                                     f'{visit}/{cmp}_gal_cat_{visit}.txt')
            else:
                fpath = os.path.join(self._output_dir, f'{visit}/{cmp}_cat_{visit}.txt')
            if not self._clear: check_file(fpath)
            handle_dict[cmp] = (fpath, open(fpath, mode='w'))
        self._handle_dict = handle_dict

        if not region:
            for p in pixels:
                self.translate_pixel(p)
        else:
            self.translate_region(region)

        # Write a line to the summary file for each cmp
        for cmp in self._object_types:
            if cmp == 'galaxy':
                continue
            else:
                self._handle_dict['summary'][1].write(f'includeobj {cmp}_cat_{self._visit}.txt.gz\n')

        for v in handle_dict.values():
            v[1].close()
            if  os.path.basename(v[0]).find('instcat') == -1:
                with open(v[0], 'rb') as f_in:
                    with gzip.open(f'{v[0]}.gz', 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

    def translate_pixel(self, pixel=9556):
        '''
        Make all instance catalog entries for the healpix pixel
        '''
        if 'star' in self._object_types:
            self._translate_pointsource_pixel(pixel, 'star')



        if 'disk' in self._object_types or 'bulge' in self._object_types:
            self._translate_galaxy_pixel(pixel)

    def translate_region(self, region):
        '''
        Make all instance catalog entries for objects within the region
        '''
        raise NotImplementedError('Instance catalog creation for objects within a region not yet supported')

    def _translate_pointsource_pixel(self, pixel, obj_type):
        # Handle stars and possibly other kinds of pointsource
        if obj_type != 'star':
            raise NotImplementedError(f'_translate_pointsource_pixel: obj_type {obj_type} not currently supported')

        star_instance = form_star_instance_columns(self._band)

        star_data_columns = [q.source_parm for q in star_instance if q.source_type == SourceType.DATA]

        star_config_columns = {q.instance_name : q.source_parm for q in star_instance if q.source_type == SourceType.CONFIG}

        #  Get columns from SkyCatalog
        collections = self._sky_cat.get_objects_by_hp(pixel,
                                                      obj_type_set=set(['star'])).get_collections()
        star_collection = collections[0]
        skydata_star = star_collection.get_native_attributes(star_data_columns)
        data_len = len(skydata_star['ra'])

        # Make a new ordered dict including everything
        star_write = OrderedDict()
        for c in star_instance:
            if c.source_type == SourceType.DATA:
                star_write[c.instance_name] = skydata_star[c.source_parm]
            elif c.source_type == SourceType.FIXED:
                star_write[c.instance_name] = np.full(data_len,
                                                      c.source_parm[0],
                                                      dtype=c.source_parm[1])
            elif c.source_type == SourceType.CONFIG:
                val = self._config.get_config_value(c.source_parm[0])
                if c.instance_name == 'galacticExtinctionModel':
                    star_write[c.instance_name] = np.full(data_len, val)
                elif c.instance_name == 'galacticRv':
                    star_write[c.instance_name] = np.full(data_len,
                                                          float(val))
                else:
                    raise ValueError(f'unknown config src {c.instance_name}')
            else:
                raise ValueError(f'unknown source type {c.source_type}')

        write_to_instance(self._handle_dict['star'][1], star_write,
                          STAR_FMT)

    def _translate_galaxy_pixel(self, pixel):
        #  Get objects from SkyCatalog.   We can use the same collection for both
        # 'disk' and 'bulge', but fetch somewhat different columns
        collections = self._sky_cat.get_objects_by_hp(pixel,
                                                      obj_type_set=set(['galaxy'])).get_collections()
        cmp_collection = collections[0]

        cmps = self._object_types.intersection({'disk', 'bulge'})
        for cmp in cmps:
            sed_root_relative = f'{cmp}_sed'
            # Following is not needed since we're not writing SED files
            # sed_root = os.path.join(self._output_dir, f'{self._visit}', sed_root_relative)
            # try:
            #     os.mkdir(sed_root)
            # except FileExistsError as e:
            #     pass
            cmp_instance = form_cmp_instance_columns(cmp, self._band)
            cmp_data_columns = [q.source_parm for q in cmp_instance if q.source_type == SourceType.DATA]
            for q in cmp_instance:
                if q.source_type == SourceType.COMPUTE:
                    cmp_data_columns.extend(q.source_parm)

            cmp_config_columns = {q.instance_name : q.source_parm for q in cmp_instance if q.source_type == SourceType.CONFIG}

            skydata_cmp = cmp_collection.get_native_attributes(cmp_data_columns)
            data_len = len(skydata_cmp['ra'])
            cmp_write = OrderedDict()
            for c in cmp_instance:
                if c.source_type == SourceType.DATA:
                    cmp_write[c.instance_name] = skydata_cmp[c.source_parm]
                elif c.source_type == SourceType.FIXED:
                    cmp_write[c.instance_name] = np.full(data_len,
                                                          c.source_parm[0],
                                                         dtype=c.source_parm[1])
                elif c.source_type == SourceType.CONFIG:     # We only have three
                    val = self._config.get_config_value(c.source_parm[0])
                    if c.instance_name == 'galacticExtinctionModel':
                        cmp_write[c.instance_name] = np.full(data_len, val)
                    elif c.instance_name == 'galacticRv' :
                        cmp_write[c.instance_name] = np.full(data_len,
                                                             float(val))
                    elif c.instance_name == 'spatialmodel' :
                        cmp_write[c.instance_name] = np.full(data_len, val)
                    else:
                        raise ValueError(f'unknown config source {c.instance_name}')
                elif c.source_type == SourceType.COMPUTE:
                    if c.instance_name == 'sedFilepath':
                        #  Comment out code having to do with SED file creation
                        # SED files will be written before or after instance
                        # catalogs.  Takes too long to do on the fly
                        # and probably only a small sample is needed.
                        # They should go somewhere independent of visit
                        # th_sed = skydata_cmp[c.source_parm[0]]
                        # redshift_hubble = skydata_cmp[c.source_parm[1]]
                        # In line below should get init parms from config
                        # magnorm_f = MagNorm()
                        # th_bins = self._config.get_tophat_parameters()

                        sed_filepaths = []
                        ####for (f_nu, r_h, g_id) in zip(th_sed, redshift_hubble, skydata_cmp['galaxy_id']):
                        for g_id in skydata_cmp['galaxy_id']:
                            fname = f'{g_id}_{cmp}_sed.txt'

                            # wv, flambda, magnorm, val_500nm = convert_tophat_sed(th_bins, f_nu, magnorm_f, r_h)

                            # Include object id in the filename
                            # sedfilepath as written in the output will be
                            # relative to some known directory
                            #####write_sed_file(os.path.join(sed_root, fname), wv, flambda)
                            sed_filepaths.append(os.path.join(sed_root_relative,
                                                              fname))
                        cmp_write[c.instance_name] = np.array(sed_filepaths)
                else:
                    raise ValueError(f'unknown column source type {c.source_type}')

            write_to_instance(self._handle_dict[cmp][1], cmp_write, CMP_FMT)
