import os
import shutil
import gzip
from collections import namedtuple, OrderedDict
import numpy as np
import pandas as pd
from desc.skycatalogs.skyCatalogs import SkyCatalog, open_catalog
from desc.skycatalogs.utils.config_utils import Config, open_config_file
from desc.skycatalogs.utils.sed_utils import MagNorm, convert_tophat_sed, write_sed_file

'''
Guts of the code to translate sky catalog to instance catalog(s)
'''

__all__ = ['Translator']

_DATA_SOURCE = 0
_CONFIG_SOURCE = 1
_FIXED = 2
_COMPUTE = 3

_STAR_FMT = '{:s} {:d} {:.14f} {:.14f} {:.8f} {:s} {:d} {:d} {:d} {:d} {:d} {:d} {:s} {:s} {:s} {:.8f} {:f}\n'

_CMP_FMT = '{:s} {:d} {:.14f} {:.14f}, {:.8f}, {:s} {:.9g} {:.9g} {:.9g} {:.9g} {:d} {:d} {:s} {:.9g} {:.9g} {:f} {:s} {:s} {:.8f} {:f}\n'

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

def _form_cmp_instance_columns(cmp, band):
    cmp_instance = [column_finder('prefix', _FIXED, ('object', np.dtype('U6'))),
                    column_finder('uniqueId', _DATA_SOURCE, 'galaxy_id'),
                    column_finder('raPhoSim', _DATA_SOURCE, 'ra'),
                    column_finder('decPhoSim', _DATA_SOURCE, 'dec'),
                    column_finder('phoSimMagNorm', _DATA_SOURCE, f'{cmp}_magnorm'),
                    column_finder('sedFilepath',_COMPUTE, [f'sed_val_{cmp}','redshift_hubble']),
                    column_finder('redshift', _DATA_SOURCE, 'redshift'),
                    column_finder('gamma1', _DATA_SOURCE, 'shear_1'),
                    column_finder('gamma2', _DATA_SOURCE, 'shear_2'),
                    column_finder('kappa', _DATA_SOURCE, 'convergence'),
                    column_finder('raOffset', _FIXED, (0, int)),
                    column_finder('decOffset', _FIXED, (0, int)),
                    column_finder('spatialmodel', _CONFIG_SOURCE,
                                  f'object_types/{cmp}_basic/spatial_model'),
                    column_finder('majorAxis', _DATA_SOURCE, f'size_{cmp}_true'),
                    column_finder('minorAxis', _DATA_SOURCE, f'size_minor_{cmp}_true'),
                    column_finder('positionAngle', _DATA_SOURCE, 'position_angle_unlensed'),
                    column_finder('internalExtinctionModel', _FIXED, ('none', np.dtype('U4'))),
                    column_finder('galacticExtinctionModel', _CONFIG_SOURCE,
                                  f'object_types/{cmp}_basic/MW_extinction'),
                    column_finder('galactivAv', _DATA_SOURCE,
                                  f'MW_av_lsst_{band}'),
                    column_finder('galacticRv', _CONFIG_SOURCE,
                                  'MW_extinction_values/r_v/value')]
    return cmp_instance

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
        if not self._clear: _check_file(summ_path)
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
            if  os.path.basename(v[0]).find('instcat') == -1:
                with open(v[0], 'rb') as f_in:
                    with gzip.open(f'{v[0]}.gz', 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

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
            star_collection = collections[0]
            skydata_star = star_collection.get_attributes(star_data_columns)
            data_len = len(skydata_star['ra'])

            # Make a new ordered dict including everything
            star_write = OrderedDict()
            for c in star_instance:
                if c.source_type == _DATA_SOURCE:
                    star_write[c.instance_name] = skydata_star[c.source_parm]
                elif c.source_type == _FIXED:
                    star_write[c.instance_name] = np.full(data_len,
                                                          c.source_parm[0],
                                                          dtype=c.source_parm[1])
                elif c.source_type == _CONFIG_SOURCE:
                    val = self._sky_cat.get_config_value(c.source_parm)
                    if c.instance_name == 'galacticExtinctionModel':
                        star_write[c.instance_name] = np.full(data_len, val)
                    elif c.instance_name == 'galacticRv':
                        star_write[c.instance_name] = np.full(data_len,
                                                              float(val))
                    else:
                        raise ValueError(f'unknown config src {c.instance_name}')
                else:
                    raise ValueError(f'unknown source type {c.source_type}')

            _write_to_instance(self._handle_dict['star'][1], star_write,
                               _STAR_FMT)

            # Also write a line to the summary file
            self._handle_dict['summary'][1].write(f'includeobj star_cat_{self._visit}.txt.gz\n')

        if 'disk' in self._object_types or 'bulge' in self._object_types:
            self._translate_galaxy_pixel(pixel)

    def _translate_galaxy_pixel(self, pixel):
        #  Get objects from SkyCatalog.   We can use the same collection for both
        # 'disk' and 'bulge', but fetch somewhat different columns
        collections = self._sky_cat.get_objects_by_hp(0, pixel,
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
            cmp_instance = _form_cmp_instance_columns(cmp, self._band)
            cmp_data_columns = [q.source_parm for q in cmp_instance if q.source_type == _DATA_SOURCE]
            for q in cmp_instance:
                if q.source_type == _COMPUTE:
                    cmp_data_columns.extend(q.source_parm)
            ####cmp_data_columns.extend([q.source_parm for q in cmp_instance if q.source_type == _COMPUTE])

            cmp_config_columns = {q.instance_name : q.source_parm for q in cmp_instance if q.source_type == _CONFIG_SOURCE}

            skydata_cmp = cmp_collection.get_attributes(cmp_data_columns)
            data_len = len(skydata_cmp['ra'])
            cmp_write = OrderedDict()
            for c in cmp_instance:
                if c.source_type == _DATA_SOURCE:
                    cmp_write[c.instance_name] = skydata_cmp[c.source_parm]
                elif c.source_type == _FIXED:
                    cmp_write[c.instance_name] = np.full(data_len,
                                                          c.source_parm[0],
                                                         dtype=c.source_parm[1])
                elif c.source_type == _CONFIG_SOURCE:     # We only have three
                    val = self._sky_cat.get_config_value(c.source_parm)
                    if c.instance_name == 'galacticExtinctionModel':
                        cmp_write[c.instance_name] = np.full(data_len, val)
                    elif c.instance_name == 'galacticRv' :
                        cmp_write[c.instance_name] = np.full(data_len,
                                                             float(val))
                    elif c.instance_name == 'spatialmodel' :
                        cmp_write[c.instance_name] = np.full(data_len, val)
                    else:
                        raise ValueError(f'unknown config source {c.instance_name}')
                elif c.source_type == _COMPUTE:
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

            _write_to_instance(self._handle_dict[cmp][1], cmp_write, _CMP_FMT)

            # Also write a line to the summary file
            self._handle_dict['summary'][1].write(f'includeobj {cmp}_cat_{self._visit}.txt.gz\n')
