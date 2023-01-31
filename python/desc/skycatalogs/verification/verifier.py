import pyarrow as pa
import pyarrow.parquet as pq
import yaml
import os
import logging

from desc.skycatalogs.utils.config_utils import open_config_file
from desc.skycatalogs.skyCatalogs import open_catalog

from desc.skycatalogs.utils.parquet_schema_utils import make_galaxy_schema, make_galaxy_flux_schema
from desc.skycatalogs.utils.parquet_schema_utils import make_star_schema, make_star_flux_schema

__all__ = ['Verifier']

_HEALPIX_PAT =  '(?P<healpix>\d+)'

class Verifier:
    def __init__(self, config_file, logname='skyCatalogs.verifier'):
        self._config_file = config_file
        #  save directory where catalog files are stored
        self._config_file = config_file
        cfg = open_config_file(config_file)
        self._cfg = cfg

        root_dir = os.getenv('SKYCATALOG_ROOT', '')
        if len(root_dir) == 0:
            root_dir = cfg.get_config_value('skycatalog_root')

        self._cat_dir = os.path.join(root_dir,
                                     cfg.get_config_value('catalog_dir'))
        self._logger = logging.getLogger(logname)

        self._healpix_pat = _HEALPIX_PAT

        if not self._logger.hasHandlers():
            loglevel = 'INFO'
            self._logger.setLevel(loglevel)
            ch = logging.StreamHandler()
            ch.setLevel(loglevel)
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            ch.setFormatter(formatter)
            self._logger.addHandler(ch)

        self._cat = None

        # Store patterns used to form file names
        # Should define user exception in case any of templates below
        # are not in config file
        g_m = self._cfg.get_config_value('object_types/galaxy/file_template')
        g_f = self._cfg.get_config_value('object_types/galaxy/flux_file_template')
        s_m = self._cfg.get_config_value('object_types/star/file_template')
        s_f = self._cfg.get_config_value('object_types/star/flux_file_template')

        patterns = {}
        patterns['galaxy_main'] = g_m.replace(self._healpix_pat, '{}')
        patterns['galaxy_flux'] = g_f.replace(self._healpix_pat, '{}')
        patterns['star_main'] = s_m.replace(self._healpix_pat, '{}')
        patterns['star_flux'] = s_f.replace(self._healpix_pat, '{}')

        self.fname_pat = patterns

    def _make_abs_path(self, hp, pattern):
        return os.path.join(self._cat_dir, pattern.format(hp))

    def verify_files_exist(self, hps):
        '''
        Check that for each hp in list all files (galaxy main, galaxy flux,
        pointsource main, pointsource flux) exist
        '''
        for hp in hps:
            i_found = 0
            for v in self.fname_pat.values():
                try:
                    f = open(self._make_abs_path(hp, v))
                except Exception as e:
                    self._logger.error(f"Failed to open file {fname}")
                    raise
                else:
                    self._logger.debug(f"Found file {self._make_abs_path(hp, v)}")
                    i_found += 1
                    f.close()
            self._logger.info(f"Found {i_found} files for hp {hp}")
        self._logger.info(f"verify_overall: All requested catalog files found")

    def verify_id_match(self, hp, object_type):
        '''
        Verify that main file and flux file for a particular healpixel
        and object type refer to the same set of objects, in the same
        order.
        Note it's not possible to use API for this since ids returned
        always come from main file
        '''

        if object_type == 'galaxy':
            id_col = 'galaxy_id'
        elif object_type == 'star':
            id_col = 'id'
        else:
            raise ValueError(f'verify_id_match: Unknown object_type {object_type}')
        m_path = self._make_abs_path(hp, self.fname_pat[object_type + '_main'])
        f_path = self._make_abs_path(hp, self.fname_pat[object_type + '_flux'])


        m_pfile = pq.ParquetFile(m_path)
        f_pfile = pq.ParquetFile(f_path)
        # For now if this fails with an exception, let the exception go off

        if m_pfile.metadata.num_rows != f_pfile.metadata.num_rows:
            self._logger.error('verify_id_match: # rows do not match')
            return False

        if m_pfile.metadata.num_row_groups != f_pfile.metadata.num_row_groups:
            self._logger.error('verify_id_match: # rows do not match')
            return False

        for i in range(m_pfile.metadata.num_row_groups):
            m_ids = m_pfile.read_row_group(i, columns=[id_col])[id_col]
            f_ids = f_pfile.read_row_group(i, columns=[id_col])[id_col]
            if m_ids != f_ids:
                return False

        return True

if __name__ == "__main__":
    #catalog_dir = 'dc2_hp'           # populated on perlmutter only
    catalog_dir = 'for_imsim_run'
    cfg_path = os.path.join(os.getenv('SCRATCH'), 'desc/skycatalogs',
                            catalog_dir, 'skyCatalog.yaml')

    ver = Verifier(cfg_path)

    ver.verify_files_exist([10066, 10067])

    for object_type in ('galaxy', 'star'):
        ok = ver.verify_id_match(10066, object_type)
        if ok:
            print(f'{object_type} ids match')
        else:
            print(f'{object_type} ids do not match!')
