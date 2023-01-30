import pyarrow as pa
import yaml
import os
import logging

from desc.skycatalogs.utils.config_utils import open_config_file
from desc.skycatalogs.skyCatalogs import open_catalog

from desc.skycatalogs.utils.parquet_schema_utils import make_galaxy_schema, make_galaxy_flux_schema
from desc.skycatalogs.utils.parquet_schema_utils import make_star_schema, make_star_flux_schema

##__all__ = ['verify_overall', 'verify_schema', 'verify_coverage', 'verify_values']
__all__ = ['Verifier']

class Verifier:
    def __init__(self, config_file, logname='skyCatalogs.verifier'):
        self._config_file = config_file
        #  save directory where catalog files are stored
        cfg = open_config_file(config_file)
        self._cfg = cfg

        root_dir = os.getenv('SKYCATALOG_ROOT', '')
        if len(root_dir) == 0:
            root_dir = cfg.get_config_value('skycatalog_root')

        self._cat_dir = os.path.join(root_dir,
                                     cfg.get_config_value('catalog_dir'))
        self._logger = logging.getLogger(logname)

        if not self._logger.hasHandlers():
            loglevel = 'INFO'
            self._logger.setLevel(loglevel)
            ch = logging.StreamHandler()
            ch.setLevel(loglevel)
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            ch.setFormatter(formatter)
            self._logger.addHandler(ch)

    def verify_overall(self, hps):
        '''
        Check that for each hp in list all files (galaxy main, galaxy flux,
        pointsource main, pointsource flux) exist
        '''

        # get the templates from self._cfg
        g_m = self._cfg.get_config_value('object_types/galaxy/file_template')
        g_f = self._cfg.get_config_value('object_types/galaxy/flux_file_template')
        s_m = self._cfg.get_config_value('object_types/star/file_template')
        s_f = self._cfg.get_config_value('object_types/star/flux_file_template')
        pat = '(?P<healpix>\d+)'
        templates = [g_m.replace(pat, '{}'), g_f.replace(pat, '{}'),
                     s_m.replace(pat, '{}'), s_f.replace(pat, '{}')]

        for hp in hps:
            for t in templates:
                fname = t.format(hp)
                # try to open os.path.join(self._cat_dir, fname)
                try:
                    f = open(os.path.join(self._cat_dir, fname))
                except Exception as e:
                    self._logger.error(f"Failed to open file {fname}")
                    raise
                else:
                    f.close()
        self._logger.info(f"verify_overall: All requested catalog files found")

if __name__ == "__main__":
    #catalog_dir = 'dc2_hp'           # populated on perlmutter only
    catalog_dir = 'for_imsim_run'
    cfg_path = os.path.join(os.getenv('SCRATCH'), 'desc/skycatalogs',
                            catalog_dir, 'skyCatalog.yaml')

    ver = Verifier(cfg_path)

    ver.verify_overall([10066, 10067])
