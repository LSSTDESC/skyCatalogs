import galsim
import numpy as np
import os
import logging
import h5py
from lsstdesc_diffsky import read_diffskypop_params
from lsstdesc_diffsky.io_utils import load_healpixel
from lsstdesc_diffsky.io_utils import load_diffsky_params
from lsstdesc_diffsky.legacy.roman_rubin_2023.dsps.data_loaders.load_ssp_data import load_ssp_templates_singlemet
from lsstdesc_diffsky.legacy.roman_rubin_2023.dsps.data_loaders.defaults import SSPDataSingleMet
from lsstdesc_diffsky.defaults import OUTER_RIM_COSMO_PARAMS
from lsstdesc_diffsky.sed.disk_bulge_sed_kernels_singlemet import calc_rest_sed_disk_bulge_knot_galpop
all_diffskypop_params = read_diffskypop_params("roman_rubin_2023")

__all__ = ['DiffskySedGenerator']


def _calculate_sed_multi(send_conn, _redshift, _mah_params, _ms_params,
                         _q_params, _fbulge_params, _fknot, _ssp_data,
                         galaxy_id, n_per):
    def _calc_seds(_redshift, _mah_params, _ms_params, _q_params,
                   _fbulge_params, _fknot, _ssp_data, l_bnd, u_bnd, n_per):
        # Documentation for calculation available here:
        # https://lsstdesc-diffsky.readthedocs.io/en/latest/demo_roman_rubin_2023_seds_singlemet.html
        args = (_redshift[l_bnd:u_bnd],
                _mah_params[l_bnd:u_bnd],
                _ms_params[l_bnd:u_bnd],
                _q_params[l_bnd:u_bnd],
                _fbulge_params[l_bnd:u_bnd],
                _fknot[l_bnd:u_bnd],
                _ssp_data,
                all_diffskypop_params,
                OUTER_RIM_COSMO_PARAMS)
        return calc_rest_sed_disk_bulge_knot_galpop(*args)

    # DSPS units
    # import astropy.units as u
    # _wave_type = u.angstrom
    # _flux_type = u.Lsun / u.Hz / u.Mpc**2
    # flux_factor = (1 * _flux_type).to(galsim.SED._fnu).value
    flux_factor = 4.0204145742268754e-16

    # Iterate over chunks of n_per
    l_bnd = 0
    u_bnd = len(galaxy_id)
    out_list = []
    lb = l_bnd
    u = min(l_bnd + n_per, u_bnd)
    for i in range((u_bnd-l_bnd)//n_per+1):
        if lb == u:
            break
        sed_info = _calc_seds(_redshift, _mah_params, _ms_params, _q_params,
                              _fbulge_params, _fknot, _ssp_data, lb, u, n_per)
        # Accumulate output chunks
        out_list.append({'galaxy_id': galaxy_id[lb:u],
                         'bulge': sed_info.rest_sed_bulge*flux_factor,
                         'disk': sed_info.rest_sed_diffuse_disk*flux_factor,
                         'knots': sed_info.rest_sed_knot*flux_factor})
        lb = u
        u = min(lb + n_per, u_bnd)

    if send_conn is not None:
        send_conn.send(out_list)
    else:
        return out_list


class DiffskySedGenerator():
    '''
    Used for evaluating and storing diffsky galaxy SEDs, which are
    represented with an adaptively thinned spectrum and stored in an hdf5 file.
    Parameters
    ----------
    logname         Where to write log output
    galaxy_truth    GCRCatalogs name for galaxy truth catalog
    output_dir      Where diffsky parquet files are
    sky_cat         To get access to skyCatalogs main files already written
                    Must be a non-null SkyCatalog object.
    skip_done       If false, overwrite existing files
    auto_loop       If true, immediately loop through pixels in 'parts'
    rel_err         Target relative tolerance for flux integral.
    wave_ang_min    Minimum wavelength to keep in SEDs.
    wave_ang_max    Maximum wavelength to keep in SEDs.
    n_per           Number of SEDs to batch calculate in diffsky.
                    Memory footprint increases nonlinearly with larger n_per
    sed_out         If SEDs are to go somewhere other than usual output_dir
    parts           Pixels for which SEDs are created (only used if auto_loop
                    is True
    '''

    def __init__(self, logname='skyCatalogs.creator', galaxy_truth=None,
                 output_dir=None, sky_cat=None, skip_done=True,
                 auto_loop=False,
                 wave_ang_min=500, wave_ang_max=100000,
                 rel_err=0.03, n_per=100000,
                 sed_out=None, parts=None):
        self._output_dir = output_dir
        self._cat = sky_cat
        self._logger = logging.getLogger(logname)
        self._skip_done = skip_done
        self._n_per = n_per
        self._sed_out = sed_out

        # Setup thinned SSP templates for evaluating SED over
        # ############### Maybe temporary ############
        from pathlib import Path
        PACKAGE_SRC_DIR = os.path.dirname(os.path.abspath(str(Path(__file__))))
        SKYCATALOGDATA_ROOT = os.path.join(PACKAGE_SRC_DIR, "data")
        SINGLE_MET = os.path.join(SKYCATALOGDATA_ROOT,
                                  "dsps_ssp_data_singlemet.h5")
        # ##############

        self._get_thinned_ssp_data(rel_err, wave_ang_min, wave_ang_max,
                                   SSP_file_name=SINGLE_MET)
        import GCRCatalogs
        gal_cat = GCRCatalogs.load_catalog(galaxy_truth)

        self._hdf5_root_dir = gal_cat.get_catalog_info()['catalog_root_dir']
        self._hdf5_name_template = gal_cat.get_catalog_info()['catalog_filename_template']

        if self._output_dir is None:
            self._output_dir = sky_cat._cat_dir
        if auto_loop:
            # Loop over necessary pixels
            for p in self._parts:
                self.generate_pixel(p)

    def _get_thinned_ssp_data(self,
                              rel_err,
                              wave_ang_min,
                              wave_ang_max,
                              SSP_file_name='dsps_ssp_data_singlemet.h5'):
        """
        Return thinned SSP templates.
        Parameters
        ----------
        rel_err         Target relative tolerance for flux integral.
        wave_ang_min    Minimum wavelength to keep in SEDs.
        wave_ang_max    Maximum wavelength to keep in SEDs.

        Side-effects
        ------------
        Saves thinned SSP data structure
        """

        # Read default SSP templates
        ssp_data = load_ssp_templates_singlemet(fn=SSP_file_name)
        ssp_wave_nm = ssp_data.ssp_wave / 10.0

        # Exclude wavelengths that will never be needed
        mask = (ssp_data.ssp_wave > wave_ang_min) & (ssp_data.ssp_wave < wave_ang_max)

        # Setup galsim SED object of sum of SSP templates
        sed_lut = galsim.LookupTable(
            x=ssp_wave_nm[mask], f=np.sum(ssp_data.ssp_flux[:, mask], axis=0)
        )
        sed = galsim.SED(sed_lut, wave_type="nm",
                         flux_type="flambda", redshift=0.0)

        # Thin template sum to target tolerance.
        sed2 = sed.thin(rel_err=rel_err, fast_search=False)
        # Create mask for which wavelengths should be kept based on thinning.
        mask2 = np.where(np.in1d(ssp_wave_nm[mask], sed2.wave_list))[0]

        # Reconstruct thinned SSP templates and save.
        thin_ssp_wave_ang = ssp_data.ssp_wave[mask][mask2]
        thin_ssp_flux = ssp_data.ssp_flux[:, mask][:, mask2]
        self.ssp_data = SSPDataSingleMet(ssp_data.ssp_lg_age_gyr,
                                         thin_ssp_wave_ang, thin_ssp_flux)

    def _combine_col(self, cnt, col1, col2, col3):
        if len(np.shape(col1)) > 1:
            tmp = np.zeros((cnt, np.shape(col1)[1]), dtype=col1.dtype)
            tmp[0:len(col1), :] = col1
            tmp[len(col1):len(col1)+len(col2), :] = col2
            tmp[len(col1)+len(col2):len(col1)+len(col2)+len(col3), :] = col3
        else:
            tmp = np.zeros(cnt, dtype=col1.dtype)
            tmp[0:len(col1)] = col1
            tmp[len(col1):len(col1)+len(col2)] = col2
            tmp[len(col1)+len(col2):len(col1)+len(col2)+len(col3)] = col3
        return tmp

    def _load_diffsky_data(self, pixel):

        hdf5_file_path = os.path.join(self._hdf5_root_dir,
                                      self._hdf5_name_template.format(0, 1,
                                                                      pixel))
        mock1, metadata = load_healpixel(hdf5_file_path)
        diffsky_param_data1 = load_diffsky_params(mock1)
        hdf5_file_path = os.path.join(self._hdf5_root_dir,
                                      self._hdf5_name_template.format(1, 2,
                                                                      pixel))
        mock2, metadata = load_healpixel(hdf5_file_path)
        diffsky_param_data2 = load_diffsky_params(mock2)
        hdf5_file_path = os.path.join(self._hdf5_root_dir,
                                      self._hdf5_name_template.format(2, 3,
                                                                      pixel))
        mock3, metadata = load_healpixel(hdf5_file_path)
        diffsky_param_data3 = load_diffsky_params(mock3)
        cnt = len(mock1['galaxy_id'])+len(mock2['galaxy_id'])+len(mock3['galaxy_id'])
        galaxy_id = self._combine_col(cnt, mock1['galaxy_id'],
                                      mock2['galaxy_id'], mock3['galaxy_id'])
        redshift = self._combine_col(cnt, mock1['redshift'],
                                     mock2['redshift'], mock3['redshift'])
        mah_params = self._combine_col(cnt, diffsky_param_data1.mah_params,
                                       diffsky_param_data2.mah_params,
                                       diffsky_param_data3.mah_params)
        ms_params = self._combine_col(cnt, diffsky_param_data1.ms_params,
                                      diffsky_param_data2.ms_params,
                                      diffsky_param_data3.ms_params)
        q_params = self._combine_col(cnt, diffsky_param_data1.q_params,
                                     diffsky_param_data2.q_params,
                                     diffsky_param_data3.q_params)
        fbulge_params = self._combine_col(cnt,
                                          diffsky_param_data1.fbulge_params,
                                          diffsky_param_data2.fbulge_params,
                                          diffsky_param_data3.fbulge_params)
        fknot = self._combine_col(cnt, diffsky_param_data1.fknot,
                                  diffsky_param_data2.fknot,
                                  diffsky_param_data3.fknot)

        return galaxy_id, redshift, mah_params, ms_params, q_params, fbulge_params, fknot

    def generate_pixel(self, pixel):
        """
        Parameters
        ----------
        pixel     Pixel to generate SEDs for.
        """

        # Check that main file for this pixel exists and has some
        # data for us  If not, issue warning and return.
        object_list = self._cat.get_object_type_by_hp(pixel, 'diffsky_galaxy')
        if len(object_list) == 0:
            self._logger.warning(f'No sky catalog data for pixel {pixel}\nCannot create SEDs')
            return

        # Load diffsky galaxy data
        diffsky_galaxy_id, redshift, mah_params, \
            ms_params, q_params, fbulge_params, fknot \
            = self._load_diffsky_data(pixel)

        # Set up output file
        output_filename = f'galaxy_sed_{pixel}.hdf5'
        if self._sed_out is not None:
            output_path = os.path.join(self._sed_out, output_filename)
        else:
            output_path = os.path.join(self._output_dir, output_filename)

        if os.path.exists(output_path):
            if not self._skip_done:
                self._logger.info(f'Overwriting file {output_path}')
            else:
                self._logger.info(f'Skipping regeneration of {output_path}')
                return

        # Set up h5 file and metadata on wavelengths
        f = h5py.File(output_path, 'w', libver='latest')
        _ = f.create_dataset('meta/wave_list',
                             maxshape=(len(self.ssp_data.ssp_wave),),
                             shape=(len(self.ssp_data.ssp_wave),),
                             dtype='f4',  # compression="gzip",
                             # compression_opts=9,
                             data=self.ssp_data.ssp_wave)
        # Set up hdf5 groups
        h5_groups = {}
        for g in np.unique(diffsky_galaxy_id//100000):
            h5_groups[g] = f.create_group('galaxy/'+str(g))

        tmp_store = np.zeros((3, len(self.ssp_data.ssp_wave)), dtype='f4')
        # Loop over object collections to do SED calculations
        # If there are multiple row groups, each is stored in a separate
        # object collection. Need to loop over them
        for object_coll in object_list.get_collections():
            galaxy_id = object_coll.get_native_attribute('galaxy_id')
            rg_written = 0
            self._logger.debug(f'Handling range 0 to {len(galaxy_id)}')

            # Limit objects to those in the matching skycatalog
            mask = np.in1d(diffsky_galaxy_id, galaxy_id)

            # Build output SED data chunks
            out_list = _calculate_sed_multi(None,
                                            redshift[mask],
                                            mah_params[mask],
                                            ms_params[mask],
                                            q_params[mask],
                                            fbulge_params[mask],
                                            fknot[mask],
                                            self.ssp_data,
                                            diffsky_galaxy_id[mask],
                                            self._n_per)
            ichunk = 0
            for chunk in out_list:
                print('chunk', ichunk)
                ichunk += 1
                for igid, gid in enumerate(chunk['galaxy_id']):
                    tmp_store[0, :] = chunk['bulge'][igid, :]
                    tmp_store[1, :] = chunk['disk'][igid, :]
                    tmp_store[2, :] = chunk['knots'][igid, :]
                    _ = h5_groups[gid//100000].\
                        create_dataset(str(gid),
                                       maxshape=(3,
                                                 len(self.ssp_data.ssp_wave),),
                                       shape=(3, len(self.ssp_data.ssp_wave),),
                                       dtype='f4',  # compression="gzip",
                                       # compression_opts=9,
                                       data=tmp_store)
            rg_written += 1

        f.close()
        self._logger.debug(f'SEDs for galaxies in {rg_written} row groups have been written')
