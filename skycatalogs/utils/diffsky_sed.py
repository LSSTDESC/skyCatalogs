"""Minimal Diffsky SED calculation used by SkyCatalogs.

This module intentionally lives outside both Diffsky and OpenCosmo.  It uses
their lower-level kernels, but returns only the disk, bulge, and knot spectra
needed by SkyCatalogs. The implementation is a candidate for eventual
upstreaming to Diffsky once its numerical equivalence and performance have
been established.

All Diffsky imports are local to callables so importing :mod:`skycatalogs`
does not make Diffsky or OpenCosmo package-wide dependencies.
"""

from collections import namedtuple
from functools import partial
import logging
import os

import numpy as np


def create_diffsky_sed_factory(config, catalog_dir, cosmology):
    """Construct the runtime factory from one Diffsky config section."""
    from skycatalogs.utils.sed_tools import DiffskySedFactory

    state_dir = config.get('sed_state_dir', 'diffsky_runtime')
    if not os.path.isabs(state_dir):
        state_dir = os.path.join(catalog_dir, state_dir)
    return DiffskySedFactory(
        state_dir, cosmology,
        object_batch_size=config.get('sed_object_batch_size', 256),
        diffsky_batch_size=config.get('sed_diffsky_batch_size', 25),
        cache_size=config.get('sed_cache_size', 8),
        rel_err=config.get('sed_rel_err', 0.03),
        wave_ang_min=config.get('sed_wave_ang_min', 500),
        wave_ang_max=config.get('sed_wave_ang_max', 100000),
        pixel_cache_size=config.get('sed_pixel_cache_size', 4),
        ssp_wave_min_micron=(
            config.get('sed_ssp_wave_min_micron') or 0.06),
        ssp_wave_max_micron=(
            config.get('sed_ssp_wave_max_micron') or 2.34),
        sed_engine=config.get('sed_engine', 'fast'),
        sed_precision=config.get('sed_precision'),
        thinning_mode=config.get('sed_thinning_mode', 'galsim'),
        emission_lines_angstrom=config.get('sed_emission_lines_angstrom'),
        emission_line_half_width_angstrom=config.get(
            'sed_emission_line_half_width_angstrom', 20.0))


_LOGGER = logging.getLogger(__name__)
_FAST_COMPILED_SHAPES = set()


def _fast_batch_shape(n_input, minimum=32):
    """Return a power-of-two JAX batch shape covering ``n_input`` rows."""
    if n_input < 1:
        raise ValueError("n_input must be positive")
    target = max(int(n_input), int(minimum))
    return 1 << (target - 1).bit_length()


def _make_fast_kernel(precision):
    """Build and return the minimal JAX disk/bulge/knot SED kernel."""
    import jax.numpy as jnp
    from diffmah import DiffmahParams
    from diffstar import DiffstarParams
    from jax import jit, tree_util

    from diffsky.experimental import mc_diffstarpop_wrappers as mcdw
    from diffsky.experimental.kernels import dbk_kernels, mc_randoms
    from diffsky.experimental.kernels import rapid_quenching as rq
    from diffsky.experimental.kernels import ssp_weight_kernels as sspwk
    from diffsky.experimental.kernels.constants import LGMET_SCATTER
    from diffsky.merging import merging_kernels, merging_model
    from diffsky.ssp_err_model import ssp_err_model

    FastSEDResult = namedtuple(
        "FastSEDResult",
        ("rest_sed_bulge", "rest_sed_disk", "rest_sed_knots"),
    )
    calculation_dtype = {"float32": jnp.float32,
                         "float64": jnp.float64}[precision]

    def cast_floating_leaf(value):
        if value is None:
            return None
        array = jnp.asarray(value)
        if jnp.issubdtype(array.dtype, jnp.floating):
            return array.astype(calculation_dtype)
        return array

    def cast_tree(value):
        return tree_util.tree_map(cast_floating_leaf, value,
                                  is_leaf=lambda item: item is None)

    @partial(jit, static_argnames=("n_t_table",))
    def kernel(
        phot_randoms,
        dbk_randoms,
        merging_randoms,
        sfh_params,
        z_obs,
        t_obs,
        mah_params,
        ssp_data,
        mzr_params,
        spspop_params,
        scatter_params,
        ssperr_params,
        merging_params,
        cosmo_params,
        fb,
        logmp_infall,
        logmhost_infall,
        t_infall,
        is_central,
        sat_weights,
        halo_indx,
        mc_merge,
        *,
        n_t_table=mcdw.N_T_TABLE,
    ):
        """Return only the three component SEDs.

        This follows Diffsky's DBK SED kernel, with two deliberate changes:

        * the total in-situ spectrum is obtained by summing the three
          component spectra instead of performing a fourth SSP contraction;
        * large diagnostic/intermediate values and merged SSP weights are not
          returned or computed when they are not needed for the component
          spectra.
        """
        upid = jnp.where(is_central == 1, -1, halo_indx).astype(int)
        lgmu_infall = logmp_infall - logmhost_infall
        gyr_since_infall = t_obs - t_infall

        t_table, sfh_table, logsm_obs, logssfr_obs = \
            mcdw.compute_diffstar_info(
                mah_params, sfh_params, t_obs, cosmo_params, fb, n_t_table)

        t_infall_from_age = t_obs - gyr_since_infall
        logmp_infall_from_ratio = lgmu_infall + logmhost_infall
        p_merge_smooth = merging_model.get_p_merge_from_merging_params(
            merging_params, logmp_infall_from_ratio, logmhost_infall,
            t_obs, t_infall_from_age, upid)

        smooth_weights = rq.get_smooth_ssp_weights_rq(
            t_table, sfh_table, logsm_obs, ssp_data, t_obs, mzr_params,
            LGMET_SCATTER, p_merge_smooth)
        burstiness = rq.get_burstiness_rq(
            phot_randoms.uran_pburst, phot_randoms.mc_is_q, logsm_obs,
            logssfr_obs, smooth_weights.age_weights,
            smooth_weights.lgmet_weights, ssp_data,
            spspop_params.burstpop_params, p_merge_smooth)

        n_gals = z_obs.size
        wave_eff = jnp.tile(ssp_data.ssp_wave, n_gals).reshape(
            (n_gals, -1))
        dust_frac_trans, _ = sspwk.compute_dust_attenuation(
            phot_randoms.uran_av, phot_randoms.uran_delta,
            phot_randoms.uran_funo, logsm_obs, logssfr_obs, ssp_data,
            z_obs, wave_eff, spspop_params.dustpop_params, scatter_params)
        frac_ssp_errors = ssp_err_model.frac_ssp_err_at_z_obs_galpop(
            ssperr_params, logsm_obs, z_obs, wave_eff)
        frac_ssp_errors = ssp_err_model.get_noisy_frac_ssp_errors(
            wave_eff, frac_ssp_errors,
            phot_randoms.delta_mag_ssp_scatter)

        age_weights = jnp.sum(burstiness.ssp_weights_mc, axis=1)
        lgmet_weights = jnp.sum(burstiness.ssp_weights_mc, axis=2)
        dbk_weights, _ = dbk_kernels._dbk_kern(
            t_obs, ssp_data, t_table, sfh_table,
            burstiness.burst_params_mc, lgmet_weights, dbk_randoms,
            logsm_obs, age_weights, p_merge_smooth)

        n_met, n_age, n_wave = ssp_data.ssp_flux.shape
        attenuation = dust_frac_trans.swapaxes(1, 2)
        component_weights = jnp.stack((
            dbk_weights.ssp_weights_bulge,
            dbk_weights.ssp_weights_disk,
            dbk_weights.ssp_weights_knots,
        ), axis=1)
        component_masses = jnp.stack((
            dbk_weights.mstar_bulge,
            dbk_weights.mstar_disk,
            dbk_weights.mstar_knots,
        ), axis=1)
        component_seds = jnp.einsum(
            "gaw,gcma,maw,gw->gcw",
            attenuation, component_weights, ssp_data.ssp_flux,
            frac_ssp_errors)
        component_seds = component_seds * component_masses[:, :, None]
        sed_bulge = component_seds[:, 0, :]
        sed_disk = component_seds[:, 1, :]
        sed_knots = component_seds[:, 2, :]

        total_in_situ = sed_bulge + sed_disk + sed_knots
        mc_p_merge = merging_kernels.get_mc_p_merge(
            merging_randoms.uran_pmerge, p_merge_smooth)
        p_merge = jnp.where(mc_merge < 1, p_merge_smooth, mc_p_merge)
        total_merged = merging_kernels.compute_x_tot_from_x_in_situ(
            total_in_situ, p_merge[:, None], sat_weights[:, None], halo_indx)

        # This is the component redistribution used by Diffsky's
        # _get_dbk_sed_kern_merging_quantities.
        scale = total_merged / total_in_situ
        return FastSEDResult(
            sed_bulge * scale, sed_disk * scale, sed_knots * scale)

    def managed(
        t_obs, mc_sfh_type, fknot, uran_fbulge, uran_av, uran_delta,
        uran_funo, uran_pburst, uran_pmerge, logm0, logtc, early_index,
        late_index, logmp_infall, logmhost_infall, central, top_host_idx,
        t_peak, lgmcrit, lgy_at_mcrit, indx_lo, indx_hi, lg_qt, qlglgdt,
        lg_drop, lg_rejuv, delta_mag_ssp_scatter, redshift_true, gal_id,
        ssp_data, param_collection, cosmology, Ob0, index=None,
    ):
        # OpenCosmo may provide the selected row index to vectorized
        # functions. It is not needed by this kernel, but accepting it keeps
        # the callable compatible across OpenCosmo versions and selection
        # sizes.
        del index
        n_input = len(redshift_true)
        if precision == "float32":
            # Casting at this boundary makes the complete calculation, rather
            # than just its stored output, use float32. Float64
            # preserves the input dtypes used by Diffsky's reference path.
            floating = cast_tree((
                t_obs, fknot, uran_fbulge, uran_av, uran_delta, uran_funo,
                uran_pburst, uran_pmerge, logm0, logtc, early_index,
                late_index, logmp_infall, logmhost_infall, t_peak, lgmcrit,
                lgy_at_mcrit, indx_lo, indx_hi, lg_qt, qlglgdt, lg_drop,
                lg_rejuv, delta_mag_ssp_scatter, redshift_true))
            (t_obs, fknot, uran_fbulge, uran_av, uran_delta, uran_funo,
             uran_pburst, uran_pmerge, logm0, logtc, early_index, late_index,
             logmp_infall, logmhost_infall, t_peak, lgmcrit, lgy_at_mcrit,
             indx_lo, indx_hi, lg_qt, qlglgdt, lg_drop, lg_rejuv,
             delta_mag_ssp_scatter, redshift_true) = floating
            ssp_data = cast_tree(ssp_data)
            param_collection = cast_tree(param_collection)
            cosmology = cast_tree(cosmology)
            Ob0 = cast_floating_leaf(Ob0)

        # The number of satellites accompanying a fixed number of central
        # systems varies from call to call. Without padding, JAX compiles and
        # retains a separate executable for every distinct galaxy-array
        # length. Power-of-two buckets trade up to 2x temporary
        # array padding for a logarithmically bounded number of compiled XLA
        # executables. 
        n_padded = _fast_batch_shape(n_input)
        shape_key = (precision, n_padded)
        if shape_key not in _FAST_COMPILED_SHAPES:
            _FAST_COMPILED_SHAPES.add(shape_key)
            _LOGGER.info(
                "Diffsky fast SED compiling JAX batch shape %d (%s)",
                n_padded, precision)
        padding = n_padded - n_input
        if padding:
            def pad_edge(value):
                value = jnp.asarray(value)
                widths = ((0, padding),) + ((0, 0),) * (value.ndim - 1)
                return jnp.pad(value, widths, mode="edge")

            per_galaxy = tuple(pad_edge(value) for value in (
                t_obs, mc_sfh_type, fknot, uran_fbulge, uran_av, uran_delta,
                uran_funo, uran_pburst, uran_pmerge, logm0, logtc,
                early_index, late_index, logmp_infall, logmhost_infall,
                t_peak, lgmcrit, lgy_at_mcrit, indx_lo, indx_hi, lg_qt,
                qlglgdt, lg_drop, lg_rejuv, delta_mag_ssp_scatter,
                redshift_true))
            (t_obs, mc_sfh_type, fknot, uran_fbulge, uran_av, uran_delta,
             uran_funo, uran_pburst, uran_pmerge, logm0, logtc, early_index,
             late_index, logmp_infall, logmhost_infall, t_peak, lgmcrit,
             lgy_at_mcrit, indx_lo, indx_hi, lg_qt, qlglgdt, lg_drop,
             lg_rejuv, delta_mag_ssp_scatter, redshift_true) = per_galaxy
            central = jnp.concatenate((
                jnp.asarray(central),
                jnp.ones(padding, dtype=jnp.asarray(central).dtype)))
            top_host_idx = jnp.concatenate((
                jnp.asarray(top_host_idx),
                jnp.arange(n_input, n_padded,
                           dtype=jnp.asarray(top_host_idx).dtype)))
        mah_params = DiffmahParams(
            logm0=logm0, logtc=logtc, early_index=early_index,
            late_index=late_index, t_peak=t_peak)
        sfh_params = DiffstarParams(
            lgmcrit=lgmcrit, lgy_at_mcrit=lgy_at_mcrit, indx_lo=indx_lo,
            indx_hi=indx_hi, lg_qt=lg_qt, qlglgdt=qlglgdt,
            lg_drop=lg_drop, lg_rejuv=lg_rejuv)
        phot_randoms = mc_randoms.PhotRandoms(
            mc_sfh_type == 0, uran_av, uran_delta, uran_funo,
            uran_pburst, delta_mag_ssp_scatter)
        dbk_randoms = mc_randoms.DBKRandoms(fknot, uran_fbulge)
        merging_randoms = mc_randoms.DiffMergeRandoms(uran_pmerge)
        result = kernel(
            phot_randoms, dbk_randoms, merging_randoms, sfh_params,
            redshift_true, t_obs, mah_params, ssp_data,
            param_collection.mzr_params, param_collection.spspop_params,
            param_collection.scatter_params, param_collection.ssperr_params,
            param_collection.merging_params, cosmology,
            Ob0 / cosmology.Om0, logmp_infall, logmhost_infall, t_peak,
            central, jnp.ones(len(redshift_true)), top_host_idx, 1)
        output_dtype = jnp.float32 if precision == "float32" else None
        return {
            "rest_sed_bulge": jnp.asarray(
                result.rest_sed_bulge[:n_input], dtype=output_dtype),
            "rest_sed_disk": jnp.asarray(
                result.rest_sed_disk[:n_input], dtype=output_dtype),
            "rest_sed_knots": jnp.asarray(
                result.rest_sed_knots[:n_input], dtype=output_dtype),
            "gal_id": gal_id,
        }

    return managed


_FAST_MANAGED = {}


def clear_fast_sed_caches():
    """Release optimized-kernel functions and JAX compilation caches."""
    _FAST_MANAGED.clear()
    _FAST_COMPILED_SHAPES.clear()
    try:
        import jax
    except ModuleNotFoundError:
        return
    jax.clear_caches()


def fast_sed_compiled_shapes():
    """Return the precision and padded size of encountered fast JAX shapes."""
    return tuple(sorted(_FAST_COMPILED_SHAPES))


def compute_dbk_seds_fast(catalog, aux_data, batch_size=25,
                          precision="float32"):
    """Compute only the component SEDs required by SkyCatalogs.

    Parameters match Diffsky's ``compute_dbk_seds_from_diffsky_mock`` except
    that this function always returns a dictionary and never inserts columns
    into the OpenCosmo catalog.

    Notes
    -----
    ``batch_size`` retains Diffsky's meaning: it is the number of complete
    central-plus-satellite systems, not necessarily the number of rows.
    """
    if precision not in ("float32", "float64"):
        raise ValueError("precision must be 'float32' or 'float64'")
    from diffsky.data_loaders.opencosmo_utils import utils

    utils.validate_batch_size(batch_size)
    if precision not in _FAST_MANAGED:
        _FAST_MANAGED[precision] = _make_fast_kernel(precision)
    managed = _FAST_MANAGED[precision]

    cosmology = utils.prep_cosmology_parameters(catalog.cosmology)
    evaluated = catalog.evaluate(
        utils.age_at_z_, vectorize=True, cosmology=cosmology, format="jax")

    input_ids = np.atleast_1d(
        evaluated.select("gal_id").get_data("numpy"))
    if len(input_ids) == 0:
        return {
            "rest_sed_bulge": np.empty((0, 0)),
            "rest_sed_disk": np.empty((0, 0)),
            "rest_sed_knots": np.empty((0, 0)),
        }

    component_names = (
        "rest_sed_bulge", "rest_sed_disk", "rest_sed_knots")
    output_dtype = np.float32 if precision == "float32" else np.float64
    n_wave = len(aux_data["ssp_data"].ssp_wave)
    result = {name: np.empty((len(input_ids), n_wave), dtype=output_dtype)
              for name in component_names}
    input_positions = {int(gal_id): index
                       for index, gal_id in enumerate(input_ids)}

    # OpenCosmo returns a scalar for a one-row selection.  Diffsky's
    # split_central_indices currently passes that scalar directly to
    # np.where, which raises on NumPy 2.x. Keep the same splitting semantics
    # while normalizing the selected column to an array first.
    central = np.atleast_1d(
        evaluated.select("central").get_data("numpy"))
    central_indices = np.flatnonzero(central)
    if batch_size == -1:
        central_batches = [central_indices]
    else:
        split_indices = np.arange(batch_size, len(central_indices), batch_size)
        central_batches = np.split(central_indices, split_indices)

    for central_rows in central_batches:
        output = evaluated.take_rows(central_rows).evaluate(
            managed,
            ssp_data=aux_data["ssp_data"],
            param_collection=aux_data["param_collection"],
            cosmology=cosmology,
            Ob0=evaluated.cosmology.Ob0,
            insert=False,
            vectorize=True,
            format="jax",
        )
        output_ids = np.atleast_1d(np.asarray(output["gal_id"]))
        destination = np.asarray(
            [input_positions[int(gal_id)] for gal_id in output_ids])
        for name in component_names:
            values = np.asarray(output[name], dtype=output_dtype)
            if values.ndim == 1:
                values = values[np.newaxis, :]
            result[name][destination] = values
        del output

    return result
