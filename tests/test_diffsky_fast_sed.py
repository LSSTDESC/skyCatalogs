"""Optional integration test for the optimized Diffsky component kernel."""

import os

import numpy as np
import pytest


@pytest.mark.parametrize(
    "n_input, expected",
    [(1, 32), (31, 32), (32, 32), (33, 64), (64, 64),
     (65, 128), (513, 1024)],
)
def test_fast_batch_shape_uses_power_of_two_buckets(n_input, expected):
    from skycatalogs.utils.diffsky_sed import _fast_batch_shape

    assert _fast_batch_shape(n_input) == expected


def test_fast_batch_shape_rejects_empty_input():
    from skycatalogs.utils.diffsky_sed import _fast_batch_shape

    with pytest.raises(ValueError, match="positive"):
        _fast_batch_shape(0)


@pytest.mark.skipif(
    not os.getenv("DIFFSKY_TEST_CATALOG"),
    reason="DIFFSKY_TEST_CATALOG is not configured",
)
def test_fast_component_seds_match_diffsky_reference():
    pytest.importorskip("diffsky")
    pytest.importorskip("opencosmo")
    from diffsky.data_loaders.opencosmo_utils import (
        compute_dbk_seds_from_diffsky_mock,
        load_diffsky_mock,
    )

    from skycatalogs.utils.diffsky_sed import compute_dbk_seds_fast

    catalog, aux_data = load_diffsky_mock(
        os.environ["DIFFSKY_TEST_CATALOG"])
    central = np.flatnonzero(
        catalog.select("central").get_data("numpy"))[:3]
    sample = catalog.take_rows(central)

    # Keep the integration test quick while exercising all physical kernels.
    ssp = aux_data["ssp_data"]
    keep = (ssp.ssp_wave >= 1450.0) & (ssp.ssp_wave <= 1900.0)
    test_aux = dict(aux_data)
    test_aux["ssp_data"] = ssp._replace(
        ssp_wave=ssp.ssp_wave[keep],
        ssp_flux=ssp.ssp_flux[..., keep],
    )

    reference = compute_dbk_seds_from_diffsky_mock(
        sample, test_aux, insert=False, batch_size=-1)
    reference_batched = compute_dbk_seds_from_diffsky_mock(
        sample, test_aux, insert=False, batch_size=2)
    optimized = compute_dbk_seds_fast(
        sample, test_aux, batch_size=-1, precision="float64")
    optimized_default = compute_dbk_seds_fast(
        sample, test_aux, batch_size=-1)
    optimized_batched = compute_dbk_seds_fast(
        sample, test_aux, batch_size=2, precision="float64")
    single = catalog.take_rows(central[:1])
    reference_single = compute_dbk_seds_from_diffsky_mock(
        single, test_aux, insert=False, batch_size=-1)
    optimized_single = compute_dbk_seds_fast(
        single, test_aux, batch_size=2, precision="float64")

    for component in (
        "rest_sed_bulge", "rest_sed_disk", "rest_sed_knots"
    ):
        np.testing.assert_allclose(
            optimized[component], reference[component],
            rtol=1.0e-11, atol=0.0, equal_nan=True)
        np.testing.assert_allclose(
            optimized_batched[component], reference_batched[component],
            rtol=1.0e-11, atol=0.0, equal_nan=True)
        np.testing.assert_allclose(
            optimized_single[component], reference_single[component],
            rtol=1.0e-11, atol=0.0, equal_nan=True)
        assert optimized_default[component].dtype == np.float32
