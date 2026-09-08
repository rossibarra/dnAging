# -*- coding: utf-8 -*-
"""p_T = 0 whenever the sample predates every mutation age the branch admits.

That is the one place the answer is certain rather than estimated, so it must
survive the table's finite age coverage and its log-age interpolation.
"""
import numpy as np
import pytest

import posterior_sample_age_infer as inf


def _tab(age, Tgrid, value=0.2, n_sample=3, n_d0=2):
    return {"table": np.full((n_d0, len(age), len(Tgrid)), value, np.float32),
            "d0": np.arange(1, n_d0 + 1), "age": np.asarray(age, float),
            "Tgrid": np.asarray(Tgrid, float), "n_sample": n_sample}


def test_branch_entirely_below_table_coverage_gives_zero():
    """Edge [0, 40] with a sample at T=50: the mutation cannot yet exist.

    The clamped fallback used to mask on age[0]=100 rather than the branch top,
    leaving T=50 unmasked and returning the raw table value 0.2.
    """
    got = inf.phi_lookup(_tab([100., 1000., 10000.], [50., 500.]), 1, 0.0, 40.0)
    assert got.tolist() == [0.0, 0.0]


def test_branch_above_the_sample_still_contributes():
    """The complementary case must not be zeroed: the edge is older than T."""
    got = inf.phi_lookup(_tab([100., 1000., 10000.], [10., 50., 500.]),
                         1, 200.0, 400.0)
    assert got[0] > 0 and got[1] > 0
    assert got[2] == 0.0                      # T=500 is above the branch top


@pytest.mark.parametrize("t_lo,t_hi", [(0.0, 40.0), (40.0, 0.0), (0.0, 0.0)])
def test_zero_holds_regardless_of_endpoint_order(t_lo, t_hi):
    got = inf.phi_lookup(_tab([100., 1000.], [50.]), 1, t_lo, t_hi)
    assert got.tolist() == [0.0]


def test_exact_age_knot_ignores_a_nan_in_the_zero_weight_neighbour():
    """0 * NaN = NaN would reject a site over a value it does not depend on."""
    tab = {"table": np.array([[[0.3, 0.3], [np.nan, np.nan], [0.5, 0.5]]], np.float32),
           "d0": np.array([1]), "age": np.array([100., 200., 1000.]),
           "Tgrid": np.array([10., 20.]), "n_sample": 2}
    got = inf.phi_lookup(tab, 1, 100.0, 100.0)
    assert np.isfinite(got).all()
    assert got.tolist() == [pytest.approx(0.3), pytest.approx(0.3)]


def test_nan_that_the_answer_does_depend_on_is_still_propagated():
    """The fix must not turn into a blanket NaN suppression."""
    tab = {"table": np.array([[[0.3, 0.3], [np.nan, np.nan], [0.5, 0.5]]], np.float32),
           "d0": np.array([1]), "age": np.array([100., 200., 1000.]),
           "Tgrid": np.array([10., 20.]), "n_sample": 2}
    got = inf.phi_lookup(tab, 1, 150.0, 150.0)      # interpolates INTO the NaN row
    assert not np.isfinite(got).all()


def test_resolver_eligibility_is_honoured():
    """Resolved-but-ineligible rows must not be returned as usable.

    PositionResolution.row_indices is request-aligned and uses -1 only for an
    UNRESOLVED coordinate, so an ineligible row still carries a real-looking
    index; eligible_mask is what decides inclusion. Reading it as `eligible`
    returned None -- no such attribute -- and skipped the filter entirely.
    """
    import sys, types
    from types import SimpleNamespace

    res = SimpleNamespace(row_indices=np.array([0, 1, 2]),
                          eligible_mask=np.array([True, False, True]),
                          resolved_mask=np.array([True, True, True]))
    pkg = types.ModuleType("normalize_tes")
    mod = types.ModuleType("normalize_tes.snp_position_resolution")
    mod.resolve_native_position_requests = lambda *a, **k: res
    pkg.snp_position_resolution = mod
    saved = {k: sys.modules.get(k) for k in
             ("normalize_tes", "normalize_tes.snp_position_resolution")}
    sys.modules["normalize_tes"] = pkg
    sys.modules["normalize_tes.snp_position_resolution"] = mod
    try:
        rows = inf._resolve_rows(object(), ["1", "1", "1"], [10, 20, 30])
    finally:
        for k, v in saved.items():
            if v is None:
                sys.modules.pop(k, None)
            else:
                sys.modules[k] = v
    assert rows.tolist() == [0, -1, 2]


def test_resolver_without_an_eligibility_mask_is_passed_through():
    """Older resolvers expose no mask; row_indices is then already the answer."""
    import sys, types
    from types import SimpleNamespace

    res = SimpleNamespace(row_indices=np.array([4, -1, 6]))
    pkg = types.ModuleType("normalize_tes")
    mod = types.ModuleType("normalize_tes.snp_position_resolution")
    mod.resolve_native_position_requests = lambda *a, **k: res
    pkg.snp_position_resolution = mod
    saved = {k: sys.modules.get(k) for k in
             ("normalize_tes", "normalize_tes.snp_position_resolution")}
    sys.modules["normalize_tes"] = pkg
    sys.modules["normalize_tes.snp_position_resolution"] = mod
    try:
        rows = inf._resolve_rows(object(), ["1", "1", "1"], [10, 20, 30])
    finally:
        for k, v in saved.items():
            if v is None:
                sys.modules.pop(k, None)
            else:
                sys.modules[k] = v
    assert rows.tolist() == [4, -1, 6]
