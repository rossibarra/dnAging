import numpy as np

import posterior_sample_age_infer as inf
import precompute_freq_trajectory_moments as pre


def test_large_tau_moments_fail_loudly():
    e1, e2 = pre.MomentEngine(26).Emoments(8, 10.0, 0.0, 1.0 / 20000)
    assert np.isnan(e1)
    assert np.isnan(e2)


def test_tau_three_remains_available_at_typical_panel_size():
    e1, e2 = pre.MomentEngine(26).Emoments(8, 3.0, 1.0, 1.0 / 20000)
    assert np.isfinite(e1)
    assert np.isfinite(e2)


def test_phi_lookup_preserves_nan():
    tab = {
        "table": np.full((1, 2, 2), np.nan),
        "d0": np.array([1]),
        "age": np.array([100.0, 1000.0]),
        "Tgrid": np.array([0.0, 10.0]),
        "n_sample": 2,
    }
    assert np.isnan(inf.phi_lookup(tab, 1, 200.0, 200.0)).all()


def test_mutation_age_max_validation():
    try:
        inf.parse_args(["--freq-table", "x", "--output", "y", "--merge", "z",
                        "--mutation-age-max", "0"])
    except SystemExit as exc:
        assert exc.code == 2
    else:
        raise AssertionError("non-positive cutoff was accepted")


def test_mutation_age_max_defaults_to_tau_three():
    args = inf.parse_args(["--freq-table", "x", "--output", "y", "--merge", "z"])
    assert args.mutation_age_max == 3.0
