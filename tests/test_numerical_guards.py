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


def test_epsilon_defaults_to_one_percent_and_is_validated():
    base = ["--freq-table", "x", "--output", "y", "--merge", "z"]
    assert inf.parse_args(base).epsilon == 0.01
    for bad in ("-0.01", "0.5", "1"):
        try:
            inf.parse_args(base + ["--epsilon", bad])
        except SystemExit as exc:
            assert exc.code == 2
        else:
            raise AssertionError(f"invalid epsilon {bad} was accepted")


def test_multiple_mapping_detection():
    assert not inf._is_multiply_mapped(np.array([0, 1, 2]))
    assert inf._is_multiply_mapped(np.array([0, 1, 1, 2]))


def test_prior_validation_accepts_valid_zero_density(tmp_path):
    prior = tmp_path / "prior.tsv"
    np.savetxt(prior, [[0, 0], [10, 1], [20, 0]])
    args = type("Args", (), {"prior_file": prior})()
    got = inf.load_prior(args, np.array([0.0, 10.0, 20.0]))
    assert np.isfinite(got).all()
    assert got[1] == 0.0


def test_prior_validation_rejects_malformed_input(tmp_path):
    bad_arrays = (
        [[0, 1]],                         # fewer than two rows
        [[0, 1, 2], [1, 1, 2]],          # wrong number of columns
        [[0, 1], [np.nan, 1]],           # non-finite
        [[10, 1], [0, 1]],               # unsorted ages
        [[0, 1], [0, 2]],                # duplicate ages
        [[0, 1], [10, -1]],              # negative density
        [[0, 0], [10, 0]],               # no positive mass
    )
    args = type("Args", (), {})()
    for i, values in enumerate(bad_arrays):
        prior = tmp_path / f"bad-{i}.tsv"
        np.savetxt(prior, values)
        args.prior_file = prior
        try:
            inf.load_prior(args, np.array([0.0, 10.0]))
        except SystemExit:
            pass
        else:
            raise AssertionError(f"malformed prior {i} was accepted")


def test_variable_panel_size_lookup_selects_exact_n():
    age = np.array([100.0, 1000.0])
    # n_panel=20 returns 0.2; n_panel=21 returns 0.8 for every d0/T.
    table = np.empty((2, 21, 2, 2), dtype=float)
    table[0].fill(0.2); table[1].fill(0.8)
    tab = {"table": table, "d0": np.arange(1, 22), "age": age,
           "Tgrid": np.array([0.0, 10.0]), "n_sample": 21,
           "n_panel": np.array([20, 21])}
    assert np.allclose(inf.phi_lookup(tab, 5, 200, 200, n_called=20), 0.2)
    assert np.allclose(inf.phi_lookup(tab, 5, 200, 200, n_called=21), 0.8)


def test_min_n_defaults_to_twenty():
    args = inf.parse_args(["--freq-table", "x", "--output", "y", "--merge", "z"])
    assert args.min_n == 20


def _write_merge_part(path, samples, grid, ll):
    path.mkdir()
    (path / "samples.txt").write_text("\n".join(samples) + "\n")
    np.save(path / "grid.npy", grid)
    np.save(path / "ll_marginal.npy", ll)


def test_merge_validates_grid_and_likelihood_shape(tmp_path):
    grid = np.array([0.0, 10.0, 20.0])
    tab = {"Tgrid": grid}
    good = tmp_path / "good"
    _write_merge_part(good, ["a", "b"], grid, np.ones((2, 3)))
    args = type("Args", (), {"merge": [good]})()
    order, got_grid, ll, _stats = inf.merge(args, tab)
    assert order == ["a", "b"]
    assert np.array_equal(got_grid, grid)
    assert np.array_equal(ll, np.ones((2, 3)))

    bad_grid = tmp_path / "bad-grid"
    _write_merge_part(bad_grid, ["a", "b"], np.array([0.0, 11.0, 20.0]),
                      np.ones((2, 3)))
    args.merge = [good, bad_grid]
    try:
        inf.merge(args, tab)
    except SystemExit as exc:
        assert "grid differs" in str(exc)
    else:
        raise AssertionError("mismatched merge grid was accepted")

    bad_shape = tmp_path / "bad-shape"
    _write_merge_part(bad_shape, ["a", "b"], grid, np.ones((1, 3)))
    args.merge = [good, bad_shape]
    try:
        inf.merge(args, tab)
    except SystemExit as exc:
        assert "likelihood shape" in str(exc)
    else:
        raise AssertionError("broadcastable likelihood shape was accepted")
