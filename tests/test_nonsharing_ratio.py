import numpy as np
import pytest
import tskit

from nonsharing_ratio import (
    cumulative_hazard,
    estimate_H,
    invert_hazard,
    lineage_counts,
    nonsharing_fraction,
    rate_from_lineage_counts,
    select_sites,
)


def test_select_sites_cutoff_and_panel_fixed():
    # 4 panel haplotypes, 5 sites.
    panel = np.array([
        [1, 1, 1, 0, 1],
        [1, 1, 0, 0, 1],
        [1, 0, 0, 0, 1],
        [1, 0, 0, 0, 1],
    ])
    targets = np.array([[0, 1, 1, 1, 1], [1, 1, 0, 0, 1]])
    ages = np.array([np.nan, 500.0, 50.0, np.nan, np.nan])
    mask = select_sites(panel, targets, ages, a0=100.0)
    # site0 panel-fixed (target 0 ancestral -> kept); site1 old; site2 young;
    # site3 absent from panel; site4 fixed in panel AND in every target -> dropped.
    assert mask.tolist() == [True, True, False, False, False]
    no_fixed = select_sites(panel, targets, ages, a0=100.0, include_panel_fixed=False)
    assert no_fixed.tolist() == [False, True, False, False, False]


def test_nonsharing_fraction_counts_pairs():
    panel = np.array([[1, 1], [1, 0], [0, 0]])  # site0 carried by 2 refs, site1 by 1
    targets = np.array([[1, 0], [0, 0], [1, 1]])
    mask = np.array([True, True])
    f, n = nonsharing_fraction(panel, targets, mask)
    assert np.all(n == 3)
    # target0 ancestral at site1 -> 1 pair of 3; target1 ancestral at both -> 3/3; target2 -> 0
    assert np.allclose(f, [1 / 3, 1.0, 0.0])


def test_missing_calls_shrink_denominator_per_target():
    panel = np.array([[1, 1], [1, 0], [0, 0]])
    targets = np.array([[0, 0], [1, 1]])
    called = np.array([[True, False], [True, True]])  # target0 missing at site1
    mask = np.array([True, True])
    f, n = nonsharing_fraction(panel, targets, mask, called)
    assert np.allclose(n, [2, 3]) and np.allclose(f, [1.0, 0.0])


def test_block_counts_sum_to_total():
    from nonsharing_ratio import nonsharing_counts
    rng = np.random.default_rng(1)
    panel = rng.integers(0, 2, (6, 40)); targets = rng.integers(0, 2, (3, 40))
    mask = rng.random(40) < 0.8
    blocks = np.repeat(np.arange(4), 10)
    num, den = nonsharing_counts(panel, targets, mask, block_ids=blocks, n_blocks=4)
    n1, d1 = nonsharing_counts(panel, targets, mask)
    assert np.allclose(num.sum(0), n1[0]) and np.allclose(den.sum(0), d1[0])


def test_nonsharing_fraction_rejects_empty_selection():
    with pytest.raises(ValueError):
        nonsharing_fraction(np.zeros((2, 3)), np.zeros((1, 3)), np.ones(3, bool))


def test_estimate_H_is_log_ratio_to_mean_control():
    assert np.isclose(estimate_H(0.3, [0.1, 0.3]), np.log(0.3 / 0.2))


def test_hazard_roundtrip_piecewise():
    edges = [0, 100, 300, np.inf]
    rates = [1e-3, 2e-4, 5e-5]
    cum = cumulative_hazard(edges, rates)
    assert np.allclose(cum[:3], [0, 0.1, 0.14])
    for T in [0.0, 40.0, 100.0, 250.0, 1000.0]:
        H = np.interp(T, edges[:-1], cum[:-1]) if T <= 300 else 0.14 + 5e-5 * (T - 300)
        if T < 100:
            H = 1e-3 * T
        elif T < 300:
            H = 0.1 + 2e-4 * (T - 100)
        assert np.isclose(invert_hazard(H, edges, rates)[0], T)
    assert invert_hazard(-0.1, edges, rates)[0] == 0.0


def test_constant_rate_inverse():
    c = 1 / 40000
    T = invert_hazard(np.log(1.1), [0, np.inf], [c])[0]
    assert np.isclose(T, np.log(1.1) / c)


def test_rate_from_lineage_counts_recovers_constant_c():
    # Deterministic coalescent: dk/dt = -C(k,2) c.
    c = 2e-5
    t = np.linspace(0, 4000, 41)
    k = np.empty_like(t)
    k[0] = 26.0
    for i in range(1, t.size):
        dt = t[i] - t[i - 1]
        # fine sub-stepping of the ODE
        kk = k[i - 1]
        for _ in range(200):
            kk -= (kk * (kk - 1) / 2) * c * dt / 200
        k[i] = kk
    est = rate_from_lineage_counts(t, k)
    assert np.allclose(est, c, rtol=0.02)


def test_lineage_counts_two_tip_tree():
    tables = tskit.TableCollection(sequence_length=10)
    a = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    b = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    r1 = tables.nodes.add_row(time=10)
    r2 = tables.nodes.add_row(time=20)
    tables.edges.add_row(0, 5, r1, a)
    tables.edges.add_row(0, 5, r1, b)
    tables.edges.add_row(5, 10, r2, a)
    tables.edges.add_row(5, 10, r2, b)
    tables.sort()
    ts = tables.tree_sequence()
    assert np.allclose(lineage_counts(ts, [0, 5, 15, 25]), [2, 2, 1, 0])


def test_edge_child_ages_uses_child_node_time():
    from nonsharing_ratio import edge_child_ages, edge_midpoint_ages
    tables = tskit.TableCollection(sequence_length=10)
    a = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    b = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    r = tables.nodes.add_row(time=1000)
    tables.edges.add_row(0, 10, r, a)
    tables.edges.add_row(0, 10, r, b)
    s = tables.sites.add_row(position=3, ancestral_state="A")
    tables.mutations.add_row(site=s, node=a, derived_state="G")
    tables.sort()
    ts = tables.tree_sequence()
    assert edge_child_ages(ts)[0] == 0.0
    assert edge_midpoint_ages(ts)[0] == 500.0
