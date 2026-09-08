# -*- coding: utf-8 -*-
"""ExactMomentEngine: the conditioning sum evaluated to enough digits.

The alternating sum sum_m binom(n-d0,m-d0)(-1)^{m-d0} M_m loses roughly 0.3*n
decimal digits, so float64 cannot express it for realistic panels. These tests
check the exact engine against the independent high-precision reference in
_reference.py (which shares no code with it), and pin the two failure modes of
the float64 path it replaces: outright NaN, and silent inaccuracy.
"""
import numpy as np
import pytest

import precompute_freq_trajectory_moments as pre
from _reference import Emoments_ref

EPS = 1.0 / 20000
POINTS = [(3.0, 1.0), (10.0, 0.0), (0.5, 0.2), (0.05, 0.01), (2.0, 1.999)]


@pytest.mark.parametrize("n", [6, 12, 26, 40])
def test_matches_independent_reference(n):
    eng = pre.ExactMomentEngine(n)
    for tau_i, tau_T in POINTS:
        for d0 in {1, 2, n // 2, n - 1, n}:
            got = eng.Emoments(d0, tau_i, tau_T, EPS)
            ref = [float(v) for v in Emoments_ref(n, d0, tau_i, tau_T, EPS)]
            for g, r in zip(got, ref):
                assert abs(g - r) <= 1e-11 * max(abs(r), 1e-12), (n, d0, tau_i, tau_T)


def test_recovers_the_value_float64_returns_as_nan():
    """tau_i=10, tau_T=0 at n=26: float64 gives NaN, the value is 0.32142857..."""
    exact = pre.ExactMomentEngine(26).Emoments(8, 10.0, 0.0, EPS)
    legacy = pre.MomentEngine(26).Emoments(8, 10.0, 0.0, EPS)
    assert np.isnan(legacy[0])
    assert exact[0] == pytest.approx(0.321428571274, abs=1e-11)
    assert exact[1] == pytest.approx(0.110837438320, abs=1e-11)


def test_is_accurate_where_float64_silently_is_not():
    exact = pre.ExactMomentEngine(26).Emoments(8, 3.0, 1.0, EPS)
    legacy = pre.MomentEngine(26).Emoments(8, 3.0, 1.0, EPS)
    ref = float(Emoments_ref(26, 8, 3.0, 1.0, EPS)[0])
    assert abs(exact[0] - ref) < 1e-11
    assert abs(legacy[0] - ref) > 1e-3          # float64 is off in the third digit


@pytest.mark.parametrize("n", [12, 26])
def test_grid_agrees_with_scalar_calls(n):
    eng = pre.ExactMomentEngine(n)
    tauT = np.linspace(0.0, 1.5, 9)
    p1, p2 = eng.grid(2.0, tauT, EPS)
    for d0 in {1, n // 2, n}:
        for it, tT in enumerate(tauT):
            a, b = eng.Emoments(d0, 2.0, float(tT), EPS)
            assert p1[d0 - 1, it] == pytest.approx(a, abs=1e-12)
            assert p2[d0 - 1, it] == pytest.approx(b, abs=1e-12)


def test_grid_is_free_of_nan_where_the_mutation_predates_the_sample():
    """Every entry with tau_T < tau_i is computable; none may come back NaN."""
    n = 40
    tauT = np.linspace(0.0, 0.9, 12)
    p1, p2 = pre.ExactMomentEngine(n).grid(1.0, tauT, EPS)
    assert np.isfinite(p1).all() and np.isfinite(p2).all()
    # the same block in float64 loses a large share of its entries
    legacy = pre.MomentEngine(n)
    lost = sum(not np.isfinite(legacy.Emoments(d0, 1.0, float(tT), EPS)[0])
               for d0 in range(1, n + 1) for tT in tauT)
    assert lost > 0.2 * n * tauT.size


@pytest.mark.parametrize("n", [12, 26])
def test_moment_constraints_hold_without_clipping(n):
    """0 <= E[p] <= 1 and E[p]^2 <= E[p^2] <= E[p] as mathematics, not by clamping."""
    eng = pre.ExactMomentEngine(n)
    p1, p2 = eng.grid(1.3, np.linspace(0.0, 1.2, 10), EPS)
    assert (p1 >= 0).all() and (p1 <= 1).all()
    assert (p2 >= p1 * p1 - 1e-12).all()
    assert (p2 <= p1 + 1e-12).all()


def test_sample_older_than_mutation_is_zero():
    eng = pre.ExactMomentEngine(12)
    assert eng.Emoments(5, 1.0, 1.0, EPS) == (0.0, 0.0)
    assert eng.Emoments(5, 1.0, 2.0, EPS) == (0.0, 0.0)
    assert eng.Efreq(5, 1.0, 2.0, EPS) == 0.0


def test_efreq_is_the_first_moment():
    eng = pre.ExactMomentEngine(12)
    for d0 in (1, 6, 12):
        assert eng.Efreq(d0, 0.7, 0.2, EPS) == eng.Emoments(d0, 0.7, 0.2, EPS)[0]


def test_precision_scales_with_panel_size():
    assert pre.ExactMomentEngine(12).dps < pre.ExactMomentEngine(40).dps
    assert pre.ExactMomentEngine(26, dps=80).dps == 80


# Converged values, established by self-convergence at dps = 100/160/260/400/500
# and cross-checked against _reference.py raised to 400 digits. A fixed 30 + n
# budget gets all three of these wrong, so they pin accuracy rather than mere
# plausibility: bounds alone cannot tell a converged value from a plausible one.
CONVERGED = {
    # (n, d0, tau_i): (E[p], E[p^2])
    (26, 21, 0.0005): (0.005836501557915401, 3.5577304407708e-05),
    (26, 13, 233.0): (0.5, 0.25862068965517243),
    (20, 3, 233.0): (0.18181818181818182, 0.03952569169960474),
}


@pytest.mark.parametrize("key", sorted(CONVERGED))
def test_matches_converged_value_where_a_fixed_budget_fails(key):
    n, d0, tau_i = key
    want1, want2 = CONVERGED[key]
    got1, got2 = pre.ExactMomentEngine(n).Emoments(d0, tau_i, 0.0, EPS)
    assert got1 == pytest.approx(want1, rel=1e-11)
    assert got2 == pytest.approx(want2, rel=1e-11)


def test_precision_grows_at_both_ends_of_the_age_range():
    """Young AND old mutation ages lose digits, for different reasons.

    Young: the partial-fraction expansion of e^{B u} cancels to order u^{m-j}.
    Old: every moment tends to the fixation probability, so the alternating
    conditioning sum cancels exactly in the limit.
    """
    eng = pre.ExactMomentEngine(26)
    mid = eng.required_dps(1.0)
    assert eng.required_dps(0.0005) > mid              # young end costs more
    assert eng.required_dps(233.0) > mid               # so does the old end
    assert eng.required_dps(5e-4) > eng.required_dps(5e-3) > eng.required_dps(5e-2)
    assert eng.required_dps(500.0) > eng.required_dps(233.0) > eng.required_dps(50.0)
    # the middle of the range is cheapest, but the 2^n binomial weights still
    # cost ~0.301*n digits everywhere, so it is not the bare 30 + n floor
    assert eng.dps <= mid < eng.required_dps(0.0005)


def test_moment_constraints_are_repaired_by_escalation_not_asserted():
    """A constraint violation means too few digits, so it must trigger a retry.

    A fixed 30+n budget returns E[p^2] = -1.8e-4 here -- negative, and 615% off.
    """
    p1, p2 = pre.ExactMomentEngine(26).Emoments(21, 0.0005, 0.0, EPS)
    assert p2 > 0
    assert p1 * p1 <= p2 <= p1


@pytest.mark.parametrize("tau_i", [5e-4, 5e-3, 0.05, 1.0, 3.0, 10.0, 50.0, 233.0])
def test_whole_age_range_is_self_consistent(tau_i):
    """Sweep the grid: every entry must satisfy the exact identities."""
    eng = pre.ExactMomentEngine(20)
    p1, p2 = eng.grid(tau_i, np.array([0.0, tau_i * 0.5]), EPS)
    assert np.isfinite(p1).all() and np.isfinite(p2).all()
    assert (p1 >= 0).all() and (p1 <= 1).all()
    assert (p2 >= p1 * p1 - 1e-12).all() and (p2 <= p1 + 1e-12).all()


def test_escalated_engines_are_reused():
    """Rebuilding per call discards the coefficient tables and the C cache."""
    pre.ExactMomentEngine._escalated.clear()
    eng = pre.ExactMomentEngine(12)
    eng.grid(1e-4, np.array([0.0]), EPS)
    n_after_first = len(pre.ExactMomentEngine._escalated)
    assert n_after_first >= 1
    eng.grid(1e-4, np.array([0.0]), EPS)
    assert len(pre.ExactMomentEngine._escalated) == n_after_first
    assert pre.ExactMomentEngine.at_precision(12, 90) is \
           pre.ExactMomentEngine.at_precision(12, 90)


def test_explicit_precision_is_a_floor_not_a_cap():
    """Accuracy is not negotiable, so --precision may only raise the budget."""
    eng = pre.ExactMomentEngine(26, dps=80)
    assert eng.dps == 80
    assert eng.required_dps(1.0) == 80              # floor respected
    assert eng.required_dps(0.0005) > 80            # but raised where needed


def test_partial_fraction_expansion_reproduces_the_matrix_exponential():
    """The closed form must equal expm(B u) -- lam_0 = lam_1 = 0 notwithstanding."""
    from scipy.linalg import expm
    for n in (6, 14):
        eng = pre.ExactMomentEngine(n)
        for u in (0.003, 0.4):
            C = eng._C(u)
            ref = expm(pre.MomentEngine(n).B * u)
            for m in range(eng.N):
                for j in range(eng.K):
                    assert float(C[m][j]) == pytest.approx(ref[m, j], abs=1e-12)
