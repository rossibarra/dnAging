# -*- coding: utf-8 -*-
"""The scientific core: the neutral WF moment recursion and its two moments."""
from __future__ import annotations

import numpy as np
import pytest
from scipy.linalg import expm

import precompute_freq_trajectory_moments as pre
from _reference import Emoments_ref

EPS = 1.0 / 2000.0          # single new copy in a diploid Ne = 1000


# ---------------------------------------------------------------------------
# generator structure
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n", [1, 4, 26])
def test_generator_is_the_closed_moment_hierarchy(n):
    """B must encode dM_k/dtau = k(k-1)/2 (M_{k-1}-M_k) on n+3 rows.

    Two invariants downstream depend on the shape: the SECOND moment needs one
    power beyond K = n+1 (E[X_T^2 X^m] reads M(u1)[j+2]), and row 0/row 1 must be
    identically zero so M_0 == 1 and M_1 is a martingale.
    """
    eng = pre.MomentEngine(n)
    assert eng.K == n + 1
    assert eng.B.shape == (n + 3, n + 3)
    for k in range(n + 3):
        lam = k * (k - 1) / 2.0
        assert eng.B[k, k] == -lam
        if k:
            assert eng.B[k, k - 1] == lam
    assert not eng.B[0].any() and not eng.B[1].any()
    assert np.allclose(np.triu(eng.B, 1), 0.0)   # lower-bidiagonal


def test_conditioning_coefficients_are_the_binomial_expansion():
    """coeff[d0][m] = C(n-d0, m-d0) (-1)^(m-d0): the alternating expansion of
    Binom(d0; n, x) that turns conditioning on a count into moments."""
    n = 7
    eng = pre.MomentEngine(n)
    assert set(eng.coeff) == set(range(1, n + 1))
    for d0 in range(1, n + 1):
        assert set(eng.coeff[d0]) == set(range(d0, n + 1))
        # sum_m coeff * x^m must equal (1-x)^(n-d0) x^d0 -- the pmf up to C(n,d0)
        for x in (0.1, 0.37, 0.9):
            got = sum(c * x ** m for m, c in eng.coeff[d0].items())
            assert got == pytest.approx(x ** d0 * (1 - x) ** (n - d0), rel=1e-12)


# ---------------------------------------------------------------------------
# Kimura's constant-Ne limit, in closed form
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("u", [0.0, 0.13, 1.0, 4.0, 20.0])
def test_kimura_constant_ne_moments_closed_form(u):
    """Under constant Ne the first three moments of the neutral diffusion are
    analytic (Kimura 1955); _moms must reproduce them exactly.

    M_0 = 1, M_1 = eps (martingale), M_2 = eps + (eps^2-eps)e^-u,
    M_3 = eps + 1.5(eps^2-eps)e^-u + (eps^3-1.5eps^2+0.5eps)e^-3u.
    A wrong sign or factor in B shows up here immediately.
    """
    eps = 0.3          # not 1/(2Ne): the closed form holds for any start point
    m = pre.MomentEngine(6)._moms(u, eps)
    e = np.exp(-u)
    assert m[0] == pytest.approx(1.0, abs=1e-14)
    assert m[1] == pytest.approx(eps, abs=1e-14)
    assert m[2] == pytest.approx(eps + (eps ** 2 - eps) * e, abs=1e-13)
    assert m[3] == pytest.approx(eps + 1.5 * (eps ** 2 - eps) * e
                                 + (eps ** 3 - 1.5 * eps ** 2 + 0.5 * eps) * e ** 3,
                                 abs=1e-13)


def test_moments_are_a_decreasing_valid_moment_sequence():
    """X in [0,1] forces M_0 >= M_1 >= ... >= 0 at every tau; a broken generator
    (or a transposed B) violates it."""
    for u in (0.05, 0.5, 5.0):
        m = pre.MomentEngine(26)._moms(u, EPS)
        assert np.all(m >= -1e-15)
        assert np.all(np.diff(m) <= 1e-15)


# ---------------------------------------------------------------------------
# the mutation-existence boundary
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("tau_T,tau_i", [(1.0, 1.0), (1.5, 1.0), (0.0, 0.0), (9.0, 0.1)])
def test_sample_older_than_mutation_is_exactly_zero(tau_T, tau_i):
    """T >= t_i: the allele does not exist yet, so BOTH moments are exactly 0.0
    -- not approximately, and not NaN. Downstream code relies on the hard zero."""
    eng = pre.MomentEngine(26)
    got = eng.Emoments(5, tau_i, tau_T, EPS)
    assert got == (0.0, 0.0)
    assert eng.Efreq(5, tau_i, tau_T, EPS) == 0.0


def test_efreq_is_the_first_component_of_emoments():
    eng = pre.MomentEngine(12)
    for d0 in (1, 6, 11):
        assert eng.Efreq(d0, 0.7, 0.2, EPS) == eng.Emoments(d0, 0.7, 0.2, EPS)[0]


# ---------------------------------------------------------------------------
# agreement with the high-precision reference
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n", [4, 8])
@pytest.mark.parametrize("tau_i,tau_T", [(0.3, 0.1), (1.0, 0.4), (3.0, 1.2)])
def test_matches_mpmath_reference_small_n(n, tau_i, tau_T):
    """Both moments must match an 80-digit closed-form reference at small n,
    where float64 cancellation is negligible.

    The reference shares no code with the engine: it builds expm(B*u) from the
    exact exponential-sum solution of the recursion, not from scipy.
    """
    eng = pre.MomentEngine(n)
    for d0 in range(1, n + 1):
        e1, e2 = eng.Emoments(d0, tau_i, tau_T, EPS)
        r1, r2 = (float(v) for v in Emoments_ref(n, d0, tau_i, tau_T, EPS))
        assert e1 == pytest.approx(r1, rel=1e-9)
        assert e2 == pytest.approx(r2, rel=1e-9)


@pytest.mark.parametrize("n", [4, 8])
def test_second_moment_index_shift_would_be_caught(n):
    """Guards the guard: the second-moment contraction is C[:, :K] @ Mu1[2:K+2].
    Re-running it with Mu1[1:K+1] (the FIRST-moment slice) or Mu1[3:K+3] must
    disagree with the reference by orders of magnitude more than the 1e-9
    tolerance above -- otherwise that test could not detect fix 5 regressing."""
    eng = pre.MomentEngine(n)
    tau_i, tau_T, K = 1.0, 0.4, n + 1
    Mu1 = eng._moms(tau_i - tau_T, EPS)
    C = expm(eng.B * tau_T)
    Mpres = eng._moms(tau_i, EPS)
    for shift in (1, 3):
        Ej = C[:, :K] @ Mu1[shift:K + shift]
        for d0 in (1, n // 2 + 1, n):
            num = den = 0.0
            for m, c in eng.coeff[d0].items():
                num += c * Ej[m]
                den += c * Mpres[m]
            ref = float(Emoments_ref(n, d0, tau_i, tau_T, EPS)[1])
            assert abs(num / den - ref) / ref > 1e-3


def test_cauchy_schwarz_holds_before_the_clamp():
    """E[p]^2 <= E[p^2] <= E[p] must hold as MATHEMATICS, not because Emoments
    clamps to those bounds.

    Checked on the unclamped high-precision reference, with a strict interior
    margin, and then the engine is shown to agree with it -- so the clamp is
    inactive here and the inequality is a real property of the recursion.
    """
    n = 26
    for tau_i, tau_T in ((0.3, 0.1), (1.0, 0.4), (3.0, 1.0)):
        for d0 in (1, 8, 13, 25):
            r1, r2 = (float(v) for v in Emoments_ref(n, d0, tau_i, tau_T, EPS))
            assert r1 * r1 < r2 < r1              # strictly interior: no clamping
            assert r2 - r1 * r1 > 1e-3 and r1 - r2 > 1e-3
            e1, e2 = pre.MomentEngine(n).Emoments(d0, tau_i, tau_T, EPS)
            assert e1 == pytest.approx(r1, rel=1e-2)
            assert e2 == pytest.approx(r2, rel=1e-2)


def test_conditional_moments_increase_with_the_present_count():
    """A larger present-day count d0 implies a larger expected past frequency;
    monotonicity in d0 is the qualitative signature of the conditioning being
    applied in the right direction (a sign error in coeff inverts it)."""
    eng = pre.MomentEngine(26)
    vals = [eng.Efreq(d0, 1.0, 0.4, EPS) for d0 in range(1, 27)]
    assert np.all(np.diff(vals) > 0)


def test_frequency_decays_toward_the_origin():
    """Walking T back toward the mutation age drives E[p_T] to eps: at tau_T ->
    tau_i the trajectory is pinned at its single-copy start."""
    eng = pre.MomentEngine(26)
    tau_i = 1.0
    vals = [eng.Efreq(13, tau_i, tt, EPS) for tt in (0.999, 0.99, 0.9, 0.5, 0.1, 0.0)]
    assert np.all(np.diff(vals) > 0)             # older sample age -> lower freq
    assert vals[0] == pytest.approx(EPS, rel=5e-2)


# ---------------------------------------------------------------------------
# forward Wright-Fisher Monte Carlo (validate_moments_vs_mc.py's approach)
# ---------------------------------------------------------------------------


def _wf_mc(n, d0_target, tau_i, tau_T, nc, reps, seed=20260903, dtf=0.005):
    """Forward WF-diffusion MC, conditioned on the panel count, as in
    validate_moments_vs_mc.py but small enough for a unit test."""
    rng = np.random.default_rng(seed)
    dt = dtf * nc
    p = np.full(reps, 1.0 / nc)
    t = tau_i * nc
    pT = None
    while t > 0:
        step = min(dt, t)
        act = (p > 0) & (p < 1)
        p[act] += np.sqrt(np.maximum(p[act] * (1 - p[act]) / nc, 0) * step) \
            * rng.standard_normal(int(act.sum()))
        p = np.clip(p, 0, 1)
        t -= step
        if pT is None and t <= tau_T * nc:
            pT = p.copy()
    surv = p > 0
    hit = rng.binomial(n, p[surv]) == d0_target
    v = pT[surv][hit]
    assert v.size > 200, "MC bin too thin to test against"
    return v.mean(), (v ** 2).mean(), v.std(ddof=1) / np.sqrt(v.size), v.size


@pytest.mark.parametrize("tau_i,tau_T", [(0.6, 0.25), (1.5, 0.6)])
def test_agrees_with_forward_wf_monte_carlo(tau_i, tau_T):
    """Independent check that the recursion is the right model, not just the
    right arithmetic: a forward WF diffusion, conditioned on the same panel
    count, must reproduce both moments.

    Tolerance = 4 * MC standard error + 0.015 absolute. The additive term covers
    the Euler-Maruyama discretisation bias of the simulator (dt = 0.005 * 2Ne);
    it shrinks with dt, but a finer step makes the test slow for no extra power.
    A wrong index in the second-moment contraction shifts E[p^2] by ~0.05-0.15
    here, well outside the band -- the exactness test above is the tight one.
    """
    n, nc = 8, 200.0
    eng = pre.MomentEngine(n)
    for d0 in (1, 2, 4):
        m1, m2, se, cnt = _wf_mc(n, d0, tau_i, tau_T, nc, 300_000)
        a1, a2 = eng.Emoments(d0, tau_i, tau_T, 1.0 / nc)
        tol = 4 * se + 0.015
        assert abs(a1 - m1) < tol, f"d0={d0} E[p] {a1} vs MC {m1} (n={cnt}, se={se})"
        assert abs(a2 - m2) < tol, f"d0={d0} E[p^2] {a2} vs MC {m2} (n={cnt}, se={se})"


def test_plug_in_squared_mean_is_materially_wrong():
    """The reason fix 5 exists: E[p]^2 is not a usable stand-in for E[p^2], so a
    diploid likelihood built from the first moment alone is biased, not merely
    imprecise. Pins a >30% relative gap at a representative table entry."""
    e1, e2 = pre.MomentEngine(26).Emoments(8, 1.0, 0.4, EPS)
    assert (e2 - e1 * e1) / e2 > 0.30
