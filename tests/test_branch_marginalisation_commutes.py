# -*- coding: utf-8 -*-
"""T1: does the within-branch age marginalisation commute with the likelihood?

Eq. (10b) averages the FREQUENCY table over the branch and only then forms the
per-site likelihood.  The estimand is the other order -- average the LIKELIHOOD
over the branch:

    prod_i  integral p(t_i | branch) * p(a_i | T, t_i) dt_i

The two agree exactly iff the per-site likelihood is affine in the tabulated
quantities, which it is: qA = eps + (1-2eps)*phi is linear in phi, and the
diploid dosage probabilities are linear in (phi, phi2) jointly.  This file
pins that, in two separate claims:

  (a) COMMUTATION, exact.  Averaging ell over a set of point ages equals ell of
      the averaged phi, when both use the same weights.  This is pure algebra,
      so it is asserted at machine precision; a failure means the production
      likelihood has acquired a term that is not affine in the table (a
      frequency-dependent epsilon, a clip that bites before the average, a
      third moment).

  (b) The production closed-form branch integral equals an independent fine
      trapezoidal average of production POINT-age lookups.  That reference
      treats the mutation-existence boundary differently -- it steps at each
      quadrature node rather than resolving it as an exact integration limit --
      so agreement is only to quadrature tolerance, and it is a genuine check
      of the knot-splitting integrator in phi_lookup rather than a restatement
      of it.

Together (a) and (b) say: marginalising the age into the frequency first, with
this integrator, gives the same site likelihood as marginalising the likelihood.

Scope.  This is about ONE site's own age.  The order that does NOT commute is
the ARG draw, which is shared across sites; that is pinned separately in
tests/test_draw_marginalization.py.
"""
from __future__ import annotations

import numpy as np
import pytest

import posterior_sample_age_infer as inf

N_PANEL = 8
EPS = 0.01

# Ages log-spaced as in a real table, and a T grid that reaches well inside the
# age range so that T >= t_i fires part-way along the branches tested below.
AGE = np.geomspace(60.0, 6000.0, 11)
TGRID = np.linspace(0.0, 3000.0, 25)


def _table():
    """A table that varies along BOTH axes, unlike conftest.make_table().

    make_table() is deliberately constant in age with every T below age[0], so
    the branch quadrature there is trivially exact and the existence mask never
    fires -- the two things this file is about.  Values need only be a valid
    (phi, phi2) pair in [0,1] with phi^2 <= phi2 <= phi; the identity under
    test does not depend on them being the real diffusion solution.
    """
    d0 = np.arange(1, N_PANEL + 1)
    a = AGE[None, :, None]
    t = TGRID[None, None, :]
    phi = (d0[:, None, None] / (N_PANEL + 1.0)) * np.sqrt(a / (a + t))
    phi2 = 0.5 * (phi ** 2 + phi)          # strictly inside [phi^2, phi]
    tab = {"table": phi.astype(np.float32), "table2": phi2.astype(np.float32),
           "d0": d0, "age": AGE, "Tgrid": TGRID, "n_sample": N_PANEL}
    # log P(d0 | n, t_i): smooth, finite, and varying in age so that the
    # weighted measure is genuinely different from the uniform one.
    tab["log_den"] = (-0.5 * np.log(AGE)[None, :] - AGE[None, :] / 4000.0
                      - 0.1 * d0[:, None]).astype(np.float64)
    return tab


def _log_den_at(tab, d0, a):
    """Mirror of phi_lookup's log_den_at: linear in table-row index, log-age."""
    la = np.log(np.clip(tab["age"], 1e-9, None))
    LD = tab["log_den"][d0 - 1]
    k = np.interp(np.log(max(a, 1e-9)), la, np.arange(len(tab["age"])))
    k0 = int(np.floor(k)); k1 = min(k0 + 1, len(tab["age"]) - 1); w = k - k0
    return float(LD[k0] if w == 0.0 else (1 - w) * LD[k0] + w * LD[k1])


def _point_rows(tab, d0, t_lo, t_hi, key, n_nodes, marginalise):
    """Production POINT-age lookups across the covered part of the branch.

    Returns (nodes, rows, weights).  The nodes start at the table's youngest
    age when the branch reaches below it, matching phi_lookup's convention of
    integrating only the covered part; the caller still normalises by the true
    branch width.
    """
    lo = max(float(t_lo), float(tab["age"][0]))
    hi = float(t_hi)
    nodes = np.linspace(lo, hi, n_nodes)
    rows = np.stack([inf.phi_lookup(tab, d0, t, t, key=key,
                                    marginalise="uniform") for t in nodes])
    if marginalise == "uniform":
        weights = np.ones_like(nodes)
    else:
        ld = np.array([_log_den_at(tab, d0, t) for t in nodes])
        weights = np.exp(ld - ld.max())
    return nodes, rows, weights


def _weighted_mean(nodes, rows, weights, width):
    """Branch average of `rows`, normalised by `width` (None -> the weight mass)."""
    num = np.trapezoid(rows * weights[:, None], nodes, axis=0)
    den = width if width is not None else np.trapezoid(weights, nodes)
    return num / den


# ---------------------------------------------------------------------------
# (a) commutation: exact, because the likelihood is affine in the table
# ---------------------------------------------------------------------------

BRANCHES = [(200.0, 2500.0),     # spans many age knots; T lands inside it
            (1000.0, 1050.0),    # narrow, inside a single knot interval
            (4000.0, 5800.0)]    # entirely older than every T on the grid


@pytest.mark.parametrize("t_lo, t_hi", BRANCHES)
@pytest.mark.parametrize("d0", [1, 3, 7])
@pytest.mark.parametrize("marginalise", ["uniform", "weighted"])
def test_haploid_likelihood_average_equals_likelihood_of_average(
        t_lo, t_hi, d0, marginalise):
    """mean_t [eps + (1-2eps) phi(t)] == eps + (1-2eps) mean_t phi(t)."""
    tab = _table()
    width = None if marginalise == "weighted" else t_hi - t_lo
    nodes, rows, w = _point_rows(tab, d0, t_lo, t_hi, "table", 2001, marginalise)

    phi_bar = _weighted_mean(nodes, rows, w, width)

    for observed_alt in (True, False):
        # ell evaluated at every age on the branch, then averaged ...
        qA_nodes = EPS + (1 - 2 * EPS) * rows
        ell_nodes = qA_nodes if observed_alt else 1.0 - qA_nodes
        ell_bar = _weighted_mean(nodes, ell_nodes, w, width)
        # ... versus ell evaluated once, at the averaged frequency.
        qA_bar = EPS + (1 - 2 * EPS) * phi_bar
        ell_of_bar = qA_bar if observed_alt else 1.0 - qA_bar
        assert np.allclose(ell_bar, ell_of_bar, rtol=1e-12, atol=1e-14)


@pytest.mark.parametrize("t_lo, t_hi", BRANCHES)
@pytest.mark.parametrize("d0", [1, 3, 7])
def test_diploid_dosage_probabilities_also_commute(t_lo, t_hi, d0):
    """P(dosage=k) is linear in (phi, phi2) jointly, so the same swap is exact.

    This is the case where plugging in a mean would be wrong if only ONE plane
    were averaged: P(k) is quadratic in r, hence nonlinear in phi alone.  It
    commutes because phi_lookup averages the second-moment plane over the same
    branch with the same weights.
    """
    tab = _table()
    width = t_hi - t_lo
    nodes, r1, w = _point_rows(tab, d0, t_lo, t_hi, "table", 2001, "uniform")
    _n2, r2, _w2 = _point_rows(tab, d0, t_lo, t_hi, "table2", 2001, "uniform")

    def dosage(phi, phi2):
        qA = EPS + (1 - 2 * EPS) * phi
        Er2 = (EPS ** 2 + 2 * EPS * (1 - 2 * EPS) * phi
               + (1 - 2 * EPS) ** 2 * phi2)
        return np.stack([1.0 - 2.0 * qA + Er2, 2.0 * (qA - Er2), Er2])

    per_node = np.stack([dosage(r1[i], r2[i]) for i in range(len(nodes))])
    averaged_likelihood = _weighted_mean(nodes, per_node.reshape(len(nodes), -1),
                                         w, width).reshape(3, -1)
    likelihood_of_averaged = dosage(_weighted_mean(nodes, r1, w, width),
                                    _weighted_mean(nodes, r2, w, width))

    assert np.allclose(averaged_likelihood, likelihood_of_averaged,
                       rtol=1e-12, atol=1e-14)
    # The dosage probabilities must still be a distribution after the swap.
    assert np.allclose(likelihood_of_averaged.sum(axis=0), 1.0, atol=1e-12)


# ---------------------------------------------------------------------------
# (b) the production integrator against an independent fine quadrature
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("t_lo, t_hi", BRANCHES)
@pytest.mark.parametrize("d0", [1, 3, 7])
@pytest.mark.parametrize("key", ["table", "table2"])
def test_closed_form_branch_integral_matches_fine_quadrature(t_lo, t_hi, d0, key):
    """phi_lookup's knot-splitting integral == trapezoid over point-age rows.

    The reference resolves the T >= t_i boundary only to node spacing, while
    phi_lookup resolves it exactly as an integration limit, so the tolerance is
    the quadrature error and not machine precision.
    """
    tab = _table()
    nodes, rows, w = _point_rows(tab, d0, t_lo, t_hi, key, 40001, "uniform")
    reference = _weighted_mean(nodes, rows, w, t_hi - t_lo)
    produced = inf.phi_lookup(tab, d0, t_lo, t_hi, key=key)
    assert np.allclose(produced, reference, rtol=2e-3, atol=2e-4)


@pytest.mark.parametrize("t_lo, t_hi", [(200.0, 2500.0), (1000.0, 1050.0)])
def test_weighted_branch_integral_matches_fine_quadrature(t_lo, t_hi):
    """Same check for --marginalise weighted, whose measure is exp(log_den)."""
    tab = _table()
    nodes, rows, w = _point_rows(tab, d0 := 3, t_lo, t_hi, "table", 40001,
                                 "weighted")
    reference = _weighted_mean(nodes, rows, w, None)
    produced = inf.phi_lookup(tab, d0, t_lo, t_hi, marginalise="weighted")
    assert np.allclose(produced, reference, rtol=2e-3, atol=2e-4)


def test_degenerate_branch_reduces_to_the_point_lookup():
    """A zero-width branch must not go through the integrator at all."""
    tab = _table()
    for t in (80.0, 800.0, 5000.0):
        assert np.array_equal(inf.phi_lookup(tab, 3, t, t),
                              inf.phi_lookup(tab, 3, t, t, key="table"))
        point = inf.phi_lookup(tab, 3, t, t)
        assert np.all(point[TGRID >= t] == 0.0)


def test_branch_reaching_below_the_table_keeps_the_true_width():
    """The low-clip convention: numerator over the covered part, denominator the
    true branch.  Exact only for T >= age[0], which is what is asserted."""
    tab = _table()
    t_lo, t_hi = 5.0, 1500.0                      # 5 < age[0] = 60
    nodes, rows, w = _point_rows(tab, 3, t_lo, t_hi, "table", 40001, "uniform")
    reference = _weighted_mean(nodes, rows, w, t_hi - t_lo)
    produced = inf.phi_lookup(tab, 3, t_lo, t_hi)
    supported = TGRID >= tab["age"][0]
    assert np.allclose(produced[supported], reference[supported],
                       rtol=2e-3, atol=2e-4)
    # Normalising by the covered width instead would inflate the site by this.
    inflation = (t_hi - t_lo) / (t_hi - tab["age"][0])
    assert inflation > 1.0
    assert not np.allclose(produced[supported],
                           reference[supported] * inflation, rtol=1e-2)
