#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
precompute_freq_trajectory_moments.py
=====================================

Build the age-conditioned frequency-trajectory table

    table[d0, i_age, i_T] = E[ p_T | present sample count d0, mutation age t_i ]

EXACTLY, via the neutral Wright-Fisher moment recursion -- no Monte Carlo, no
diffusion PDE, no rare-bin coverage problem. This is the established route
(Griffiths 2003; Song & Steinrucken 2012; used e.g. in the allele-age /
selection literature): the neutral moment hierarchy is closed,

    dM_k/dtau = (k(k-1)/2) (M_{k-1} - M_k),      M_k(tau) = E[X(tau)^k],

so for a sample of n chromosomes we only need moments up to order n+2 (n+1 for the
first moment, one more for the second), one matrix exponential of an (n+3) x (n+3)
bidiagonal generator, and the Binomial(n, x) sampling does the conditioning on the
observed count d0. Validated against a forward WF-diffusion Monte Carlo (agreement
to MC noise for both moments, rare bins included; validate_moments_vs_mc.py).

Time-varying Ne is handled EXACTLY by the diffusion-time change (the neutral
diffusion has no drift, so Ne enters only through the clock):

    tau(t) = \\int_0^t dt' / (2 Ne(t')),

a cumulative sum over your Ne windows. Everything downstream is Ne-free in tau.

Model / quantities
------------------
Let X(u) be the derived-allele population frequency, u = diffusion time since the
mutation arose (u=0 origin, u=tau_i present). A new mutation starts at a single
copy, X(0) = 1/(2 Ne(t_i)). With sample age T (before present) at diffusion time
tau_T, define u1 = tau_i - tau_T (origin -> sample age). Using the closed moment
map C(dt) = expm(B*dt) and moments M(u) = expm(B*u) m0,

    E[X_T * X_pres^m] = sum_j C(tau_T)[m, j] * M(u1)[j+1]
    Binom(d0; n, x) proportional to sum_{m>=d0} C(n-d0, m-d0) (-1)^(m-d0) x^m
    E[p_T | d0, t_i] = ( sum_m coeff_{d0}[m] E[X_T X_pres^m] )
                       / ( sum_m coeff_{d0}[m] M_pres[m] )

The SECOND moment comes from the same contraction shifted one index (each extra
factor of X_T raises the power of x in E[X_pres^m | X_T = x] by one), so it costs
nothing but one more moment order:

    E[X_T^2 * X_pres^m] = sum_j C(tau_T)[m, j] * M(u1)[j+2]
    E[p_T^2 | d0, t_i] = ( sum_m coeff_{d0}[m] E[X_T^2 X_pres^m] )
                         / ( sum_m coeff_{d0}[m] M_pres[m] )      # same denominator

It is required for DIPLOID genotype likelihoods, which are nonlinear in the latent
frequency (E[X^2] != E[X]^2); the first moment alone suffices for one haploid
Bernoulli observation.

For T >= t_i the allele does not yet exist, so p_T = 0 (and p_T^2 = 0).

Demography input
----------------
The coalescence-Ne TSV from RILAB/argtest scripts/coalescence_ne_plots_from_ts.py
(columns include series, time_left, time_right, effective_population_size; Ne =
1/(2*rate), a diploid effective size, so the population carries 2*Ne gene copies
and a single new copy has frequency 1/(2 Ne)). Uses --ne-series (default
posterior_mean) as a step function.

Output (.npz, --output)
-----------------------
    table   float32 (n_n, n_d0, n_age, n_T)   E[p_T | n,d0,t_i]
    table2  float32 (n_n, n_d0, n_age, n_T)   E[p_T^2 | n,d0,t_i]
    n_panel int     (n_n,)                 called-panel sizes min_n..n_sample
    d0      int     (n_d0,)               present counts 1..n
    age     float   (n_age,)              mutation ages t_i (generations)
    Tgrid   float   (n_T,)                sample ages T (generations)
    n_sample int                           n chromosomes (= 26)
    meta    (json)  parameters / provenance

Numerics note: the alternating binomial sums can lose all float64 precision for
large diffusion times. Entries with fewer than two estimated significant decimal
digits are written as NaN rather than clipped into apparently valid probabilities.
"""

from __future__ import annotations

import argparse
import json
import sys
from math import comb
from pathlib import Path

import numpy as np
from scipy.linalg import expm


# ---------------------------------------------------------------------------
# Demography: tau(t) and Ne(t) from the coalescence-Ne TSV
# ---------------------------------------------------------------------------


def load_demography(path, series="posterior_mean"):
    """Return (tau_of_t, ne_of_t) callables from the Ne windows.

    tau(t) = cumulative \\int dt/(2 Ne); ne(t) is the step-function Ne. Both
    extrapolate the last window's Ne beyond the tabulated range.
    """
    lefts, rights, nes = [], [], []
    header = None
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if header is None:
            header = parts
            ci = {c: k for k, c in enumerate(header)}
            for need in ("time_left", "time_right", "effective_population_size"):
                if need not in ci:
                    raise SystemExit(f"Ne TSV missing column {need!r}; got {header}")
            si = ci.get("series")
            continue
        if si is not None and series is not None and parts[si] != series:
            continue
        try:
            lo = float(parts[ci["time_left"]]); hi = float(parts[ci["time_right"]])
            ne = float(parts[ci["effective_population_size"]])
        except (ValueError, IndexError):
            continue
        if np.isfinite(ne) and ne > 0 and hi > lo:
            lefts.append(lo); rights.append(hi); nes.append(ne)
    if not nes:
        raise SystemExit(f"No usable rows for series={series!r} in {path}")
    o = np.argsort(lefts)
    L = np.array(lefts)[o]; R = np.array(rights)[o]; NE = np.array(nes)[o]
    # the cumsum below is tau at a window's LEFT edge only if the windows tile the
    # axis: a gap or an overlap mis-scales tau for every later window (and a t in a
    # gap would get a NEGATIVE elapsed term), so refuse rather than assume a demography
    tol = 1e-6 * np.maximum(np.abs(R[:-1]), np.maximum(np.abs(L[1:]), 1.0))
    bad = np.flatnonzero(np.abs(R[:-1] - L[1:]) > tol)
    if bad.size:
        i = int(bad[0]); kind = "gap" if L[i + 1] > R[i] else "overlap"
        raise SystemExit(f"Ne windows are not contiguous for series={series!r} in {path}: "
                         f"{kind} between window {i} [{L[i]:.10g},{R[i]:.10g}] and window "
                         f"{i+1} [{L[i+1]:.10g},{R[i+1]:.10g}] "
                         f"(time_right={R[i]:.10g} != time_left={L[i+1]:.10g})")
    # cumulative tau at each window's right edge
    seg = (R - L) / (2.0 * NE)
    cum_at_R = np.concatenate([[0.0], np.cumsum(seg)])   # cum_at_R[i] = tau(L[i])

    def ne_of_t(t):
        t = np.atleast_1d(np.asarray(t, float))
        idx = np.clip(np.searchsorted(R, t, side="right"), 0, len(NE) - 1)
        return NE[idx]

    def tau_of_t(t):
        t = float(t)
        if t <= L[0]:
            return t / (2.0 * NE[0])
        i = int(np.clip(np.searchsorted(R, t, side="right"), 0, len(NE) - 1))
        base = cum_at_R[i]                    # tau at left edge of window i
        # max(...,0): tau must never run backwards even if a window is degenerate
        return base + max(t - L[i], 0.0) / (2.0 * NE[i])

    return tau_of_t, ne_of_t, (L, R, NE)


# ---------------------------------------------------------------------------
# Neutral WF moment machinery
# ---------------------------------------------------------------------------


class MomentEngine:
    """Closed neutral moment recursion for a sample of n chromosomes."""

    def __init__(self, n):
        self.n = n
        self.K = n + 1                        # highest moment order needed
        K = self.K
        # one power BEYOND K: E[X_T^2 X_pres^m] needs M(u1)_{j+2} (see Emoments).
        # B is lower-bidiagonal, hence expm(B) lower-triangular, so the extra row
        # leaves every entry the first-moment contraction uses untouched.
        B = np.zeros((K + 2, K + 2))
        for k in range(1, K + 2):
            B[k, k] = -k * (k - 1) / 2.0
            B[k, k - 1] = k * (k - 1) / 2.0   # dM_k/dtau = k(k-1)/2 (M_{k-1}-M_k)
        self.B = B
        # binomial-conditioning coefficients coeff[d0][m] for m=d0..n
        self.coeff = {}
        for d0 in range(1, n + 1):
            self.coeff[d0] = {m: comb(n - d0, m - d0) * (-1) ** (m - d0)
                              for m in range(d0, n + 1)}

        # For a sum S of floating-point terms, eps * sum(abs(term)) / abs(S)
        # estimates its relative roundoff amplification. This cutoff corresponds
        # to only roughly 1--2 trustworthy decimal digits in float64.
        self.max_cancellation = 1.25e14

    def _moms(self, u, eps):
        m0 = np.array([eps ** k for k in range(self.K + 2)], dtype=np.float64)
        m0[0] = 1.0
        return expm(self.B * u) @ m0

    def Emoments(self, d0, tau_i, tau_T, eps):
        """(E[p_T | d0, t_i], E[p_T^2 | d0, t_i]); tau_i=tau(t_i), tau_T=tau(T),
        eps=1/(2Ne(t_i)). The SECOND moment is needed for a diploid genotype
        likelihood, which is nonlinear in the latent frequency."""
        if tau_T >= tau_i:                    # sample older than the mutation
            return 0.0, 0.0
        u1 = tau_i - tau_T
        Mu1 = self._moms(u1, eps)
        C = expm(self.B * tau_T)
        Mpres = self._moms(tau_i, eps)
        # E[X_T^k X_pres^m] = sum_j C[m,j] Mu1[j+k]: C[m,j] multiplies x^j in
        # E[X_pres^m | X_T=x], so each extra factor of X_T shifts the index by one
        EjX = C[:, :self.K] @ Mu1[1:self.K + 1]
        EjX2 = C[:, :self.K] @ Mu1[2:self.K + 2]
        num = num2 = den = 0.0
        abs_num = abs_num2 = abs_den = 0.0
        for m, c in self.coeff[d0].items():
            td = c * Mpres[m]; t1 = c * EjX[m]; t2 = c * EjX2[m]
            den += td; num += t1; num2 += t2
            abs_den += abs(td); abs_num += abs(t1); abs_num2 += abs(t2)

        def unreliable(total, absolute_total):
            return (not np.isfinite(total) or total == 0.0 or
                    absolute_total / abs(total) > self.max_cancellation)

        if unreliable(den, abs_den) or unreliable(num, abs_num):
            return np.nan, np.nan
        p1 = float(num / den)
        if unreliable(num2, abs_num2):
            p2 = np.nan
        else:
            p2 = float(num2 / den)
        # These are exact moment constraints. A violation is evidence of numerical
        # failure, not something clipping can repair.
        tol = 100.0 * np.finfo(np.float64).eps
        if p1 < -tol or p1 > 1.0 + tol:
            return np.nan, np.nan
        p1 = float(np.clip(p1, 0.0, 1.0))
        if np.isnan(p2) or p2 < p1 * p1 - tol or p2 > p1 + tol:
            return p1, np.nan
        return p1, float(np.clip(p2, p1 * p1, p1))

    def Efreq(self, d0, tau_i, tau_T, eps):
        """E[p_T | d0, t_i]; the first moment alone (haploid observations)."""
        return self.Emoments(d0, tau_i, tau_T, eps)[0]


# ---------------------------------------------------------------------------
# Table build
# ---------------------------------------------------------------------------


def build_table(args):
    tau_of_t, ne_of_t, windows = load_demography(args.ne, series=args.ne_series)
    panel_sizes = np.arange(args.min_n, args.n_sample + 1, dtype=np.int64)
    Tgrid = (np.loadtxt(args.t_grid) if args.t_grid else
             np.linspace(args.t_min, args.t_max, args.n_t)).astype(np.float64)
    age = np.geomspace(max(args.age_min, 1.0), args.age_max, args.n_age)

    tauT = np.array([tau_of_t(T) for T in Tgrid])
    shape = (len(panel_sizes), args.n_sample, args.n_age, len(Tgrid))
    table = np.full(shape, np.nan, dtype=np.float32)
    table2 = np.full(shape, np.nan, dtype=np.float32)
    for inx, n in enumerate(panel_sizes):
        eng = MomentEngine(int(n))
        for ia, t_i in enumerate(age):
            tau_i = tau_of_t(t_i)
            eps = 1.0 / (2.0 * float(ne_of_t(t_i)[0]))
            for id0, d0 in enumerate(range(1, n + 1)):
                row = np.array([eng.Emoments(d0, tau_i, tt, eps) for tt in tauT])
                table[inx, id0, ia] = row[:, 0]
                table2[inx, id0, ia] = row[:, 1]
            if not args.quiet:
                print(f"[n={n} age {ia+1}/{args.n_age}] t_i={t_i:.3g} "
                      f"tau_i={tau_i:.3g}", file=sys.stderr)
    age_tau = np.array([tau_of_t(t_i) for t_i in age], dtype=np.float64)
    if age_tau[-1] <= 3.0:
        raise SystemExit(f"--age-max={args.age_max:g} reaches only tau={age_tau[-1]:.6g}; "
                         "increase --age-max so the table extends beyond the default "
                         "inference cutoff tau=3")
    return (table, table2, np.arange(1, args.n_sample + 1), panel_sizes, age,
            age_tau, Tgrid, windows)


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="precompute_freq_trajectory_moments.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Exact E[p_T | d0, t_i] via the neutral WF moment recursion "
                    "under a coalescence-Ne(t) curve (time-change to diffusion "
                    "time). No Monte Carlo, no rare-bin gaps.")
    p.add_argument("--ne", type=Path, required=True,
                   help="coalescence-ne-estimates.tsv (RILAB/argtest).")
    p.add_argument("--ne-series", default="posterior_mean")
    p.add_argument("--n-sample", type=int, default=26,
                   help="sample chromosomes (the ARG panel) [26].")
    p.add_argument("--min-n", type=int, default=20,
                   help="smallest called-panel size to precompute [20].")
    p.add_argument("--t-min", type=float, default=0.0)
    p.add_argument("--t-max", type=float, default=30000.0,
                   help="max sample age T (generations) [30000].")
    p.add_argument("--n-t", type=int, default=300)
    p.add_argument("--t-grid", type=Path, default=None,
                   help="explicit sample-age grid file (overrides t-min/max/n-t).")
    p.add_argument("--age-min", type=float, default=10.0)
    p.add_argument("--age-max", type=float, default=4e7,
                   help="max mutation age t_i (generations) [4e7].")
    p.add_argument("--n-age", type=int, default=100,
                   help="log-spaced mutation-age grid points [100].")
    p.add_argument("--output", type=Path, required=True, help="output .npz")
    p.add_argument("--quiet", action="store_true")
    args = p.parse_args(argv)
    if not 2 <= args.min_n <= args.n_sample:
        p.error("--min-n must satisfy 2 <= min-n <= n-sample")

    table, table2, d0, n_panel, age, age_tau, Tgrid, windows = build_table(args)
    meta = {"n_sample": args.n_sample, "ne_file": str(args.ne),
            "ne_series": args.ne_series, "method": "neutral WF moment recursion",
            "ne_windows": int(len(windows[0])),
            # format 4 adds an n_panel axis for partially called panel sites
            "format_version": 4, "planes": ["table (E[p_T])", "table2 (E[p_T^2])"]}
    np.savez_compressed(args.output, table=table, table2=table2, d0=d0, age=age,
                        n_panel=n_panel, age_tau=age_tau, Tgrid=Tgrid,
                        n_sample=args.n_sample, min_n=args.min_n,
                        meta=json.dumps(meta))
    if not args.quiet:
        print(f"[precompute-moments] wrote {args.output}  shape={table.shape}",
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
