#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
precompute_freq_trajectory_moments.py
=====================================

Build the age-conditioned frequency-trajectory table

    table[i_n, d0, i_age, i_T]
        = E[ p_T | called panel size n, present count d0, mutation age t_i ]

via the neutral Wright-Fisher moment recursion -- no Monte Carlo, no diffusion
PDE, no rare-bin coverage problem. The recursion and the binomial conditioning are
EXACT within the neutral diffusion; the table itself is a numerical tabulation of
them (float32 entries on a finite log-spaced age grid, read back by interpolation),
so "exact" here always means exact as a statement about the model, never about the
arithmetic -- see the numerics note at the end. This is the established route
(Griffiths 2003; Song & Steinrucken 2012; used e.g. in the allele-age /
selection literature): the neutral moment hierarchy is closed,

    dM_k/dtau = (k(k-1)/2) (M_{k-1} - M_k),      M_k(tau) = E[X(tau)^k],

so for a sample of n chromosomes we only need moments up to order n+2 (n+1 for the
first moment, one more for the second) from matrix exponentials of an (n+3) x (n+3)
bidiagonal generator, and the Binomial(n, x) sampling does the conditioning on the
observed count d0. Emoments() takes THREE such exponentials per (n, d0, t_i, T)
entry -- see the cost note at the end. Validated against a forward WF-diffusion
Monte Carlo (agreement to MC noise for both moments, rare bins included;
validate_moments_vs_mc.py).

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

Numerics note: the conditioning step sums binom(n-d0,m-d0)(-1)^{m-d0} M_m with
alternating signs, which loses roughly 0.3 decimal digits per chromosome. float64
cannot carry that for realistic panels, so the table is built by
ExactMomentEngine, which works at 30 + n digits and reaches the value exactly (it
agrees with the independent reference in tests/_reference.py to the last digit
tested). The float64 engine remains available as --float64 for comparison against
older tables, but it fails in two ways and should not be used to build one: at
n=40 it returns 59% of entries as NaN, and -- less obviously -- entries that pass
its max_cancellation guard can already be wrong in the third significant digit
(n=26, d0=8, tau_i=3, tau_T=1: 0.470396 against the true 0.466832). Only the
legacy path writes NaN for numerical failure; the exact path writes NaN solely for
d0 above a given panel size, and 0 where T >= t_i.

Cost note: one exponential of this size is trivial, but build_table() loops over
every panel size n, every d0 in 1..n, every mutation age and every sample age.
Emoments() used to recompute all three of expm(B*u1), expm(B*tau_T) and
expm(B*tau_i) per call, ~7.8e6 exponentials of dimension 23..29 for the default
grid (n = 20..26, 100 log-spaced ages, 300 sample ages). All three depend only on
(n, tau_i, tau_T) and NOT on d0, and ExactMomentEngine.grid() exploits that: it
evaluates a whole (d0, T) block per mutation age, caching C = e^{B tau_T} per
tau_T and M(tau_i) per age, so each entry costs O(n^2) high-precision multiplies
instead of an exponential. It also never forms a general matrix exponential at
all -- B is lower-bidiagonal with eigenvalues -k(k-1)/2, so e^{B u} has an exact
partial-fraction expansion whose rational coefficients are computed once per panel
size. Net effect: the default build is ~36 min on one core at 56 digits, i.e. the
higher precision costs less than the arrangement it replaced.
"""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

import mpmath as mp
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


class ExactMomentEngine:
    """MomentEngine's conditional evaluated to `dps` digits instead of in float64.

    The conditioning step needs the alternating sum

        sum_{m=d0}^{n} binom(n-d0, m-d0) (-1)^{m-d0} M_m ,

    which loses roughly 0.3*n decimal digits to cancellation. float64 therefore
    runs out of precision as the panel grows: at n=40, 59% of table entries trip
    MomentEngine.max_cancellation and are returned as NaN, and -- worse -- entries
    that survive the guard can already be wrong in the third significant digit
    (n=26, d0=8, tau_i=3, tau_T=1: float64 gives 0.470396, the true value is
    0.466832). Carrying `dps` digits removes the loss; the guard then never fires.

    Doing that naively would mean an mpmath matrix exponential per (d0, T) entry.
    Two structural facts make it affordable instead:

    * B is lower-bidiagonal with eigenvalues lam_k = -k(k-1)/2, so e^{B u} has an
      exact partial-fraction expansion in the e^{lam_k u}, with rational
      coefficients computed once per panel size (`_coeffs`). No general-purpose
      matrix exponential is ever formed.
    * Of the three exponentials Emoments needs, C = e^{B tau_T} depends only on
      tau_T and M(tau_i) only on tau_i -- neither depends on d0. `grid()`
      evaluates a whole (d0, T) block per mutation age, caching C per tau_T, so
      each entry costs O(n^2) high-precision multiplies rather than an O(n^3)
      exponential.

    This is the same closed form the independent reference in tests/_reference.py
    uses; the two implementations are kept separate on purpose so the reference
    stays an independent check.
    """

    def __init__(self, n, dps=None):
        self.n = int(n)
        self.K = self.n + 1
        self.N = self.K + 2          # E[X_T^2 X_pres^m] reaches M_{j+2}
        # ~0.3 digits lost per chromosome, plus float64's 16 and a safety margin
        self.dps = int(dps) if dps else 30 + self.n
        self.coeff = {d0: {m: comb(self.n - d0, m - d0) * (-1) ** (m - d0)
                           for m in range(d0, self.n + 1)}
                      for d0 in range(1, self.n + 1)}
        self._lam_int = [-(k * (k - 1) // 2) for k in range(self.N)]
        with mp.workdps(self.dps):
            self._lam = [mp.mpf(l) for l in self._lam_int]
            # rational -> mpf once; these are hit O(n^2) times per table entry
            self._Ccoef = {ji: [(k, mp.mpf(v.numerator) / v.denominator)
                                for k, v in row.items()]
                           for ji, row in self._coeffs().items()}
        self._Ccache = {}

    def _coeffs(self):
        """D[(j, i)][k] with (e^{B u})_{ij} = sum_k D[(j,i)][k] e^{lam_k u}.

        For a lower-bidiagonal generator the (i, j) entry is a single cascade
        j -> i, giving the standard partial-fraction form

            (prod_{m=j+1}^{i} c_m) * sum_{k=j}^{i} e^{lam_k u}
                                     / prod_{l != k} (lam_k - lam_l),  c_m = -lam_m.

        lam_0 = lam_1 = 0 are the only repeated eigenvalues, and c_1 = 0 kills the
        prefactor of every block whose range spans both, so no repeated root ever
        reaches a denominator (asserted below).
        """
        lam = [Fraction(l) for l in self._lam_int]
        c = [-l for l in lam]
        D = {}
        for j in range(self.N):
            pref = Fraction(1)
            for i in range(j, self.N):
                if i > j:
                    pref *= c[i]
                row = {}
                if pref != 0:
                    for k in range(j, i + 1):
                        den = Fraction(1)
                        for l in range(j, i + 1):
                            if l != k:
                                assert lam[k] != lam[l], (
                                    "repeated eigenvalue reached a partial-fraction "
                                    f"denominator at (j={j}, i={i}, k={k}, l={l})")
                                den *= lam[k] - lam[l]
                        row[k] = pref / den
                D[(j, i)] = row
        return D

    def _C(self, tau_T):
        """C[m][j] = coefficient of y^j in E[X(tau_T)^m | X(0)=y]. Cached per tau_T."""
        key = float(tau_T)
        hit = self._Ccache.get(key)
        if hit is not None:
            return hit
        with mp.workdps(self.dps):
            t = mp.mpf(float(tau_T))
            e = [mp.e ** (l * t) for l in self._lam]
            C = [[mp.mpf(0)] * self.K for _ in range(self.N)]
            for j in range(self.K):
                for i in range(j, self.N):
                    row = self._Ccoef[(j, i)]
                    if row:
                        C[i][j] = mp.fsum(v * e[k] for k, v in row)
        self._Ccache[key] = C
        return C

    def moms(self, u, eps):
        """M_k(u) = E[X(u)^k] for a mutation entering at frequency eps."""
        with mp.workdps(self.dps):
            t = mp.mpf(float(u))
            e = [mp.e ** (l * t) for l in self._lam]
            x0 = mp.mpf(float(eps))
            m0 = [mp.mpf(1)] + [x0 ** k for k in range(1, self.N)]
            out = []
            for i in range(self.N):
                s = mp.mpf(0)
                for j in range(i + 1):
                    row = self._Ccoef[(j, i)]
                    if row and m0[j]:
                        s += m0[j] * mp.fsum(v * e[k] for k, v in row)
                out.append(s)
            return out

    def grid(self, tau_i, tauT, eps):
        """(E[p_T|d0,t_i], E[p_T^2|d0,t_i]) for every (d0, T), shape (n, len(tauT)).

        Entries with tau_T >= tau_i (sample older than the mutation) are 0, matching
        MomentEngine.Emoments.
        """
        tauT = np.asarray(tauT, dtype=np.float64)
        p1 = np.zeros((self.n, tauT.size))
        p2 = np.zeros((self.n, tauT.size))
        with mp.workdps(self.dps):
            Mpres = self.moms(tau_i, eps)
            den = {d0: mp.fsum(c * Mpres[m] for m, c in self.coeff[d0].items())
                   for d0 in range(1, self.n + 1)}
            for it, tT in enumerate(tauT):
                if tT >= tau_i:
                    continue
                C = self._C(tT)
                Mu1 = self.moms(tau_i - tT, eps)
                EjX = [mp.fsum(C[m][j] * Mu1[j + 1] for j in range(self.K))
                       for m in range(self.N)]
                EjX2 = [mp.fsum(C[m][j] * Mu1[j + 2] for j in range(self.K))
                        for m in range(self.N)]
                for d0, cs in self.coeff.items():
                    d = den[d0]
                    if d == 0:
                        continue
                    a = mp.fsum(c * EjX[m] for m, c in cs.items()) / d
                    b = mp.fsum(c * EjX2[m] for m, c in cs.items()) / d
                    p1[d0 - 1, it] = float(a)
                    p2[d0 - 1, it] = float(b)
        # exact moment constraints; with dps digits these hold, so assert rather
        # than clip -- a violation now means a real bug, not roundoff
        bad = (p1 < -1e-12) | (p1 > 1 + 1e-12) | (p2 < p1 * p1 - 1e-12) | (p2 > p1 + 1e-12)
        if bad.any():
            i = np.argwhere(bad)[0]
            raise AssertionError(
                f"moment constraints violated at d0={i[0]+1}, tau_T={tauT[i[1]]:g}, "
                f"tau_i={tau_i:g}: E[p]={p1[tuple(i)]!r}, E[p^2]={p2[tuple(i)]!r}")
        return np.clip(p1, 0.0, 1.0), np.clip(p2, p1 * p1, p1)

    def Emoments(self, d0, tau_i, tau_T, eps):
        """Scalar form, for parity with MomentEngine.Emoments."""
        if tau_T >= tau_i:
            return 0.0, 0.0
        a, b = self.grid(tau_i, [tau_T], eps)
        return float(a[d0 - 1, 0]), float(b[d0 - 1, 0])

    def Efreq(self, d0, tau_i, tau_T, eps):
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
        eng = (MomentEngine(int(n)) if args.float64
               else ExactMomentEngine(int(n), dps=args.precision))
        for ia, t_i in enumerate(age):
            tau_i = tau_of_t(t_i)
            eps = 1.0 / (2.0 * float(ne_of_t(t_i)[0]))
            if args.float64:
                for id0, d0 in enumerate(range(1, n + 1)):
                    row = np.array([eng.Emoments(d0, tau_i, tt, eps) for tt in tauT])
                    table[inx, id0, ia] = row[:, 0]
                    table2[inx, id0, ia] = row[:, 1]
            else:
                # one (d0, T) block per age: C = e^{B tau_T} and M(tau_i) are
                # shared across d0, so they are computed once instead of n times
                p1, p2 = eng.grid(tau_i, tauT, eps)
                table[inx, :n, ia] = p1
                table2[inx, :n, ia] = p2
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
    p.add_argument("--precision", type=int, default=None,
                   help="digits carried through the conditioning sum "
                        "[30 + n_sample]. The sum loses ~0.3 digits per "
                        "chromosome to cancellation, so float64 is not enough.")
    p.add_argument("--float64", action="store_true",
                   help="use the legacy float64 engine. Loses ~0.3*n digits to "
                        "cancellation: at n=40 it drops 59%% of entries as NaN, "
                        "and surviving entries can be wrong in the third digit. "
                        "For comparison against old tables only.")
    p.add_argument("--quiet", action="store_true")
    args = p.parse_args(argv)
    if not 2 <= args.min_n <= args.n_sample:
        p.error("--min-n must satisfy 2 <= min-n <= n-sample")
    if args.precision is not None and args.precision < 20:
        p.error("--precision below 20 digits defeats the purpose; omit it for the "
                "panel-size default")

    table, table2, d0, n_panel, age, age_tau, Tgrid, windows = build_table(args)
    meta = {"n_sample": args.n_sample, "ne_file": str(args.ne),
            "ne_series": args.ne_series, "method": "neutral WF moment recursion",
            "ne_windows": int(len(windows[0])),
            # format 4 adds an n_panel axis for partially called panel sites
            "format_version": 4, "planes": ["table (E[p_T])", "table2 (E[p_T^2])"],
            # provenance: which arithmetic produced this table
            "arithmetic": ("float64 (legacy, cancellation-limited)" if args.float64
                           else f"mpmath dps={args.precision or 30 + args.n_sample}")}
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
