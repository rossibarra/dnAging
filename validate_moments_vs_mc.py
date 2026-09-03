#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
validate_moments_vs_mc.py
=========================

Reproduces the validation behind this pipeline: the exact neutral WF moment
recursion for E[p_T | d0, t_i] AND E[p_T^2 | d0, t_i] (the engine in
precompute_freq_trajectory_moments.py) is checked against a forward Wright-Fisher-
diffusion Monte Carlo, under a CONSTANT Ne (where the MC is cheap). They should
agree to MC noise, including rare d0. The second moment is what makes the diploid
genotype likelihood correct, and it comes from the same contraction shifted one
index (E[X_T^2 X_pres^m] = sum_j C[m,j] M(u1)[j+2]).

    python validate_moments_vs_mc.py            # prints a comparison table

This is a self-contained sanity check (numpy + scipy only); it does not touch the
repo, the store, or any VCF.
"""
from __future__ import annotations
import numpy as np
from math import comb
from scipy.linalg import expm

NE = 1000.0            # diploid Ne
NC = 2 * NE            # gene copies; single new copy freq = 1/NC
N = 26                 # sample haplotypes


def moment_engine(n):
    K = n + 1
    B = np.zeros((K + 2, K + 2))          # one power beyond K: E[X_T^2 X^m] needs M_{j+2}
    for k in range(1, K + 2):
        B[k, k] = -k * (k - 1) / 2.0
        B[k, k - 1] = k * (k - 1) / 2.0
    return B, K


def analytic(d0, t_i, T):
    """(E[p_T | d0, t_i], E[p_T^2 | d0, t_i])."""
    B, K = moment_engine(N)
    eps = 1.0 / NC
    tau_i, tau_T = t_i / NC, T / NC
    if tau_T >= tau_i:
        return 0.0, 0.0
    def moms(u):
        m0 = np.array([eps ** k for k in range(K + 2)]); m0[0] = 1.0
        return expm(B * u) @ m0
    Mu1 = moms(tau_i - tau_T); C = expm(B * tau_T); Mpres = moms(tau_i)
    EjX = C[:, :K] @ Mu1[1:K + 1]         # E[X_T   X_pres^m]
    EjX2 = C[:, :K] @ Mu1[2:K + 2]        # E[X_T^2 X_pres^m]
    num = num2 = den = 0.0
    for m in range(d0, N + 1):
        c = comb(N - d0, m - d0) * (-1) ** (m - d0)
        den += c * Mpres[m]; num += c * EjX[m]; num2 += c * EjX2[m]
    return num / den, num2 / den


def mc(d0_target, t_i, T, M=4_000_000, seed=0):
    """(MC E[p_T|d0], MC E[p_T^2|d0], n_replicates) from a forward WF diffusion."""
    rng = np.random.default_rng(seed)
    dt = 0.005 * NC; p = np.full(M, 1.0 / NC); t = t_i; pT = None
    while t > 0:
        step = min(dt, t); act = (p > 0) & (p < 1)
        p[act] += np.sqrt(np.maximum(p[act] * (1 - p[act]) / NC, 0) * step) * rng.standard_normal(act.sum())
        p = np.clip(p, 0, 1); t -= step
        if pT is None and t <= T:
            pT = p.copy()
    x0 = p; surv = x0 > 0
    d0 = rng.binomial(N, x0[surv]); pv = pT[surv]; mm = d0 == d0_target
    if mm.sum() <= 50:
        return float("nan"), float("nan"), int(mm.sum())
    return float(pv[mm].mean()), float((pv[mm] ** 2).mean()), int(mm.sum())


if __name__ == "__main__":
    print(f"constant Ne={NE:.0f} (Nc={NC:.0f}), n={N}\n")
    print(f"{'d0':>3} {'t_i':>6} {'T':>6} | {'E[p]':>7} {'MC':>7} | {'E[p^2]':>7} "
          f"{'MC':>7} | {'E[p]^2':>7} {'MC_n':>6}")
    for d0, t_i, T in [(8, 4000, 1000), (5, 3000, 800), (12, 5000, 1500),
                       (2, 2000, 500), (2, 4000, 1000), (1, 1500, 300)]:
        a, a2 = analytic(d0, t_i, T); mv, mv2, cnt = mc(d0, t_i, T)
        # E[p]^2 is printed to show how far the plug-in (squared mean) sits from the
        # true second moment -- the error the diploid likelihood used to make
        print(f"{d0:>3} {t_i:>6} {T:>6} | {a:>7.3f} {mv:>7.3f} | {a2:>7.3f} "
              f"{mv2:>7.3f} | {a*a:>7.3f} {cnt:>6}")
