#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
validate_moments_vs_mc.py
=========================

Reproduces the validation behind this pipeline: the exact neutral WF moment
recursion for E[p_T | d0, t_i] (the engine in precompute_freq_trajectory_moments.py)
is checked against a forward Wright-Fisher-diffusion Monte Carlo, under a CONSTANT
Ne (where the MC is cheap). They should agree to MC noise, including rare d0.

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
    B = np.zeros((K + 1, K + 1))
    for k in range(1, K + 1):
        B[k, k] = -k * (k - 1) / 2.0
        B[k, k - 1] = k * (k - 1) / 2.0
    return B, K


def analytic(d0, t_i, T):
    B, K = moment_engine(N)
    eps = 1.0 / NC
    tau_i, tau_T = t_i / NC, T / NC
    if tau_T >= tau_i:
        return 0.0
    def moms(u):
        m0 = np.array([eps ** k for k in range(K + 1)]); m0[0] = 1.0
        return expm(B * u) @ m0
    Mu1 = moms(tau_i - tau_T); C = expm(B * tau_T); Mpres = moms(tau_i)
    EjX = C[:, :K] @ Mu1[1:K + 1]
    num = den = 0.0
    for m in range(d0, N + 1):
        c = comb(N - d0, m - d0) * (-1) ** (m - d0)
        den += c * Mpres[m]; num += c * EjX[m]
    return num / den


def mc(d0_target, t_i, T, M=4_000_000, seed=0):
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
    return (float(pv[mm].mean()), int(mm.sum())) if mm.sum() > 50 else (float("nan"), int(mm.sum()))


if __name__ == "__main__":
    print(f"constant Ne={NE:.0f} (Nc={NC:.0f}), n={N}\n")
    print(f"{'d0':>3} {'t_i':>6} {'T':>6} | {'moment':>7} {'MC':>7} {'MC_n':>6}")
    for d0, t_i, T in [(8, 4000, 1000), (5, 3000, 800), (12, 5000, 1500),
                       (2, 2000, 500), (2, 4000, 1000), (1, 1500, 300)]:
        a = analytic(d0, t_i, T); mv, cnt = mc(d0, t_i, T)
        print(f"{d0:>3} {t_i:>6} {T:>6} | {a:>7.3f} {mv:>7.3f} {cnt:>6}")
