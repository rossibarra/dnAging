# -*- coding: utf-8 -*-
"""Independent high-precision reference for the neutral moment recursion.

Deliberately shares NO code with precompute_freq_trajectory_moments: the moment
map C(u) = expm(B u) is built from the CLOSED-FORM solution of
dM_k/du = lambda_k (M_{k-1} - M_k), lambda_k = k(k-1)/2, rather than from a
general-purpose matrix exponential. Each C[m,j] is an exact rational combination
of exp(-lambda_l u), evaluated in mpmath, so the reference is limited only by
`dps` -- it is not subject to the float64 cancellation the engine suffers.
"""
from __future__ import annotations

from fractions import Fraction
from math import comb

import mpmath as mp

DPS = 80


def _lam(k):
    return k * (k - 1) // 2


def _expsum_C(sz, ):
    """f[m][j] = {lambda: exact rational coeff} with C(u)[m,j] = sum coeff e^{-lam u}.

    P_m(u,y) = E[X(u)^m | X(0)=y] = sum_j C[m,j] y^j obeys the same recursion in m,
    and it decouples across powers of y, so each column j is one scalar cascade
    with initial condition delta_{mj}.
    """
    f = [[None] * sz for _ in range(sz)]
    for j in range(sz):
        for m in range(sz):
            if m < j:
                f[m][j] = {}
            elif m == j:
                f[m][j] = {_lam(j): Fraction(1)}
            elif m == 1:                       # lambda_1 = 0: dM/du = 0
                f[m][j] = {}
            else:
                lm = _lam(m)
                prev = f[m - 1][j]
                out = {}
                acc = Fraction(0)
                for l, a in prev.items():
                    assert l != lm, "lambda collision breaks the partial fractions"
                    c = a * Fraction(lm, lm - l)
                    out[l] = out.get(l, Fraction(0)) + c
                    acc += c
                out[lm] = out.get(lm, Fraction(0)) - acc   # enforce f(0) = 0
                f[m][j] = {l: c for l, c in out.items() if c != 0}
    return f


_CACHE = {}


def C_matrix(n, u):
    """Exact expm(B*u) for the (n+3)x(n+3) neutral moment generator, in mpmath."""
    sz = n + 3
    f = _CACHE.setdefault(sz, _expsum_C(sz))
    with mp.workdps(DPS):
        u = mp.mpf(u)
        ex = {}
        out = [[mp.mpf(0)] * sz for _ in range(sz)]
        for m in range(sz):
            for j in range(sz):
                s = mp.mpf(0)
                for l, c in f[m][j].items():
                    if l not in ex:
                        ex[l] = mp.e ** (-mp.mpf(l) * u)
                    s += mp.mpf(c.numerator) / mp.mpf(c.denominator) * ex[l]
                out[m][j] = s
        return out


def moments(n, u, eps):
    """M(u) with M_k(0) = eps^k, i.e. C(u) @ [eps^j]."""
    C = C_matrix(n, u)
    sz = n + 3
    with mp.workdps(DPS):
        e = mp.mpf(eps)
        pw = [e ** j for j in range(sz)]
        return [sum(C[m][j] * pw[j] for j in range(sz)) for m in range(sz)]


def Emoments_ref(n, d0, tau_i, tau_T, eps):
    """(E[p_T|d0,t_i], E[p_T^2|d0,t_i]) at DPS digits -- unclamped, unclipped."""
    if tau_T >= tau_i:
        return mp.mpf(0), mp.mpf(0)
    K = n + 1
    with mp.workdps(DPS):
        Mu1 = moments(n, tau_i - tau_T, eps)
        C = C_matrix(n, tau_T)
        Mpres = moments(n, tau_i, eps)
        num = num2 = den = mp.mpf(0)
        for m in range(d0, n + 1):
            c = mp.mpf(comb(n - d0, m - d0) * (-1) ** (m - d0))
            EjX = sum(C[m][j] * Mu1[j + 1] for j in range(K))
            EjX2 = sum(C[m][j] * Mu1[j + 2] for j in range(K))
            den += c * Mpres[m]
            num += c * EjX
            num2 += c * EjX2
        return num / den, num2 / den
