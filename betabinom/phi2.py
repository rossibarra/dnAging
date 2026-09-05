"""phi2(n_T, nleaf, u1, tauT) for STRADDLING edges: condition on both

    * exactly 1 of n_T panel lineages carries it at T
    * nleaf of n modern samples carry it today

  num = sum_{i,j} C(n_T-1,i)(-1)^i C(n-nleaf,j)(-1)^j E[X_T^{2+i} X_pres^{nleaf+j}]
  den = same with X_T^{1+i}
  E[X_T^a X_pres^m] = sum_l K(tauT)[m,l] M(u1)[a+l]

u1 = diffusion time from mutation origin to T; tauT = T to present.
Needs moments to order n_T+1+n, so it is built in mpmath.
"""
import os, sys, time, numpy as np, mpmath as mp
from math import comb
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "tests"))
from _reference import C_matrix

def phi2(nT, nleaf, u1, tauT, Ne, n=26, dps=80):
    mp.mp.dps=dps
    eps=mp.mpf(1)/(2*Ne)
    amax=nT+1; mmax=n
    need=amax+mmax                        # highest moment index of M(u1)
    Cm=C_matrix(max(mmax,2)-3+3, mp.mpf(tauT))   # K(tauT), rows m=0..mmax
    Mv_C=C_matrix(need-3+3, mp.mpf(u1))
    sz=need+1
    m0=[mp.mpf(1)]+[eps**k for k in range(1,sz)]
    M=[mp.fsum(Mv_C[r][c]*m0[c] for c in range(sz)) for r in range(sz)]
    def EX(a,m):
        return mp.fsum(Cm[m][l]*M[a+l] for l in range(mmax+1))
    num=den=mp.mpf(0)
    for i in range(nT):
        ci=mp.mpf(comb(nT-1,i))*(-1)**i
        for j in range(n-nleaf+1):
            cj=mp.mpf(comb(n-nleaf,j))*(-1)**j
            w=ci*cj; m=nleaf+j
            num+=w*EX(2+i,m); den+=w*EX(1+i,m)
    if den<=0 or num<=0: return float('nan')
    v=num/den
    return float(v) if 0<v<1 else float('nan')

if __name__=="__main__":
    NE=100_000; tau=lambda t: t/(2*NE)
    t0=time.time()
    v=phi2(24,1,tau(500.),tau(2500.),NE)
    print(f"one evaluation: {time.time()-t0:.1f}s   phi2(nT=24,nleaf=1)={v:.6g}")
    from phi import Phi
    P=Phi(np.geomspace(0.5,2e6,70),26,NE)
    print(f"\nphi (ignores nleaf) for nT=24, age=500: {P(1,24,500.):.6g}")
    print(f"\n{'nleaf':>7}{'phi2':>14}{'ratio to phi':>15}")
    base=P(1,24,500.)
    for nl in (1,2,3,5,8):
        t0=time.time(); v=phi2(24,nl,tau(500.),tau(2500.),NE)
        print(f"{nl:>7}{v:>14.6g}{v/base:>15.3f}   ({time.time()-t0:.1f}s)")
