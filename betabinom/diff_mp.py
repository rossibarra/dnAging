"""High-precision version of MATH.md's diffusion conditional E[p_T | d_0, t_i].

The float64 implementation loses the whole quantity to cancellation for large
panels (59% of table cells at n=40).  B is lower-bidiagonal with eigenvalues
lam_k = -k(k-1)/2, so e^{B tau} has an exact partial-fraction expansion whose
coefficients are rationals computed once.  Substituting it for all three matrix
exponentials collapses the conditional to a bilinear form

    num_{d0}(tau_T, tau_i) = a^T F[d0] b,   a_k = e^{lam_k tau_T},
                                            b_k = e^{lam_k (tau_i - tau_T)}
    den_{d0}(tau_i)        = G[d0] . e^{lam tau_i}

with F, G exact.  Only the exponentials and the final sums need mpmath, so the
per-cell cost is O(n K^2) high-precision multiplies instead of a matrix
exponential.
"""
import numpy as np, mpmath as mp
from fractions import Fraction as Fr
from math import comb

def _expB_coeffs(N):
    c=[Fr(k*(k-1),2) for k in range(N)]; lam=[-x for x in c]
    D={}
    for j in range(N):
        pref=Fr(1)
        for i in range(j,N):
            if i>j: pref*=c[i]
            d={}
            if pref!=0:
                for k in range(j,i+1):
                    den=Fr(1)
                    for l in range(j,i+1):
                        if l!=k: den*=(lam[k]-lam[l])
                    d[k]=pref/den
            D[(j,i)]=d
    return D,lam

def build(n, Ne, dps=50):
    """Exact F, G and lam for a panel of n chromosomes with x0 = 1/(2Ne)."""
    mp.mp.dps=dps
    K=n+1; N=K+2
    D,lamF=_expB_coeffs(N)
    lam=[mp.mpf(int(x.numerator))/int(x.denominator) for x in lamF]
    x0=Fr(1,2*Ne)
    m0=[Fr(1)]+[x0**k for k in range(1,N)]
    # A2[i][k]: (e^{B t} m0)_i = sum_k A2[i][k] e^{lam_k t}
    A2=[[Fr(0)]*N for _ in range(N)]
    for j in range(N):
        for i in range(j,N):
            for k,cf in D[(j,i)].items(): A2[i][k]+=cf*m0[j]
    A={d0:{m:comb(n-d0,m-d0)*(-1)**(m-d0) for m in range(d0,n+1)}
       for d0 in range(1,n+1)}
    F=[]; G=[]
    for d0 in range(1,n+1):
        # P[j'][k] = sum_m A[d0][m] * D[(j',m)][k]
        P=[[Fr(0)]*N for _ in range(K)]
        for m,cA in A[d0].items():
            for jp in range(min(K,m+1)):     # C is lower-triangular: C[m,jp]=0 for jp>m
                for k,cf in D[(jp,m)].items(): P[jp][k]+=cA*cf
        Fd=[[Fr(0)]*N for _ in range(N)]
        for jp in range(K):
            row=A2[jp+1]
            for k in range(N):
                p=P[jp][k]
                if p==0: continue
                for kp in range(N):
                    if row[kp]: Fd[k][kp]+=p*row[kp]
        nz=[(k,kp,mp.mpf(int(v.numerator))/int(v.denominator))
            for k in range(N) for kp in range(N) for v in (Fd[k][kp],) if v!=0]
        F.append(nz)
        g=[Fr(0)]*N
        for m,cA in A[d0].items():
            for k in range(N): g[k]+=cA*A2[m][k]
        G.append([mp.mpf(int(v.numerator))/int(v.denominator) for v in g])
    return lam, F, G

def table(n, Ne, tig, grid, dps=50, floor=mp.mpf(10)**-40):
    """tab[d0-1, i_ti, i_T] = E[p_T | d0, t_i], evaluated at `dps` digits."""
    mp.mp.dps=dps
    lam,F,G=build(n,Ne,dps)
    N=len(lam); two=mp.mpf(2*Ne)
    tab=np.zeros((n,len(tig),len(grid)))
    for ia,ti in enumerate(tig):
        tt=mp.mpf(float(ti))/two
        ei=[mp.e**(l*tt) for l in lam]
        den=[mp.fdot(zip(G[d],ei)) for d in range(n)]
        for iT,T in enumerate(grid):
            if T>=ti: continue
            tT=mp.mpf(float(T))/two
            a=[mp.e**(l*tT) for l in lam]
            b=[mp.e**(l*(tt-tT)) for l in lam]
            ab=[[a[k]*b[kp] for kp in range(N)] for k in range(N)]
            for d in range(n):
                dv=den[d]
                if dv<=floor: continue
                num=mp.fsum(v*ab[k][kp] for k,kp,v in F[d])
                if num<=floor: continue
                v=num/dv
                if 0<v<1: tab[d,ia,iT]=float(v)
    return tab
