"""phi(k, n_T, a) = E[X(T) | k of n_T panel lineages carry it at T,
                            mutation age a generations at T]

    phi = sum_j C(n_T-k,j)(-1)^j M_{k+1+j}(tau_a)
        / sum_j C(n_T-k,j)(-1)^j M_{k+j}(tau_a)

M_m(tau) are the neutral-diffusion moments from a single new copy.  The
alternating sums cancel catastrophically, so the table is built in mpmath
(exact moment map, dps=60) and only then rounded to float64.
"""
import os, sys, numpy as np, mpmath as mp
from math import comb
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "tests"))
from _reference import C_matrix

def build(ages, max_n, Ne, dps=60):
    mp.mp.dps=dps
    eps=mp.mpf(1)/(2*Ne); two_ne=mp.mpf(2*Ne)
    sz=max_n+3
    tab=np.full((max_n+2, max_n+2, len(ages)), np.nan)   # [k, n_T, age]
    for ia,a in enumerate(ages):
        u=mp.mpf(a)/two_ne
        C=C_matrix(max_n,u)
        m0=[mp.mpf(1)]+[eps**m for m in range(1,sz)]
        M=[mp.fsum(C[m][j]*m0[j] for j in range(sz)) for m in range(sz)]
        for nT in range(2,max_n+1):
            for k in range(1,nT):
                num=mp.fsum(mp.mpf(comb(nT-k,j))*(-1)**j*M[k+1+j] for j in range(nT-k+1))
                den=mp.fsum(mp.mpf(comb(nT-k,j))*(-1)**j*M[k+j]   for j in range(nT-k+1))
                if den>0 and num>0:
                    v=num/den
                    if 0<v<1: tab[k,nT,ia]=float(v)
    return tab

class Phi:
    def __init__(self, ages, max_n, Ne):
        self.ages=np.asarray(ages,float); self.la=np.log(self.ages)
        self.tab=build(self.ages,max_n,Ne); self.max_n=max_n
    def __call__(self, k, nT, a):
        """a may be an array; log-interpolated in age, clamped at the ends."""
        k=int(np.clip(k,1,nT-1)); nT=int(np.clip(nT,2,self.max_n))
        row=self.tab[k,nT]
        good=np.isfinite(row)
        return np.interp(np.log(np.clip(a,self.ages[0],self.ages[-1])),
                         self.la[good], row[good])

if __name__=="__main__":
    NE=100_000
    ages=np.geomspace(0.5, 2e6, 70)
    P=Phi(ages,26,NE)
    print(f"table built: {P.tab.shape}, finite entries "
          f"{np.isfinite(P.tab).sum():,}/{P.tab.size:,}")
    print(f"\nNe={NE:,}  1/(2Ne)={1/(2*NE):.2e}")
    print(f"{'age':>10}" + "".join(f"{f'k={k},n=24':>13}" for k in (1,3,8)))
    for a in (1,10,100,1000,1e4,1e5,1e6):
        print(f"{a:>10,.0f}" + "".join(f"{P(k,24,a):>13.6f}" for k in (1,3,8)))
    print(f"{'k/(nT+1)':>10}" + "".join(f"{k/25:>13.6f}" for k in (1,3,8)))
    np.savez("phi_table.npz", tab=P.tab, ages=P.ages, Ne=NE)
    print("\nsaved phi_table.npz")
