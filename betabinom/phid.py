"""phi with its denominator exposed, so the age integral can be done correctly.

  den(k,n_T,a) = sum_j C(n_T-k,j)(-1)^j M_{k+j}(tau_a)     ( ∝ P(k of n_T | a) )
  num(k,n_T,a) = sum_j C(n_T-k,j)(-1)^j M_{k+1+j}(tau_a)
  phi = num/den

Marginalising the mutation age over an edge must weight by den, not uniformly:
  p = sum_q w_q num(a_q) / sum_q w_q den(a_q)
"""
import os, sys, numpy as np, mpmath as mp
from math import comb
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "tests"))
from _reference import C_matrix

class PhiD:
    def __init__(self, ages, max_n, Ne, dps=60):
        mp.mp.dps=dps
        self.ages=np.asarray(ages,float); self.la=np.log(self.ages); self.max_n=max_n
        eps=mp.mpf(1)/(2*Ne); sz=max_n+3
        self.num=np.full((max_n+2,max_n+2,len(ages)),np.nan)
        self.den=np.full((max_n+2,max_n+2,len(ages)),np.nan)
        for ia,a in enumerate(self.ages):
            C=C_matrix(max_n, mp.mpf(a)/(2*Ne))
            m0=[mp.mpf(1)]+[eps**m for m in range(1,sz)]
            M=[mp.fsum(C[m][j]*m0[j] for j in range(sz)) for m in range(sz)]
            for nT in range(2,max_n+1):
                for k in range(1,nT):
                    nu=mp.fsum(mp.mpf(comb(nT-k,j))*(-1)**j*M[k+1+j] for j in range(nT-k+1))
                    de=mp.fsum(mp.mpf(comb(nT-k,j))*(-1)**j*M[k+j]   for j in range(nT-k+1))
                    if de>0 and nu>0 and 0<nu/de<1:
                        self.num[k,nT,ia]=float(nu); self.den[k,nT,ia]=float(de)
    def _row(self,arr,k,nT):
        k=int(np.clip(k,1,nT-1)); nT=int(np.clip(nT,2,self.max_n))
        r=arr[k,nT]; g=np.isfinite(r); return self.la[g], r[g]
    def nd(self,k,nT,a):
        """(num, den) at ages a, interpolated in log-age and log-value."""
        la=np.log(np.clip(a,self.ages[0],self.ages[-1]))
        x1,y1=self._row(self.num,k,nT); x2,y2=self._row(self.den,k,nT)
        return (np.exp(np.interp(la,x1,np.log(y1))),
                np.exp(np.interp(la,x2,np.log(y2))))
    def integrate(self, k, nT, ages, w):
        """den-weighted marginal over the mutation age: sum(w*num)/sum(w*den).

        `ages` may be 1-D (one edge) or 2-D (rows = grid points, cols = quadrature
        nodes); a 2-D call returns one value per row, which keeps the hot loop
        vectorised."""
        nu,de=self.nd(k,nT,ages)
        num=(w*nu).sum(-1); den=(w*de).sum(-1)
        return np.where(den>0, num/np.where(den>0,den,1.0), 0.0)

if __name__=="__main__":
    NE=50_000
    P=PhiD(np.geomspace(0.5,4e6,80),26,NE)
    W=np.full(12,1.); W[0]=W[-1]=.5
    print(f"Ne={NE:,}   edge from T to T+A, k=1, n_T=24\n")
    print(f"{'edge span A':>13}{'uniform avg':>14}{'den-weighted':>15}{'ratio':>8}")
    for A in (200.,1000.,5000.,20000.,100000.):
        ages=np.linspace(0.5,A,12)
        nu,de=P.nd(1,24,ages)
        unif=float((W*(nu/de)).sum()/W.sum())
        corr=float((W*nu).sum()/(W*de).sum())
        print(f"{A:>13,.0f}{unif:>14.6f}{corr:>15.6f}{corr/unif:>8.2f}")
