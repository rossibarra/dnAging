"""Bias-corrected age posterior for the bigsims 100 Mb sample."""
import sys, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

d=np.load(sys.argv[1])
grid,ll,boot=d["grid"],d["loglik"],d["boot"]
nsites,ncar,Ne=int(d["nsites"]),int(d["ncarried"]),int(d["Ne"])
TRUTH=float(sys.argv[2]) if len(sys.argv)>2 else None
MAP=grid[ll.argmax()]

# pooled calibration from the 20 true-ARG simulations: MAP = sl*T + ic
T1=np.array([4840,9057,7316,2178,2142,454,8363,1995,8216,9293],float)
M1=np.array([5642,9420,7405,4130,2115,1108,10428,2619,10176,9672],float)
T2=np.array([5604,7776,5999,1799,5134,8921,7777,8918,1265,5456],float)
M2=np.array([6649,9672,7153,3627,6649,10428,9924,11435,4634,6649],float)
T=np.concatenate([T1,T2]); M=np.concatenate([M1,M2])
X=np.column_stack([np.ones(20),T]); beta,*_=np.linalg.lstsq(X,M,rcond=None)
ic,sl=beta; s2=((M-X@beta)**2).sum()/18; cov=s2*np.linalg.inv(X.T@X)
That=(MAP-ic)/sl
g=np.array([-1/sl,-(MAP-ic)/sl**2])
se=np.sqrt((boot.std(ddof=1)/sl)**2 + g@cov@g)

fig,ax=plt.subplots(figsize=(7.6,4.6))
xs=np.linspace(That-4.2*se, That+4.2*se, 900)
dens=np.exp(-0.5*((xs-That)/se)**2)/(se*np.sqrt(2*np.pi))
ax.fill_between(xs,0,dens,color="#4c72b0",alpha=.30,zorder=2)
ax.plot(xs,dens,color="#4c72b0",lw=2.1,zorder=3)
lo,hi=That-1.96*se, That+1.96*se
m=(xs>=lo)&(xs<=hi)
ax.fill_between(xs[m],0,dens[m],color="#4c72b0",alpha=.28,zorder=2)
ax.axvline(That,color="#4c72b0",lw=1.4,ls="-",zorder=4)
if TRUTH is not None:
    ax.axvline(TRUTH,color="black",lw=1.6,ls=":",zorder=5,
               label=f"true age {TRUTH:,.0f}")
ax.plot([],[],color="#4c72b0",lw=2.1,
        label=f"estimate {That:,.0f}   95% {lo:,.0f}-{hi:,.0f}")
ax.set_xlim(xs[0],xs[-1]); ax.set_ylim(0,dens.max()*1.16)
ax.set_xlabel("age of the ancient sample (generations before present)")
ax.set_ylabel("posterior density")
ax.set_title(f"bigsims simulation_01: bias-corrected age posterior\n"
             f"100 Mb, {nsites:,} sites, $N_e$={Ne:,}, 26 modern haplotypes, true ARG",
             fontsize=10)
ax.legend(loc="upper left",fontsize=9,framealpha=.95)
ax.grid(alpha=.25,lw=.6)
fig.tight_layout(); fig.savefig("bigsim_posterior.png",dpi=170)
print(f"wrote betabinom/bigsim_posterior.png")
print(f"  raw MAP        {MAP:,.0f}")
print(f"  corrected      {That:,.0f}  SE {se:.0f}   95% {lo:,.0f} - {hi:,.0f}")
if TRUTH is not None:
    print(f"  truth          {TRUTH:,.0f}   error {That-TRUTH:+,.0f} "
          f"({abs(That-TRUTH)/TRUTH*100:.1f}%)   z = {(That-TRUTH)/se:+.2f}")
