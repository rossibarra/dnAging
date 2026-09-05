"""Estimated vs true age for the ten supplied simulations, both eps settings."""
import re, sys, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

def parse(path):
    rows=[]
    for ln in open(path):
        m=re.match(r'\s*(simulation_\d+)\s+(\d+)\s+(\d+)\s+(\d+)-(\d+)', ln)
        if m: rows.append(tuple(map(float,m.groups()[1:])))
    return np.array(rows)          # T_true, MAP, blo, bhi

runs=[(sys.argv[1], r"$\varepsilon=10^{-3}$", "#c44e52", "o"),
      (sys.argv[2], r"$\varepsilon=10^{-6}$", "#4c72b0", "s")]
fig,ax=plt.subplots(figsize=(6.4,6.0))
lim=[0,14000]
ax.plot(lim,lim,color="0.35",lw=1.2,ls="--",zorder=1,label="truth ($y=x$)")
for path,lab,col,mk in runs:
    A=parse(path)
    if not len(A): continue
    T,M,lo,hi=A.T
    off = -90 if mk=="o" else 90
    ax.errorbar(T+off,M,yerr=[M-lo,hi-M],fmt=mk,ms=6,color=col,ecolor=col,
                elinewidth=1.3,capsize=3,alpha=.85,zorder=3,
                label=f"{lab}   bias {np.mean(M-T):+.0f}")
    sl,ic=np.polyfit(T,M,1)
    xs=np.array(lim,float)
    ax.plot(xs,sl*xs+ic,color=col,lw=1.4,alpha=.55,zorder=2)
ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
ax.set_xlabel("true age of ancient sample (generations)")
ax.set_ylabel("estimated age, MAP (generations)")
ax.set_title("betabinom likelihood on 10 true-ARG simulations\n"
             r"$N_e$=50,000, 10 Mb, 26 modern haplotypes; bars are block-bootstrap 95%",
             fontsize=10)
ax.legend(loc="upper left",fontsize=9,framealpha=.95)
ax.grid(alpha=.25,lw=.6)
fig.tight_layout(); fig.savefig("age_vs_truth.png",dpi=170)
print("wrote betabinom/age_vs_truth.png")
for path,lab,_,_ in runs:
    A=parse(path)
    if not len(A): continue
    T,M,lo,hi=A.T; sl,ic=np.polyfit(T,M,1)
    e=M-T
    print(f"  {lab:<22} n={len(T)}  RMSE {np.sqrt((e**2).mean()):>5.0f}  MAE {np.abs(e).mean():>5.0f}  "
          f"bias {e.mean():>+6.0f}  slope {sl:.3f}  int {ic:+.0f}  "
          f"r={np.corrcoef(T,M)[0,1]:.3f}  cov {int(((lo<=T)&(T<=hi)).sum())}/{len(T)}")
