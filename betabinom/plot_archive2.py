"""Blind Archive2 predictions vs truth, raw and Archive1-corrected."""
import re,sys,numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
T=np.array([5604,7776,5999,1799,5134,8921,7777,8918,1265,5456],float)
rows=[]
for ln in open(sys.argv[1]):
    m=re.match(r'\s*simulation_\d+\s+[\d,]+\s+[\d,]+\s+[\d,]+\s+(\d+)\s+(\d+)-(\d+)',ln)
    if m: rows.append(tuple(map(float,m.groups())))
A=np.array(rows); raw,lo,hi=A.T
SL,IC=1.026,748.0
cor,clo,chi=(raw-IC)/SL,(lo-IC)/SL,(hi-IC)/SL
fig,ax=plt.subplots(figsize=(6.4,6.0)); lim=[0,13000]
ax.plot(lim,lim,color="0.35",lw=1.2,ls="--",zorder=1,label="truth ($y=x$)")
for M,l,h,lab,col,mk,off in ((raw,lo,hi,"raw MAP","#c44e52","o",-80),
                             (cor,clo,chi,"Archive1-corrected","#4c72b0","s",80)):
    e=M-T
    ax.errorbar(T+off,M,yerr=[M-l,h-M],fmt=mk,ms=6,color=col,ecolor=col,elinewidth=1.3,
                capsize=3,alpha=.85,zorder=3,
                label=f"{lab}\n  RMSE {np.sqrt((e**2).mean()):.0f}, bias {e.mean():+.0f}, "
                      f"cov {int(((l<=T)&(T<=h)).sum())}/10")
    sl,ic=np.polyfit(T,M,1); xs=np.array(lim,float)
    ax.plot(xs,sl*xs+ic,color=col,lw=1.4,alpha=.5,zorder=2)
ax.set_xlim(lim);ax.set_ylim(lim);ax.set_aspect("equal")
ax.set_xlabel("true age (generations)");ax.set_ylabel("estimated age, MAP (generations)")
ax.set_title("Archive2: blind predictions vs truth\n"
             r"$N_e$=50,000, 10 Mb, panel-only ARG, haploid ancient; $\varepsilon=10^{-6}$",fontsize=10)
ax.legend(loc="upper left",fontsize=8.5,framealpha=.95);ax.grid(alpha=.25,lw=.6)
fig.tight_layout();fig.savefig("archive2_blind.png",dpi=170)
print("wrote betabinom/archive2_blind.png")
