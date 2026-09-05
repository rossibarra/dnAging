"""Estimated vs true age with +/-1 bootstrap SD error bars (Archive2, blind)."""
import re,sys,numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
T=np.array([5604,7776,5999,1799,5134,8921,7777,8918,1265,5456],float)
rows=[]
for ln in open(sys.argv[1]):
    m=re.match(r'\s*simulation_\d+\s+[\d,]+\s+[\d,]+\s+[\d,]+\s+(\d+)\s+(\d+)\s+(\d+)-(\d+)',ln)
    if m: rows.append(tuple(map(float,m.groups())))
A=np.array(rows)
if len(A)<len(T):
    sys.exit(f"only {len(A)} of {len(T)} simulations present yet")
raw,sd,lo,hi=A.T
SL,IC=1.026,748.0
cor,csd=(raw-IC)/SL, sd/SL
fig,ax=plt.subplots(figsize=(6.6,6.1)); lim=[0,13000]
ax.plot(lim,lim,color="0.35",lw=1.2,ls="--",zorder=1,label="truth ($y=x$)")
for M,S,lab,col,mk,off in ((raw,sd,"raw MAP","#c44e52","o",-80),
                           (cor,csd,"Archive1-corrected","#4c72b0","s",80)):
    e=M-T
    ax.errorbar(T+off,M,yerr=S,fmt=mk,ms=6.5,color=col,ecolor=col,elinewidth=1.5,
                capsize=4,alpha=.88,zorder=3,
                label=f"{lab}\n  RMSE {np.sqrt((e**2).mean()):.0f}   bias {e.mean():+.0f}"
                      f"   median SD {np.median(S):.0f}")
    sl,ic=np.polyfit(T,M,1); xs=np.array(lim,float)
    ax.plot(xs,sl*xs+ic,color=col,lw=1.4,alpha=.5,zorder=2)
ax.set_xlim(lim);ax.set_ylim(lim);ax.set_aspect("equal")
ax.set_xlabel("true age of ancient sample (generations)")
ax.set_ylabel("estimated age, MAP (generations)")
ax.set_title("Archive2 blind predictions, error bars are $\\pm$1 block-bootstrap SD\n"
             r"$N_e$=50,000, 10 Mb, 26 modern haplotypes, haploid ancient, $\varepsilon=10^{-6}$",
             fontsize=9.5)
ax.legend(loc="upper left",fontsize=8.5,framealpha=.95);ax.grid(alpha=.25,lw=.6)
fig.tight_layout();fig.savefig("archive2_sd.png",dpi=170)
print("wrote betabinom/archive2_sd.png")
for M,S,lab in ((raw,sd,"raw"),(cor,csd,"corrected")):
    e=M-T; z=e/S
    print(f"  {lab:<11} RMSE {np.sqrt((e**2).mean()):>5.0f}  bias {e.mean():>+6.0f}  "
          f"scatter SD {e.std(ddof=1):>5.0f}  median bootstrap SD {np.median(S):>5.0f}  "
          f"|z|>2: {(abs(z)>2).sum()}/10")
