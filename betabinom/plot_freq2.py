import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
A=np.load("freq_compare2.npy"); tr,bb,df=A[:,0].astype(float),A[:,1].astype(float),A[:,2].astype(float)
rng=np.random.default_rng(0); sub=rng.choice(len(A),12000,replace=False)
fig,axes=plt.subplots(1,2,figsize=(12.4,6.1),gridspec_kw=dict(width_ratios=[1,1]))

ax=axes[0]
ax.plot([0,1],[0,1],color="0.25",lw=1.3,ls="--",zorder=5,label="truth ($y=x$)")
ax.scatter(tr[sub],df[sub],s=3,color="#c44e52",alpha=.10,lw=0,zorder=2)
ax.scatter(tr[sub],bb[sub],s=3,color="#4c72b0",alpha=.10,lw=0,zorder=3)
edges=np.linspace(0,1,31); ctr=(edges[:-1]+edges[1:])/2
for y,col,lab in ((bb,"#4c72b0","betabinom  $\\varphi(k,n_T,a)$"),
                  (df,"#c44e52","diffusion  $E[p_T\\mid d_0,t_i]$")):
    mu=np.array([y[(tr>=a_)&(tr<b_)].mean() if ((tr>=a_)&(tr<b_)).sum()>50 else np.nan
                 for a_,b_ in zip(edges[:-1],edges[1:])])
    sd=np.array([y[(tr>=a_)&(tr<b_)].std() if ((tr>=a_)&(tr<b_)).sum()>50 else np.nan
                 for a_,b_ in zip(edges[:-1],edges[1:])])
    ax.plot(ctr,mu,color=col,lw=2.4,zorder=6); ax.fill_between(ctr,mu-sd,mu+sd,color=col,alpha=.16,lw=0,zorder=1)
    e=y-tr
    ax.plot([],[],color=col,lw=2.4,label=f"{lab}\n   bias {e.mean():+.4f}  "
            f"RMSE {np.sqrt((e**2).mean()):.4f}  r {np.corrcoef(tr,y)[0,1]:.3f}")
ax.set_xlim(0,1);ax.set_ylim(0,1);ax.set_aspect("equal")
ax.set_xlabel("true derived frequency at $T$  (400 reference chromosomes at $T$)")
ax.set_ylabel("estimated derived frequency at $T$")
ax.set_title("full range",fontsize=10)
ax.legend(loc="upper left",fontsize=8,framealpha=.95); ax.grid(alpha=.22,lw=.6)

ax=axes[1]
lim=0.15
ax.plot([0,lim],[0,lim],color="0.25",lw=1.3,ls="--",zorder=5)
e2=np.linspace(0,lim,16); c2=(e2[:-1]+e2[1:])/2
for y,col,lab in ((bb,"#4c72b0","betabinom"),(df,"#c44e52","diffusion")):
    mu=np.array([y[(tr>=a_)&(tr<b_)].mean() if ((tr>=a_)&(tr<b_)).sum()>50 else np.nan
                 for a_,b_ in zip(e2[:-1],e2[1:])])
    se=np.array([y[(tr>=a_)&(tr<b_)].std()/np.sqrt(max(((tr>=a_)&(tr<b_)).sum(),1))
                 for a_,b_ in zip(e2[:-1],e2[1:])])
    ax.plot(c2,mu,color=col,lw=2.4,marker='o',ms=4,zorder=6)
    ax.fill_between(c2,mu-1.96*se,mu+1.96*se,color=col,alpha=.30,lw=0,zorder=1)
    m=tr<0.1; ee=(y-tr)[m]
    ax.plot([],[],color=col,lw=2.4,label=f"{lab}   true<0.1:\n   bias {ee.mean():+.4f}  "
            f"RMSE {np.sqrt((ee**2).mean()):.4f}")
ax.set_xlim(0,lim);ax.set_ylim(0,lim);ax.set_aspect("equal")
ax.set_xlabel("true derived frequency at $T$"); ax.set_ylabel("binned mean estimate")
ax.set_title("rare end, where the age signal lives\n(bands are 95% CI on the binned mean)",fontsize=10)
ax.legend(loc="upper left",fontsize=8.5,framealpha=.95); ax.grid(alpha=.22,lw=.6)

fig.suptitle(f"Estimated vs true allele frequency at time $T$: {len(A):,} (SNP,$T$) pairs, "
             f"100 x 1 Mb sims, $N_e$=10,000, 20-haplotype panel, 5 times each",fontsize=10.5)
fig.tight_layout(rect=[0,0,1,0.955]); fig.savefig("freq_compare2.png",dpi=165)
print("wrote betabinom/freq_compare2.png")
for lo,hi in ((0,.02),(.02,.05),(.05,.1),(.1,.3),(.3,.7),(.7,1.01)):
    m=(tr>=lo)&(tr<hi)
    print(f"  true [{lo:.2f},{hi:.2f})  n={m.sum():>7,}   "
          f"bb bias {(bb-tr)[m].mean():+.4f}   df bias {(df-tr)[m].mean():+.4f}")
