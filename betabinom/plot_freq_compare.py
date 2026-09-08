import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
A=np.load("freq_compare.npy"); tr,bb,df=A[:,0],A[:,1],A[:,2]
rng=np.random.default_rng(0)
sub=rng.choice(len(A),min(9000,len(A)),replace=False)

fig,ax=plt.subplots(figsize=(7.4,7.0))
ax.plot([0,1],[0,1],color="0.25",lw=1.3,ls="--",zorder=5,label="truth ($y=x$)")
ax.scatter(tr[sub],df[sub],s=5,color="#c44e52",alpha=.16,lw=0,zorder=2)
ax.scatter(tr[sub],bb[sub],s=5,color="#4c72b0",alpha=.16,lw=0,zorder=3)

edges=np.linspace(0,1,26); ctr=(edges[:-1]+edges[1:])/2
for y,col,lab in ((bb,"#4c72b0","betabinom  $\\varphi(k,n_T,a)$"),
                  (df,"#c44e52","diffusion  $E[p_T\\mid d_0,t_i]$")):
    mu=np.array([y[(tr>=a)&(tr<b)].mean() if ((tr>=a)&(tr<b)).sum()>20 else np.nan
                 for a,b in zip(edges[:-1],edges[1:])])
    sd=np.array([y[(tr>=a)&(tr<b)].std() if ((tr>=a)&(tr<b)).sum()>20 else np.nan
                 for a,b in zip(edges[:-1],edges[1:])])
    ax.plot(ctr,mu,color=col,lw=2.4,zorder=6)
    ax.fill_between(ctr,mu-sd,mu+sd,color=col,alpha=.18,lw=0,zorder=1)
    e=y-tr
    ax.plot([],[],color=col,lw=2.4,
            label=f"{lab}\n   bias {e.mean():+.4f}  RMSE {np.sqrt((e**2).mean()):.4f}"
                  f"  r {np.corrcoef(tr,y)[0,1]:.3f}")
ax.set_xlim(0,1); ax.set_ylim(0,1); ax.set_aspect("equal")
ax.set_xlabel("true derived frequency at $T$  (200 reference chromosomes sampled at $T$)")
ax.set_ylabel("estimated derived frequency at $T$")
ax.set_title("Estimated vs true allele frequency at time $T$\n"
             f"{len(A):,} (SNP, $T$) pairs: 296 SNPs x 100 random times, "
             f"$N_e$=50,000, 26-haplotype panel", fontsize=10)
ax.legend(loc="upper left",fontsize=8.5,framealpha=.95)
ax.grid(alpha=.22,lw=.6)
fig.tight_layout(); fig.savefig("freq_compare.png",dpi=170)
print("wrote betabinom/freq_compare.png")
lo=tr<0.1
for y,lab in ((bb,"betabinom"),(df,"diffusion")):
    e=y-tr
    print(f"  {lab:<11} all: bias {e.mean():+.4f} RMSE {np.sqrt((e**2).mean()):.4f}"
          f"   |  true<0.1: bias {e[lo].mean():+.4f} RMSE {np.sqrt((e[lo]**2).mean()):.4f}")
