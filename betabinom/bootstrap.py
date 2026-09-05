"""Block bootstrap over genomic windows: honest interval for the age posterior.

The composite likelihood treats linked sites as independent, so its curvature
overstates precision.  Resampling contiguous blocks preserves LD within a block
and gives the sampling distribution of the estimate.
"""
import numpy as np, msprime, tskit, sys
from phi import Phi
NE,NMOD,L,RR,MU,EPS = 100_000,26,2e6,1e-8,1e-8,1e-3
NBLOCK,NBOOT = 100, 400
grid=np.linspace(200.,9000.,60); NQ=24
P=Phi(np.geomspace(0.5,2e6,70),26,NE)
W=np.full(NQ,1.); W[0]=W[-1]=.5; WS=W.sum()

def per_site(T_TRUE,SEED):
    ts=msprime.sim_ancestry(samples=[msprime.SampleSet(NMOD,time=0,ploidy=1),
        msprime.SampleSet(1,time=T_TRUE,ploidy=1)],population_size=NE,sequence_length=L,
        recombination_rate=RR,ploidy=2,random_seed=SEED)
    ts=msprime.sim_mutations(ts,rate=MU,random_seed=SEED+7,discrete_genome=False)
    anc=ts.num_samples-1; ar=list(ts.samples()).index(anc)
    gt={v.site.position:int(v.genotypes[ar]!=0) for v in ts.variants()}
    tm=ts.simplify(samples=[s for s in ts.samples() if s!=anc])
    LL=[]; POS=[]
    for tree in tm.trees():
        internal=np.sort([tree.time(u) for u in tree.nodes() if tree.is_internal(u)])
        nTv=NMOD-np.searchsorted(internal,grid,side='right'); cache={}
        for site in tree.sites():
            if len(site.mutations)!=1: continue
            m=site.mutations[0]; par=tree.parent(m.node)
            if par==tskit.NULL: continue
            nl=tree.num_samples(m.node)
            if not (0<nl<NMOD): continue
            tc,tp=tree.time(m.node),tree.time(par)
            if m.node not in cache:
                cache[m.node]=np.sort([tree.time(u) for u in tree.nodes(root=m.node)
                                       if tree.is_internal(u)])
            kv=nl-np.searchsorted(cache[m.node],grid,side='right')
            p=np.zeros(len(grid))
            iA=np.flatnonzero((grid<tc)&(nTv>=2)&(kv>=1)&(kv<nTv))
            if len(iA):
                nd=np.linspace(tc,tp,NQ)
                for iT in iA: p[iT]=(P(kv[iT],nTv[iT],nd-grid[iT])*W).sum()/WS
            for iT in np.flatnonzero((grid>=tc)&(grid<tp)&(nTv>=2)):
                nd=np.linspace(grid[iT],tp,NQ)
                p[iT]=(P(1,nTv[iT],nd-grid[iT])*W).sum()/WS*(tp-grid[iT])/(tp-tc)
            r=np.clip(EPS+(1-2*EPS)*p,1e-15,1-1e-15)
            LL.append((np.log(r) if gt.get(site.position,0) else np.log1p(-r)).astype(np.float32))
            POS.append(site.position)
    return np.array(LL), np.array(POS)

def ci(ll):
    post=np.exp(ll-ll.max()); post/=np.trapezoid(post,grid)
    cdf=np.concatenate([[0],np.cumsum((post[:-1]+post[1:])/2*np.diff(grid))]); cdf/=cdf[-1]
    return float(np.interp(.025,cdf,grid)), float(np.interp(.975,cdf,grid))

print(f"{'seed':>5}{'T_true':>8}{'MAP':>7}   {'composite 95%':>18}{'cov':>5}   "
      f"{'block-bootstrap 95%':>22}{'cov':>5}")
rng=np.random.default_rng(0)
for T_TRUE,SEED in ((2500.,11),(2500.,22),(2500.,33),(1000.,44),(6000.,55)):
    LL,POS=per_site(T_TRUE,SEED)
    tot=LL.sum(0); mapv=grid[tot.argmax()]
    clo,chi=ci(tot)
    blk=np.minimum((POS/ (L/NBLOCK)).astype(int), NBLOCK-1)
    sums=np.zeros((NBLOCK,len(grid)))
    for b in range(NBLOCK):
        s=blk==b
        if s.any(): sums[b]=LL[s].sum(0)
    maps=np.empty(NBOOT)
    for j in range(NBOOT):
        idx=rng.integers(0,NBLOCK,NBLOCK)
        maps[j]=grid[sums[idx].sum(0).argmax()]
    blo,bhi=np.percentile(maps,[2.5,97.5])
    print(f"{SEED:>5}{T_TRUE:>8.0f}{mapv:>7.0f}   {clo:>8.0f}-{chi:<9.0f}"
          f"{('YES' if clo<=T_TRUE<=chi else 'no'):>5}   "
          f"{blo:>10.0f}-{bhi:<11.0f}{('YES' if blo<=T_TRUE<=bhi else 'no'):>5}")
