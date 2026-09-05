import numpy as np, msprime, tskit, sys
from phi import Phi
NE,NMOD,L,RR,MU,EPS = 100_000,26,1e6,1e-8,1e-8,1e-3
grid=np.linspace(200.,9000.,60); NQ=24
P=Phi(np.geomspace(0.5,2e6,70),26,NE)
W=np.full(NQ,1.0); W[0]=W[-1]=0.5; WS=W.sum()

def run(T_TRUE,SEED):
    ts=msprime.sim_ancestry(samples=[msprime.SampleSet(NMOD,time=0,ploidy=1),
        msprime.SampleSet(1,time=T_TRUE,ploidy=1)],population_size=NE,sequence_length=L,
        recombination_rate=RR,ploidy=2,random_seed=SEED)
    ts=msprime.sim_mutations(ts,rate=MU,random_seed=SEED+7,discrete_genome=False)
    anc=ts.num_samples-1; ar=list(ts.samples()).index(anc)
    gt={v.site.position:int(v.genotypes[ar]!=0) for v in ts.variants()}
    tm=ts.simplify(samples=[s for s in ts.samples() if s!=anc])
    ll=np.zeros(len(grid)); n=0
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
            # edges above T: nodes fixed on (tc,tp), k varies with T
            iA=np.flatnonzero((grid<tc)&(nTv>=2))
            if len(iA):
                nodes=np.linspace(tc,tp,NQ)
                for k in np.unique(kv[iA]):
                    sel=iA[(kv[iA]==k)&(kv[iA]<nTv[iA])]
                    if not len(sel): continue
                    ages=nodes[None,:]-grid[sel][:,None]
                    for r,iT in enumerate(sel):
                        p[iT]=(P(k,nTv[iT],ages[r])*W).sum()/WS
            # straddling: k=1, nodes on (T,tp)
            for iT in np.flatnonzero((grid>=tc)&(grid<tp)&(nTv>=2)):
                nd=np.linspace(grid[iT],tp,NQ)
                p[iT]=(P(1,nTv[iT],nd-grid[iT])*W).sum()/WS*(tp-grid[iT])/(tp-tc)
            r=np.clip(EPS+(1-2*EPS)*p,1e-15,1-1e-15)
            ll+= np.log(r) if gt.get(site.position,0) else np.log1p(-r)
            n+=1
    post=np.exp(ll-ll.max()); post/=np.trapezoid(post,grid)
    cdf=np.concatenate([[0],np.cumsum((post[:-1]+post[1:])/2*np.diff(grid))]); cdf/=cdf[-1]
    q=lambda a: float(np.interp(a,cdf,grid))
    lo,hi=q(.025),q(.975)
    return n,grid[ll.argmax()],lo,hi,('YES' if lo<=T_TRUE<=hi else 'no'),ll,post

print(f"{'seed':>5}{'T_true':>8}{'sites':>8}{'MAP':>7}{'95% CI':>15}{'covers':>8}")
res=None
for T_TRUE,SEED in ((2500.,11),(2500.,22),(2500.,33),(1000.,44),(6000.,55)):
    n,mapv,lo,hi,cov,ll,post=run(T_TRUE,SEED)
    if SEED==11: res=(grid,post)
    print(f"{SEED:>5}{T_TRUE:>8.0f}{n:>8,}{mapv:>7.0f}{lo:>7.0f}-{hi:<7.0f}{cov:>8}")
g,post=res
sel=post>post.max()*1e-3; gs,ps=g[sel],post[sel]
print("\nposterior, seed 11 (truth 2500)")
for row in range(12,-1,-1):
    thr=ps.max()*row/12
    print("  |"+''.join('#' if v>=thr else ' ' for v in
        np.interp(np.linspace(gs[0],gs[-1],72),gs,ps)))
print("  +"+"-"*72)
print(f"   {gs[0]:<24.0f}{(gs[0]+gs[-1])/2:^24.0f}{gs[-1]:>24.0f}  generations")
