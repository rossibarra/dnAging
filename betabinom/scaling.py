"""Does more data fix it?  Same model, increasing sequence length."""
import importlib.util, numpy as np, msprime, tskit
spec=importlib.util.spec_from_file_location("pre",
    "/Users/jeffreyross-ibarra/src/dnAging/precompute_freq_trajectory_moments.py")
pre=importlib.util.module_from_spec(spec); spec.loader.exec_module(pre)
NE,NMOD,RR,MU,EPS,SEED,T_TRUE=100_000,26,1e-8,1e-8,1e-3,11,2500.
grid=np.linspace(200.,9000.,60)
tau=lambda t:t/(2.0*NE); eps0=1.0/(2.0*NE)
tig=np.geomspace(30.,6e5,36); eng=pre.MomentEngine(NMOD)
tab=np.zeros((NMOD,len(tig),len(grid)))
for id0,d0 in enumerate(range(1,NMOD+1)):
    for ia,ti in enumerate(tig):
        for iT,T in enumerate(grid):
            if T<ti: tab[id0,ia,iT]=eng.Emoments(d0,tau(ti),tau(T),eps0)[0]
tab=np.nan_to_num(tab); ltig=np.log(tig)
def dp(d0,tm_,iT):
    x=np.interp(np.log(np.clip(tm_,tig[0],tig[-1])),ltig,np.arange(len(tig)))
    i0=int(np.clip(x,0,len(tig)-2)); w=x-i0
    return (1-w)*tab[d0-1,i0,iT]+w*tab[d0-1,i0+1,iT]

print(f"truth = {T_TRUE:g}\n{'seq len':>9}{'sites':>9}{'MAP':>8}{'95% CI':>16}{'CI width':>10}{'covers':>8}")
for L in (0.5e6, 1e6, 2e6, 4e6):
    ts=msprime.sim_ancestry(samples=[msprime.SampleSet(NMOD,time=0,ploidy=1),
         msprime.SampleSet(1,time=T_TRUE,ploidy=1)],population_size=NE,sequence_length=L,
         recombination_rate=RR,ploidy=2,random_seed=SEED)
    ts=msprime.sim_mutations(ts,rate=MU,random_seed=SEED+7,discrete_genome=False)
    anc=ts.num_samples-1; ar=list(ts.samples()).index(anc)
    gt={v.site.position:int(v.genotypes[ar]!=0) for v in ts.variants()}
    tm=ts.simplify(samples=[s for s in ts.samples() if s!=anc])
    P=[];G=[]
    for tree in tm.trees():
        internal=np.sort([tree.time(u) for u in tree.nodes() if tree.is_internal(u)])
        nT=NMOD-np.searchsorted(internal,grid,side='right'); cache={}
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
            k=nl-np.searchsorted(cache[m.node],grid,side='right')
            p=np.where(grid<tc,k/(nT+1.0),0.0)
            for iT in np.flatnonzero((grid>=tc)&(grid<tp)):
                nq=6; nd=np.linspace(grid[iT],tp,nq); w=np.full(nq,1.);w[0]=w[-1]=.5
                p[iT]=(np.array([dp(nl,t_,iT) for t_ in nd])*w).sum()/w.sum()*(tp-grid[iT])/(tp-tc)
            P.append(p.astype(np.float32)); G.append(gt.get(site.position,0))
    P=np.array(P); G=np.array(G,bool)
    r=np.clip(EPS+(1-2*EPS)*P,1e-15,1-1e-15)
    ll=np.log(r[G]).sum(0)+np.log1p(-r[~G]).sum(0)
    post=np.exp(ll-ll.max()); post/=np.trapezoid(post,grid)
    cdf=np.concatenate([[0],np.cumsum((post[:-1]+post[1:])/2*np.diff(grid))]); cdf/=cdf[-1]
    q=lambda a: float(np.interp(a,cdf,grid))
    lo,hi=q(.025),q(.975)
    print(f"{L/1e6:>7.1f}Mb{P.shape[0]:>9,}{grid[ll.argmax()]:>8.0f}"
          f"{lo:>8.0f}-{hi:<7.0f}{hi-lo:>10.0f}"
          f"{('YES' if lo<=T_TRUE<=hi else 'no'):>8}")
