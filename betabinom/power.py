"""If the model were EXACTLY right, how much could the data distinguish T?
   Expected loglik penalty for using T instead of T_true:
       D(T) = sum_i KL( Bern(p_i(T_true)) || Bern(p_i(T)) )
   D(T)=2 is roughly a 2-unit loglik drop, i.e. the edge of a 95% interval."""
import importlib.util, numpy as np, msprime, tskit
spec=importlib.util.spec_from_file_location("pre",
    "/Users/jeffreyross-ibarra/src/dnAging/precompute_freq_trajectory_moments.py")
pre=importlib.util.module_from_spec(spec); spec.loader.exec_module(pre)
NMOD,L,RR,MU,EPS,SEED,T_TRUE=26,1e6,1e-8,1e-8,1e-3,11,2500.
grid=np.linspace(200.,9000.,45)

def run(NE):
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
    ts=msprime.sim_ancestry(samples=[msprime.SampleSet(NMOD,time=0,ploidy=1),
         msprime.SampleSet(1,time=T_TRUE,ploidy=1)],population_size=NE,sequence_length=L,
         recombination_rate=RR,ploidy=2,random_seed=SEED)
    ts=msprime.sim_mutations(ts,rate=MU,random_seed=SEED+7,discrete_genome=False)
    anc=ts.num_samples-1
    tm=ts.simplify(samples=[s for s in ts.samples() if s!=anc])
    P=[]
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
            P.append(p.astype(np.float64))
    P=np.array(P)
    R=np.clip(EPS+(1-2*EPS)*P,1e-12,1-1e-12)
    i0=int(np.argmin(abs(grid-T_TRUE))); q=R[:,i0][:,None]
    KL=(q*np.log(q/R)+(1-q)*np.log((1-q)/(1-R))).sum(0)
    return grid,KL,P.shape[0]

for NE in (100_000, 10_000):
    g,KL,n=run(NE)
    print(f"\nNe={NE:,}   {n:,} sites   tau at T_true = {T_TRUE/(2*NE):.4f}")
    print(f"  expected loglik penalty D(T) for the WRONG T, if the model were exact:")
    for T in (500,1000,1500,2000,3000,4000,6000,9000):
        j=int(np.argmin(abs(g-T)))
        print(f"     T={T:>5}   D = {KL[j]:>9.1f}")
    band=g[KL<2.0]
    if len(band): print(f"  ==> T values indistinguishable from truth (D<2): "
                        f"{band.min():.0f} - {band.max():.0f}")
