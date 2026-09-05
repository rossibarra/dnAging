"""Is the small-p over-prediction model error or numerical error?
Compare the production settings (NQ=24, 70-age table, log-interp) against a
much finer quadrature and a much finer table, on the same sites."""
import numpy as np, msprime, tskit
from phi import Phi
NE,NMOD,L,RR,MU = 100_000,26,5e5,1e-8,1e-8
TIMES=[1000.,2500.]
coarse=Phi(np.geomspace(0.5,2e6,70),26,NE)
fine  =Phi(np.geomspace(0.05,2e6,240),26,NE)
def integ(Pt,k,nT,lo,tp,tc,T,NQ):
    nd=np.linspace(lo,tp,NQ); w=np.full(NQ,1.); w[0]=w[-1]=.5
    return float((Pt(k,nT,nd-T)*w).sum()/w.sum()*(tp-lo)/(tp-tc))
ts=msprime.sim_ancestry(samples=[msprime.SampleSet(NMOD,time=0,ploidy=1)]+
    [msprime.SampleSet(1,time=t,ploidy=1) for t in TIMES],population_size=NE,
    sequence_length=L,recombination_rate=RR,ploidy=2,random_seed=77)
ts=msprime.sim_mutations(ts,rate=MU,random_seed=78,discrete_genome=False)
tm=ts.simplify(samples=list(ts.samples())[:NMOD])
rec=[]
for tree in tm.trees():
    internal=np.sort([tree.time(u) for u in tree.nodes() if tree.is_internal(u)])
    cache={}
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
        for T in TIMES:
            nT=NMOD-np.searchsorted(internal,T,side='right')
            if nT<2 or T>=tp: continue
            if T>=tc: k,lo=1,T
            else:
                k=nl-np.searchsorted(cache[m.node],T,side='right'); lo=tc
                if k<1 or k>=nT: continue
            a=integ(coarse,k,nT,lo,tp,tc,T,24)      # production
            b=integ(fine,  k,nT,lo,tp,tc,T,400)     # reference
            rec.append((a,b))
R=np.array(rec); a,b=R.T
print(f"{len(R):,} site-T evaluations\n")
print(f"{'p range (reference)':<22}{'n':>8}{'sum coarse':>12}{'sum fine':>11}{'coarse/fine':>13}")
for lo,hi in [(0,1e-3),(1e-3,3e-3),(3e-3,1e-2),(1e-2,3e-2),(3e-2,0.1),(0.1,0.4),(0.4,1.01)]:
    s=(b>=lo)&(b<hi)
    if s.sum()<5: continue
    print(f"[{lo:g},{hi:g})".ljust(22)+f"{s.sum():>8,}{a[s].sum():>12.2f}"
          f"{b[s].sum():>11.2f}{a[s].sum()/b[s].sum():>13.4f}")
print(f"\nTOTAL{'':17}{len(R):>8,}{a.sum():>12.1f}{b.sum():>11.1f}{a.sum()/b.sum():>13.4f}")
