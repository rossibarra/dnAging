"""Does the mutant clade's shape explain the residual?

At fixed (k, n_T) the model gives one p.  But k lineages at T descending to
`nleaf` leaves have coalesced fast if nleaf/k is large -- under the structured
coalescent the mutant class coalesces at rate ~ 1/X, so fast coalescence implies
a SMALL frequency.  If that is real, carriage should fall with nleaf/k at fixed
(k, n_T).
"""
import numpy as np, msprime, tskit
from collections import defaultdict
from phi import Phi
NE,NMOD,L,RR,MU = 100_000,26,2e6,1e-8,1e-8
TIMES=[500.,1000.,2000.,2500.,4000.]
NQ=24
P=Phi(np.geomspace(0.5,2e6,70),26,NE)
W=np.full(NQ,1.); W[0]=W[-1]=.5; WS=W.sum()

rows=[]
for rep in range(8):
    ss=[msprime.SampleSet(NMOD,time=0,ploidy=1)]+[msprime.SampleSet(1,time=t,ploidy=1) for t in TIMES]
    ts=msprime.sim_ancestry(samples=ss,population_size=NE,sequence_length=L,
        recombination_rate=RR,ploidy=2,random_seed=1200+rep)
    ts=msprime.sim_mutations(ts,rate=MU,random_seed=6100+rep,discrete_genome=False)
    smp=list(ts.samples()); r_={t:smp.index(smp[NMOD+i]) for i,t in enumerate(TIMES)}
    gts={v.site.position:v.genotypes for v in ts.variants()}
    tm=ts.simplify(samples=smp[:NMOD])
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
            g=gts.get(site.position)
            if g is None: continue
            if m.node not in cache:
                cache[m.node]=np.sort([tree.time(u) for u in tree.nodes(root=m.node)
                                       if tree.is_internal(u)])
            sub=cache[m.node]
            for T in TIMES:
                if T>=tc: continue                       # non-straddling only
                nT=NMOD-np.searchsorted(internal,T,side='right')
                k=nl-np.searchsorted(sub,T,side='right')
                if k<1 or k>=nT: continue
                nd=np.linspace(tc,tp,NQ)
                p=float((P(k,nT,nd-T)*W).sum()/WS)
                rows.append((k,nT,nl,p,int(g[r_[T]]!=0)))
A=np.array(rows,float)
print(f"{len(A):,} (site,T) observations, edges above T only\n")

# within (k,n_T) strata with enough data, split by leaves-per-lineage nleaf/k
print("residual vs clade compression, pooled over (k,n_T) strata")
print(f"  {'nleaf/k':<12}{'n':>9}{'observed':>10}{'predicted':>11}{'obs/pred':>10}")
k_,nT_,nl_,p_,g_=A.T
ratio=nl_/k_
edges=[1.0,1.15,1.4,1.8,2.5,100.]
# stratify: compare within (k,nT) cells so the comparison is like-for-like
tot=defaultdict(lambda:[0,0.,0])
for kk in range(1,20):
    for nn in range(5,27):
        s=(k_==kk)&(nT_==nn)
        if s.sum()<400: continue
        for b in range(len(edges)-1):
            q=s&(ratio>=edges[b])&(ratio<edges[b+1])
            if q.sum()<30: continue
            d=tot[b]; d[0]+=g_[q].sum(); d[1]+=p_[q].sum(); d[2]+=q.sum()
for b in sorted(tot):
    c,p,n=tot[b]
    print(f"  [{edges[b]:.2f},{edges[b+1]:.2f})  {n:>9,.0f}{c:>10,.0f}{p:>11,.1f}{c/p:>10.4f}")
