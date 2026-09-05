"""For straddling edges phi(1,n_T,a) ignores nleaf -- the size of the clade the
mutant lineage founds below T.  Does carriage depend on it?"""
import numpy as np, msprime, tskit
from phi import Phi
NE,NMOD,L,RR,MU = 100_000,26,2e6,1e-8,1e-8
TIMES=[500.,1000.,2000.,2500.,4000.]; NQ=24
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
        for site in tree.sites():
            if len(site.mutations)!=1: continue
            m=site.mutations[0]; par=tree.parent(m.node)
            if par==tskit.NULL: continue
            nl=tree.num_samples(m.node)
            if not (0<nl<NMOD): continue
            tc,tp=tree.time(m.node),tree.time(par)
            g=gts.get(site.position)
            if g is None: continue
            for T in TIMES:
                if not (tc<=T<tp): continue          # STRADDLING only
                nT=NMOD-np.searchsorted(internal,T,side='right')
                if nT<2: continue
                nd=np.linspace(T,tp,NQ)
                p=float((P(1,nT,nd-T)*W).sum()/WS*(tp-T)/(tp-tc))
                rows.append((nl,p,int(g[r_[T]]!=0)))
A=np.array(rows); nl_,p_,g_=A.T
print(f"{len(A):,} straddling observations\n")
print(f"{'nleaf (clade size)':<20}{'n':>9}{'observed':>10}{'predicted':>11}{'obs/pred':>10}")
for lo,hi in [(1,2),(2,3),(3,5),(5,9),(9,27)]:
    s=(nl_>=lo)&(nl_<hi)
    if s.sum()<50: continue
    o,pr=g_[s].sum(),p_[s].sum()
    lbl=f"{lo}" if hi==lo+1 else f"{lo}-{hi-1}"
    print(f"{lbl:<20}{s.sum():>9,}{o:>10,.0f}{pr:>11,.1f}{(o/pr if pr>0 else 0):>10.4f}")
o,pr=g_.sum(),p_.sum()
print(f"{'ALL':<20}{len(A):>9,}{o:>10,.0f}{pr:>11,.1f}{o/pr:>10.4f}")
