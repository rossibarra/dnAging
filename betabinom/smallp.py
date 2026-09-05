"""Calibration on a log scale, where the likelihood actually has leverage."""
import numpy as np, msprime, tskit
from phi import Phi
NE,NMOD,L,RR,MU,EPS = 100_000,26,2e6,1e-8,1e-8,1e-3
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
                nT=NMOD-np.searchsorted(internal,T,side='right')
                if nT<2 or T>=tp: continue
                if T>=tc:
                    nd=np.linspace(T,tp,NQ); p=(P(1,nT,nd-T)*W).sum()/WS*(tp-T)/(tp-tc)
                else:
                    k=nl-np.searchsorted(sub,T,side='right')
                    if k<1 or k>=nT: continue
                    nd=np.linspace(tc,tp,NQ); p=(P(k,nT,nd-T)*W).sum()/WS
                rows.append((p,int(g[r_[T]]!=0)))
A=np.array(rows); p_,g_=A.T
print(f"{len(A):,} observations\n")
print(f"{'p range':<20}{'n':>9}{'observed':>10}{'predicted':>11}{'obs/pred':>10}{'lev share':>11}")
edges=[0,1e-4,3e-4,1e-3,3e-3,1e-2,3e-2,0.1,0.2,0.4,1.01]
lev_tot=0.
levs=[]
for i in range(len(edges)-1):
    s=(p_>=edges[i])&(p_<edges[i+1])
    if s.sum()==0: levs.append(0); continue
    r=np.clip(EPS+(1-2*EPS)*p_[s],1e-12,1-1e-12)
    lev=(p_[s]**2/(r*(1-r))).sum()     # Fisher information weight for a scale change
    levs.append(lev); lev_tot+=lev
for i in range(len(edges)-1):
    s=(p_>=edges[i])&(p_<edges[i+1])
    if s.sum()==0: continue
    o,pr=g_[s].sum(),p_[s].sum()
    lbl=f"[{edges[i]:g},{edges[i+1]:g})"
    print(f"{lbl:<20}{s.sum():>9,}{o:>10,.0f}{pr:>11,.1f}"
          f"{(o/pr if pr>0 else 0):>10.4f}{levs[i]/lev_tot*100:>10.1f}%")
