"""Calibrate phi2 on straddling sites (cached on a coarse grid for speed)."""
import numpy as np, msprime, tskit
from phi import Phi
from phi2 import phi2
NE,NMOD,L,RR,MU = 100_000,26,2e6,1e-8,1e-8
T=2500.; NQ=8
tau=lambda t: t/(2*NE)
P=Phi(np.geomspace(0.5,2e6,70),26,NE)
W=np.full(NQ,1.); W[0]=W[-1]=.5; WS=W.sum()
cache={}
def p2(nT,nleaf,age):
    key=(nT,min(nleaf,8),int(np.log(max(age,1.))*3))
    if key not in cache:
        a=np.exp(key[2]/3.)
        cache[key]=phi2(nT,min(nleaf,8),tau(a),tau(T),NE)
    return cache[key]
obs=pr1=pr2=0.; bins={}
for rep in range(4):
    ts=msprime.sim_ancestry(samples=[msprime.SampleSet(NMOD,time=0,ploidy=1),
        msprime.SampleSet(1,time=T,ploidy=1)],population_size=NE,sequence_length=L,
        recombination_rate=RR,ploidy=2,random_seed=1200+rep)
    ts=msprime.sim_mutations(ts,rate=MU,random_seed=6100+rep,discrete_genome=False)
    anc=ts.num_samples-1; ar=list(ts.samples()).index(anc)
    gt={v.site.position:int(v.genotypes[ar]!=0) for v in ts.variants()}
    tm=ts.simplify(samples=[s for s in ts.samples() if s!=anc])
    for tree in tm.trees():
        internal=np.sort([tree.time(u) for u in tree.nodes() if tree.is_internal(u)])
        nT=NMOD-np.searchsorted(internal,T,side='right')
        if nT<2: continue
        for site in tree.sites():
            if len(site.mutations)!=1: continue
            m=site.mutations[0]; par=tree.parent(m.node)
            if par==tskit.NULL: continue
            nl=tree.num_samples(m.node)
            if not (0<nl<NMOD): continue
            tc,tp=tree.time(m.node),tree.time(par)
            if not (tc<=T<tp): continue
            nd=np.linspace(T,tp,NQ); f=(tp-T)/(tp-tc)
            a1=float((P(1,nT,nd-T)*W).sum()/WS*f)
            v=np.array([p2(nT,nl,x) for x in (nd-T)])
            a2=float(np.nansum(v*W)/WS*f)
            g=gt.get(site.position,0)
            obs+=g; pr1+=a1; pr2+=a2
            b=min(nl,4); d=bins.setdefault(b,[0,0.,0.,0]); d[0]+=g;d[1]+=a1;d[2]+=a2;d[3]+=1
print(f"straddling sites at T={T:g}, {sum(d[3] for d in bins.values()):,} sites"
      f"   ({len(cache)} phi2 evaluations cached)\n")
print(f"{'nleaf':>7}{'n':>9}{'obs':>8}{'phi pred':>11}{'obs/phi':>9}{'phi2 pred':>12}{'obs/phi2':>10}")
for b in sorted(bins):
    o,x1,x2,n=bins[b]
    print(f"{b:>7}{n:>9,}{o:>8,.0f}{x1:>11.1f}{o/x1:>9.3f}{x2:>12.1f}{o/x2:>10.3f}")
print(f"{'ALL':>7}{sum(d[3] for d in bins.values()):>9,}{obs:>8,.0f}"
      f"{pr1:>11.1f}{obs/pr1:>9.3f}{pr2:>12.1f}{obs/pr2:>10.3f}")
