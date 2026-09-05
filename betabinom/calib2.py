"""Calibration gate for the unified phi(k,n_T,age) likelihood."""
import numpy as np, msprime, tskit
from phi import Phi
NE,NMOD,L,RR,MU = 100_000,26,2e6,1e-8,1e-8
TIMES=[500.,1000.,2000.,2500.,4000.]
NQ=32
P=Phi(np.geomspace(0.5,2e6,70),26,NE)

def p_site(k, nT, tc, tp, T):
    """uniform mutation placement on (tc,tp); below T contributes 0."""
    lo=max(T,tc)
    if tp<=lo: return 0.0
    nodes=np.linspace(lo,tp,NQ); w=np.full(NQ,1.0); w[0]=w[-1]=0.5
    vals=P(k,nT,nodes-T)
    return float((vals*w).sum()/w.sum()*(tp-lo)/(tp-tc))

old=[0,0.]; strad=[0,0.]; yng=[0,0.]; sb={}; ob={}
for rep in range(6):
    ss=[msprime.SampleSet(NMOD,time=0,ploidy=1)]+[msprime.SampleSet(1,time=t,ploidy=1) for t in TIMES]
    ts=msprime.sim_ancestry(samples=ss,population_size=NE,sequence_length=L,
        recombination_rate=RR,ploidy=2,random_seed=900+rep)
    ts=msprime.sim_mutations(ts,rate=MU,random_seed=4000+rep,discrete_genome=False)
    smp=list(ts.samples()); rows={t:smp.index(smp[NMOD+i]) for i,t in enumerate(TIMES)}
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
                if nT<2: continue
                c=int(g[rows[T]]!=0)
                if T>=tp:
                    yng[0]+=c; continue
                if T>=tc:
                    k=1; pr=p_site(1,nT,tc,tp,T); f=(tp-T)/(tp-tc)
                    strad[0]+=c; strad[1]+=pr
                    b=min(int(f*5),4); d=sb.setdefault(b,[0,0.,0]); d[0]+=c;d[1]+=pr;d[2]+=1
                else:
                    k=nl-np.searchsorted(sub,T,side='right')
                    if k<1 or k>=nT: continue
                    pr=p_site(k,nT,tc,tp,T)
                    old[0]+=c; old[1]+=pr
                    b=min(int(pr*5),4); d=ob.setdefault(b,[0,0.,0]); d[0]+=c;d[1]+=pr;d[2]+=1

print(f"{'case':<26}{'observed':>10}{'predicted':>11}{'obs/pred':>10}")
for nm,v in (("edge older than T",old),("edge straddles T",strad)):
    print(f"{nm:<26}{v[0]:>10,}{v[1]:>11,.0f}{v[0]/v[1]:>10.4f}")
print(f"{'edge younger than T (p=0)':<26}{yng[0]:>10,}{'0':>11}{'--':>10}")
to=old[0]+strad[0]+yng[0]; tp_=old[1]+strad[1]
print(f"{'TOTAL':<26}{to:>10,}{tp_:>11,.0f}{to/tp_:>10.4f}")
print(f"\nstraddling, by fraction of edge above T   (was 0.00/0.09/0.17/0.31/1.07)")
print(f"  {'frac':<12}{'n':>9}{'observed':>10}{'predicted':>11}{'obs/pred':>10}")
for b in sorted(sb):
    c,p,n=sb[b]; print(f"  [{b/5:.1f},{(b+1)/5:.1f})   {n:>9,}{c:>10,}{p:>11,.1f}{(c/p if p else 0):>10.4f}")
print(f"\nedge older than T, by predicted p")
print(f"  {'p':<12}{'n':>9}{'observed':>10}{'predicted':>11}{'obs/pred':>10}")
for b in sorted(ob):
    c,p,n=ob[b]; print(f"  [{b/5:.1f},{(b+1)/5:.1f})   {n:>9,}{c:>10,}{p:>11,.1f}{(c/p if p else 0):>10.4f}")
