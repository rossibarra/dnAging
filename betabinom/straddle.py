"""Calibrate the three cases separately."""
import msprime, tskit, numpy as np
NE,NMOD,L,RR,MU = 100_000,26,2e6,1e-8,1e-8
TIMES=[500.,1000.,2000.,2500.,4000.]
older=[0,0.]; strad=[0,0.]; post=[0,0.]     # [observed carried, predicted sum]
sbins={}
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
                c=int(g[rows[T]]!=0)
                if T>=tp:
                    post[0]+=c; post[1]+=0.0
                elif T>=tc:
                    f=(tp-T)/(tp-tc); pr=f/(nT+1.0)
                    strad[0]+=c; strad[1]+=pr
                    b=min(int(f*5),4); d=sbins.setdefault(b,[0,0.,0])
                    d[0]+=c; d[1]+=pr; d[2]+=1
                else:
                    k=nl-np.searchsorted(sub,T,side='right')
                    older[0]+=c; older[1]+=k/(nT+1.0)
print(f"{'case':<34}{'observed':>10}{'predicted':>11}{'obs/pred':>10}")
print(f"{'edge older than T  (k/(n_T+1))':<34}{older[0]:>10,}{older[1]:>11,.0f}{older[0]/older[1]:>10.4f}")
print(f"{'edge straddles T   (frac/(n_T+1))':<34}{strad[0]:>10,}{strad[1]:>11,.0f}{strad[0]/strad[1]:>10.4f}")
print(f"{'edge younger than T  (p=0)':<34}{post[0]:>10,}{'0':>11}{'--':>10}")
tot_o=older[0]+strad[0]+post[0]; tot_p=older[1]+strad[1]
print(f"{'TOTAL':<34}{tot_o:>10,}{tot_p:>11,.0f}{tot_o/tot_p:>10.4f}")
print(f"\nstraddling, split by fraction of the edge above T:")
print(f"  {'frac bin':<12}{'n':>9}{'observed':>10}{'predicted':>11}{'obs/pred':>10}")
for b in sorted(sbins):
    c,p,n=sbins[b]
    print(f"  [{b/5:.1f},{(b+1)/5:.1f})   {n:>9,}{c:>10,}{p:>11,.0f}{(c/p if p else 0):>10.4f}")
