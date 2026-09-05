import numpy as np, msprime, tskit
from phi import Phi
NE,NMOD,L,RR,MU,EPS,SEED,T_TRUE=100_000,26,1e6,1e-8,1e-8,1e-3,11,2500.
grid=np.array([1000.,2500.,4000.,7658.]); NQ=24
P=Phi(np.geomspace(0.5,2e6,70),26,NE)
W=np.full(NQ,1.); W[0]=W[-1]=.5; WS=W.sum()
ts=msprime.sim_ancestry(samples=[msprime.SampleSet(NMOD,time=0,ploidy=1),
    msprime.SampleSet(1,time=T_TRUE,ploidy=1)],population_size=NE,sequence_length=L,
    recombination_rate=RR,ploidy=2,random_seed=SEED)
ts=msprime.sim_mutations(ts,rate=MU,random_seed=SEED+7,discrete_genome=False)
anc=ts.num_samples-1; ar=list(ts.samples()).index(anc)
gt={v.site.position:int(v.genotypes[ar]!=0) for v in ts.variants()}
tm=ts.simplify(samples=[s for s in ts.samples() if s!=anc])
S={i:{'old':[0,0.,0],'str':[0,0.,0],'yng':[0,0.,0]} for i in range(len(grid))}
LL={i:0. for i in range(len(grid))}
for tree in tm.trees():
    internal=np.sort([tree.time(u) for u in tree.nodes() if tree.is_internal(u)])
    nTv=NMOD-np.searchsorted(internal,grid,side='right'); cache={}
    for site in tree.sites():
        if len(site.mutations)!=1: continue
        m=site.mutations[0]; par=tree.parent(m.node)
        if par==tskit.NULL: continue
        nl=tree.num_samples(m.node)
        if not (0<nl<NMOD): continue
        tc,tp=tree.time(m.node),tree.time(par); g=gt.get(site.position,0)
        if m.node not in cache:
            cache[m.node]=np.sort([tree.time(u) for u in tree.nodes(root=m.node)
                                   if tree.is_internal(u)])
        for i,T in enumerate(grid):
            nT=nTv[i]
            if nT<2: continue
            if T>=tp: cat='yng'; p=0.
            elif T>=tc:
                cat='str'; nd=np.linspace(T,tp,NQ)
                p=(P(1,nT,nd-T)*W).sum()/WS*(tp-T)/(tp-tc)
            else:
                cat='old'; k=nl-np.searchsorted(cache[m.node],T,side='right')
                if k<1 or k>=nT: continue
                nd=np.linspace(tc,tp,NQ); p=(P(k,nT,nd-T)*W).sum()/WS
            d=S[i][cat]; d[0]+=g; d[1]+=p; d[2]+=1
            r=min(max(EPS+(1-2*EPS)*p,1e-15),1-1e-15)
            LL[i]+= np.log(r) if g else np.log1p(-r)
print(f"{'T':>7} {'case':>5} {'sites':>8} {'obs':>7} {'pred':>9} {'obs/pred':>9}")
for i,T in enumerate(grid):
    for cat in ('old','str','yng'):
        c,p,n=S[i][cat]
        if n==0 and cat=='yng': continue
        print(f"{T:>7.0f} {cat:>5} {n:>8,} {c:>7,} {p:>9.1f} {(f'{c/p:.4f}' if p>0 else '--'):>9}")
    to=sum(S[i][x][0] for x in S[i]); tp_=sum(S[i][x][1] for x in S[i])
    print(f"{'':>7} {'ALL':>5} {'':>8} {to:>7,} {tp_:>9.1f} {to/tp_:>9.4f}   loglik {LL[i]:.1f}\n")
