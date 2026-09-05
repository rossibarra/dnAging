"""Does pseudohaploid sampling change anything?

Archive1's ancient individuals are diploid, so we can compare:
  (a) one haplotype used throughout      -- what was validated
  (b) pseudohaploid: at each site pick one of the two at random

Marginally these are the same Bernoulli(X(T)) at every site, so the MAP should
agree.  Pseudohaploid switches chromosome between sites, so it carries less LD and
its block-bootstrap SD should if anything be smaller.
"""
import sys, os, glob, json, collections
import numpy as np, tskit
from phid import PhiD
ARCH=sys.argv[1]; EPS=float(sys.argv[2]) if len(sys.argv)>2 else 1e-6
NQ=12; NBLOCK=100; NBOOT=300
grid=np.linspace(100.,20000.,80)
W=np.full(NQ,1.); W[0]=W[-1]=.5
_c={}
def get(Ne):
    if Ne not in _c: _c[Ne]=PhiD(np.geomspace(0.5,4e6,80),26,Ne)
    return _c[Ne]

def run(path, mode, seed):
    ts=tskit.load(path)
    Ne=int(json.loads(ts.provenance(0).record)["parameters"]["population_size"])
    P=get(Ne)
    times=collections.Counter(ts.node(s).time for s in ts.samples()); T_TRUE=max(times)
    anc=[s for s in ts.samples() if ts.node(s).time==T_TRUE]
    mod=[s for s in ts.samples() if ts.node(s).time==0.0]; NMOD=len(mod)
    smp=list(ts.samples()); rows=[smp.index(a) for a in anc]
    rng=np.random.default_rng(seed)
    gt={}
    for v in ts.variants():
        g=v.genotypes
        if mode=="hap":  gt[v.site.position]=int(g[rows[0]]!=0)
        else:            gt[v.site.position]=int(g[rows[rng.integers(0,2)]]!=0)
    tm=ts.simplify(samples=mod)
    LL=[];POS=[]
    for tree in tm.trees():
        internal=np.sort([tree.time(u) for u in tree.nodes() if tree.is_internal(u)])
        nTv=NMOD-np.searchsorted(internal,grid,side='right'); cache={}
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
            kv=nl-np.searchsorted(cache[m.node],grid,side='right')
            p=np.zeros(len(grid)); nd=np.linspace(tc,tp,NQ)
            iA=np.flatnonzero((grid<tc)&(nTv>=2)&(kv>=1)&(kv<nTv))
            if len(iA):
                key=kv[iA]*100+nTv[iA]
                for u in np.unique(key):
                    sel=iA[key==u]; k,nT=int(u//100),int(u%100)
                    p[sel]=P.integrate(k,nT,nd[None,:]-grid[sel][:,None],W)
            for iT in np.flatnonzero((grid>=tc)&(grid<tp)&(nTv>=2)):
                q=np.linspace(grid[iT],tp,NQ)
                p[iT]=float(P.integrate(1,nTv[iT],q-grid[iT],W))*(tp-grid[iT])/(tp-tc)
            r=np.clip(EPS+(1-2*EPS)*p,1e-15,1-1e-15)
            LL.append((np.log(r) if gt.get(site.position,0) else np.log1p(-r)).astype(np.float32))
            POS.append(site.position)
    LL=np.array(LL);POS=np.array(POS)
    blk=np.minimum((POS/(ts.sequence_length/NBLOCK)).astype(int),NBLOCK-1)
    sums=np.zeros((NBLOCK,len(grid)))
    for b in range(NBLOCK):
        s=blk==b
        if s.any(): sums[b]=LL[s].sum(0)
    r2=np.random.default_rng(0)
    bm=np.array([grid[sums[r2.integers(0,NBLOCK,NBLOCK)].sum(0).argmax()] for _ in range(NBOOT)])
    return T_TRUE, grid[LL.sum(0).argmax()], bm.std(ddof=1)

print(f"{'sim':<15}{'T_true':>8}{'hap MAP':>9}{'hap SD':>8}{'pseudo MAP':>12}{'pseudo SD':>11}",flush=True)
Th=[];Mh=[];Sh=[];Mp=[];Sp=[]
for d in sorted(glob.glob(f"{ARCH}/simulation_*")):
    f=os.path.join(d,os.path.basename(d)+".trees")
    T,mh,sh=run(f,"hap",0); _,mp,sp=run(f,"pseudo",12345)
    Th.append(T);Mh.append(mh);Sh.append(sh);Mp.append(mp);Sp.append(sp)
    print(f"{os.path.basename(d):<15}{T:>8.0f}{mh:>9.0f}{sh:>8.0f}{mp:>12.0f}{sp:>11.0f}",flush=True)
Th,Mh,Sh,Mp,Sp=map(np.array,(Th,Mh,Sh,Mp,Sp))
for lab,M,S in (("one haplotype",Mh,Sh),("pseudohaploid",Mp,Sp)):
    e=M-Th
    print(f"\n{lab:<16} RMSE {np.sqrt((e**2).mean()):>6.0f}  bias {e.mean():>+7.0f}  "
          f"median bootstrap SD {np.median(S):>6.0f}")
print(f"\nMAP agreement: mean |diff| {np.abs(Mp-Mh).mean():.0f} generations "
      f"({np.abs(Mp-Mh).mean()/Th.mean()*100:.1f}% of mean true age)")
print(f"SD ratio pseudo/hap: {np.median(Sp)/np.median(Sh):.3f}  (<1 if LD is reduced)")
