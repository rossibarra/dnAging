"""Same as run_archive.py but marginalising the mutation age correctly:
   p = sum(w*num)/sum(w*den)   instead of   mean(w*num/den)."""
import sys, os, glob, json, collections
import numpy as np, tskit
from phid import PhiD
ARCH=sys.argv[1]
EPS=float(sys.argv[2]) if len(sys.argv)>2 else 1e-3
EPS_=EPS; NQ=12; NBLOCK=100; NBOOT=300
grid=np.linspace(100.,20000.,80)
W=np.full(NQ,1.); W[0]=W[-1]=.5
_c={}
def get(Ne):
    if Ne not in _c: _c[Ne]=PhiD(np.geomspace(0.5,4e6,80),26,Ne)
    return _c[Ne]

def run(path):
    ts=tskit.load(path)
    Ne=json.loads(ts.provenance(0).record)["parameters"]["population_size"]
    P=get(int(Ne))
    times=collections.Counter(ts.node(s).time for s in ts.samples()); T_TRUE=max(times)
    anc=[s for s in ts.samples() if ts.node(s).time==T_TRUE]
    mod=[s for s in ts.samples() if ts.node(s).time==0.0]; NMOD=len(mod)
    row=list(ts.samples()).index(anc[0])
    gt={v.site.position:int(v.genotypes[row]!=0) for v in ts.variants()}
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
            p=np.zeros(len(grid))
            nd_full=np.linspace(tc,tp,NQ)
            iA=np.flatnonzero((grid<tc)&(nTv>=2)&(kv>=1)&(kv<nTv))
            if len(iA):
                key=kv[iA]*100+nTv[iA]
                for u in np.unique(key):                  # vectorised per (k,n_T)
                    sel=iA[key==u]; k,nT=int(u//100),int(u%100)
                    p[sel]=P.integrate(k,nT,nd_full[None,:]-grid[sel][:,None],W)
            for iT in np.flatnonzero((grid>=tc)&(grid<tp)&(nTv>=2)):
                nd=np.linspace(grid[iT],tp,NQ)
                # the below-T part of the edge contributes 0 to the numerator and
                # nothing to the k=1-at-T conditioning, so it only rescales by frac
                p[iT]=float(P.integrate(1,nTv[iT],nd-grid[iT],W))*(tp-grid[iT])/(tp-tc)
            r=np.clip(EPS+(1-2*EPS)*p,1e-15,1-1e-15)
            LL.append((np.log(r) if gt.get(site.position,0) else np.log1p(-r)).astype(np.float32))
            POS.append(site.position)
    LL=np.array(LL);POS=np.array(POS);tot=LL.sum(0)
    blk=np.minimum((POS/(ts.sequence_length/NBLOCK)).astype(int),NBLOCK-1)
    sums=np.zeros((NBLOCK,len(grid)))
    for b in range(NBLOCK):
        s=blk==b
        if s.any(): sums[b]=LL[s].sum(0)
    rng=np.random.default_rng(0)
    bm=np.array([grid[sums[rng.integers(0,NBLOCK,NBLOCK)].sum(0).argmax()] for _ in range(NBOOT)])
    return T_TRUE, grid[tot.argmax()], *np.percentile(bm,[2.5,97.5])

print(f"eps={EPS:g}\n{'sim':<15}{'T_true':>8}{'MAP':>8}{'bootstrap 95%':>20}{'cov':>5}",flush=True)
T_=[];M_=[];cov=0
for d in sorted(glob.glob(f"{ARCH}/simulation_*")):
    f=os.path.join(d,os.path.basename(d)+".trees")
    T,mv,blo,bhi=run(f)
    c='Y' if blo<=T<=bhi else '.'; cov+= c=='Y'; T_.append(T);M_.append(mv)
    print(f"{os.path.basename(d):<15}{T:>8.0f}{mv:>8.0f}{blo:>10.0f}-{bhi:<9.0f}{c:>5}",flush=True)
T_=np.array(T_);M_=np.array(M_);sl,ic=np.polyfit(T_,M_,1)
print(f"\nCORRECTED  MAP = {sl:.3f}*T {ic:+.0f}   r={np.corrcoef(T_,M_)[0,1]:.3f}"
      f"   mean bias {(M_-T_).mean():+.0f}   coverage {cov}/10")
print(f"reference (eps=1e-3, corrected): MAP = 1.023*T +1117  r=0.962  bias +1239  cov 8/10")
