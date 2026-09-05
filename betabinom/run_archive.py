"""Run the phi likelihood on the supplied true-ARG simulations."""
import sys, os, glob, json, time, collections
import numpy as np, tskit
from phi import Phi
ARCH=sys.argv[1] if len(sys.argv)>1 else "arch"
EPS=1e-3; NQ=12; NBLOCK=100; NBOOT=300
grid=np.linspace(100.,20000.,80)
W=np.full(NQ,1.); W[0]=W[-1]=.5; WS=W.sum()
_cache={}
def get_phi(Ne):
    if Ne not in _cache:
        _cache[Ne]=Phi(np.geomspace(0.5,4e6,80),26,Ne)
    return _cache[Ne]

def run(path):
    ts=tskit.load(path)
    Ne=json.loads(ts.provenance(0).record)["parameters"]["population_size"]
    P=get_phi(int(Ne))
    times=collections.Counter(ts.node(s).time for s in ts.samples())
    T_TRUE=max(times)
    anc=[s for s in ts.samples() if ts.node(s).time==T_TRUE]
    mod=[s for s in ts.samples() if ts.node(s).time==0.0]
    NMOD=len(mod)
    row=list(ts.samples()).index(anc[0])                 # ONE ancient haplotype
    gt={v.site.position:int(v.genotypes[row]!=0) for v in ts.variants()}
    tm=ts.simplify(samples=mod)
    LL=[]; POS=[]
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
            iA=np.flatnonzero((grid<tc)&(nTv>=2)&(kv>=1)&(kv<nTv))
            if len(iA):
                nd=np.linspace(tc,tp,NQ)
                key=(kv[iA]*100+nTv[iA])
                for u in np.unique(key):                  # vectorised per (k,n_T)
                    sel=iA[key==u]; k,nT=int(u//100),int(u%100)
                    ages=nd[None,:]-grid[sel][:,None]
                    p[sel]=(P(k,nT,ages)*W).sum(1)/WS
            for iT in np.flatnonzero((grid>=tc)&(grid<tp)&(nTv>=2)):
                nd=np.linspace(grid[iT],tp,NQ)
                p[iT]=(P(1,nTv[iT],nd-grid[iT])*W).sum()/WS*(tp-grid[iT])/(tp-tc)
            r=np.clip(EPS+(1-2*EPS)*p,1e-15,1-1e-15)
            LL.append((np.log(r) if gt.get(site.position,0) else np.log1p(-r)).astype(np.float32))
            POS.append(site.position)
    LL=np.array(LL); POS=np.array(POS)
    tot=LL.sum(0); mapv=grid[tot.argmax()]
    post=np.exp(tot-tot.max()); post/=np.trapezoid(post,grid)
    cdf=np.concatenate([[0],np.cumsum((post[:-1]+post[1:])/2*np.diff(grid))]); cdf/=cdf[-1]
    clo,chi=(float(np.interp(a,cdf,grid)) for a in (.025,.975))
    blk=np.minimum((POS/(ts.sequence_length/NBLOCK)).astype(int),NBLOCK-1)
    sums=np.zeros((NBLOCK,len(grid)))
    for b in range(NBLOCK):
        s=blk==b
        if s.any(): sums[b]=LL[s].sum(0)
    rng=np.random.default_rng(0)
    bm=np.array([grid[sums[rng.integers(0,NBLOCK,NBLOCK)].sum(0).argmax()] for _ in range(NBOOT)])
    blo,bhi=np.percentile(bm,[2.5,97.5])
    return T_TRUE,LL.shape[0],mapv,clo,chi,blo,bhi

print(f"{'sim':<15}{'T_true':>8}{'sites':>8}{'MAP':>8}{'composite 95%':>18}{'c':>3}"
      f"{'bootstrap 95%':>20}{'b':>3}{'sec':>6}")
cc=bb=0; n=0
for d in sorted(glob.glob(f"{ARCH}/simulation_*")):
    f=os.path.join(d,os.path.basename(d)+".trees")
    t0=time.time(); T,ns,mv,clo,chi,blo,bhi=run(f)
    c='Y' if clo<=T<=chi else '.'; b='Y' if blo<=T<=bhi else '.'
    cc+= c=='Y'; bb+= b=='Y'; n+=1
    print(f"{os.path.basename(d):<15}{T:>8.0f}{ns:>8,}{mv:>8.0f}"
          f"{clo:>9.0f}-{chi:<8.0f}{c:>3}{blo:>10.0f}-{bhi:<9.0f}{b:>3}{time.time()-t0:>6.0f}")
print(f"\ncoverage: composite {cc}/{n}   block-bootstrap {bb}/{n}")
