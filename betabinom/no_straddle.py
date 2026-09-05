"""Is the +892 intercept coming from the young-mutation (straddling) term?

Retain only sites whose edge never straddles ANY T in the grid -- i.e. the edge
lies entirely above grid.max() or entirely below grid.min().  The retained set is
then the same at every T, so the likelihood stays a comparison over one fixed set
of observations, and all remaining T-dependence comes from k(T) and n_T(T).
"""
import sys, os, glob, json, time, collections
import numpy as np, tskit
from phi import Phi
ARCH=sys.argv[1]
EPS=1e-3; NQ=12; NBLOCK=100; NBOOT=300
grid=np.linspace(100.,20000.,80)
W=np.full(NQ,1.); W[0]=W[-1]=.5; WS=W.sum()
_c={}
def get_phi(Ne):
    if Ne not in _c: _c[Ne]=Phi(np.geomspace(0.5,4e6,80),26,Ne)
    return _c[Ne]

def run(path, drop_straddle):
    ts=tskit.load(path)
    Ne=json.loads(ts.provenance(0).record)["parameters"]["population_size"]
    P=get_phi(int(Ne))
    times=collections.Counter(ts.node(s).time for s in ts.samples())
    T_TRUE=max(times)
    anc=[s for s in ts.samples() if ts.node(s).time==T_TRUE]
    mod=[s for s in ts.samples() if ts.node(s).time==0.0]; NMOD=len(mod)
    row=list(ts.samples()).index(anc[0])
    gt={v.site.position:int(v.genotypes[row]!=0) for v in ts.variants()}
    tm=ts.simplify(samples=mod)
    LL=[]; POS=[]; nstr=0
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
            straddles = (tc < grid[-1]) and (tp > grid[0])
            if drop_straddle and straddles:
                nstr+=1; continue
            if m.node not in cache:
                cache[m.node]=np.sort([tree.time(u) for u in tree.nodes(root=m.node)
                                       if tree.is_internal(u)])
            kv=nl-np.searchsorted(cache[m.node],grid,side='right')
            p=np.zeros(len(grid))
            iA=np.flatnonzero((grid<tc)&(nTv>=2)&(kv>=1)&(kv<nTv))
            if len(iA):
                nd=np.linspace(tc,tp,NQ); key=kv[iA]*100+nTv[iA]
                for u in np.unique(key):
                    sel=iA[key==u]; k,nT=int(u//100),int(u%100)
                    p[sel]=(P(k,nT,nd[None,:]-grid[sel][:,None])*W).sum(1)/WS
            if not drop_straddle:
                for iT in np.flatnonzero((grid>=tc)&(grid<tp)&(nTv>=2)):
                    nd=np.linspace(grid[iT],tp,NQ)
                    p[iT]=(P(1,nTv[iT],nd-grid[iT])*W).sum()/WS*(tp-grid[iT])/(tp-tc)
            r=np.clip(EPS+(1-2*EPS)*p,1e-15,1-1e-15)
            LL.append((np.log(r) if gt.get(site.position,0) else np.log1p(-r)).astype(np.float32))
            POS.append(site.position)
    LL=np.array(LL); POS=np.array(POS); tot=LL.sum(0)
    blk=np.minimum((POS/(ts.sequence_length/NBLOCK)).astype(int),NBLOCK-1)
    sums=np.zeros((NBLOCK,len(grid)))
    for b in range(NBLOCK):
        s=blk==b
        if s.any(): sums[b]=LL[s].sum(0)
    rng=np.random.default_rng(0)
    bm=np.array([grid[sums[rng.integers(0,NBLOCK,NBLOCK)].sum(0).argmax()] for _ in range(NBOOT)])
    blo,bhi=np.percentile(bm,[2.5,97.5])
    return T_TRUE, LL.shape[0], nstr, grid[tot.argmax()], blo, bhi

for drop in (False, True):
    tag="straddling DROPPED" if drop else "full model"
    print(f"\n=== {tag} ===")
    print(f"{'sim':<15}{'T_true':>8}{'sites':>8}{'dropped':>9}{'MAP':>8}{'bootstrap 95%':>20}{'cov':>5}")
    T_=[];M_=[];cov=0
    for d in sorted(glob.glob(f"{ARCH}/simulation_*")):
        f=os.path.join(d,os.path.basename(d)+".trees")
        T,ns,nstr,mv,blo,bhi=run(f,drop)
        c='Y' if blo<=T<=bhi else '.'; cov+= c=='Y'
        T_.append(T); M_.append(mv)
        print(f"{os.path.basename(d):<15}{T:>8.0f}{ns:>8,}{nstr:>9,}{mv:>8.0f}"
              f"{blo:>10.0f}-{bhi:<9.0f}{c:>5}")
    T_=np.array(T_,float); M_=np.array(M_,float); sl,ic=np.polyfit(T_,M_,1)
    print(f"  regression MAP = {sl:.3f}*T + {ic:+.0f}   r={np.corrcoef(T_,M_)[0,1]:.3f}"
          f"   mean bias {(M_-T_).mean():+.0f}   coverage {cov}/10")
