"""Blind run on Archive2: panel-only .trees + haploid ancient VCF.

Differences from run_archive.py:
  * the ancient sample is NOT in the tree sequence -- it comes from the VCF
  * VCF POS = int(tree position) + 1  (tskit 1-based convention)
  * the ancient VCF lists only carried sites; absence means not carried
  * Ne is not recorded, so it is estimated per simulation by Watterson
"""
import sys, os, glob, collections
import numpy as np, tskit
from phid import PhiD
ARCH=sys.argv[1]; EPS=float(sys.argv[2]) if len(sys.argv)>2 else 1e-6
MU=1e-8; NQ=12; NBLOCK=100; NBOOT=300
NE_FIXED=float(os.environ.get("NE_FIXED","0")) or None   # override the Watterson estimate
grid=np.linspace(100.,20000.,80)
W=np.full(NQ,1.); W[0]=W[-1]=.5
_c={}
def get(Ne):
    key=int(round(Ne/1000.)*1000)
    if key not in _c: _c[key]=PhiD(np.geomspace(0.5,4e6,80),26,key)
    return _c[key],key

def run(d):
    b=os.path.basename(d); ts=tskit.load(f"{d}/{b}.trees")
    NMOD=ts.num_samples
    a_n=sum(1.0/i for i in range(1,NMOD))
    Ne_hat=NE_FIXED or (ts.num_sites/a_n/ts.sequence_length)/(4*MU)
    P,Ne=get(Ne_hat)
    carried=set()
    for ln in open(f"{d}/{b}_ancient.vcf"):
        if ln.startswith('#'): continue
        f=ln.split('\t')
        if f[9].strip().startswith('1'): carried.add(int(f[1])-1)
    LL=[];POS=[]
    for tree in ts.trees():
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
            p=np.zeros(len(grid)); nd_full=np.linspace(tc,tp,NQ)
            iA=np.flatnonzero((grid<tc)&(nTv>=2)&(kv>=1)&(kv<nTv))
            if len(iA):
                key=kv[iA]*100+nTv[iA]
                for u in np.unique(key):
                    sel=iA[key==u]; k,nT=int(u//100),int(u%100)
                    p[sel]=P.integrate(k,nT,nd_full[None,:]-grid[sel][:,None],W)
            for iT in np.flatnonzero((grid>=tc)&(grid<tp)&(nTv>=2)):
                nd=np.linspace(grid[iT],tp,NQ)
                p[iT]=float(P.integrate(1,nTv[iT],nd-grid[iT],W))*(tp-grid[iT])/(tp-tc)
            r=np.clip(EPS+(1-2*EPS)*p,1e-15,1-1e-15)
            g=int(site.position) in carried
            LL.append((np.log(r) if g else np.log1p(-r)).astype(np.float32))
            POS.append(site.position)
    LL=np.array(LL);POS=np.array(POS);tot=LL.sum(0)
    blk=np.minimum((POS/(ts.sequence_length/NBLOCK)).astype(int),NBLOCK-1)
    sums=np.zeros((NBLOCK,len(grid)))
    for bi in range(NBLOCK):
        s=blk==bi
        if s.any(): sums[bi]=LL[s].sum(0)
    rng=np.random.default_rng(0)
    bm=np.array([grid[sums[rng.integers(0,NBLOCK,NBLOCK)].sum(0).argmax()] for _ in range(NBOOT)])
    ncar=sum(1 for x in POS if int(x) in carried)
    return Ne, LL.shape[0], ncar, grid[tot.argmax()], *np.percentile(bm,[2.5,97.5])

print(f"eps={EPS:g}   (BLIND - no true ages available in these files)",flush=True)
print(f"{'sim':<15}{'Ne_hat':>8}{'sites':>8}{'carried':>9}{'MAP':>8}{'bootstrap 95%':>20}",flush=True)
for d in sorted(glob.glob(f"{ARCH}/simulation_*")):
    Ne,n,nc,mv,lo,hi=run(d)
    print(f"{os.path.basename(d):<15}{Ne:>8,}{n:>8,}{nc:>9,}{mv:>8.0f}{lo:>10.0f}-{hi:<9.0f}",flush=True)
