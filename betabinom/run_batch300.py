"""betabinom age estimates across the 300-simulation batch.

Each simulation has its own Ne, mu, recombination ratio, panel size and true age,
so the phi table is built per simulation at that Ne and panel size.  The point is
to see whether the additive offset found on the 20 fixed-parameter simulations
holds when Ne, mu and panel size vary.
"""
import os, sys, csv, time, argparse
import multiprocessing as mp
import numpy as np, tskit
from phid import PhiD

EPS=1e-6; NQ=12; NBLOCK=100; NBOOT=200
GRID=np.linspace(50., 15000., 100)
W=np.full(NQ,1.); W[0]=W[-1]=.5

def one(args):
    row, simdir = args
    sim=row["sim"]; Ne=int(row["Ne"]); nmod=int(row["n_modern"]); age=float(row["age"])
    t0=time.time()
    try:
        ts=tskit.load(os.path.join(simdir, f"{sim}.trees"))
        if ts.num_samples != nmod:
            raise RuntimeError(f"panel size {ts.num_samples} != {nmod}")
        # tskit's VCF POS convention differs between versions (POS=pos in 1.0.3,
        # POS=pos+1 in 1.0.0), so detect the offset against the tree's own sites
        # rather than assuming one.
        raw=[]
        for ln in open(os.path.join(simdir, f"{sim}_ancient.vcf")):
            if ln.startswith('#'): continue
            f=ln.split('\t'); raw.append((int(f[1]), f[9].strip().startswith('1')))
        tree_pos=set(int(s_.position) for s_ in ts.sites())
        allp={p for p,_ in raw}
        off=max((0,-1,1), key=lambda o: len({p+o for p in allp} & tree_pos))
        hit=len({p+off for p in allp} & tree_pos)
        if hit < 0.5*len(tree_pos):
            raise RuntimeError(f"VCF positions do not align: best offset {off} "
                               f"matches only {hit}/{len(tree_pos)} tree sites")
        carried={p+off for p,c in raw if c}
        P=PhiD(np.geomspace(0.5, 4e6, 60), nmod, Ne, dps=45)
        LL=[]; POS=[]
        for tree in ts.trees():
            internal=np.sort([tree.time(u) for u in tree.nodes() if tree.is_internal(u)])
            nTv=nmod-np.searchsorted(internal, GRID, side='right'); cache={}
            for site in tree.sites():
                if len(site.mutations)!=1: continue
                m=site.mutations[0]; par=tree.parent(m.node)
                if par==tskit.NULL: continue
                nl=tree.num_samples(m.node)
                if not (0<nl<nmod): continue
                tc,tp=tree.time(m.node),tree.time(par)
                if m.node not in cache:
                    cache[m.node]=np.sort([tree.time(u) for u in tree.nodes(root=m.node)
                                           if tree.is_internal(u)])
                kv=nl-np.searchsorted(cache[m.node], GRID, side='right')
                p=np.zeros(len(GRID)); nd=np.linspace(tc,tp,NQ)
                iA=np.flatnonzero((GRID<tc)&(nTv>=2)&(kv>=1)&(kv<nTv))
                if len(iA):
                    key=kv[iA]*100+nTv[iA]
                    for u in np.unique(key):
                        sel=iA[key==u]; k,nT=int(u//100),int(u%100)
                        p[sel]=P.integrate(k,nT,nd[None,:]-GRID[sel][:,None],W)
                for iT in np.flatnonzero((GRID>=tc)&(GRID<tp)&(nTv>=2)):
                    q=np.linspace(GRID[iT],tp,NQ)
                    p[iT]=float(P.integrate(1,nTv[iT],q-GRID[iT],W))*(tp-GRID[iT])/(tp-tc)
                r=np.clip(EPS+(1-2*EPS)*p,1e-15,1-1e-15)
                LL.append((np.log(r) if int(site.position) in carried
                           else np.log1p(-r)).astype(np.float32))
                POS.append(site.position)
        LL=np.array(LL); POS=np.array(POS); tot=LL.sum(0)
        blk=np.minimum((POS/(ts.sequence_length/NBLOCK)).astype(int), NBLOCK-1)
        sums=np.zeros((NBLOCK,len(GRID)))
        for b in range(NBLOCK):
            s=blk==b
            if s.any(): sums[b]=LL[s].sum(0)
        rng=np.random.default_rng(0)
        bm=np.array([GRID[sums[rng.integers(0,NBLOCK,NBLOCK)].sum(0).argmax()]
                     for _ in range(NBOOT)])
        MAP=float(GRID[tot.argmax()])
        out=dict(row, MAP=MAP, boot_sd=float(bm.std(ddof=1)),
                 lo=float(np.percentile(bm,2.5)), hi=float(np.percentile(bm,97.5)),
                 n_used=LL.shape[0], err=MAP-age, wall=round(time.time()-t0,1))
        print(f"  {sim}  Ne={Ne:>6,} n={nmod:>2} age={age:>7.0f}  MAP={MAP:>7.0f}  "
              f"err={MAP-age:>+7.0f}  {out['wall']:>6.1f}s", flush=True)
        return out
    except Exception as e:
        print(f"  !! {sim} FAILED {type(e).__name__}: {e}", flush=True)
        return dict(row, error=f"{type(e).__name__}: {e}")

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--simdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--procs", type=int, default=12)
    ap.add_argument("--limit", type=int, default=0)
    a=ap.parse_args()
    rows=list(csv.DictReader(open(os.path.join(a.simdir,"parameters.tsv")), delimiter='\t'))
    rows=[r for r in rows if not r.get("error")]
    if a.limit: rows=rows[:a.limit]
    print(f"[est] {len(rows)} simulations, {a.procs} processes", flush=True)
    t0=time.time()
    with mp.Pool(a.procs, maxtasksperchild=2) as pool:
        res=pool.map(one, [(r,a.simdir) for r in rows], chunksize=1)
    cols=list(rows[0].keys())+["MAP","boot_sd","lo","hi","n_used","err","wall","error"]
    with open(a.out,"w") as f:
        f.write("\t".join(cols)+"\n")
        for r in res: f.write("\t".join(str(r.get(c,"")) for c in cols)+"\n")
    ok=[r for r in res if "error" not in r]
    print(f"[est] {len(ok)}/{len(res)} in {(time.time()-t0)/60:.1f} min -> {a.out}", flush=True)
