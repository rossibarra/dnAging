"""Same 300 simulations, but with MATH.md's DIFFUSION likelihood: p is
E[p_T | d_0, t_i] marginalised over the mutation's edge, conditioning on the
present-day panel count rather than on the tree at T.

The three matrix exponentials in Emoments depend only on (n, tau_i, tau_T) and NOT
on d_0, so they are hoisted out of the d_0 loop -- the optimisation MATH.md notes
is unimplemented.  Without it this would be ~40x slower.
"""
import os, sys, csv, time, argparse, importlib.util
import diff_mp
import multiprocessing as mp
import numpy as np, tskit
from math import comb
from scipy.linalg import expm
spec=importlib.util.spec_from_file_location("pre",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir,
                 "precompute_freq_trajectory_moments.py"))
pre=importlib.util.module_from_spec(spec); spec.loader.exec_module(pre)

EPS=1e-6; NQ=12; NBLOCK=100; NBOOT=200
GRID=np.linspace(50., 15000., 100)
W=np.full(NQ,1.); W[0]=W[-1]=.5
MAXCANC=1.25e14

def diffusion_table(n, Ne, tig, grid):
    """table[d0-1, i_ti, i_T] = E[p_T | d0, t_i], exponentials hoisted over d0."""
    eng=pre.MomentEngine(n); B=eng.B; K=eng.K; eps=1.0/(2.0*Ne)
    tau=lambda t: t/(2.0*Ne)
    coeff={d0:{m:comb(n-d0,m-d0)*(-1)**(m-d0) for m in range(d0,n+1)}
           for d0 in range(1,n+1)}
    m0=np.array([eps**k for k in range(K+2)]); m0[0]=1.0
    tab=np.zeros((n,len(tig),len(grid)))
    for ia,ti in enumerate(tig):
        ti_t=tau(ti)
        Mpres=expm(B*ti_t)@m0
        for iT,T in enumerate(grid):
            if T>=ti: continue
            tT=tau(T)
            Mu1=expm(B*(ti_t-tT))@m0
            C=expm(B*tT)
            EjX=C[:, :K]@Mu1[1:K+1]
            for d0,cs in coeff.items():
                num=den=anum=aden=0.0
                for m,c in cs.items():
                    t1=c*EjX[m]; td=c*Mpres[m]
                    num+=t1; den+=td; anum+=abs(t1); aden+=abs(td)
                if den<=0 or num<=0: continue
                if aden/abs(den)>MAXCANC or anum/abs(num)>MAXCANC: continue
                v=num/den
                if 0.0<v<1.0: tab[d0-1,ia,iT]=v
    return tab

def one(args):
    row, simdir = args
    sim=row["sim"]; Ne=int(row["Ne"]); nmod=int(row["n_modern"]); age=float(row["age"])
    t0=time.time()
    try:
        ts=tskit.load(os.path.join(simdir, f"{sim}.trees"))
        raw=[]
        for ln in open(os.path.join(simdir, f"{sim}_ancient.vcf")):
            if ln.startswith('#'): continue
            f=ln.split('\t'); raw.append((int(f[1]), f[9].strip().startswith('1')))
        tree_pos=set(int(s_.position) for s_ in ts.sites())
        allp={p for p,_ in raw}
        off=max((0,-1,1), key=lambda o: len({p+o for p in allp} & tree_pos))
        carried={p+off for p,c in raw if c}
        tig=np.geomspace(30., 4e6, 40)
        tab=diff_mp.table(nmod, Ne, tig, GRID, dps=50); ltig=np.log(tig)
        LL=[]; POS=[]
        for tree in ts.trees():
            for site in tree.sites():
                if len(site.mutations)!=1: continue
                m=site.mutations[0]; par=tree.parent(m.node)
                if par==tskit.NULL: continue
                d0=tree.num_samples(m.node)
                if not (0<d0<nmod): continue
                tc,tp=tree.time(m.node),tree.time(par)
                p=np.zeros(len(GRID))
                for iT in np.flatnonzero(GRID<tp):
                    lo=max(GRID[iT],tc); nd=np.linspace(lo,tp,NQ)
                    x=np.interp(np.log(np.clip(nd,tig[0],tig[-1])),ltig,np.arange(len(tig)))
                    i0=np.clip(x.astype(int),0,len(tig)-2); fr=x-i0
                    vals=(1-fr)*tab[d0-1,i0,iT]+fr*tab[d0-1,i0+1,iT]
                    p[iT]=float((vals*W).sum()/W.sum())*(tp-lo)/(tp-tc)
                r=np.clip(EPS+(1-2*EPS)*p,1e-15,1-1e-15)
                LL.append((np.log(r) if int(site.position) in carried
                           else np.log1p(-r)).astype(np.float32))
                POS.append(site.position)
        LL=np.array(LL); POS=np.array(POS); tot=LL.sum(0)
        blk=np.minimum((POS/(ts.sequence_length/NBLOCK)).astype(int),NBLOCK-1)
        sums=np.zeros((NBLOCK,len(GRID)))
        for b in range(NBLOCK):
            s=blk==b
            if s.any(): sums[b]=LL[s].sum(0)
        rng=np.random.default_rng(0)
        bm=np.array([GRID[sums[rng.integers(0,NBLOCK,NBLOCK)].sum(0).argmax()]
                     for _ in range(NBOOT)])
        MAP=float(GRID[tot.argmax()])
        print(f"  {sim}  Ne={Ne:>6,} n={nmod:>2} age={age:>7.0f}  MAP={MAP:>7.0f}  "
              f"err={MAP-age:>+7.0f}  {time.time()-t0:>6.1f}s", flush=True)
        return dict(row, MAP_mp=MAP, boot_sd_mp=float(bm.std(ddof=1)),
                    err_mp=MAP-age, n_used_mp=LL.shape[0])
    except Exception as e:
        print(f"  !! {sim} {type(e).__name__}: {e}", flush=True)
        return dict(row, fail=f"{type(e).__name__}: {e}")

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--simdir",required=True); ap.add_argument("--out",required=True)
    ap.add_argument("--procs",type=int,default=12); ap.add_argument("--limit",type=int,default=0)
    a=ap.parse_args()
    rows=[r for r in csv.DictReader(open(os.path.join(a.simdir,"parameters.tsv")),
                                    delimiter='\t') if not r.get("error")]
    if a.limit: rows=rows[:a.limit]
    print(f"[diff-mp] {len(rows)} sims, {a.procs} procs", flush=True)
    t0=time.time()
    with mp.Pool(a.procs, maxtasksperchild=2) as pool:
        res=pool.map(one, [(r,a.simdir) for r in rows], chunksize=1)
    cols=list(rows[0].keys())+["MAP_mp","boot_sd_mp","err_mp","n_used_mp","fail"]
    with open(a.out,"w") as f:
        f.write("\t".join(cols)+"\n")
        for r in res: f.write("\t".join(str(r.get(c,"")) for c in cols)+"\n")
    print(f"[diff-mp] done in {(time.time()-t0)/60:.1f} min -> {a.out}", flush=True)
