"""100 simple 1 Mb simulations, sample of 20, betabinom vs diffusion vs truth
at 5 random times per tree.

The population frequency at time T is not recoverable from a time-0 sample, so
each simulation also draws NREF reference chromosomes AT each of its 5 times; the
derived proportion among those is the truth on the x axis.  The two estimators see
only the 20-haplotype panel tree.
"""
import os, sys, csv, time, importlib.util, argparse
import multiprocessing as mp
import numpy as np, msprime, tskit
from phid import PhiD
spec=importlib.util.spec_from_file_location("pre",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir,
                 "precompute_freq_trajectory_moments.py"))
pre=importlib.util.module_from_spec(spec); spec.loader.exec_module(pre)

NE, MU, RR, L = 10_000, 1e-8, 1e-8, 1e6
NPANEL, NREF, NTIME = 20, 400, 5
NQ=12; W=np.full(NQ,1.); W[0]=W[-1]=.5
tau=lambda t: t/(2.0*NE)

def one(i):
    rng=np.random.default_rng(70000+i)
    times=np.sort(rng.uniform(100, 10_000, NTIME))
    sets=[msprime.SampleSet(NPANEL, time=0, ploidy=1)] + \
         [msprime.SampleSet(NREF, time=float(t), ploidy=1) for t in times]
    ts=msprime.sim_ancestry(samples=sets, population_size=NE, sequence_length=L,
                            recombination_rate=RR, ploidy=2,
                            random_seed=int(rng.integers(1,2**31-1)))
    ts=msprime.sim_mutations(ts, rate=MU, random_seed=int(rng.integers(1,2**31-1)))
    ts=ts.delete_sites([s.id for s in ts.sites() if len(s.mutations)!=1])
    smp=list(ts.samples()); nt=np.array([ts.node(s).time for s in smp])
    mod=np.flatnonzero(nt==0.0)
    ref={float(t): np.flatnonzero(nt==t) for t in times}
    truth={}
    for v in ts.variants():
        g=v.genotypes
        truth[v.site.position]=np.array([(g[ref[float(t)]]!=0).mean() for t in times])
    ts_mod=ts.simplify(samples=[smp[j] for j in mod])
    P=PhiD(np.geomspace(0.5,4e6,60), NPANEL, NE, dps=45)
    tig=np.geomspace(30.,4e6,40); eng=pre.MomentEngine(NPANEL); eps0=1.0/(2*NE)
    dtab=np.zeros((NPANEL,len(tig),NTIME))
    for id0,d0 in enumerate(range(1,NPANEL+1)):
        for ia,ti in enumerate(tig):
            for iT,T in enumerate(times):
                if T<ti: dtab[id0,ia,iT]=eng.Emoments(d0,tau(ti),tau(T),eps0)[0]
    dtab=np.nan_to_num(dtab); ltig=np.log(tig)
    out=[]
    for tree in ts_mod.trees():
        internal=np.sort([tree.time(u) for u in tree.nodes() if tree.is_internal(u)])
        nTv=NPANEL-np.searchsorted(internal,times,side='right')
        for site in tree.sites():
            if len(site.mutations)!=1: continue
            m=site.mutations[0]; par=tree.parent(m.node)
            if par==tskit.NULL: continue
            d0=tree.num_samples(m.node)
            if not (0<d0<NPANEL): continue
            tc,tp=tree.time(m.node),tree.time(par)
            sub=np.sort([tree.time(u) for u in tree.nodes(root=m.node) if tree.is_internal(u)])
            kv=d0-np.searchsorted(sub,times,side='right')
            tr=truth.get(site.position)
            if tr is None: continue
            for iT,T in enumerate(times):
                nT=nTv[iT]
                if nT<2: continue
                if T>=tp: bb=df=0.0
                else:
                    lo=max(T,tc); w=(tp-lo)/(tp-tc); nd=np.linspace(lo,tp,NQ)
                    k=1 if T>=tc else kv[iT]
                    if k<1 or k>=nT: continue
                    bb=float(P.integrate(k,nT,nd-T,W))*w
                    x=np.interp(np.log(np.clip(nd,tig[0],tig[-1])),ltig,np.arange(len(tig)))
                    i0=np.clip(x.astype(int),0,len(tig)-2); fr=x-i0
                    vals=(1-fr)*dtab[d0-1,i0,iT]+fr*dtab[d0-1,i0+1,iT]
                    df=float((vals*W).sum()/W.sum())*w
                out.append((tr[iT],bb,df,T,d0,nT))
    return np.array(out,dtype=np.float32)

if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("--n",type=int,default=100)
    ap.add_argument("--procs",type=int,default=4); a=ap.parse_args()
    t0=time.time()
    with mp.Pool(a.procs, maxtasksperchild=4) as pool:
        parts=[]
        for j,r in enumerate(pool.imap_unordered(one, range(a.n)), 1):
            parts.append(r)
            if j%10==0: print(f"  {j}/{a.n} sims, {sum(len(p) for p in parts):,} pairs, "
                              f"{(time.time()-t0)/60:.1f} min", flush=True)
    A=np.concatenate(parts); np.save("freq_compare2.npy", A)
    print(f"[done] {len(A):,} (SNP,T) pairs from {a.n} sims -> freq_compare2.npy")
    for j,lab in ((1,"betabinom"),(2,"diffusion")):
        e=A[:,j]-A[:,0]
        print(f"  {lab:<11} bias {e.mean():+.4f}  RMSE {np.sqrt((e**2).mean()):.4f}  "
              f"r {np.corrcoef(A[:,0],A[:,j])[0,1]:.3f}")
