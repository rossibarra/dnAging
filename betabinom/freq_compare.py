"""True vs estimated allele frequency at time T, for both estimators.

The panel-only tree (26 modern haplotypes) is all either estimator may use.  The
TRUE population frequency at T is measured from NREF reference chromosomes sampled
at T, which is why this needs a purpose-built simulation rather than an existing
tree: the existing trees contain only time-0 samples.

  x         true derived frequency among the reference samples at T
  blue      betabinom  phi(k, n_T, a) marginalised over the mutation's edge
  red       diffusion  E[p_T | d_0, t_i] marginalised over the same edge
"""
import os, sys, importlib.util, numpy as np, msprime, tskit
from phid import PhiD
spec=importlib.util.spec_from_file_location("pre",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir,
                 "precompute_freq_trajectory_moments.py"))
pre=importlib.util.module_from_spec(spec); spec.loader.exec_module(pre)

NE, MU, RR, L = 50_000, 1e-8, 1e-8, 5e5
NMOD, NREF, NTIME, NSITE = 26, 200, 100, 300
SEED = 4242
tau = lambda t: t/(2.0*NE)

rng = np.random.default_rng(SEED)
times = np.sort(rng.uniform(100, 10_000, NTIME))
sets = [msprime.SampleSet(NMOD, time=0, ploidy=1)] + \
       [msprime.SampleSet(NREF, time=float(t), ploidy=1) for t in times]
print(f"[sim] {NMOD} modern + {NTIME}x{NREF} reference samples, {L/1e6:g} Mb", flush=True)
ts = msprime.sim_ancestry(samples=sets, population_size=NE, sequence_length=L,
                          recombination_rate=RR, ploidy=2, random_seed=SEED)
ts = msprime.sim_mutations(ts, rate=MU, random_seed=SEED+1)
ts = ts.delete_sites([s.id for s in ts.sites() if len(s.mutations)!=1])
print(f"[sim] {ts.num_samples:,} samples, {ts.num_sites:,} biallelic sites", flush=True)

smp = list(ts.samples())
node_time = np.array([ts.node(s).time for s in smp])
mod_rows = np.flatnonzero(node_time == 0.0)
ref_rows = {float(t): np.flatnonzero(node_time == t) for t in times}

# true frequency at each T, per site
truth = {}
for v in ts.variants():
    g = v.genotypes
    truth[v.site.position] = np.array([(g[ref_rows[float(t)]] != 0).mean() for t in times])
print(f"[sim] truth measured at {NTIME} times", flush=True)

ts_mod = ts.simplify(samples=[smp[i] for i in mod_rows])
P = PhiD(np.geomspace(0.5, 4e6, 80), NMOD, NE)

# diffusion table E[p_T | d0, t_i] on the same time grid
tig = np.geomspace(30., 4e6, 44); eng = pre.MomentEngine(NMOD); eps0 = 1.0/(2*NE)
dtab = np.zeros((NMOD, len(tig), NTIME))
for id0, d0 in enumerate(range(1, NMOD+1)):
    for ia, ti in enumerate(tig):
        for iT, T in enumerate(times):
            if T < ti: dtab[id0, ia, iT] = eng.Emoments(d0, tau(ti), tau(T), eps0)[0]
dtab = np.nan_to_num(dtab); ltig = np.log(tig)
print(f"[sim] diffusion table built", flush=True)

NQ = 12; W = np.full(NQ, 1.); W[0] = W[-1] = .5
rows = []
sites = list(ts_mod.sites())
pick = set(rng.choice(len(sites), min(NSITE, len(sites)), replace=False))
si = -1
for tree in ts_mod.trees():
    internal = np.sort([tree.time(u) for u in tree.nodes() if tree.is_internal(u)])
    nTv = NMOD - np.searchsorted(internal, times, side='right')
    for site in tree.sites():
        si += 1
        if si not in pick: continue
        if len(site.mutations) != 1: continue
        m = site.mutations[0]; par = tree.parent(m.node)
        if par == tskit.NULL: continue
        d0 = tree.num_samples(m.node)
        if not (0 < d0 < NMOD): continue
        tc, tp = tree.time(m.node), tree.time(par)
        sub = np.sort([tree.time(u) for u in tree.nodes(root=m.node) if tree.is_internal(u)])
        kv = d0 - np.searchsorted(sub, times, side='right')
        tr = truth.get(site.position)
        if tr is None: continue
        for iT, T in enumerate(times):
            nT = nTv[iT]
            if nT < 2: continue
            if T >= tp:                                   # mutation postdates T
                bb = df = 0.0
            else:
                lo = max(T, tc); w = (tp-lo)/(tp-tc)
                nd = np.linspace(lo, tp, NQ)
                k = 1 if T >= tc else kv[iT]
                if k < 1 or k >= nT: continue
                bb = float(P.integrate(k, nT, nd-T, W))*w
                x = np.interp(np.log(np.clip(nd, tig[0], tig[-1])), ltig, np.arange(len(tig)))
                i0 = np.clip(x.astype(int), 0, len(tig)-2); fr = x-i0
                vals = (1-fr)*dtab[d0-1, i0, iT] + fr*dtab[d0-1, i0+1, iT]
                df = float((vals*W).sum()/W.sum())*w
            rows.append((tr[iT], bb, df, T, k if T < tp else 0, nT))
A = np.array(rows)
np.save("freq_compare.npy", A)
print(f"[sim] {len(A):,} (SNP, T) pairs -> freq_compare.npy", flush=True)
for j, lab in ((1, "betabinom"), (2, "diffusion")):
    e = A[:, j]-A[:, 0]
    print(f"  {lab:<11} bias {e.mean():+.4f}  RMSE {np.sqrt((e**2).mean()):.4f}  "
          f"r {np.corrcoef(A[:,0], A[:,j])[0,1]:.3f}")
