"""Test the proposed ARG-based likelihood against coalescent simulation.

Proposal: for an ancient pseudohaploid sampled at time T, and a mutation that
exists at T, P(carries) is a beta-binomial predictive given (k, n_T) read off
the modern-only tree at time T:
    k   = ancestral lineages at T that descend from the mutation
    n_T = ancestral lineages at T in total
Candidates: k/n_T  vs  k/(n_T+1).
"""
import msprime, tskit, numpy as np
from collections import defaultdict

NE, T, NMOD = 10_000, 2_000.0, 26
REPS, L, R, MU = 60, 2e6, 1e-8, 1e-8

def lineages_at(tree, T):
    """nodes whose edge spans time T (i.e. the ancestral lineages at T)."""
    out=[]
    for u in tree.nodes():
        p = tree.parent(u)
        if p == tskit.NULL: continue
        if tree.time(u) <= T < tree.time(p): out.append(u)
    return out

hits = defaultdict(lambda: [0,0])       # (k,n) -> [carried, total]
zero_below = [0,0]                      # mutations younger than T
for rep in range(REPS):
    ts = msprime.sim_ancestry(
        samples=[msprime.SampleSet(NMOD, time=0, ploidy=1),
                 msprime.SampleSet(1, time=T, ploidy=1)],
        population_size=NE, sequence_length=L, recombination_rate=R,
        ploidy=2, random_seed=1000+rep)
    ts = msprime.sim_mutations(ts, rate=MU, random_seed=5000+rep, discrete_genome=False)
    anc = ts.num_samples-1                       # the ancient sample id
    mod = [s for s in ts.samples() if s != anc]
    ts_mod = ts.simplify(samples=mod, filter_sites=False)

    gt = {v.site.position: v.genotypes for v in ts.variants()}
    anc_idx = list(ts.samples()).index(anc)

    for tree_m, site_list in zip(ts_mod.trees(), None or [None]*ts_mod.num_trees):
        pass
    # walk sites in the modern-only ts, pairing with the full-ts genotypes
    tm = ts_mod.first()
    lin_cache = {}
    for site in ts_mod.sites():
        while tm.interval.right <= site.position: tm.next()
        if len(site.mutations) != 1: continue
        mut = site.mutations[0]
        t_mut = mut.time
        g = gt.get(site.position)
        if g is None: continue
        carried = int(g[anc_idx] != 0)
        if t_mut < T:                     # mutation younger than the sample
            zero_below[1]+=1; zero_below[0]+=carried; continue
        key=(tm.index,)
        if key not in lin_cache: lin_cache[key]=lineages_at(tm, T)
        lin = lin_cache[key]
        n_T = len(lin)
        if n_T == 0: continue
        c = mut.node
        if tm.time(c) > T:                # T is below the mutation's child node
            k = sum(1 for u in lin if u==c or (c in list(tm.ancestors(u))))
        else:                             # edge straddles T -> single lineage
            k = 1
        if k==0 or k>n_T: continue
        h=hits[(k,n_T)]; h[1]+=1; h[0]+=carried

print(f"mutations younger than T: carried {zero_below[0]} / {zero_below[1]}"
      f"  (model predicts 0)")
print(f"\n{'k':>3} {'n_T':>4} {'trials':>7} {'observed':>9} {'k/n_T':>8} {'k/(n_T+1)':>10}")
tot_n=tot_n1=tot_obs=0
for (k,n),(c,t) in sorted(hits.items()):
    if t < 400: continue
    obs=c/t
    print(f"{k:>3} {n:>4} {t:>7} {obs:>9.4f} {k/n:>8.4f} {k/(n+1):>10.4f}")
    tot_obs+=c; tot_n+=t*k/n; tot_n1+=t*k/(n+1)
print(f"\npooled over shown bins: observed {tot_obs:.0f}"
      f" | k/n_T predicts {tot_n:.0f} | k/(n_T+1) predicts {tot_n1:.0f}")
