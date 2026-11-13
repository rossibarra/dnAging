#!/usr/bin/env python3
"""msprime validation of the non-sharing-ratio estimator (pipeline/nonsharing_ratio.py).

Each replicate simulates `n_panel` modern haploid panel haplotypes (these build the
ARG in the real pipeline), `n_ctrl` held-out modern controls, and one ancient
haploid at each requested age, all under one coalescent history.

Site selection uses either the true mutation age (`true`) or the age of the
midpoint of the mutation's edge in the tree sequence simplified to the panel
(`edge`, a stand-in for a dated ARG).  Inferred ARGs (SINGER + POLEGON) are the
next step: pass their per-site dated ages through `estimate_ages` the same way.

Run from the repo root, through hpc_run, e.g.
  HPC_MEM=8G ~/.claude/bin/hpc_run 'python validation/nonsharing_ratio_validation.py --out /tmp/x.tsv'
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import msprime
import numpy as np
import tskit

sys.path[:0] = [str(Path(__file__).resolve().parents[1] / "pipeline")]
from nonsharing_ratio import (  # noqa: E402
    estimate_ages, lineage_counts, rate_from_lineage_counts,
)


def simulate(ne, length, mu, rho, n_panel, n_ctrl, ages, seed):
    """Return ts and the sample-column indices of panel, controls, ancient."""
    # ploidy=1 with population_size=2*Ne gives pairwise rate c = 1/(2 Ne).
    sets = [msprime.SampleSet(n_panel, time=0, ploidy=1),
            msprime.SampleSet(n_ctrl, time=0, ploidy=1)]
    sets += [msprime.SampleSet(1, time=float(a), ploidy=1) for a in ages]
    ts = msprime.sim_ancestry(
        samples=sets, ploidy=1, population_size=2 * ne,
        sequence_length=length, recombination_rate=rho,
        discrete_genome=True, random_seed=seed)
    ts = msprime.sim_mutations(
        ts, rate=mu, model=msprime.BinaryMutationModel(),
        discrete_genome=True, random_seed=seed + 1)
    return ts


def site_ages(ts, panel_ids):
    """True (oldest mutation) and edge-midpoint ages per site, NaN where absent."""
    n = ts.num_sites
    true_age = np.full(n, np.nan)
    for m in ts.mutations():
        t = m.time
        if np.isnan(true_age[m.site]) or t > true_age[m.site]:
            true_age[m.site] = t
    sub = ts.simplify(panel_ids, filter_sites=False)
    edge_mid = np.full(n, np.nan)
    times = sub.tables.nodes.time
    for tree in sub.trees():
        for site in tree.sites():
            for m in site.mutations:
                p = tree.parent(m.node)
                if p == tskit.NULL:
                    continue
                mid = 0.5 * (times[m.node] + times[p])
                if np.isnan(edge_mid[site.id]) or mid > edge_mid[site.id]:
                    edge_mid[site.id] = mid
    return true_age, edge_mid


def run_replicate(args, rep):
    ages = np.array(args.ages)
    seed = args.seed + 1000 * rep
    ts = simulate(args.ne, args.length, args.mu, args.rho,
                  args.n_panel, args.n_ctrl, ages, seed)
    geno = ts.genotype_matrix()  # (sites, samples), 0/1
    p = args.n_panel
    panel = geno[:, :p].T
    ctrl = geno[:, p:p + args.n_ctrl].T
    anc = geno[:, p + args.n_ctrl:].T
    true_age, edge_mid = site_ages(ts, list(range(p)))

    grid = np.linspace(0, args.grid_max, args.grid_n + 1)
    c_true = np.full(args.grid_n, 1.0 / (2 * args.ne))
    # c(t) estimated from the panel's own lineage counts
    panel_ts = ts.simplify(list(range(p)))
    kappa = lineage_counts(panel_ts, grid)
    c_est = rate_from_lineage_counts(grid, kappa)
    edges = np.concatenate([grid, [np.inf]])
    rows = []
    for sel_name, sel_ages in (("true", true_age), ("edge", edge_mid)):
        for c_name, c in (("c_true", c_true), ("c_arg", c_est)):
            rates = np.concatenate([c, [c[-1] if c[-1] > 0 else 1.0 / (2 * args.ne)]])
            out = estimate_ages(panel, anc, ctrl, sel_ages, args.a0, edges, rates)
            for a, th, h in zip(ages, out["T_hat"], out["H"]):
                rows.append((rep, sel_name, c_name, a, th, h, out["n_sites"]))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--ne", type=float, default=20000)
    ap.add_argument("--length", type=float, default=20e6)
    ap.add_argument("--mu", type=float, default=1.25e-8)
    ap.add_argument("--rho", type=float, default=1e-8)
    ap.add_argument("--n-panel", type=int, default=26)
    ap.add_argument("--n-ctrl", type=int, default=8)
    ap.add_argument("--ages", type=float, nargs="+",
                    default=[250, 500, 1000, 2000, 3000, 4000, 6000, 8000])
    ap.add_argument("--a0", type=float, default=12000)
    ap.add_argument("--grid-max", type=float, default=20000)
    ap.add_argument("--grid-n", type=int, default=100)
    args = ap.parse_args()
    rows = []
    for rep in range(args.reps):
        rows += run_replicate(args, rep)
    with open(args.out, "w") as fh:
        fh.write("rep\tselection\tc\ttrue_age\tT_hat\tH\tn_sites\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    # compact summary: mean bias per (selection, c, age)
    import collections
    agg = collections.defaultdict(list)
    for rep, s, c, a, th, h, ns in rows:
        agg[(s, c, a)].append(th - a)
    print("selection\tc\ttrue_age\tmean_bias\tse\tn")
    for (s, c, a), v in sorted(agg.items()):
        v = np.array(v)
        se = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else float("nan")
        print(f"{s}\t{c}\t{a:g}\t{v.mean():.0f}\t{se:.0f}\t{len(v)}")


if __name__ == "__main__":
    main()
