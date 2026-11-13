#!/usr/bin/env python3
"""Does including panel-fixed sites remove the young bias of the non-sharing ratio?

Design rule 2 of the report: sites derived in every panel haplotype are older than the
panel's local common ancestor, hence older than any cutoff A0, and must be included
undated.  Using only panel-segregating sites (what the CLI sees if the ARG's site list
lacks the fixed sites) drops the oldest part of the allele-age distribution, which biases
old ancient samples young.

Arms (same simulated replicates, same estimator, c(t) = true 1/(2 Ne)):
  with_fixed     panel-fixed derived sites included (rule 2)
  no_fixed       only panel-segregating sites older than A0
  wrong_pol_q    with_fixed plus a fraction q of mislabelled sites: panel-monomorphic
                 *ancestral* sites at which some target carries a new derived allele are
                 treated as 'panel fixed derived' (what a wrong polarity fallback does)

Run from the repo root through hpc_run, e.g.
  HPC_MEM=8G HPC_TIME=01:00:00 ~/.claude/bin/hpc_run \
    'python validation/nonsharing_ratio_fixedsites_validation.py --out results/nonsharing_ratio/fixedsites.tsv'
"""
from __future__ import annotations

import argparse
import collections
import sys
from pathlib import Path

import numpy as np

sys.path[:0] = [str(Path(__file__).resolve().parents[1] / "pipeline"),
                str(Path(__file__).resolve().parent)]
from nonsharing_ratio import estimate_ages  # noqa: E402
from nonsharing_ratio_validation import simulate, site_ages  # noqa: E402


def replicate(args, rep):
    ages = np.array(args.ages)
    ts = simulate(args.ne, args.length, args.mu, args.rho, args.n_panel, args.n_ctrl,
                  ages, args.seed + 1000 * rep)
    geno = ts.genotype_matrix()
    p = args.n_panel
    panel = geno[:, :p].T
    ctrl = geno[:, p:p + args.n_ctrl].T
    anc = geno[:, p + args.n_ctrl:].T
    true_age, edge_mid = site_ages(ts, list(range(p)))
    sel_ages = true_age if args.selection == "true" else edge_mid
    edges = np.array([0.0, args.grid_max, np.inf])
    rates = np.array([1.0 / (2 * args.ne)] * 2)
    count = panel.sum(axis=0)
    n_fixed = int((count == p).sum())
    n_seg_old = int(((count > 0) & (count < p) & (sel_ages > args.a0)).sum())
    rng = np.random.default_rng(args.seed + rep)
    rows = []

    def add(arm, out):
        for a, th, h in zip(ages, out["T_hat"], out["H"]):
            rows.append((rep, arm, a, th, h, out["n_sites"], n_fixed, n_seg_old))

    add("with_fixed", estimate_ages(panel, anc, ctrl, sel_ages, args.a0, edges, rates, True))
    add("no_fixed", estimate_ages(panel, anc, ctrl, sel_ages, args.a0, edges, rates, False))
    # wrong polarity: panel-monomorphic-ancestral sites that some target carries a new allele at.
    # Under a wrong label the panel's (ancestral) allele is 'derived' in the panel and a target
    # carrying the new mutation looks ancestral, i.e. non-sharing.
    cand = np.nonzero((count == 0) & ((ctrl.sum(0) + anc.sum(0)) > 0))[0]
    for q in args.polarity_error:
        pick = cand[rng.random(cand.size) < q]
        if pick.size == 0:
            continue
        panel2 = np.concatenate([panel, np.ones((p, pick.size), dtype=panel.dtype)], axis=1)
        ctrl2 = np.concatenate([ctrl, 1 - ctrl[:, pick]], axis=1)
        anc2 = np.concatenate([anc, 1 - anc[:, pick]], axis=1)
        ages2 = np.concatenate([sel_ages, np.full(pick.size, np.inf)])
        add(f"wrong_pol_{q:g}", estimate_ages(panel2, anc2, ctrl2, ages2, args.a0, edges, rates, True))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--seed", type=int, default=777)
    ap.add_argument("--ne", type=float, default=20000)
    ap.add_argument("--length", type=float, default=20e6)
    ap.add_argument("--mu", type=float, default=1.25e-8)
    ap.add_argument("--rho", type=float, default=1e-8)
    ap.add_argument("--n-panel", type=int, default=26)
    ap.add_argument("--n-ctrl", type=int, default=8)
    ap.add_argument("--ages", type=float, nargs="+", default=[250, 1000, 2000, 4000, 6000, 8000])
    ap.add_argument("--a0", type=float, default=12000)
    ap.add_argument("--grid-max", type=float, default=20000)
    ap.add_argument("--selection", choices=["true", "edge"], default="edge")
    ap.add_argument("--polarity-error", type=float, nargs="*", default=[0.01, 0.1])
    args = ap.parse_args()
    rows = []
    for rep in range(args.reps):
        rows += replicate(args, rep)
    with open(args.out, "w") as fh:
        fh.write("rep\tarm\ttrue_age\tT_hat\tH\tn_sites\tn_panel_fixed\tn_seg_old\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    agg = collections.defaultdict(list)
    for rep, arm, a, th, h, *_ in rows:
        agg[(arm, a)].append(th - a)
    print("arm\ttrue_age\tmean_bias\tse\tn")
    for (arm, a), v in sorted(agg.items(), key=lambda kv: (kv[0][0], kv[0][1])):
        v = np.array(v)
        se = v.std(ddof=1) / np.sqrt(len(v)) if len(v) > 1 else float("nan")
        print(f"{arm}\t{a:g}\t{v.mean():.0f}\t{se:.0f}\t{len(v)}")


if __name__ == "__main__":
    main()
