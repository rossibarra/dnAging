#!/usr/bin/env python3
"""How much of the age-estimate variance is the ancient lineage draw itself?

Block bootstrap takes coverage from 0.46 to 0.88, but the residual sits almost
entirely in old samples (0.82-0.85 above 4,000 generations against 0.94-0.95
below).  Every block of one genome shares the same realised ancient lineage and
the same modern panel, so any within-genome resampling scheme is blind to
variance carried by those single draws -- which is also why 0.2 cM and 2 cM
blocks gave identical coverage.

This isolates that component.  Each simulation places N ancient haploids at the
*same* age in the *same* tree sequence, so across those N estimates the tree,
the panel, the demography and the true age are all held fixed and the only thing
that varies is which lineage was drawn.  The spread of T-hat across them is
therefore the lineage-draw standard deviation, directly measured.

Read it against two numbers already in hand: the across-replicate error SD from
batch300 (~1,200 generations in the 4,000-10,000 band), and the block-bootstrap
interval half-width (~1,340).  If the lineage spread is comparable to the former,
the missing variance is the lineage draw and no site-level or block-level
resampling can recover it.  If it is small at every age, the hypothesis is dead.

**One traversal serves every lineage.**  `derived_probability_uniform_edge_grid`
depends on the tree, the mutation's edge and the candidate age -- not on which
ancient sample is being scored.  So the per-site probability vector is computed
once and accumulated into N separate log-likelihoods according to each lineage's
genotype, which is N times cheaper than calling the production driver N times and
guarantees every lineage sees an identical likelihood surface.

Per-block curves are retained so each lineage also gets its own block-bootstrap
interval, making within-tree spread and per-lineage interval width directly
comparable.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import tskit

sys.path.insert(0, str(Path(__file__).resolve().parent))
from insertion_likelihood import (PiecewiseConstantNe,
                                  derived_probability_uniform_edge_grid)
from msprime_insertion_validation import vcf_call_matrix
from posterior_sample_age_infer import summarize


def run(args):
    directory = args.replicate_dir
    metadata = json.loads((directory / "metadata.json").read_text())
    n_modern = int(metadata["n_modern_haploid"])
    ts = tskit.load(directory / "known_modern_arg.trees")
    names, matrix = vcf_call_matrix(directory / "all_samples.vcf.gz", n_modern)
    n_anc = len(names)
    # Ne source, in order: an explicit --ne-table; the scalar some simulators
    # record in metadata; otherwise the replicate's own epoch table, which
    # simulate_variable_ne_error.py always writes (it records `epochs`, not a
    # scalar `Ne`, so the scalar lookup alone would fail on its output).
    if args.ne_table:
        ne = PiecewiseConstantNe.from_tsv(args.ne_table, args.ne_series, False)
    elif "Ne" in metadata:
        ne = float(metadata["Ne"])
    else:
        epochs = directory / "constant_ne_epochs.tsv"
        if not epochs.exists():
            raise SystemExit(f"no Ne in metadata and no {epochs}")
        import csv as _csv
        with epochs.open() as h:
            series = {r.get("series") for r in _csv.DictReader(h, delimiter="\t")}
        series.discard(None)
        ne = PiecewiseConstantNe.from_tsv(
            epochs, series.pop() if len(series) == 1 else None, False)
    grid = np.arange(args.age_min, args.age_max + 0.5 * args.age_step,
                     args.age_step, dtype=float)

    edges = np.linspace(0.0, float(ts.sequence_length), args.n_blocks + 1)
    # (lineage, block, age): one traversal fills all of it.
    block_ll = np.zeros((n_anc, args.n_blocks, len(grid)), dtype=float)
    used = 0

    for site in ts.sites():
        if len(site.mutations) != 1:
            continue
        mutation = site.mutations[0]
        d0, genotypes = matrix[site.id]
        tree = ts.at(site.position)
        if int(tree.num_samples(mutation.node)) != d0:
            continue
        parent = tree.parent(mutation.node)
        if parent == tskit.NULL:
            continue
        probability = derived_probability_uniform_edge_grid(
            tree, mutation.node, grid,
            float(ts.node(mutation.node).time),
            float(ts.node(parent).time), ne)
        if not np.all(np.isfinite(probability)):
            continue
        probability = np.clip(probability, 1e-300, 1 - 1e-15)
        observed = np.clip(args.epsilon + (1 - 2 * args.epsilon) * probability,
                           1e-300, 1 - 1e-15)
        carry = np.log(observed)
        absent = np.log1p(-observed)
        b = min(int(np.searchsorted(edges, site.position, side="right") - 1),
                args.n_blocks - 1)
        g = np.asarray(genotypes, dtype=bool)
        block_ll[g, b] += carry
        block_ll[~g, b] += absent
        used += 1

    rng = np.random.default_rng(args.seed)
    truths = metadata["true_ancient_ages"]
    rows = []
    for i, name in enumerate(names):
        total = block_ll[i].sum(axis=0)
        summary, _ = summarize(grid, total)
        w = rng.multinomial(args.n_blocks,
                            np.full(args.n_blocks, 1.0 / args.n_blocks),
                            size=args.n_boot).astype(float)
        boots = grid[np.argmax(w @ block_ll[i], axis=1)]
        lo, hi = (float(x) for x in np.percentile(boots, [2.5, 97.5]))
        truth = float(truths[i])
        rows.append({"lineage": name, "true_T": truth,
                     "map_T": summary["map_T"], "mean_T": summary["mean_T"],
                     "boot_lo": lo, "boot_hi": hi,
                     "boot_contains_true": int(lo <= truth <= hi),
                     "boot_width": hi - lo})

    maps = np.array([r["map_T"] for r in rows])
    truth0 = float(truths[0])
    out = {"replicate_dir": str(directory), "n_lineages": n_anc,
           "sites_used": used, "true_age": truth0,
           "all_ages_equal": bool(len(set(round(t, 6) for t in truths)) == 1),
           "lineage_map_estimates": maps.tolist(),
           "lineage_spread_SD": float(maps.std(ddof=1)),
           "lineage_spread_IQR": float(np.subtract(*np.percentile(maps, [75, 25]))),
           "lineage_mean_error": float(maps.mean() - truth0),
           "median_bootstrap_width": float(np.median([r["boot_width"] for r in rows])),
           "bootstrap_coverage_within_tree": float(
               np.mean([r["boot_contains_true"] for r in rows])),
           "per_lineage": rows}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    np.savez_compressed(args.output.with_suffix(".blocks.npz"),
                        Tgrid=grid, block_log_likelihood=block_ll)
    print(json.dumps({k: v for k, v in out.items() if k != "per_lineage"},
                     indent=2), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--replicate-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--n-blocks", type=int, default=50)
    p.add_argument("--n-boot", type=int, default=1000)
    p.add_argument("--epsilon", type=float, default=0.0)
    p.add_argument("--age-min", type=float, default=0.0)
    p.add_argument("--age-max", type=float, default=15000.0)
    p.add_argument("--age-step", type=float, default=50.0)
    p.add_argument("--ne-table", type=Path, default=None)
    p.add_argument("--ne-series", default=None)
    p.add_argument("--seed", type=int, default=17)
    run(p.parse_args())
