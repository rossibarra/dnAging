#!/usr/bin/env python3
"""Godambe (sandwich) interval calibration for the insertion composite likelihood.

The insertion estimator's point bias is solved, but its intervals are far too
narrow: nominal 95% coverage is 139/300 on the batch300 benchmark.  That is the
signature of a composite likelihood -- linked sites are multiplied as if
independent, so the likelihood is over-peaked by roughly the ratio of the true
information to the naive curvature.

For a scalar T with composite score U(T), sensitivity H = -E[dU/dT] and
variability J = Var(U), the Godambe variance is J/H^2 while the naive curvature
gives 1/H.  The inflation is J/H, and under a correctly specified likelihood the
second Bartlett identity gives H = J so it collapses to the usual answer.

Our intervals are equal-tailed from the normalised exp(logCL) rather than read
off curvature, so the matching fix is a **magnitude adjustment**: rescale

    logCL(T)  ->  c * logCL(T),      c = H / J

Two properties make this the right form here.  Scaling a log-likelihood cannot
move its argmax, so the MAP -- and therefore the point-bias result -- is exactly
preserved.  And it produces a whole calibrated surface rather than a single
standard deviation, which matters because these posteriors can be asymmetric.

**Estimating J needs blocks.**  Sites within a block are linked; blocks far apart
in genetic distance are nearly independent.  This accumulates the per-site log
likelihood into `--n-blocks` contiguous genomic blocks, so the total is exactly
the quantity the production run computes (block sums add back to it), and

    H_hat = -d2/dT2 [ sum_b logCL_b ]          at T_hat
    J_hat = (B/(B-1)) * sum_b (U_b - Ubar)^2   at T_hat,   U_b = d/dT logCL_b

Both derivatives are taken numerically on the age grid, which is dense and
regular, so central differences are adequate.

**This corrects for linkage only.**  It treats the ARG as given, so if the ARG
posterior is itself overconfident no site-level sandwich can see it.  With a true
ARG and exact mutation times the exact-time diffusion run still had 0.68
coverage, implying J/H ~ 4 from linkage alone, against ~10 implied here -- so a
residual gap after correction would point at the ARG term rather than at this
one.

Mirrors `msprime_insertion_validation.infer` site-for-site; the only change is
where the per-site contribution is added.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import tskit

sys.path[:0] = [str(Path(__file__).resolve().parents[1] / d) for d in ("pipeline", "validation")]
from insertion_likelihood import (PiecewiseConstantNe,
                                  derived_probability_uniform_edge_grid)
from msprime_insertion_validation import replicate_dir, vcf_call_matrix
from posterior_sample_age_infer import summarize


def block_loglikelihoods(args):
    """Per-block log-likelihood curves; their sum is the production total."""
    directory = replicate_dir(args.input_root, args.replicate)
    metadata = json.loads((directory / "metadata.json").read_text())
    ts = tskit.load(directory / "known_modern_arg.trees")
    ancient_names, matrix = vcf_call_matrix(
        directory / "all_samples.vcf.gz", args.n_modern)
    if args.ancient_sample not in ancient_names:
        raise ValueError(f"ancient sample {args.ancient_sample!r} not in VCF")
    ancient_index = ancient_names.index(args.ancient_sample)
    calls = {site: (d0, int(g[ancient_index]))
             for site, (d0, g) in matrix.items()}
    ne = (PiecewiseConstantNe.from_tsv(args.ne_table, args.ne_series, False)
          if args.ne_table else float(metadata["Ne"]))
    grid = np.arange(args.age_min, args.age_max + 0.5 * args.age_step,
                     args.age_step, dtype=float)

    edges = np.linspace(0.0, float(ts.sequence_length), args.n_blocks + 1)
    block_ll = np.zeros((args.n_blocks, len(grid)), dtype=float)
    used = np.zeros(args.n_blocks, dtype=np.int64)
    dropped = 0

    for site in ts.sites():
        if len(site.mutations) != 1:
            dropped += 1
            continue
        mutation = site.mutations[0]
        d0, ancient_genotype = calls[site.id]
        tree = ts.at(site.position)
        if int(tree.num_samples(mutation.node)) != d0:
            dropped += 1
            continue
        parent = tree.parent(mutation.node)
        if parent == tskit.NULL:
            dropped += 1
            continue
        probability = derived_probability_uniform_edge_grid(
            tree, mutation.node, grid,
            float(ts.node(mutation.node).time),
            float(ts.node(parent).time), ne)
        if not np.all(np.isfinite(probability)):
            dropped += 1
            continue
        probability = np.clip(probability, 1e-300, 1 - 1e-15)
        observed = np.clip(args.epsilon + (1 - 2 * args.epsilon) * probability,
                           1e-300, 1 - 1e-15)
        term = (np.log(observed) if ancient_genotype
                else np.log1p(-observed))
        b = min(int(np.searchsorted(edges, site.position, side="right") - 1),
                args.n_blocks - 1)
        block_ll[b] += term
        used[b] += 1

    truth = (float(metadata["true_ancient_ages"][ancient_index])
             if "true_ancient_ages" in metadata
             else float(metadata["true_ancient_age"]))
    return grid, block_ll, used, dropped, truth


def godambe_scale(grid, block_ll):
    """Return (c, H, J, T_hat) from the block decomposition.

    H and J are evaluated at the composite MAP. A non-positive or non-finite H
    means the surface is not locally concave there, in which case no scalar
    rescaling is meaningful and c is returned as NaN rather than silently
    clamped into plausibility.
    """
    total = block_ll.sum(axis=0)
    k = int(np.argmax(total))
    step = float(grid[1] - grid[0])
    # Keep one grid point either side of the MAP for a central difference.
    k = min(max(k, 1), len(grid) - 2)
    T_hat = float(grid[k])

    H = -(total[k + 1] - 2.0 * total[k] + total[k - 1]) / step ** 2
    scores = (block_ll[:, k + 1] - block_ll[:, k - 1]) / (2.0 * step)
    B = block_ll.shape[0]
    J = float(np.sum((scores - scores.mean()) ** 2) * B / (B - 1))

    if not np.isfinite(H) or H <= 0 or not np.isfinite(J) or J <= 0:
        return float("nan"), float(H), float(J), T_hat
    return float(H / J), float(H), float(J), T_hat


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-root", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--replicate", type=int, required=True)
    p.add_argument("--n-modern", type=int, required=True)
    p.add_argument("--n-blocks", type=int, default=50)
    p.add_argument("--ancient-sample", default="ancient_01")
    p.add_argument("--epsilon", type=float, default=0.0)
    p.add_argument("--age-min", type=float, default=0.0)
    p.add_argument("--age-max", type=float, default=15000.0)
    p.add_argument("--age-step", type=float, default=50.0)
    p.add_argument("--ne-table", type=Path, default=None)
    p.add_argument("--ne-series", default=None)
    args = p.parse_args()

    grid, block_ll, used, dropped, truth = block_loglikelihoods(args)
    total = block_ll.sum(axis=0)
    c, H, J, T_hat = godambe_scale(grid, block_ll)

    naive, _ = summarize(grid, total)
    out = {"replicate": args.replicate, "true_T": truth,
           "n_blocks": args.n_blocks, "sites_used": int(used.sum()),
           "sites_dropped": int(dropped),
           "empty_blocks": int((used == 0).sum()),
           "T_hat_for_derivatives": T_hat,
           "H": H, "J": J, "godambe_c": c,
           "implied_variance_inflation": (float(J / H) if np.isfinite(c) and c > 0
                                          else float("nan")),
           "naive": naive,
           "naive_contains_true": int(naive["ci95_lower_T"] <= truth
                                      <= naive["ci95_upper_T"])}
    if np.isfinite(c) and c > 0:
        adjusted, _ = summarize(grid, c * total)
        out["adjusted"] = adjusted
        out["adjusted_contains_true"] = int(adjusted["ci95_lower_T"] <= truth
                                            <= adjusted["ci95_upper_T"])
        # A positive rescaling cannot move the argmax; assert rather than
        # assume, since a violation would mean the point result had moved.
        out["map_unchanged"] = bool(
            abs(adjusted["map_T"] - naive["map_T"]) < 1e-9)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        args.output_dir / f"replicate_{args.replicate:03d}.blocks.npz",
        Tgrid=grid, block_log_likelihood=block_ll, sites_per_block=used)
    (args.output_dir / f"replicate_{args.replicate:03d}.json").write_text(
        json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2), flush=True)


if __name__ == "__main__":
    main()
