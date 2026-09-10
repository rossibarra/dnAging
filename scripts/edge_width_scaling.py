#!/usr/bin/env python3
"""Does the point-to-edge bias scale with the width of the edge?

The paired msprime run showed that replacing an exact mutation time by its true
containing edge moves MAP bias from +10 to +364 generations at Ne=50K.  Two very
different mechanisms produce that:

  (i) the edge is simply WIDE.  The uniform-on-edge prior spreads the mutation
      age over a range far larger than the age resolution being sought, and the
      damage would then vanish as edges narrow; or
  (ii) the edge REPRESENTATION discards information -- eq. (4) keeps only the
      descendant count and the endpoints -- in which case narrowing the edge
      need not help, because what is lost is not the width.

This stratifies by edge width to separate them.  The obvious version of that
test is confounded: edge width is tightly coupled to allele-frequency class
(at Ne=50K, median width runs 17K generations for singletons to 137K for
d0 >= 13), so a bare width split also splits on d0 and would attribute a
frequency-class effect to width.

The design here removes that.  Within each width stratum BOTH arms are
evaluated -- exact point age and edge interval -- on the *identical* site set,
and the reported quantity is the paired difference.  Any effect of d0, mutation
age or site count is common to the two arms of a stratum and cancels in the
difference.  It also fixes a smaller confound in the original paired run: exact
mode drops mutations younger than the table's youngest age while edge mode keeps
them, so the two arms there used 7,354,754 and 7,631,501 sites rather than one
shared set.  Here a site enters only if it passes both filters.

Usage mirrors msprime_exact_time_validation.py: `infer` one replicate, then
`merge`.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import tskit

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from posterior_sample_age_infer import load_table, phi_lookup, summarize

# Edge width in generations.  Chosen from the observed distribution at Ne=50K
# (5th percentile 4.1K, median 50.3K, 95th 258K) so that every stratum is
# populated and the narrowest is comparable to the ancient ages being estimated.
WIDTH_EDGES = (0.0, 10_000.0, 30_000.0, 100_000.0, float("inf"))


def stratum_labels():
    out = []
    for lo, hi in zip(WIDTH_EDGES[:-1], WIDTH_EDGES[1:]):
        out.append(f"{lo/1000:.0f}k-{'inf' if hi == float('inf') else f'{hi/1000:.0f}k'}")
    return out


def vcf_calls(path: Path, n_modern: int):
    """Reuse the validation script's reader so the call parsing is identical."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from msprime_exact_time_validation import vcf_calls as reader
    return reader(path, n_modern)


def infer(args):
    directory = args.input_root / f"replicate_{args.replicate:03d}"
    metadata = json.loads((directory / "metadata.json").read_text())
    ts = tskit.load(directory / "known_modern_arg.trees")
    if ts.num_samples != args.n_modern:
        raise ValueError(f"modern ARG has {ts.num_samples} samples")
    calls = vcf_calls(directory / "all_samples.vcf.gz", args.n_modern)
    tab = load_table(args.frequency_table)
    grid = np.asarray(tab["Tgrid"], dtype=np.float64)
    cutoff = args.mutation_age_max_tau * 2.0 * float(metadata["Ne"])
    age_floor = float(tab["age"][0])

    n_strata = len(WIDTH_EDGES) - 1
    ll_exact = np.zeros((n_strata, len(grid)))
    ll_edge = np.zeros((n_strata, len(grid)))
    tally = [{"sites": 0, "sum_width": 0.0, "sum_d0": 0.0, "sum_age": 0.0,
              "carried": 0} for _ in range(n_strata)]
    dropped = {"multi_mutation": 0, "count_mismatch": 0, "root": 0,
               "exact_filtered": 0, "edge_filtered": 0, "nonfinite": 0}

    for tree in ts.trees():
        for site in tree.sites():
            if len(site.mutations) != 1:
                dropped["multi_mutation"] += 1
                continue
            mutation = site.mutations[0]
            if tskit.is_unknown_time(mutation.time):
                raise ValueError(f"unknown time for mutation {mutation.id}")
            d0, ancient_genotype = calls[site.id]
            if d0 != int(tree.num_samples(mutation.node)):
                dropped["count_mismatch"] += 1
                continue
            parent = tree.parent(mutation.node)
            if parent == tskit.NULL:
                dropped["root"] += 1
                continue
            t_exact = float(mutation.time)
            t_lo = float(ts.node(mutation.node).time)
            t_hi = min(float(ts.node(parent).time), cutoff)
            # A site must be admissible to BOTH arms or it enters neither.
            if t_exact >= cutoff or t_exact < age_floor:
                dropped["exact_filtered"] += 1
                continue
            if t_lo >= cutoff or t_hi <= t_lo:
                dropped["edge_filtered"] += 1
                continue

            phi_e = phi_lookup(tab, d0, t_exact, t_exact, n_called=args.n_modern)
            phi_b = phi_lookup(tab, d0, t_lo, t_hi, n_called=args.n_modern)
            if (phi_e is None or phi_b is None
                    or not np.all(np.isfinite(phi_e))
                    or not np.all(np.isfinite(phi_b))):
                dropped["nonfinite"] += 1
                continue

            k = int(np.searchsorted(WIDTH_EDGES, t_hi - t_lo, side="right") - 1)
            k = min(max(k, 0), n_strata - 1)
            for phi, acc in ((phi_e, ll_exact), (phi_b, ll_edge)):
                q = np.clip((1.0 - args.epsilon) * phi
                            + args.epsilon * (1.0 - phi), 1e-300, 1.0)
                acc[k] += (np.log(q) if ancient_genotype
                           else np.log(np.clip(1.0 - q, 1e-300, 1.0)))
            tally[k]["sites"] += 1
            tally[k]["sum_width"] += t_hi - t_lo
            tally[k]["sum_d0"] += d0
            tally[k]["sum_age"] += t_exact
            tally[k]["carried"] += int(bool(ancient_genotype))

    true_T = float(metadata["true_ancient_age"])
    rows = []
    for k, label in enumerate(stratum_labels()):
        n = tally[k]["sites"]
        row = {"replicate": args.replicate, "true_T": true_T,
               "stratum": label, "sites": n}
        if n:
            row.update(mean_width=tally[k]["sum_width"] / n,
                       mean_d0=tally[k]["sum_d0"] / n,
                       mean_mutation_age=tally[k]["sum_age"] / n,
                       carried_fraction=tally[k]["carried"] / n)
            for arm, acc in (("exact", ll_exact), ("edge", ll_edge)):
                summary, _density = summarize(grid, acc[k])
                row[f"map_T_{arm}"] = summary["map_T"]
                row[f"mean_T_{arm}"] = summary["mean_T"]
                row[f"ci95_width_{arm}"] = (summary["ci95_upper_T"]
                                            - summary["ci95_lower_T"])
        rows.append(row)

    # Pooled over strata: the whole-replicate paired contrast on shared sites.
    # This is the original paired test with its site-set confound removed, so it
    # is the number the stratum rows should be read against.
    total = sum(t["sites"] for t in tally)
    pooled = {"replicate": args.replicate, "true_T": true_T, "stratum": "all",
              "sites": total}
    if total:
        pooled.update(
            mean_width=sum(t["sum_width"] for t in tally) / total,
            mean_d0=sum(t["sum_d0"] for t in tally) / total,
            mean_mutation_age=sum(t["sum_age"] for t in tally) / total,
            carried_fraction=sum(t["carried"] for t in tally) / total)
        for arm, acc in (("exact", ll_exact), ("edge", ll_edge)):
            summary, _density = summarize(grid, acc.sum(axis=0))
            pooled[f"map_T_{arm}"] = summary["map_T"]
            pooled[f"mean_T_{arm}"] = summary["mean_T"]
            pooled[f"ci95_width_{arm}"] = (summary["ci95_upper_T"]
                                           - summary["ci95_lower_T"])
    rows.append(pooled)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / f"replicate_{args.replicate:03d}.json").write_text(
        json.dumps({"rows": rows, "dropped": dropped,
                    "width_edges": list(WIDTH_EDGES[:-1]) + ["inf"],
                    "Ne": metadata["Ne"], "epsilon": args.epsilon}, indent=2) + "\n")
    print(json.dumps({"replicate": args.replicate, "dropped": dropped,
                      "rows": rows}, indent=2), flush=True)


def merge(args):
    rows = []
    for path in sorted(args.input_dir.glob("replicate_*.json")):
        rows.extend(json.loads(path.read_text())["rows"])
    if not rows:
        raise SystemExit(f"no replicate json under {args.input_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    cols = ["replicate", "true_T", "stratum", "sites", "mean_width", "mean_d0",
            "mean_mutation_age", "carried_fraction", "map_T_exact", "map_T_edge",
            "mean_T_exact", "mean_T_edge", "ci95_width_exact", "ci95_width_edge"]
    with (args.output_dir / "per_replicate_strata.tsv").open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    summary = []
    for label in stratum_labels() + ["all"]:
        sel = [r for r in rows if r["stratum"] == label and "map_T_edge" in r]
        if not sel:
            continue
        true = np.array([r["true_T"] for r in sel])
        exact = np.array([r["map_T_exact"] for r in sel])
        edge = np.array([r["map_T_edge"] for r in sel])
        paired = edge - exact
        n = len(sel)
        summary.append({
            "stratum": label, "replicates": n,
            "mean_sites_per_replicate": float(np.mean([r["sites"] for r in sel])),
            "mean_width": float(np.mean([r.get("mean_width", np.nan) for r in sel])),
            "mean_d0": float(np.mean([r.get("mean_d0", np.nan) for r in sel])),
            "bias_exact": float(np.mean(exact - true)),
            "bias_edge": float(np.mean(edge - true)),
            "paired_shift_edge_minus_exact": float(np.mean(paired)),
            "paired_shift_se": float(np.std(paired, ddof=1) / np.sqrt(n)),
            "rmse_exact": float(np.sqrt(np.mean((exact - true) ** 2))),
            "rmse_edge": float(np.sqrt(np.mean((edge - true) ** 2))),
        })
    with (args.output_dir / "width_scaling_summary.tsv").open("w") as fh:
        cols = list(summary[0])
        fh.write("\t".join(cols) + "\n")
        for r in summary:
            fh.write("\t".join(f"{r[c]:.4g}" if isinstance(r[c], float) else str(r[c])
                               for c in cols) + "\n")
    print(json.dumps(summary, indent=2))
    print(f"\nwrote {args.output_dir}/width_scaling_summary.tsv and "
          f"per_replicate_strata.tsv")


def build_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)
    i = sub.add_parser("infer")
    i.add_argument("--input-root", type=Path, required=True)
    i.add_argument("--frequency-table", type=Path, required=True)
    i.add_argument("--output-dir", type=Path, required=True)
    i.add_argument("--replicate", type=int, required=True)
    i.add_argument("--n-modern", type=int, default=26)
    i.add_argument("--epsilon", type=float, default=0.0)
    i.add_argument("--mutation-age-max-tau", type=float, default=3.0)
    i.set_defaults(func=infer)
    m = sub.add_parser("merge")
    m.add_argument("--input-dir", type=Path, required=True)
    m.add_argument("--output-dir", type=Path, required=True)
    m.set_defaults(func=merge)
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    args.func(args)
