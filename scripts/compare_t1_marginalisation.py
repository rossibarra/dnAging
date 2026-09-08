#!/usr/bin/env python3
"""Compare uniform and denominator-weighted branch-age marginalisation.

This harness operates directly on the true tree sequences from the matched Ne
simulation sweep.  It groups mutations sharing the same descendant count and
branch-time interval, then evaluates exactly the production ``phi_lookup`` for
both T1 alternatives on the same sites.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import tskit

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from posterior_sample_age_infer import load_table, phi_lookup, summarize


def ancient_alt_ids(path: Path) -> tuple[set[int], int]:
    """Return site IDs carrying ALT; legacy missing calls mean ancestral here."""
    alt = set()
    rows = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            site_id = int(fields[2])
            gt_keys = fields[8].split(":")
            gt = fields[9].split(":")[gt_keys.index("GT")]
            if gt == "1":
                alt.add(site_id)
            elif gt not in {"0", "."}:
                raise ValueError(f"unexpected ancient haploid GT {gt!r} at {path}:{rows + 1}")
            rows += 1
    return alt, rows


def site_groups(ts: tskit.TreeSequence, alt_ids: set[int], cutoff: float):
    groups = defaultdict(lambda: [0, 0])
    skipped = defaultdict(int)
    for site in ts.sites():
        if len(site.mutations) != 1:
            skipped["multiple_mutations"] += 1
            continue
        mut = site.mutations[0]
        tree = ts.at(site.position)
        parent = tree.parent(mut.node)
        if parent == tskit.NULL:
            skipped["root_mutation"] += 1
            continue
        d0 = tree.num_samples(mut.node)
        if not 1 <= d0 < ts.num_samples:
            skipped["monomorphic"] += 1
            continue
        lo = float(ts.node(mut.node).time)
        hi = min(float(ts.node(parent).time), cutoff)
        if lo >= cutoff or hi <= lo:
            skipped["older_than_cutoff"] += 1
            continue
        counts = groups[(int(d0), lo, hi)]
        counts[1 if site.id in alt_ids else 0] += 1
    return groups, skipped


def infer_one(sim_dir: Path, table_path: Path, ne: int, truth: float):
    stem = sim_dir.name
    ts = tskit.load(sim_dir / f"{stem}.trees")
    alt_ids, n_vcf = ancient_alt_ids(sim_dir / f"{stem}_ancient.vcf")
    if n_vcf != ts.num_sites:
        raise ValueError(f"{stem}: VCF has {n_vcf} rows but tree sequence has {ts.num_sites} sites")
    tab = load_table(table_path)
    cutoff = 6.0 * ne
    groups, skipped = site_groups(ts, alt_ids, cutoff)
    rows = []
    for method in ("uniform", "weighted"):
        ll = np.zeros_like(tab["Tgrid"], dtype=float)
        for (d0, lo, hi), (n_ref, n_alt) in groups.items():
            p = phi_lookup(tab, d0, lo, hi, n_called=ts.num_samples,
                           marginalise=method)
            if p is None or not np.all(np.isfinite(p)):
                raise ValueError(f"{stem}: invalid {method} lookup for {(d0, lo, hi)}")
            p = np.clip(p, 1e-300, 1.0 - 1e-15)
            ll += n_alt * np.log(p) + n_ref * np.log1p(-p)
        summary, _ = summarize(tab["Tgrid"], ll)
        rows.append({"ne": ne, "simulation": stem, "true_T": truth,
                     "method": method, "sites_used": sum(sum(x) for x in groups.values()),
                     "unique_intervals": len(groups), **summary})
    return rows, dict(skipped)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root-template", default="/tmp/dnAging-ne{ne}-infinite-10mb")
    parser.add_argument("--table-template", default="/tmp/dnAging-ne{ne}-t1-format5-table.npz")
    parser.add_argument("--ne", nargs="+", type=int, default=[10000, 50000, 100000, 200000])
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    jobs = []
    for ne in args.ne:
        root = Path(args.root_template.format(ne=ne))
        table = Path(args.table_template.format(ne=ne))
        manifest = json.loads((root / "run.json").read_text())
        truth = {str(row["simulation"]): float(row["true_age"])
                 for row in manifest["replicates"]}
        for sim_dir in sorted(root.glob("simulation_*")):
            keys = (sim_dir.name, sim_dir.name.rsplit("_", 1)[-1], str(int(sim_dir.name.rsplit("_", 1)[-1])))
            true_t = next((truth[k] for k in keys if k in truth), None)
            if true_t is None:
                raise KeyError(f"no truth value for {sim_dir.name}")
            jobs.append((sim_dir, table, ne, true_t))
    out = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(infer_one, *job): job for job in jobs}
        for future in as_completed(futures):
            sim_dir, _table, ne, _truth = futures[future]
            rows, skipped = future.result()
            out.extend(rows)
            print(ne, sim_dir.name,
                  {r["method"]: round(r["mean_T"], 1) for r in rows}, skipped,
                  flush=True)
    out.sort(key=lambda row: (row["ne"], row["simulation"], row["method"]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(out[0])
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(out)


if __name__ == "__main__":
    main()
