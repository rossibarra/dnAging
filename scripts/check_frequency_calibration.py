#!/usr/bin/env python3
"""Check ancient-carriage calibration, binning by predicted probability."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import tskit

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from posterior_sample_age_infer import load_table, phi_lookup


def alt_ids(path):
    carried = set()
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            gt = fields[9].split(":")[fields[8].split(":").index("GT")]
            if gt == "1":
                carried.add(int(fields[2]))
            elif gt not in {"0", "."}:
                raise ValueError(f"unexpected haploid GT {gt!r}")
    return carried


def narrow_table(tab, true_t):
    grid = tab["Tgrid"]
    right = int(np.searchsorted(grid, true_t, side="right"))
    ids = np.unique(np.clip([right - 1, right], 0, len(grid) - 1))
    out = dict(tab)
    out["Tgrid"] = grid[ids]
    out["table"] = tab["table"][..., ids]
    if "table2" in tab:
        out["table2"] = tab["table2"][..., ids]
    return out


def one_sim(sim_dir, table_path, ne, true_t):
    stem = sim_dir.name
    ts = tskit.load(sim_dir / f"{stem}.trees")
    carried = alt_ids(sim_dir / f"{stem}_ancient.vcf")
    tab = narrow_table(load_table(table_path), true_t)
    cutoff = 6.0 * ne
    rows = []
    for site in ts.sites():
        if len(site.mutations) != 1:
            continue
        mut = site.mutations[0]
        tree = ts.at(site.position)
        parent = tree.parent(mut.node)
        if parent == tskit.NULL:
            continue
        d0 = int(tree.num_samples(mut.node))
        if not 1 <= d0 < ts.num_samples:
            continue
        lo = float(ts.node(mut.node).time)
        hi = min(float(ts.node(parent).time), cutoff)
        if lo >= cutoff or hi <= lo:
            continue
        phi = phi_lookup(tab, d0, lo, hi, n_called=ts.num_samples,
                         marginalise="uniform")
        if phi is None or not np.all(np.isfinite(phi)):
            continue
        pred = float(np.interp(true_t, tab["Tgrid"], phi))
        rows.append((ne, stem, d0, pred, int(site.id in carried)))
    return rows


def binned(rows, bins):
    output = []
    categories = (("d0=1", lambda d: d == 1),
                  ("d0=2-3", lambda d: 2 <= d <= 3),
                  ("d0>=4", lambda d: d >= 4),
                  ("all", lambda d: True))
    for ne in sorted({r[0] for r in rows}):
        ne_rows = [r for r in rows if r[0] == ne]
        for label, keep in categories:
            selected = [r for r in ne_rows if keep(r[2])]
            if not selected:
                continue
            selected.sort(key=lambda r: r[3])
            for ib, idx in enumerate(np.array_split(np.arange(len(selected)), bins), 1):
                block = [selected[i] for i in idx]
                pred = np.array([r[3] for r in block])
                obs = np.array([r[4] for r in block])
                output.append({"ne": ne, "d0_group": label, "bin": ib,
                               "n": len(block), "predicted": pred.mean(),
                               "observed": obs.mean(),
                               "observed_se": np.sqrt(obs.mean() * (1 - obs.mean()) / len(obs)),
                               "difference": obs.mean() - pred.mean()})
    return output


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ne", nargs="+", type=int, default=[10000, 50000, 100000, 200000])
    p.add_argument("--root-template", default="/tmp/dnAging-ne{ne}-infinite-10mb")
    p.add_argument("--table-template", default="/tmp/dnAging-ne{ne}-t1-format5-table.npz")
    p.add_argument("--bins", type=int, default=10)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args()
    jobs = []
    for ne in args.ne:
        root = Path(args.root_template.format(ne=ne))
        manifest = json.loads((root / "run.json").read_text())
        truth = {f"simulation_{r['simulation']:02d}": r["true_age"]
                 for r in manifest["replicates"]}
        for sim_dir in sorted(root.glob("simulation_*")):
            jobs.append((sim_dir, Path(args.table_template.format(ne=ne)),
                         ne, truth[sim_dir.name]))
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(one_sim, *job) for job in jobs]
        for future in as_completed(futures):
            block = future.result()
            rows.extend(block)
            print(f"completed {block[0][0]} {block[0][1]}: {len(block)} sites", flush=True)
    summary = binned(rows, args.bins)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(summary)


if __name__ == "__main__":
    main()
