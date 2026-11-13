#!/usr/bin/env python3
"""Collect per-sample simulation results (constant size, bottleneck) into one small CSV.

Reads the two batches' per_sample_results.tsv (gitignored by the repo-wide *.tsv rule) and writes
results/nonsharing_ratio/sim_inferred_vs_true.csv (tracked) so the README figure can be rebuilt from a clean checkout.
Columns: demography, rep, sample, true_age, then for each A0: T_arg (mean over 5 POLEGON draws, c(t) from the ARG)
and T_truetree (true tree sequence, c(t) from that tree).  Simulated data only; no real samples.
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BATCHES = {
    "constant": ROOT / "batches/nonsharing_ratio_realsims_20261112/per_sample_results.tsv",
    "bottleneck": ROOT / "batches/nonsharing_ratio_realsims_bottleneck_20261112/per_sample_results.tsv",
}
A0S = ("10000", "20000", "50000")
out = ROOT / "results/nonsharing_ratio/sim_inferred_vs_true.csv"
with open(out, "w", newline="") as fo:
    w = csv.writer(fo)
    w.writerow(["demography", "rep", "sample", "true_age"] + [f"{k}_{a}" for a in A0S for k in ("T_arg", "T_truetree")])
    n = 0
    for dem, p in BATCHES.items():
        for r in csv.DictReader(open(p), delimiter="\t"):
            w.writerow([dem, r["rep"], r["sample"], r["true_age"]]
                       + [r[f"{c}_{a}"] for a in A0S for c in ("T", "Ttrue")])
            n += 1
print("wrote", out, n, "rows")
