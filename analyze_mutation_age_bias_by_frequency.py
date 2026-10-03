#!/usr/bin/env python3
"""Compare recorded ARG mutation ages with truth by modern allele count."""
import csv
import gzip
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import tszip


ROOT = Path("msprime_variable_ne_error")
BINS = [(1, 1, "1"), (2, 2, "2"), (3, 4, "3–4"), (5, 8, "5–8"),
        (9, 13, "9–13"), (14, 18, "14–18"), (19, 22, "19–22"),
        (23, 25, "23–25")]


def load_truth(label):
    ans = {}
    with gzip.open(ROOT / "simulations" / label / "mutation_truth.tsv.gz", "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            ans[row["site_id"]] = (float(row["mutation_time"]), int(row["modern_derived_count"]))
    return ans


def main():
    out = ROOT / "mutation_age_frequency"
    out.mkdir(parents=True, exist_ok=False)
    records = []
    for rep in range(1, 11):
        label = f"replicate_{rep:03d}"
        truth = load_truth(label)
        site_by_singer = {}
        with (ROOT / "singer" / label / "position_map.tsv").open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                site_by_singer[int(row["singer_pos"])] = row["site_id"]
        eligible = set(np.load(ROOT / "insertion_posterior" / "prepared_frequency" /
                               label / "calls" / "chrchr1.npz")["position"].tolist())

        # The known ARG is an ordinary tskit .trees file.
        import tskit
        true_ts = tskit.load(ROOT / "simulations" / label / "known_modern_arg.trees")
        true_age = {str(site.id): float(site.mutations[0].time) for site in true_ts.sites()
                    if len(site.mutations) == 1}

        by_site = {sid: [] for pos, sid in site_by_singer.items() if pos in eligible}
        for draw in range(50, 100):
            ts = tszip.decompress(ROOT / "argtest" / label / "out" / "combined" /
                                  f"run.combined.{draw}.tsz")
            for site in ts.sites():
                sid = site_by_singer.get(int(site.position))
                if sid in by_site and len(site.mutations) == 1 and np.isfinite(site.mutations[0].time):
                    by_site[sid].append(float(site.mutations[0].time))

        for sid, ages in by_site.items():
            if not ages or sid not in truth or sid not in true_age:
                continue
            known, count = truth[sid]
            records.append((rep, sid, count, known, true_age[sid], float(np.mean(ages)),
                            float(np.median(ages)), len(ages)))

    header = ["replicate", "site_id", "modern_derived_count", "known_mutation_age",
              "true_arg_mutation_age", "estimated_arg_mean_mutation_age",
              "estimated_arg_median_mutation_age", "n_arg_draws"]
    with (out / "per_site.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header); writer.writerows(records)

    rows = []
    for lo, hi, label in BINS:
        subset = [r for r in records if lo <= r[2] <= hi]
        for name, col in [("true_ARG", 4), ("estimated_ARG", 5)]:
            err = np.array([r[col] - r[3] for r in subset])
            ratio = np.array([r[col] / r[3] for r in subset])
            rows.append({"frequency_bin": label, "arg": name, "sites": len(subset),
                         "mean_bias": float(np.mean(err)), "median_bias": float(np.median(err)),
                         "mae": float(np.mean(np.abs(err))), "rmse": float(np.sqrt(np.mean(err**2))),
                         "median_log2_ratio": float(np.median(np.log2(ratio))),
                         "mean_true_age": float(np.mean([r[3] for r in subset]))})
    with (out / "summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)

    x = np.arange(len(BINS)); fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))
    for name, color, marker in [("true_ARG", "#377eb8", "o"), ("estimated_ARG", "#e41a1c", "s")]:
        rr = [r for r in rows if r["arg"] == name]
        axes[0].plot(x, [r["mean_bias"] for r in rr], marker=marker, color=color, label=name.replace("_", " "))
        axes[1].plot(x, [r["mae"] for r in rr], marker=marker, color=color)
        axes[2].plot(x, [r["median_log2_ratio"] for r in rr], marker=marker, color=color)
    for ax in axes:
        ax.axhline(0, color="0.5", lw=0.8); ax.set_xticks(x, [b[2] for b in BINS], rotation=35)
        ax.set_xlabel("Modern derived-allele count (n=26)")
    axes[0].set_ylabel("Mean mutation-age bias (generations)"); axes[0].legend(frameon=False)
    axes[1].set_ylabel("Mean absolute error (generations)")
    axes[2].set_ylabel("Median log2(recorded age / true age)")
    fig.tight_layout(); fig.savefig(out / "mutation_age_bias_by_frequency.png", dpi=220)
    (out / "run.json").write_text(json.dumps({"replicates": 10, "estimated_arg_draws": list(range(50, 100)),
        "sites": len(records), "frequency_definition": "observed modern derived/ALT count",
        "estimated_arg_site_summary": "mean explicit mutation time across available posterior ARG draws"}, indent=2) + "\n")


if __name__ == "__main__":
    main()
