#!/usr/bin/env python3
"""Quantify polarity switches among high-frequency mutation-age errors."""
import csv
import gzip
import json
from pathlib import Path

import numpy as np
import tszip


R = Path("msprime_variable_ne_error")


def main():
    rows = []
    for rep in range(1, 11):
        label = f"replicate_{rep:03d}"
        truth = {}
        with gzip.open(R / "simulations" / label / "mutation_truth.tsv.gz", "rt") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                d = int(row["modern_derived_count"])
                if 23 <= d <= 25:
                    truth[row["site_id"]] = (float(row["mutation_time"]), d)
        sid_by_pos = {}
        with (R / "singer" / label / "position_map.tsv").open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["site_id"] in truth:
                    sid_by_pos[int(row["singer_pos"])] = row["site_id"]
        eligible = set(np.load(R / "insertion_posterior" / "prepared_frequency" / label /
                               "calls" / "chrchr1.npz")["position"].tolist())
        for draw in range(50, 100):
            ts = tszip.decompress(R / "argtest" / label / "out" / "combined" /
                                  f"run.combined.{draw}.tsz")
            for site in ts.sites():
                pos = int(site.position)
                sid = sid_by_pos.get(pos)
                if pos not in eligible or sid is None or len(site.mutations) != 1:
                    continue
                mut = site.mutations[0]
                if not np.isfinite(mut.time):
                    continue
                tree = ts.at(site.position)
                inferred_d = tree.num_samples(mut.node)
                true_age, d = truth[sid]
                # Exact count labels avoid treating unrelated count errors as switches.
                if inferred_d == d:
                    state = "correct_polarity_and_count"
                elif inferred_d == 26 - d:
                    state = "complementary_count_polarity_switch"
                else:
                    state = "other_count_mismatch"
                rows.append((rep, draw, sid, d, inferred_d, state,
                             float(mut.time) - true_age, float(mut.time) / true_age))

    out = R / "mutation_age_frequency" / "polarity_23_25"
    out.mkdir(parents=True, exist_ok=False)
    fields = ["replicate", "draw", "site_id", "true_derived_count", "inferred_descendant_count",
              "classification", "age_error", "age_ratio"]
    with (out / "draw_site.tsv").open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t"); w.writerow(fields); w.writerows(rows)
    result = {"draw_site_assignments": len(rows), "classes": {}}
    errors = np.array([r[6] for r in rows])
    negative_total = -errors[errors < 0].sum()
    for state in ("correct_polarity_and_count", "complementary_count_polarity_switch",
                  "other_count_mismatch"):
        rr = [r for r in rows if r[5] == state]
        e = np.array([r[6] for r in rr])
        result["classes"][state] = {
            "n": len(rr), "fraction": len(rr) / len(rows),
            "mean_error": float(e.mean()), "median_error": float(np.median(e)),
            "median_ratio": float(np.median([r[7] for r in rr])),
            "fraction_underestimated": float(np.mean(e < 0)),
            "share_of_total_negative_error_magnitude": float(-e[e < 0].sum() / negative_total),
        }
    # Also describe the most severe underestimates.
    for cutoff in (10000, 25000, 50000):
        severe = [r for r in rows if r[6] <= -cutoff]
        result[f"error_at_most_minus_{cutoff}"] = {
            "n": len(severe),
            "fraction_switched": float(np.mean([r[5] == "complementary_count_polarity_switch"
                                                 for r in severe])) if severe else None,
            "fraction_other_mismatch": float(np.mean([r[5] == "other_count_mismatch"
                                                       for r in severe])) if severe else None,
        }
    (out / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
