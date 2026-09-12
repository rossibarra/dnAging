#!/usr/bin/env python3
"""Merge variable-Ne insertion results and analyse demographic effects."""
import argparse, csv, json
from pathlib import Path
import numpy as np

def exposure(epochs, age):
    vals, weights = [], []
    for e in epochs:
        left, right = float(e["time_left"]), float(e["time_right"])
        width = max(0.0, min(float(age), right) - left)
        if width:
            vals.append(float(e["effective_population_size"])); weights.append(width)
    vals, weights = np.asarray(vals), np.asarray(weights)
    logs = np.log(vals); mean = np.average(logs, weights=weights)
    return float(vals.min()), float(np.exp(mean)), float(np.sqrt(np.average((logs-mean)**2, weights=weights)))

def bins(x, y):
    edges = np.quantile(x, [0, .25, .5, .75, 1]); out = []
    for j in range(4):
        take = (x >= edges[j]) & ((x <= edges[j+1]) if j == 3 else (x < edges[j+1]))
        out.append({"lower": float(edges[j]), "upper": float(edges[j+1]),
                    "n": int(take.sum()), "mean_abs_error": float(y[take].mean())})
    return out

def main(a):
    rows = []
    for rep in range(1, a.n_replicates + 1):
        meta = json.loads((a.root/"simulations"/f"replicate_{rep:03d}"/"metadata.json").read_text())
        sizes = np.asarray([e["effective_population_size"] for e in meta["epochs"]])
        for j in range(1, a.n_ancient + 1):
            sample = f"ancient_{j:02d}"
            result = json.loads((a.root/"insertion_inference"/sample/f"replicate_{rep:03d}.json").read_text())
            mn, gm, sd = exposure(meta["epochs"], result["true_T"])
            error = float(result["map_T"] - result["true_T"])
            rows.append({"replicate": rep, "ancient_sample": sample,
                **{k: result[k] for k in ("true_T","map_T","mean_T","ci95_lower_T","ci95_upper_T","contains_true_95")},
                "error": error, "abs_error": abs(error), "min_ne_to_T": mn,
                "geomean_ne_to_T": gm, "sd_log_ne_to_T": sd,
                "global_min_ne": float(sizes.min()),
                "global_log_ne_range": float(np.log(sizes.max()/sizes.min())),
                "realised_epsilon": meta["realised_flip_rate"]})
    col = lambda key: np.asarray([r[key] for r in rows], float)
    err, ae, truth = col("error"), col("abs_error"), col("true_T")
    raw = np.column_stack([np.log(col("min_ne_to_T")), col("sd_log_ne_to_T"), truth])
    z = (raw-raw.mean(0))/raw.std(0); X = np.column_stack([np.ones(len(rows)), z])
    summary = {"n_samples": len(rows), "map_bias": float(err.mean()),
        "map_mae": float(ae.mean()), "map_rmse": float(np.sqrt(np.mean(err**2))),
        "ci95_coverage": float(col("contains_true_95").mean()),
        "standardized_predictors": ["log_min_ne_to_T","sd_log_ne_to_T","true_T"],
        "absolute_error_coefficients": np.linalg.lstsq(X, ae, rcond=None)[0][1:].tolist(),
        "signed_error_coefficients": np.linalg.lstsq(X, err, rcond=None)[0][1:].tolist(),
        "min_ne_quartiles": bins(col("min_ne_to_T"), ae),
        "ne_variation_quartiles": bins(col("sd_log_ne_to_T"), ae)}
    a.output_dir.mkdir(parents=True, exist_ok=True)
    (a.output_dir/"run.json").write_text(json.dumps(summary, indent=2)+"\n")
    with (a.output_dir/"age_demography_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    axes[0].scatter(truth, col("map_T"), s=8, alpha=.5); axes[0].plot([0,10000],[0,10000],"k--")
    axes[0].set(xlabel="True age", ylabel="MAP age")
    axes[1].scatter(col("min_ne_to_T"), ae, s=8, alpha=.5); axes[1].set(xscale="log", xlabel="Minimum Ne from 0 to T", ylabel="Absolute error")
    axes[2].scatter(col("sd_log_ne_to_T"), ae, s=8, alpha=.5); axes[2].set(xlabel="Time-weighted SD(log Ne) from 0 to T", ylabel="Absolute error")
    fig.tight_layout(); fig.savefig(a.output_dir/"accuracy_vs_demography.png", dpi=180)
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    p=argparse.ArgumentParser(); p.add_argument("--root",type=Path,required=True); p.add_argument("--output-dir",type=Path,required=True); p.add_argument("--n-replicates",type=int,default=100); p.add_argument("--n-ancient",type=int,default=10); main(p.parse_args())
