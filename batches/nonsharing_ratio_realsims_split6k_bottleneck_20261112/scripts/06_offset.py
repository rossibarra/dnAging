#!/usr/bin/env python3
"""Offset ratio (median ancient f / mean control f) per (rep, A0), from the ARG draws and the true tree, and the
controls-as-pseudo-targets summary.  Run via hpc_run from the repo root AFTER 05_analyze.py (controls_pseudotargets.tsv)."""
from pathlib import Path
import numpy as np
B = Path("batches/nonsharing_ratio_realsims_split6k_bottleneck_20261112")
DRAWS = [75, 93, 112, 131, 149]; A0S = [10000, 20000, 50000]; REPS = [1, 2, 3, 4, 5]
def fs(p):
    z = np.load(p); ic = z["is_control"]; f = z["num"].sum(0) / z["den"].sum(0)
    return np.median(f[~ic]), f[ic].mean(), np.mean(np.log(f[~ic] / f[ic].mean()))
with open(B / "offset_by_rep.tsv", "w") as fh, open(B / "offset_summary.tsv", "w") as fs_:
    fh.write("rep\ta0\tmedian_ancient_f\tmean_control_f\tratio_arg_meanOverDraws\tmean_H_arg\tratio_truetree\tmean_H_truetree\n")
    fs_.write("a0\tmean_ratio_arg\tmin_ratio_arg\tmax_ratio_arg\tmean_ratio_truetree\tmean_median_ancient_f\tmean_control_f\n")
    for a0 in A0S:
        o = []
        for rep in REPS:
            fr = [fs(B / f"draws/rep{rep}/draw{d}_a0{a0}.npz") for d in DRAWS]
            ft = fs(B / f"draws/rep{rep}/true_a0{a0}.npz")
            row = (np.mean([x[0] for x in fr]), np.mean([x[1] for x in fr]), np.mean([x[0] / x[1] for x in fr]), np.mean([x[2] for x in fr]), ft[0] / ft[1], ft[2])
            o.append(row); fh.write(f"{rep}\t{a0}\t" + "\t".join(f"{v:.5f}" for v in row) + "\n")
        o = np.array(o)
        fs_.write(f"{a0}\t{o[:,2].mean():.4f}\t{o[:,2].min():.4f}\t{o[:,2].max():.4f}\t{o[:,4].mean():.4f}\t{o[:,0].mean():.5f}\t{o[:,1].mean():.5f}\n")
print(open(B / "offset_summary.tsv").read())
# pseudo-target summary (leave-one-out controls vs other controls)
import collections
d = collections.defaultdict(list)
for l in (B / "controls_pseudotargets.tsv").read_text().split("\n")[1:]:
    if l:
        src, rep, a0, dr, c, H, T = l.split("\t"); d[(src, int(a0))].append((float(H), float(T)))
print("source a0 n mean_H mean_T median_T")
for (src, a0), v in sorted(d.items()):
    v = np.array(v); print(src, a0, len(v), f"{v[:,0].mean():.4f}", f"{v[:,1].mean():.0f}", f"{np.median(v[:,1]):.0f}")
