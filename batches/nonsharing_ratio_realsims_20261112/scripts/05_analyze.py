#!/usr/bin/env python3
"""Collect per-sample results, summary tables, sanity checks and the figure. Run via hpc_run from the repo root."""
import json, sys, glob
from pathlib import Path
import numpy as np
sys.path[:0] = ["pipeline"]
import nonsharing_ratio as nr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

B = Path("batches/nonsharing_ratio_realsims_20261112")
A0S = [10000, 20000, 50000]
REPS = [1, 2, 3, 4, 5]
BINS = [(0, 1500), (1500, 3000), (3000, 4500), (4500, 6000)]
DRAWS = [75, 93, 112, 131, 149]

def meanT(rep, a0):
    """across-draw mean/sd of T_hat per ancient sample, flags, plus selection stats from the draw NPZs."""
    Ts, fl, sel, usable, nm = [], [], [], [], None
    for d in DRAWS:
        p = B / f"draws/rep{rep}/draw{d}_a0{a0}.npz"
        if not p.exists(): continue
        z = np.load(p)
        ic = z["is_control"]; nm = [str(x) for x in z["targets"][~ic]]
        T, H, f = nr.ages_from_counts(z["num"], z["den"], ic, z["edges"], z["rates"])
        Ts.append(T); fl.append([str(x) for x in z["flags"]])
        m = json.loads(str(z["meta"])); sel.append(m["n_selected_sites"]); usable.append(m["n_usable_sites"])
    return nm, np.array(Ts), fl, sel, usable

def pseudo_controls(npz):
    """held-out controls as pseudo-targets, leave-one-out reference = other controls."""
    z = np.load(npz); ic = z["is_control"]
    num = z["num"].sum(0); den = z["den"].sum(0); f = num / den
    fc = f[ic]; out = []
    for j in range(fc.size):
        H = np.log(fc[j] / np.delete(fc, j).mean())
        out.append((H, float(nr.invert_hazard(np.array([H]), z["edges"], z["rates"])[0])))
    return out

rows = []; sel_rows = []
for rep in REPS:
    truth = [l.split("\t") for l in (B / f"sim/rep{rep}/truth.tsv").read_text().split("\n")[1:] if l]
    ta = {n: float(a) for _, n, a in truth}
    res = {}
    for a0 in A0S:
        nm, Ts, fl, sel, usable = meanT(rep, a0)
        if not len(Ts): continue
        res[a0] = (nm, Ts, fl)
        sel_rows.append((rep, a0, len(Ts), np.mean(sel), np.mean(usable)))
    for j, (_, n, a) in enumerate(truth):
        r = dict(rep=rep, sample=n, true_age=float(a))
        for a0 in A0S:
            if a0 not in res: continue
            nm, Ts, fl = res[a0]; k = nm.index(n)
            r[f"T_{a0}"] = Ts[:, k].mean(); r[f"sd_{a0}"] = Ts[:, k].std(ddof=1) if len(Ts) > 1 else np.nan
            r[f"flag_{a0}"] = ";".join(sorted({f[k] for f in fl if f[k]}))
        # true-ARG reference
        for a0 in A0S:
            p = B / f"draws/rep{rep}/true_a0{a0}.npz"
            if p.exists():
                z = np.load(p); ic = z["is_control"]; T, H, f = nr.ages_from_counts(z["num"], z["den"], ic, z["edges"], z["rates"])
                r[f"Ttrue_{a0}"] = T[[str(x) for x in z["targets"][~ic]].index(n)]
        rows.append(r)

cols = ["rep", "sample", "true_age"] + [f"{p}_{a}" for a in A0S for p in ("T", "sd", "flag", "Ttrue")]
with open(B / "per_sample_results.tsv", "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join("" if r.get(c) is None else (f"{r[c]:.1f}" if isinstance(r[c], float) else str(r[c])) for c in cols) + "\n")

def stats(x, y):
    ok = np.isfinite(y); x, y = x[ok], y[ok]; e = y - x
    sl = np.polyfit(x, y, 1)[0] if len(x) > 2 else np.nan
    return dict(n=len(x), mean_bias=e.mean(), median_bias=np.median(e), rmse=np.sqrt((e ** 2).mean()),
                r=np.corrcoef(x, y)[0, 1], slope=sl)

out = []
x = np.array([r["true_age"] for r in rows])
for label, key in (("ARG(SINGER+POLEGON)", "T"), ("true tree", "Ttrue")):
    for a0 in A0S:
        y = np.array([r.get(f"{key}_{a0}", np.nan) for r in rows], float)
        if not np.isfinite(y).any(): continue
        s = stats(x, y)
        bb = []
        for lo, hi in BINS:
            m = (x >= lo) & (x < hi) & np.isfinite(y); bb.append((y[m] - x[m]).mean() if m.any() else np.nan)
        nflag = sum(1 for r in rows if r.get(f"flag_{a0}")) if key == "T" else ""
        out.append((label, a0, s, bb, nflag))
with open(B / "summary_by_a0.tsv", "w") as fh:
    fh.write("source\ta0\tn\tmean_bias\tmedian_bias\trmse\tpearson_r\tslope\tbias_0-1500\tbias_1500-3000\tbias_3000-4500\tbias_4500-6000\tn_targets_flagged_any_draw\n")
    for label, a0, s, bb, nf in out:
        fh.write(f"{label}\t{a0}\t{s['n']}\t{s['mean_bias']:.0f}\t{s['median_bias']:.0f}\t{s['rmse']:.0f}\t{s['r']:.3f}\t{s['slope']:.3f}\t"
                 + "\t".join(f"{b:.0f}" for b in bb) + f"\t{nf}\n")
with open(B / "site_selection.tsv", "w") as fh:
    fh.write("rep\ta0\tn_draws\tmean_selected_sites\tmean_usable_sites\tfraction\n")
    for rep, a0, nd, s, u in sel_rows: fh.write(f"{rep}\t{a0}\t{nd}\t{s:.0f}\t{u:.0f}\t{s/u:.4f}\n")

# controls as pseudo-targets
with open(B / "controls_pseudotargets.tsv", "w") as fh:
    fh.write("source\trep\ta0\tdraw\tcontrol\tH\tT_hat\n")
    for rep in REPS:
        for a0 in A0S:
            for d in DRAWS + ["true"]:
                p = B / (f"draws/rep{rep}/true_a0{a0}.npz" if d == "true" else f"draws/rep{rep}/draw{d}_a0{a0}.npz")
                if not p.exists(): continue
                for j, (H, T) in enumerate(pseudo_controls(p)):
                    fh.write(f"{'true' if d=='true' else 'arg'}\t{rep}\t{a0}\t{d}\t{j+1}\t{H:.5f}\t{T:.1f}\n")

# figure
cmap = plt.get_cmap("tab10"); fig, axs = plt.subplots(1, 3, figsize=(15, 5), sharex=True, sharey=True)
for ax, a0 in zip(axs, A0S):
    for rep in REPS:
        rr = [r for r in rows if r["rep"] == rep and f"T_{a0}" in r]
        if not rr: continue
        ax.errorbar([r["true_age"] for r in rr], [r[f"T_{a0}"] for r in rr], yerr=[r[f"sd_{a0}"] for r in rr],
                    fmt="o", ms=5, color=cmap(rep - 1), alpha=0.85, capsize=2, lw=0.8, label=f"rep {rep}")
    ax.plot([0, 6500], [0, 6500], "k--", lw=1)
    ax.set_title(f"A0 = {a0:,} generations (age-mode child)"); ax.set_xlabel("true age (generations)")
    ax.set_xlim(-200, 6500)
axs[0].set_ylabel("inferred age T_hat (generations)\n(mean over draws, bars = across-draw SD)"); axs[0].legend(title="replicate", loc="upper left")
fig.tight_layout(); fig.savefig(B / "inferred_vs_true.png", dpi=150); fig.savefig(B / "inferred_vs_true.svg")
print(open(B / "summary_by_a0.tsv").read()); print(open(B / "site_selection.tsv").read())
