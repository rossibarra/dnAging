#!/usr/bin/env python3
"""Collect per-sample results, summary tables, sanity checks and the figure. Run via hpc_run from the repo root."""
import json, sys, glob
from pathlib import Path
import numpy as np
sys.path[:0] = ["pipeline"]
import nonsharing_ratio as nr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ARMS = {"A": "", "B": "B"}   # arm -> file prefix infix (A: controls from A; B: controls from B)
B = Path("batches/nonsharing_ratio_realsims_split6k_Bcontrols_20261113")
A0S = [10000, 20000, 50000]
REPS = [1, 2, 3, 4, 5]
BINS = [(0, 1500), (1500, 3000), (3000, 4500), (4500, 6000)]
DRAWS = [75, 93, 112, 131, 149]

def meanT(rep, a0, tag="", X=""):
    """across-draw mean/sd of T_hat per ancient sample, flags, plus selection stats from the draw NPZs."""
    Ts, fl, sel, usable, nm = [], [], [], [], None
    for d in DRAWS:
        p = B / f"draws/rep{rep}/draw{X}{d}{tag}_a0{a0}.npz"
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


def fstats(npz):
    """median f over ancient targets, mean f over controls (the quantities whose ratio is exp(H) at the sample-median)."""
    z = np.load(npz); ic = z["is_control"]
    f = z["num"].sum(0) / z["den"].sum(0)
    return float(np.median(f[~ic])), float(f[ic].mean()), float(np.mean(np.log(f[~ic] / f[ic].mean())))

def hrow(npz):
    z = np.load(npz); ic = z["is_control"]; f = z["num"].sum(0) / z["den"].sum(0)
    return f, ic

rows = []; sel_rows = []; off_rows = []
for arm, X in ARMS.items():
  for rep in REPS:
    truth = [l.split("\t") for l in (B / f"sim/rep{rep}/truth.tsv").read_text().split("\n")[1:] if l]
    res = {}; resd = {}
    for a0 in A0S:
        nm, Ts, fl, sel, usable = meanT(rep, a0, "", X)
        if not len(Ts): continue
        res[a0] = (nm, Ts, fl)
        sel_rows.append((arm, rep, a0, len(Ts), np.mean(sel), np.mean(usable)))
        nmd, Tsd, fld, _, _ = meanT(rep, a0, "_demes", X)
        if len(Tsd): resd[a0] = (nmd, Tsd, fld)
        # offset: per draw median ancient f / mean control f, then mean over draws; true tree separately
        fr = [fstats(B / f"draws/rep{rep}/draw{X}{d}_a0{a0}.npz") for d in DRAWS if (B / f"draws/rep{rep}/draw{X}{d}_a0{a0}.npz").exists()]
        pt = B / f"draws/rep{rep}/true{X}_a0{a0}.npz"
        ft = fstats(pt) if pt.exists() else (np.nan,) * 3
        off_rows.append((arm, rep, a0, np.mean([x[0] for x in fr]), np.mean([x[1] for x in fr]),
                         np.mean([x[0] / x[1] for x in fr]), np.mean([x[2] for x in fr]), ft[0] / ft[1], ft[2]))
    for j, (_, n, a) in enumerate(truth):
        r = dict(arm=arm, rep=rep, sample=n, true_age=float(a))
        for a0 in A0S:
            if a0 not in res: continue
            nm, Ts, fl = res[a0]; k = nm.index(n)
            r[f"T_{a0}"] = Ts[:, k].mean(); r[f"sd_{a0}"] = Ts[:, k].std(ddof=1) if len(Ts) > 1 else np.nan
            r[f"flag_{a0}"] = ";".join(sorted({f[k] for f in fl if f[k]}))
            if a0 in resd:
                nm, Ts, fl = resd[a0]; k = nm.index(n)
                r[f"Td_{a0}"] = Ts[:, k].mean(); r[f"sdd_{a0}"] = Ts[:, k].std(ddof=1) if len(Ts) > 1 else np.nan
                r[f"flagd_{a0}"] = ";".join(sorted({f[k] for f in fl if f[k]}))
            for key, fn in (("Ttrue", "true"), ("Ttruedemes", "truedemes")):
                p = B / f"draws/rep{rep}/{fn}{X}_a0{a0}.npz"
                if p.exists():
                    z = np.load(p); ic = z["is_control"]; T, H, f = nr.ages_from_counts(z["num"], z["den"], ic, z["edges"], z["rates"])
                    r[f"{key}_{a0}"] = T[[str(x) for x in z["targets"][~ic]].index(n)]
        rows.append(r)

cols = ["arm", "rep", "sample", "true_age"] + [f"{p}_{a}" for a in A0S for p in ("T", "sd", "flag", "Td", "sdd", "flagd", "Ttrue", "Ttruedemes")]
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
for arm in ARMS:
    rr_ = [r for r in rows if r["arm"] == arm]
    x = np.array([r["true_age"] for r in rr_])
    for label, key in (("ARG (c(t) from ARG)", "T"), ("ARG (c(t) from true demography)", "Td"), ("true tree (c(t) from tree)", "Ttrue"), ("true tree (c(t) from true demography)", "Ttruedemes")):
        for a0 in A0S:
            y = np.array([r.get(f"{key}_{a0}", np.nan) for r in rr_], float)
            if not np.isfinite(y).any(): continue
            s = stats(x, y)
            bb = []
            for lo, hi in BINS:
                m = (x >= lo) & (x < hi) & np.isfinite(y); bb.append((y[m] - x[m]).mean() if m.any() else np.nan)
            nflag = sum(1 for r in rr_ if r.get({"T": f"flag_{a0}", "Td": f"flagd_{a0}"}.get(key, "x"))) if key in ("T", "Td") else ""
            out.append((arm, label, a0, s, bb, nflag))
with open(B / "summary_by_a0.tsv", "w") as fh:
    fh.write("arm\tsource\ta0\tn\tmean_bias\tmedian_bias\trmse\tpearson_r\tslope\tbias_0-1500\tbias_1500-3000\tbias_3000-4500\tbias_4500-6000\tn_targets_flagged_any_draw\n")
    for arm, label, a0, s, bb, nf in out:
        fh.write(f"{arm}\t{label}\t{a0}\t{s['n']}\t{s['mean_bias']:.0f}\t{s['median_bias']:.0f}\t{s['rmse']:.0f}\t{s['r']:.3f}\t{s['slope']:.3f}\t"
                 + "\t".join(f"{b:.0f}" for b in bb) + f"\t{nf}\n")
with open(B / "offset_by_rep.tsv", "w") as fh:
    fh.write("arm\trep\ta0\tmedian_ancient_f\tmean_control_f\tratio_arg_meanOverDraws\tmean_H_arg\tratio_truetree\tmean_H_truetree\n")
    for o in off_rows: fh.write("\t".join([o[0], str(o[1]), str(o[2])] + [f"{v:.5f}" for v in o[3:]]) + "\n")
with open(B / "offset_summary.tsv", "w") as fh:
    fh.write("arm\ta0\tmean_ratio_arg\tmin_ratio_arg\tmax_ratio_arg\tmean_ratio_truetree\tmean_median_ancient_f\tmean_control_f\n")
    for arm in ARMS:
        for a0 in A0S:
            o = [x for x in off_rows if x[0] == arm and x[2] == a0]
            if not o: continue
            fh.write(f"{arm}\t{a0}\t{np.mean([x[5] for x in o]):.4f}\t{min(x[5] for x in o):.4f}\t{max(x[5] for x in o):.4f}\t{np.mean([x[7] for x in o]):.4f}\t{np.mean([x[3] for x in o]):.5f}\t{np.mean([x[4] for x in o]):.5f}\n")
with open(B / "site_selection.tsv", "w") as fh:
    fh.write("arm\trep\ta0\tn_draws\tmean_selected_sites\tmean_usable_sites\tfraction\n")
    for arm, rep, a0, nd, s, u in sel_rows: fh.write(f"{arm}\t{rep}\t{a0}\t{nd}\t{s:.0f}\t{u:.0f}\t{s/u:.4f}\n")

# controls as pseudo-targets.  Leave-one-out within each arm's own control set (arm B: B-controls vs other B-controls, ~0 by
# construction), PLUS cross pseudo-targets: B-controls dated against the A-control mean and A-controls against the B-control mean
# (f values taken from the arm-A and arm-B runs of the same draw/A0; these are separate run-draw calls over the same panel tree).
def ctrl_f(npz):
    z = np.load(npz); ic = z["is_control"]; return (z["num"].sum(0) / z["den"].sum(0))[ic], z["edges"], z["rates"]
with open(B / "controls_pseudotargets.tsv", "w") as fh:
    fh.write("kind\tsource\trep\ta0\tdraw\tcontrol\tH\tT_hat\n")
    for rep in REPS:
        for a0 in A0S:
            for src, d in [("arg", d) for d in DRAWS] + [("true", "true")]:
                pa = B / (f"draws/rep{rep}/true_a0{a0}.npz" if d == "true" else f"draws/rep{rep}/draw{d}_a0{a0}.npz")
                pb = B / (f"draws/rep{rep}/trueB_a0{a0}.npz" if d == "true" else f"draws/rep{rep}/drawB{d}_a0{a0}.npz")
                if not (pa.exists() and pb.exists()): continue
                fa, edges, rates = ctrl_f(pa); fb, _, _ = ctrl_f(pb)
                for kind, tgt, ref in (("LOO_Acontrols", fa, None), ("LOO_Bcontrols", fb, None), ("B_vs_Aref", fb, fa), ("A_vs_Bref", fa, fb)):
                    for j in range(tgt.size):
                        refm = np.delete(tgt, j).mean() if ref is None else ref.mean()
                        H = np.log(tgt[j] / refm)
                        T = float(nr.invert_hazard(np.array([H]), edges, rates)[0])
                        fh.write(f"{kind}\t{src}\t{rep}\t{a0}\t{d}\t{j+1}\t{H:.5f}\t{T:.1f}\n")

# figures: rows = arms, cols = A0
def make_fig(key, sdkey, title, name):
    cmap = plt.get_cmap("tab10"); fig, axs = plt.subplots(2, 3, figsize=(15, 9.5), sharex=True, sharey=True)
    for i, (arm, lab) in enumerate((("A", "controls from A (panel population)"), ("B", "controls from B (ancient population)"))):
        for ax, a0 in zip(axs[i], A0S):
            for rep in REPS:
                rr = [r for r in rows if r["arm"] == arm and r["rep"] == rep and f"{key}_{a0}" in r]
                if not rr: continue
                ax.errorbar([r["true_age"] for r in rr], [r[f"{key}_{a0}"] for r in rr], yerr=[r[f"{sdkey}_{a0}"] for r in rr],
                            fmt="o", ms=5, color=cmap(rep - 1), alpha=0.85, capsize=2, lw=0.8, label=f"rep {rep}")
            ax.plot([0, 6500], [0, 6500], "k--", lw=1)
            ax.set_title(f"{lab}; A0 = {a0:,}", fontsize=10); ax.set_xlim(-200, 6500)
        axs[i][0].set_ylabel("inferred age T_hat (generations)\n(mean over draws, bars = across-draw SD)")
    for ax in axs[1]: ax.set_xlabel("true age (generations)")
    axs[0][0].legend(title="replicate", loc="upper left")
    fig.suptitle(title); fig.tight_layout(); fig.savefig(B / f"{name}.png", dpi=150); fig.savefig(B / f"{name}.svg")
make_fig("T", "sd", "Split-6000 + bottleneck; ancient from B, panel from A; c(t) from the ARG lineage counts; top: A controls, bottom: B controls", "inferred_vs_true")
make_fig("Td", "sdd", "Same, c(t) from the TRUE (deme A) demography", "inferred_vs_true_truedemography_ct")
print(open(B / "summary_by_a0.tsv").read()); print(open(B / "offset_summary.tsv").read()); print(open(B / "site_selection.tsv").read())
