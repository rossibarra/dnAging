"""Chunk/pooling check on a run-draw NPZ (per-block num/den). No dates in this file.
usage: python -I block_check.py NPZ OUTDIR [--ages FILE]   (FILE: one value per ancient column, 'nan' if undated;
kept OUTSIDE the repo; only aggregate statistics are written)."""
import sys, json, numpy as np
from scipy import stats
sys.path.insert(0, "/quobyte/jrigrp/beil/logan_collab/dnAging/pipeline")
import nonsharing_ratio as nr
npz, out = sys.argv[1], sys.argv[2]
ages = np.loadtxt(sys.argv[4]) if len(sys.argv) > 4 else None
z = np.load(npz, allow_pickle=True)
num, den, ic = z["num"], z["den"], z["is_control"]; edges, rates = z["edges"], z["rates"]
nb = num.shape[0]; NCH = 10
chunks = np.array_split(np.arange(nb), NCH)
L = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); L.append(s)
def est(idx):
    T, H, f = nr.ages_from_counts(num[idx], den[idx], ic, edges, rates)
    fc = (num[idx].sum(0) / den[idx].sum(0))[ic]
    return T, H, f, fc
def corr(T, m):
    ok = m & np.isfinite(T)
    if ok.sum() < 5 or np.ptp(T[ok]) == 0: return np.nan, np.nan
    return stats.pearsonr(T[ok], ages[ok])[0], stats.spearmanr(T[ok], ages[ok])[0]
m = np.isfinite(ages) if ages is not None else None
rows = []
P("chunk blocks first_Mb last_Mb ctrl_f anc_med_f offset(anc/ctrl) sd_f_anc medT frac_flag r_pearson r_spearman")
for ci, idx in enumerate(chunks):
    T, H, f, fc = est(idx)
    r, rho = corr(T, m) if m is not None else (np.nan, np.nan)
    flag = np.mean(H <= 0)
    rows.append((ci, np.median(f), fc.mean(), r, rho))
    P(ci, f"{idx[0]}-{idx[-1]}", idx[0] * 5, (idx[-1] + 1) * 5, f"{fc.mean():.4f}", f"{np.median(f):.4f}",
      f"{np.median(f) / fc.mean():.3f}", f"{f.std(ddof=1):.4f}", f"{np.median(T):.0f}", f"{flag:.2f}",
      f"{r:.3f}", f"{rho:.3f}")
rows = np.array(rows)
P(f"\nacross 10 chunks: mean offset {np.mean(rows[:,1]/rows[:,2]):.3f} SD {np.std(rows[:,1]/rows[:,2],ddof=1):.3f}; "
  f"mean r {np.nanmean(rows[:,3]):.3f} SD {np.nanstd(rows[:,3],ddof=1):.3f} min {np.nanmin(rows[:,3]):.3f} max {np.nanmax(rows[:,3]):.3f}")
P("offset vs chunk index spearman: %.2f (p=%.2f)" % stats.spearmanr(rows[:, 0], rows[:, 1] / rows[:, 2]))
# per-sample f across chunks: SD across chunks of per-sample f, relative; between-sample SD of f using all
Fch = np.array([est(i)[2] for i in chunks])  # chunk x sample
Fall = est(np.arange(nb))[2]
P(f"\nper-sample f: mean over samples of SD across chunks {Fch.std(0,ddof=1).mean():.4f} "
  f"(mean f {Fch.mean():.4f}, CV {Fch.std(0,ddof=1).mean()/Fch.mean():.2f}); between-sample SD of all-chunk f {Fall.std(ddof=1):.4f}")
# split-half reliability of per-sample f (odd vs even chunks), 200 random splits
rng = np.random.default_rng(1)
rel = []
for _ in range(200):
    p = rng.permutation(NCH); a, b = np.concatenate([chunks[i] for i in p[:5]]), np.concatenate([chunks[i] for i in p[5:]])
    rel.append(stats.pearsonr(est(a)[2], est(b)[2])[0])
P(f"split-half (5 vs 5 chunks) correlation of per-sample f across 50 ancient samples: mean {np.mean(rel):.3f} SD {np.std(rel):.3f}")
# pooling curve
P("\npooling curve: k chunks chosen at random (500 draws; k=10 is one set)")
P("k  mean_r  sd_r  frac_r>0.3  mean_rho  sd_f_between_samples(mean)  sd_of_T_median_over_draws  mean_offset")
curve = {}
for k in (1, 2, 3, 5, 10):
    rs, rhos, sds, offs = [], [], [], []
    for _ in range(500 if k < 10 else 1):
        sel = np.sort(rng.choice(NCH, k, replace=False)); idx = np.concatenate([chunks[i] for i in sel])
        T, H, f, fc = est(idx); r, rho = corr(T, m) if m is not None else (np.nan, np.nan)
        rs.append(r); rhos.append(rho); sds.append(f.std(ddof=1) / f.mean()); offs.append(np.median(f) / fc.mean())
    rs = np.array(rs)
    P(k, f"{np.nanmean(rs):.3f} {np.nanstd(rs):.3f} {np.nanmean(rs>0.3):.3f} {np.nanmean(rhos):.3f} {np.mean(sds):.3f} - {np.mean(offs):.3f}")
    curve[k] = rs
# per-sample noise of f vs k: SD over random k-subsets of a sample's f / mean f
P("\nper-sample relative SD of f across random k-chunk subsets (mean over samples; CV~1/sqrt(k) expected if noise-dominated)")
for k in (1, 2, 5):
    Fs = np.array([est(np.concatenate([chunks[i] for i in rng.choice(NCH, k, replace=False)]))[2] for _ in range(300)])
    P(k, f"{(Fs.std(0, ddof=1) / Fs.mean(0)).mean():.3f}")
# per-block r too
P("\nper-block (5 Mb) r_pearson:")
if m is not None:
    rb = [corr(est(np.array([b]))[0], m)[0] for b in range(nb)]
    P(" ".join(f"{x:.2f}" for x in rb)); P(f"mean {np.nanmean(rb):.3f} SD {np.nanstd(rb,ddof=1):.3f}; n(r>0)={np.sum(np.array(rb)>0)}/{nb}")
open(out + "/block_check_output.txt", "w").write("\n".join(L) + "\n")
