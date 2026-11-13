"""Aggregate over draws. f = num.sum(0)/den.sum(0) per target (checked against T_hat sanity below).
Usage: python -I aggregate.py out_dir [dates.txt]  ; dates file (scratchpad only, one value per ancient sample, 'NA' allowed)
prints aggregates only."""
import sys, glob, re, numpy as np
from scipy import stats
files = sorted(glob.glob(sys.argv[1] + "/draw_*.deid.npz"), key=lambda s: int(re.search(r"draw_(\d+)", s).group(1)))
draws = [int(re.search(r"draw_(\d+)", s).group(1)) for s in files]
Z = [np.load(f, allow_pickle=False) for f in files]
ic = Z[0]["is_control"]; nA = int((~ic).sum())
assert all((z["is_control"] == ic).all() and (z["targets"] == Z[0]["targets"]).all() for z in Z)
f = np.array([z["num"].sum(0) / z["den"].sum(0) for z in Z])      # draws x targets
T = np.array([z["T_hat"] for z in Z])                              # draws x ancient
print("draws", draws); print("n ancient", nA, "n controls", int(ic.sum()))
fa, fc = f[:, ~ic], f[:, ic]
print("\nPer-draw: ancient median f, control mean f, offset, ancient median T_hat, #T_hat==0")
for d, a, c, t in zip(draws, fa, fc, T): print(d, "%.4f %.4f %.4f %.0f %d" % (np.median(a), c.mean(), np.median(a)-c.mean(), np.median(t), (t==0).sum()))
def icc(X):  # one-way ICC(1): rows=draws (k), cols=samples (n)
    k, n = X.shape; gm = X.mean(); sm = X.mean(0)
    msb = k * ((sm - gm) ** 2).sum() / (n - 1); msw = ((X - sm) ** 2).sum() / (n * (k - 1))
    return (msb - msw) / (msb + (k - 1) * msw)
for name, X in [("f ancient", fa), ("T_hat ancient", T), ("f controls", fc)]:
    within = X.std(0, ddof=1).mean(); between = X.mean(0).std(ddof=1)
    print("\n%s: mean within-sample SD across draws %.4g; SD across sample means %.4g; ratio within/between %.2f; ICC(1) %.3f" % (name, within, between, within/between, icc(X)))
k = len(draws); iu = np.triu_indices(k, 1)
for name, X in [("f ancient", fa), ("T_hat ancient", T)]:
    P = np.corrcoef(X); S = stats.spearmanr(X.T).correlation
    print("%s between-draw correlation over %d ancient samples: Pearson mean %.3f (min %.3f max %.3f); Spearman mean %.3f (min %.3f max %.3f)" % (name, nA, P[iu].mean(), P[iu].min(), P[iu].max(), S[iu].mean(), S[iu].min(), S[iu].max()))
print("\nOffset across draws: ancient median f mean %.4f SD %.4f; control mean f mean %.4f SD %.4f; offset mean %.4f SD %.4f min %.4f max %.4f" % (np.median(fa,1).mean(), np.median(fa,1).std(ddof=1), fc.mean(1).mean(), fc.mean(1).std(ddof=1), (np.median(fa,1)-fc.mean(1)).mean(), (np.median(fa,1)-fc.mean(1)).std(ddof=1), (np.median(fa,1)-fc.mean(1)).min(), (np.median(fa,1)-fc.mean(1)).max()))
print("per-control f range (min of draw means .. max):", np.round(fc.mean(0),4))
print("fraction of ancient samples with f > max control f, per draw:", [round(float((a > c.max()).mean()),2) for a,c in zip(fa,fc)])
# averaged over draws using pooled num/den too
fpool = np.array([z["num"].sum(0) for z in Z]).sum(0) / np.array([z["den"].sum(0) for z in Z]).sum(0)
Tm = T.mean(0)
print("\nT_hat mean over draws: median %.0f, min %.0f, max %.0f, SD across samples %.1f" % (np.median(Tm), Tm.min(), Tm.max(), Tm.std(ddof=1)))
if len(sys.argv) > 2:
    y = np.array([np.nan if v == "NA" else float(v) for v in open(sys.argv[3-1]).read().split()]); m = ~np.isnan(y); yy = y[m]
    def rep(lbl, x):
        p = stats.pearsonr(x[m], yy); s = stats.spearmanr(x[m], yy)
        return "%s: Pearson %.3f (p=%.2g), Spearman %.3f (p=%.2g)" % (lbl, p[0], p[1], s.correlation, s.pvalue)
    print("\nn dated", int(m.sum()))
    for d, t in zip(draws, T): print(rep("draw %d T_hat" % d, t))
    print(rep("mean T_hat over draws", Tm)); print(rep("mean f over draws", fa.mean(0)))
    r1 = [stats.pearsonr(t[m], yy)[0] for t in T]; print("single-draw Pearson: mean %.3f min %.3f max %.3f" % (np.mean(r1), min(r1), max(r1)))
    rng = np.random.default_rng(1)
    for kk in (1, 2, 5, k):
        rs = [stats.pearsonr(T[rng.choice(k, kk, replace=False)].mean(0)[m], yy)[0] for _ in range(200)]
        print("subsets of %d draws: mean Pearson %.3f (range %.3f..%.3f)" % (kk, np.mean(rs), min(rs), max(rs)))
    # does f (not T_hat) carry rank info among ancient
    print(rep("per-sample f residual after removing draw mean; sanity f-pooled", fpool[~ic]))
