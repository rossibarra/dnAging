"""Three-way: betabinom, diffusion (float64, as MATH.md's pipeline computes it),
and diffusion at 50 digits (the same math without the cancellation loss)."""
import csv, numpy as np
D="/Users/jeffreyross-ibarra/Projects/mutrates/simbatch300/"
L=lambda f: {r["sim"]:r for r in csv.DictReader(open(D+f),delimiter='\t')}
bb,d64,dmp=L("estimates.tsv"),L("estimates_diff.tsv"),L("estimates_diffmp.tsv")
sims=sorted(s for s in bb if s in d64 and s in dmp
            and d64[s].get("MAP_diff") and dmp[s].get("MAP_mp"))
G=lambda d,k: np.array([float(d[s][k]) for s in sims])
age=G(bb,"age"); Ne=G(bb,"Ne"); nm=G(bb,"n_modern"); n=len(sims)
M={"betabinom":G(bb,"MAP"),"diffusion f64":G(d64,"MAP_diff"),
   "diffusion 50dp":G(dmp,"MAP_mp")}
S={"betabinom":G(bb,"boot_sd"),"diffusion f64":G(d64,"boot_sd_diff"),
   "diffusion 50dp":G(dmp,"boot_sd_mp")}
print(f"n = {n} simulations, all three estimators\n")
print("=== RAW ===")
print(f"{'':<16}{'bias':>8}{'med err':>9}{'SD':>8}{'RMSE':>8}{'slope':>8}{'r':>7}{'95%cov':>8}")
for k,m in M.items():
    e=m-age; sl=np.linalg.lstsq(np.c_[np.ones(n),age],m,rcond=None)[0][1]
    print(f"{k:<16}{e.mean():>+8.0f}{np.median(e):>+9.0f}{e.std(ddof=1):>8.0f}"
          f"{np.sqrt((e**2).mean()):>8.0f}{sl:>8.3f}{np.corrcoef(age,m)[0,1]:>7.3f}"
          f"{np.mean(abs(e)<1.96*S[k]):>8.2f}")

print("\n=== after a leave-one-out LINEAR calibration in diffusion time (+1/n_panel) ===")
corr={}
for k,m in M.items():
    tm=m/(2*Ne); tt=age/(2*Ne); X=np.c_[np.ones(n),tm,1.0/nm]
    pred=np.array([X[i]@np.linalg.lstsq(X[np.arange(n)!=i],tt[np.arange(n)!=i],
                                        rcond=None)[0] for i in range(n)])
    e=pred*2*Ne-age; corr[k]=e
    b=np.linalg.lstsq(X,tt,rcond=None)[0]; s=S[k]*b[1]
    Z=np.c_[np.ones(n),Ne/1e4,nm,age/1e3]
    bq=np.linalg.lstsq(Z,e,rcond=None)[0]; r=e-Z@bq
    se=np.sqrt(np.diag(np.linalg.inv(Z.T@Z))*(r@r)/(n-4))
    print(f"{k:<16}bias {e.mean():>+6.0f}  med|err| {np.median(abs(e)):>5.0f}  "
          f"RMSE {np.sqrt((e**2).mean()):>5.0f}  95%cov {np.mean(abs(e)<1.96*s):.2f}   "
          f"resid t: Ne {bq[1]/se[1]:>+5.2f}  n_mod {bq[2]/se[2]:>+6.2f}  "
          f"age {bq[3]/se[3]:>+5.2f}")

print("\n=== median |err| by panel size (calibrated) ===")
print(f"  {'n_modern':<14}{'n':>4}" + "".join(f"{k:>17}" for k in M))
for lo,hi in ((10,18),(18,26),(26,34),(34,41)):
    k=(nm>=lo)&(nm<hi)
    print(f"  [{lo},{hi}){'':<7}{k.sum():>4}" +
          "".join(f"{np.median(abs(corr[q][k])):>17.0f}" for q in M))

print("\n=== paired head-to-head (calibrated |err|) ===")
a=corr["betabinom"]
for k in ("diffusion f64","diffusion 50dp"):
    b=corr[k]; d=abs(a)-abs(b); se=d.std(ddof=1)/np.sqrt(n)
    print(f"  betabinom vs {k:<16} wins {int((abs(a)<abs(b)).sum()):>3}/{n}   "
          f"mean |err| diff {d.mean():>+6.0f} +- {se:.0f}   t={d.mean()/se:>+6.2f}")
d=abs(corr["diffusion f64"])-abs(corr["diffusion 50dp"]); se=d.std(ddof=1)/np.sqrt(n)
print(f"  f64 vs 50dp (same math)        {'':>12} "
      f"mean |err| diff {d.mean():>+6.0f} +- {se:.0f}   t={d.mean()/se:>+6.2f}")
print(f"\n  MAP correlation f64 vs 50dp: r={np.corrcoef(M['diffusion f64'],M['diffusion 50dp'])[0,1]:.3f}"
      f"   mean |MAP diff| {np.abs(M['diffusion f64']-M['diffusion 50dp']).mean():.0f} gen")
