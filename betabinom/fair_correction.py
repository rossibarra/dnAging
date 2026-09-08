"""Give each estimator a leave-one-out LINEAR calibration in diffusion time,
with panel size as a covariate, so neither is penalised for a bias that a
user could calibrate away."""
import csv, numpy as np
D="/Users/jeffreyross-ibarra/Projects/mutrates/simbatch300/"
bb={r["sim"]:r for r in csv.DictReader(open(D+"estimates.tsv"),delimiter='\t')}
df={r["sim"]:r for r in csv.DictReader(open(D+"estimates_diff.tsv"),delimiter='\t')}
sims=sorted(s for s in bb if s in df and df[s].get("MAP_diff"))
G=lambda d,k: np.array([float(d[s][k]) for s in sims])
age=G(bb,"age"); Ne=G(bb,"Ne"); nm=G(bb,"n_modern")
M={"betabinom":G(bb,"MAP"), "diffusion":G(df,"MAP_diff")}
S={"betabinom":G(bb,"boot_sd"), "diffusion":G(df,"boot_sd_diff")}
n=len(sims); print(f"n={n}\n")
out={}
for name,m in M.items():
    tau_m=m/(2*Ne); tau_t=age/(2*Ne)
    X=np.c_[np.ones(n), tau_m, 1.0/nm]          # calibrate tau_true ~ f(tau_hat, panel)
    pred=np.empty(n)
    for i in range(n):
        k=np.arange(n)!=i
        b=np.linalg.lstsq(X[k],tau_t[k],rcond=None)[0]
        pred[i]=X[i]@b
    e=pred*2*Ne-age
    out[name]=e
    # propagate the SD through the same linear map (slope on tau cancels units)
    b_all=np.linalg.lstsq(X,tau_t,rcond=None)[0]
    s=S[name]*b_all[1]
    print(f"{name:<10} LOO-calibrated:  bias {e.mean():>+6.0f}  med|err| "
          f"{np.median(abs(e)):>5.0f}  RMSE {np.sqrt((e**2).mean()):>5.0f}  "
          f"SD {e.std(ddof=1):>5.0f}  95%cov {np.mean(abs(e)<1.96*s):.2f}")
    Z=np.c_[np.ones(n),Ne/1e4,nm,age/1e3]
    bq,*_=np.linalg.lstsq(Z,e,rcond=None); r=e-Z@bq
    se=np.sqrt(np.diag(np.linalg.inv(Z.T@Z))*(r@r)/(n-4))
    print(f"{'':<10} residual dependence:  " + "  ".join(
        f"{q}: t={bq[j]/se[j]:>+6.2f}" for j,q in
        enumerate(["","Ne/10k","n_mod","age/1k"]) if j)
        + f"   R2={1-(r@r)/((e-e.mean())@(e-e.mean())):.3f}")
    print()
a,b=out["betabinom"],out["diffusion"]
d=abs(a)-abs(b)
print(f"betabinom wins {int((abs(a)<abs(b)).sum())}/{n};  paired mean |err| diff "
      f"{d.mean():+.0f} +- {d.std(ddof=1)/np.sqrt(n):.0f}  t={d.mean()/(d.std(ddof=1)/np.sqrt(n)):+.2f}")
print("\nby panel size, median |err| after calibration:")
for lo,hi in ((10,18),(18,26),(26,34),(34,41)):
    k=(nm>=lo)&(nm<hi)
    print(f"  n_modern [{lo},{hi})  n={k.sum():>3}   betabinom {np.median(abs(a[k])):>5.0f}"
          f"   diffusion {np.median(abs(b[k])):>5.0f}")
