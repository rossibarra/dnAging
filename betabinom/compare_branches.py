"""#2 Head-to-head: betabinom (conditioning on the ARG at T) vs diffusion
(conditioning on the present-day panel count) over the same 300 simulations."""
import csv, numpy as np, sys

D="/Users/jeffreyross-ibarra/Projects/mutrates/simbatch300/"
bb={r["sim"]:r for r in csv.DictReader(open(D+"estimates.tsv"),delimiter='\t')}
df={r["sim"]:r for r in csv.DictReader(open(D+"estimates_diff.tsv"),delimiter='\t')}
sims=[s for s in bb if s in df and not df[s].get("fail") and df[s].get("MAP_diff")
      and not bb[s].get("fail")]
sims.sort()
g=lambda d,k: np.array([float(d[s][k]) for s in sims])
age=g(bb,"age"); Ne=g(bb,"Ne"); nmod=g(bb,"n_modern")
mb=np.array([float(bb[s].get("MAP") or bb[s].get("MAP_bb")) for s in sims])
md=g(df,"MAP_diff")
sb=np.array([float(bb[s].get("boot_sd") or bb[s].get("boot_sd_bb")) for s in sims])
sd_=g(df,"boot_sd_diff")
print(f"n = {len(sims)} simulations with both estimators\n")

def block(name, m, s):
    e=m-age; tau=e/(2*Ne)
    A=np.c_[np.ones_like(age),age]; sl,ic=np.linalg.lstsq(A,m,rcond=None)[0][::-1]
    r=np.corrcoef(age,m)[0,1]
    cov=np.mean(np.abs(e)<1.96*s)
    # LOO constant correction, in tau
    lo=np.array([np.delete(tau,i).mean() for i in range(len(e))])
    el=m-2*Ne*lo-age
    print(f"{name:<12}  bias {e.mean():>+7.0f}  med {np.median(e):>+7.0f}  "
          f"SD {e.std(ddof=1):>6.0f}  RMSE {np.sqrt((e**2).mean()):>6.0f}")
    print(f"{'':<12}  slope {sl:.3f}  intercept {ic:>+7.0f}  r {r:.3f}  "
          f"95% cov {cov:.2f}  med SD {np.median(s):>5.0f}")
    print(f"{'':<12}  tau-bias {tau.mean():.4f} +- {tau.std(ddof=1):.4f}   "
          f"after LOO tau-correction: RMSE {np.sqrt((el**2).mean()):>6.0f} "
          f"bias {el.mean():>+5.0f} med|err| {np.median(abs(el)):>5.0f}")
    return e, tau, el

print("=== RAW ACCURACY ===")
eb,tb,elb=block("betabinom",mb,sb); print()
ed,td,eld=block("diffusion",md,sd_)

print("\n=== PARAMETER DEPENDENCE of the error (t-stats) ===")
def reg(y,label):
    X=np.c_[np.ones_like(age),Ne/1e4,nmod,age/1e3]
    b,res,*_=np.linalg.lstsq(X,y,rcond=None)
    resid=y-X@b; dof=len(y)-X.shape[1]
    se=np.sqrt(np.diag(np.linalg.inv(X.T@X))*(resid@resid)/dof)
    R2=1-(resid@resid)/((y-y.mean())@(y-y.mean()))
    print(f"  {label:<26} " + "  ".join(
        f"{n}: t={b[i]/se[i]:>+6.2f}" for i,n in
        enumerate(["const","Ne/10k","n_mod","age/1k"]) if i) + f"   R2={R2:.3f}")
reg(eb,"betabinom err (gen)"); reg(tb,"betabinom err (tau)")
reg(ed,"diffusion err (gen)"); reg(td,"diffusion err (tau)")

print("\n=== HEAD TO HEAD, after each gets its own LOO tau-correction ===")
wb=np.abs(elb)<np.abs(eld)
print(f"  betabinom closer to truth in {wb.sum()}/{len(sims)} sims ({100*wb.mean():.0f}%)")
print(f"  median |err|   betabinom {np.median(abs(elb)):>6.0f}   "
      f"diffusion {np.median(abs(eld)):>6.0f}")
print(f"  RMSE           betabinom {np.sqrt((elb**2).mean()):>6.0f}   "
      f"diffusion {np.sqrt((eld**2).mean()):>6.0f}")
d=np.abs(elb)-np.abs(eld)
print(f"  paired diff in |err| (bb - diff): {d.mean():>+6.0f} +- "
      f"{d.std(ddof=1)/np.sqrt(len(d)):.0f} (SE)   t={d.mean()/(d.std(ddof=1)/np.sqrt(len(d))):>+.2f}")
print(f"  correlation of the two errors: r={np.corrcoef(elb,eld)[0,1]:.3f}")

print("\n=== by Ne quartile: bias in generations / in tau ===")
qs=np.quantile(Ne,[0,.25,.5,.75,1.]); qs[-1]+=1
print(f"  {'Ne bin':<22} {'n':>4} {'bb gen':>8} {'df gen':>8} {'bb tau':>8} {'df tau':>8}")
for i in range(4):
    s=(Ne>=qs[i])&(Ne<qs[i+1])
    print(f"  [{qs[i]:.3g}, {qs[i+1]:.3g}){'':<4} {s.sum():>4} {eb[s].mean():>8.0f} "
          f"{ed[s].mean():>8.0f} {tb[s].mean():>8.4f} {td[s].mean():>8.4f}")

for nm,e,t in (("betabinom",eb,tb),("diffusion",ed,td)):
    bg=[e[(Ne>=qs[i])&(Ne<qs[i+1])].mean() for i in range(4)]
    bt=[t[(Ne>=qs[i])&(Ne<qs[i+1])].mean() for i in range(4)]
    print(f"  {nm} spread across Ne:  {max(bg)/min(bg):.2f}x in gen, "
          f"{max(bt)/min(bt):.2f}x in tau")
