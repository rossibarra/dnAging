"""Is the estimator biased, or just noisy?  Many seeds at one true age."""
import numpy as np, sys
import bootstrap as B          # reuse per_site / grid / ci
T_TRUE=2500.
maps=[]; cov_c=0; cov_b=0; widths=[]
rng=np.random.default_rng(1)
print(f"{'seed':>5}{'MAP':>7}{'composite 95%':>19}{'bootstrap 95%':>21}")
for SEED in range(101,113):
    LL,POS=B.per_site(T_TRUE,SEED)
    tot=LL.sum(0); mapv=B.grid[tot.argmax()]; maps.append(mapv)
    clo,chi=B.ci(tot)
    blk=np.minimum((POS/(B.L/B.NBLOCK)).astype(int),B.NBLOCK-1)
    sums=np.zeros((B.NBLOCK,len(B.grid)))
    for b in range(B.NBLOCK):
        s=blk==b
        if s.any(): sums[b]=LL[s].sum(0)
    bm=np.array([B.grid[sums[rng.integers(0,B.NBLOCK,B.NBLOCK)].sum(0).argmax()]
                 for _ in range(300)])
    blo,bhi=np.percentile(bm,[2.5,97.5]); widths.append(bhi-blo)
    cov_c+= clo<=T_TRUE<=chi; cov_b+= blo<=T_TRUE<=bhi
    print(f"{SEED:>5}{mapv:>7.0f}{clo:>9.0f}-{chi:<9.0f}{blo:>10.0f}-{bhi:<10.0f}")
m=np.array(maps)
print(f"\ntruth {T_TRUE:.0f}   n={len(m)} seeds, L={B.L/1e6:g}Mb")
print(f"  MAP  median {np.median(m):.0f}   mean {m.mean():.0f}   "
      f"IQR {np.percentile(m,25):.0f}-{np.percentile(m,75):.0f}")
print(f"  fraction of MAPs above truth: {(m>T_TRUE).mean():.2f}   (0.5 if unbiased)")
print(f"  coverage: composite {cov_c}/{len(m)}   block-bootstrap {cov_b}/{len(m)}")
print(f"  median bootstrap CI width: {np.median(widths):.0f} generations")
