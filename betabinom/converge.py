"""Does the estimate converge to the truth as sequence length grows?
If the residual is variance, MAPs tighten ON the truth; if bias, they tighten
somewhere else."""
import numpy as np
import bootstrap as B
T_TRUE=2500.
print(f"truth {T_TRUE:.0f}\n{'L':>7}{'seeds':>7}{'MAPs':>34}{'median':>9}{'>truth':>8}")
for L in (2e6, 8e6):
    B.L=L
    maps=[]
    for SEED in (201,202,203,204):
        LL,POS=B.per_site(T_TRUE,SEED)
        maps.append(B.grid[LL.sum(0).argmax()])
    m=np.array(maps)
    print(f"{L/1e6:>5.0f}Mb{len(m):>7}{str([int(x) for x in m]):>34}"
          f"{np.median(m):>9.0f}{(m>T_TRUE).sum():>5}/{len(m)}")
