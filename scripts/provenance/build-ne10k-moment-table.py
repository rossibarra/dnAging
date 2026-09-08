import json
import os
import sys

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, "/Users/jeffreyross-ibarra/src/dnAging")
from precompute_freq_trajectory_moments import MomentEngine

NE = float(os.environ.get("TABLE_NE", "10000"))
AGE = np.geomspace(10.0, 10 * NE, 100)
TGRID = np.arange(0.0, 15_000.0 + 10.0, 20.0)
OUT = os.environ.get("TABLE_OUT", "/tmp/dnAging-ne10000-infinite-100mb/frequency_table_ne10k.npz")
eng = MomentEngine(26)
B = eng.B
x0 = 1 / (2 * NE)
Cgrid = [expm(B * (T / (2 * NE))) for T in TGRID]
table = np.full((1, 26, len(AGE), len(TGRID)), np.nan, dtype=np.float32)

for ia, age in enumerate(AGE):
    tau_i = age / (2 * NE)
    Mpres = eng._moms(tau_i, x0)
    den = np.full(26, np.nan)
    abs_den = np.full(26, np.nan)
    for d0 in range(1, 27):
        terms = np.array([c * Mpres[m] for m, c in eng.coeff[d0].items()])
        den[d0 - 1] = terms.sum()
        abs_den[d0 - 1] = np.abs(terms).sum()
    for it, T in enumerate(TGRID):
        tau_T = T / (2 * NE)
        if tau_T >= tau_i:
            table[0, :, ia, it] = 0.0
            continue
        Mu1 = eng._moms(tau_i - tau_T, x0)
        EjX = Cgrid[it][:, :eng.K] @ Mu1[1:eng.K + 1]
        for d0 in range(1, 27):
            terms = np.array([c * EjX[m] for m, c in eng.coeff[d0].items()])
            num = terms.sum()
            reliable = (np.isfinite(num) and np.isfinite(den[d0 - 1]) and
                        num != 0 and den[d0 - 1] != 0 and
                        np.abs(terms).sum() / abs(num) <= eng.max_cancellation and
                        abs_den[d0 - 1] / abs(den[d0 - 1]) <= eng.max_cancellation)
            if reliable:
                p = num / den[d0 - 1]
                if -1e-12 <= p <= 1 + 1e-12:
                    table[0, d0 - 1, ia, it] = np.clip(p, 0, 1)
    print(f"age {ia+1}/100: {age:.3f}", flush=True)

meta = {
    "n_sample": 26, "min_n": 26, "Ne": NE,
    "method": "optimized neutral WF moment recursion", "format_version": 4,
    "planes": ["table (E[p_T])"],
}
np.savez_compressed(OUT, table=table, d0=np.arange(1, 27), age=AGE,
                    n_panel=np.array([26]), age_tau=AGE / (2 * NE), Tgrid=TGRID,
                    n_sample=26, min_n=26, meta=json.dumps(meta))
print(f"wrote {OUT} {table.shape}; NaNs={np.isnan(table).sum()}", flush=True)
