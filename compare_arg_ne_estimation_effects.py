#!/usr/bin/env python3
"""Compare true/estimated ARG and Ne on the same first ten simulations."""
import csv, json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path("msprime_variable_ne_error")
def read(path):
    with path.open() as h: return list(csv.DictReader(h,delimiter="\t"))
base=[r for r in read(ROOT/"insertion_results"/"age_demography_summary.tsv") if int(r["replicate"]) <= 10]
true=read(ROOT/"insertion_posterior"/"true_ne"/"summary"/"ages.tsv")
estimated=read(ROOT/"insertion_posterior"/"estimated_ne"/"summary"/"ages.tsv")
sets={"true_ARG_true_Ne":base,"estimated_ARG_true_Ne":true,"estimated_ARG_estimated_Ne":estimated}
stats={}
for name,rows in sets.items():
    truth=np.array([float(r["true_T"]) for r in rows]); maps=np.array([float(r["map_T"]) for r in rows])
    lo=np.array([float(r["ci95_lower_T"]) for r in rows]); hi=np.array([float(r["ci95_upper_T"]) for r in rows])
    err=maps-truth
    stats[name]={"n":len(rows),"bias":float(err.mean()),"mae":float(abs(err).mean()),
                 "rmse":float(np.sqrt(np.mean(err**2))),"ci95_coverage":float(np.mean((lo<=truth)&(truth<=hi)))}
stats["incremental_effects"]={
    "ARG_estimation_delta_MAE":stats["estimated_ARG_true_Ne"]["mae"]-stats["true_ARG_true_Ne"]["mae"],
    "Ne_estimation_delta_MAE":stats["estimated_ARG_estimated_Ne"]["mae"]-stats["estimated_ARG_true_Ne"]["mae"]}
out=ROOT/"insertion_posterior"/"comparison"; out.mkdir(parents=True,exist_ok=False)
(out/"run.json").write_text(json.dumps(stats,indent=2)+"\n")
fig,axes=plt.subplots(1,3,figsize=(15,4.8),sharex=True,sharey=True)
for ax,(name,rows) in zip(axes,sets.items()):
    truth=np.array([float(r["true_T"]) for r in rows]); maps=np.array([float(r["map_T"]) for r in rows])
    ax.scatter(truth,maps,s=12,alpha=.6); ax.plot([0,15000],[0,15000],"k--",lw=1)
    ax.set(title=name.replace("_","\n"),xlabel="True T",xlim=(0,15000),ylim=(0,15000))
axes[0].set_ylabel("MAP T"); fig.tight_layout(); fig.savefig(out/"three_way_comparison.png",dpi=220)
print(json.dumps(stats,indent=2))
