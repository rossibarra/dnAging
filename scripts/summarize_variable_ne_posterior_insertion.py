#!/usr/bin/env python3
"""Aggregate posterior-ARG insertion results across variable-Ne simulations; prove."""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=Path("msprime_variable_ne_error"))
    p.add_argument("--n-replicates",type=int,default=10)
    p.add_argument("--scenario",choices=("true_ne","estimated_ne"),required=True)
    p.add_argument("--output",type=Path,default=Path("msprime_variable_ne_error/insertion_posterior/summary"))
    a=p.parse_args(); rows=[]
    for rep in range(1,a.n_replicates+1):
        label=f"replicate_{rep:03d}"
        truth=json.loads((a.root/"simulations"/label/"metadata.json").read_text())["true_ancient_ages"]
        with (a.root/"insertion_posterior"/a.scenario/"merged"/label/"ages_table.tsv").open() as h:
            estimates=list(csv.DictReader(h,delimiter="\t"))
        if len(estimates) != len(truth): raise ValueError(f"{label}: {len(estimates)} estimates != {len(truth)} truths")
        for i,(est,t) in enumerate(zip(estimates,truth),1):
            row={"replicate":rep,"sample":est["sample"],"true_T":float(t)}
            row.update({k:float(est[k]) for k in ("map_T","mean_T","median_T","ci95_lower_T","ci95_upper_T")})
            row["contains_true_95"]=int(row["ci95_lower_T"] <= row["true_T"] <= row["ci95_upper_T"])
            rows.append(row)
    true=np.array([r["true_T"] for r in rows]); maps=np.array([r["map_T"] for r in rows])
    means=np.array([r["mean_T"] for r in rows]); covered=np.array([r["contains_true_95"] for r in rows],bool)
    ne_label=("true replicate-specific trajectory" if a.scenario=="true_ne" else
              "replicate-specific ARGtest posterior-mean trajectory")
    result={"n_replicates":a.n_replicates,"n_samples":len(rows),"arg_draws_per_replicate":50,
            "scenario":a.scenario,"ne":ne_label,
            "map_bias_estimated_minus_true":float(np.mean(maps-true)),"map_mae":float(np.mean(abs(maps-true))),
            "map_rmse":float(np.sqrt(np.mean((maps-true)**2))),
            "mean_bias_estimated_minus_true":float(np.mean(means-true)),
            "ci95_coverage_count":int(covered.sum()),"ci95_coverage_fraction":float(covered.mean())}
    a.output.mkdir(parents=True,exist_ok=True)
    with (a.output/"ages.tsv").open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(rows)
    (a.output/"run.json").write_text(json.dumps(result,indent=2)+"\n")
    fig,ax=plt.subplots(figsize=(6,6))
    ax.vlines(true, [r["ci95_lower_T"] for r in rows],
              [r["ci95_upper_T"] for r in rows], color="0.7", lw=.6)
    ax.scatter(true,maps,s=10,alpha=.55)
    lim=max(15000,float(max(true.max(),maps.max()))); ax.plot([0,lim],[0,lim],"k--",lw=1)
    ax.set(xlabel="True sample age (generations)",ylabel="Estimated MAP age (generations)",
           title=f"Insertion: 50 posterior ARG draws + {a.scenario.replace('_',' ')}",xlim=(0,lim),ylim=(0,lim))
    fig.tight_layout(); fig.savefig(a.output/"map_vs_true.png",dpi=220)
    print(json.dumps(result,indent=2))
if __name__=="__main__": main()
