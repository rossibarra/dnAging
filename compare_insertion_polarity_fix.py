#!/usr/bin/env python3
"""Compare buggy and polarity-fixed estimated-ARG insertion benchmarks."""
import csv,json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from betabinom_real_data import summarize

R=Path("msprime_variable_ne_error")
def read(path):
    with path.open() as h: return list(csv.DictReader(h,delimiter="\t"))
def merged_rows(tag):
    rows=[]
    for rep in range(1,11):
        lab=f"replicate_{rep:03d}"; truth=json.loads((R/"simulations"/lab/"metadata.json").read_text())["true_ancient_ages"]
        est=read(R/"insertion_posterior"/tag/"merged"/lab/"ages_table.tsv")
        for x,t in zip(est,truth): rows.append({"replicate":rep,"sample":x["sample"],"true_T":float(t),**{k:float(x[k]) for k in ("map_T","ci95_lower_T","ci95_upper_T")}})
    return rows
def metrics(rows):
    t=np.array([float(x["true_T"]) for x in rows]); m=np.array([float(x["map_T"]) for x in rows]); e=m-t
    return {"n":len(rows),"bias":float(e.mean()),"mae":float(abs(e).mean()),"rmse":float(np.sqrt(np.mean(e**2))),
            "ci95_coverage":float(np.mean([(x["ci95_lower_T"]<=x["true_T"]<=x["ci95_upper_T"]) for x in rows]))}
def main():
    old_true=read(R/"insertion_posterior"/"true_ne"/"summary"/"ages.tsv")
    old_est=read(R/"insertion_posterior"/"estimated_ne"/"summary"/"ages.tsv")
    sets={"buggy_ARG + true Ne":old_true,"fixed_ARG + true Ne":merged_rows("true_ne_polarity_fixed_merged"),
          "buggy_ARG + estimated Ne":old_est,"fixed_ARG + estimated Ne":merged_rows("estimated_ne_polarity_fixed_merged")}
    stats={k:metrics(v) for k,v in sets.items()}
    out=R/"insertion_posterior"/"polarity_fix_comparison"; out.mkdir(parents=True,exist_ok=False)
    (out/"run.json").write_text(json.dumps(stats,indent=2)+"\n")
    fig,axes=plt.subplots(2,2,figsize=(10,9),sharex=True,sharey=True)
    for ax,(name,rows) in zip(axes.flat,sets.items()):
        t=np.array([float(x["true_T"]) for x in rows]); m=np.array([float(x["map_T"]) for x in rows])
        ax.scatter(t,m,s=13,alpha=.65); ax.plot([0,15000],[0,15000],"k--",lw=1); ax.set_title(name); ax.set(xlim=(0,15000),ylim=(0,15000),xlabel="True T",ylabel="MAP T")
    fig.tight_layout(); fig.savefig(out/"before_after_map_vs_true.png",dpi=220)

    bins=[(1,1),(2,2),(3,4),(5,8),(9,13),(14,18),(19,22),(23,25)]; freq=[]
    for model,tag in (("buggy estimated ARG","true_ne_frequency_alt_d0"),("fixed estimated ARG","true_ne_frequency_alt_d0_polarity_fixed")):
        for lo,hi in bins:
            errors=[]
            for rep in range(1,11):
                lab=f"replicate_{rep:03d}"; truth=json.loads((R/"simulations"/lab/"metadata.json").read_text())["true_ancient_ages"]
                base=R/"insertion_posterior"/tag/"merged"/lab; grid=np.load(base/"grid.npy"); ll=np.load(base/"ll_marginal_by_d0.npy")[lo:hi+1].sum(0)
                errors.extend(summarize(grid,x)[0]-t for x,t in zip(ll,truth))
            e=np.asarray(errors); freq.append({"model":model,"d0_lo":lo,"d0_hi":hi,"bias":float(e.mean()),"mae":float(abs(e).mean()),"rmse":float(np.sqrt(np.mean(e**2)))})
    (out/"frequency.json").write_text(json.dumps(freq,indent=2)+"\n")
    labels=[str(a) if a==b else f"{a}-{b}" for a,b in bins]; fig,ax=plt.subplots(figsize=(7,5)); ax.axhline(0,color="k",lw=1)
    for model in ("buggy estimated ARG","fixed estimated ARG"):
        z=[x for x in freq if x["model"]==model]; ax.plot(labels,[x["bias"] for x in z],"o-",label=model)
    ax.set(xlabel="Modern derived count",ylabel="Mean MAP bias (generations)"); ax.legend(frameon=False); fig.tight_layout(); fig.savefig(out/"frequency_bias_before_after.png",dpi=220)
    print(json.dumps(stats,indent=2))
if __name__=="__main__": main()
