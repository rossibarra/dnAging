#!/usr/bin/env python3
"""Compare insertion age bias across present-day derived-count bins."""
import csv,json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path[:0] = [str(Path(__file__).resolve().parents[1] / d) for d in ("pipeline", "validation")]
from betabinom_real_data import summarize

R=Path("msprime_variable_ne_error")
BINS=[(1,1),(2,2),(3,4),(5,8),(9,13),(14,18),(19,22),(23,25)]
rows=[]
for model in ("true_ARG","estimated_ARG"):
  for rep in range(1,11):
    label=f"replicate_{rep:03d}"; truth=json.loads((R/"simulations"/label/"metadata.json").read_text())["true_ancient_ages"]
    if model=="true_ARG":
      base=R/"insertion_frequency"/"true_arg"/label; grid=np.load(base/"grid.npy"); ll=np.load(base/"ll_by_d0.npy"); sites=np.load(base/"sites_by_d0.npy")
    else:
      base=R/"insertion_posterior"/"true_ne_frequency_alt_d0"/"merged"/label; grid=np.load(base/"grid.npy"); ll=np.load(base/"ll_marginal_by_d0.npy"); sites=np.load(base/"sites_by_d0.npy")
    for lo,hi in BINS:
      combined=ll[lo:hi+1].sum(axis=0); n=int(sites[lo:hi+1].sum())
      for j,(lp,t) in enumerate(zip(combined,truth),1):
        est=summarize(grid,lp); rows.append({"model":model,"replicate":rep,"sample":f"ancient_{j:02d}","d0_lo":lo,"d0_hi":hi,"sites":n,"true_T":t,"map_T":est[0],"bias":est[0]-t,"abs_error":abs(est[0]-t)})
out=R/"insertion_frequency"/"summary"; out.mkdir(parents=True,exist_ok=False)
with (out/"age_by_frequency.tsv").open("w",newline="") as h:
 w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(rows)
summary=[]
for model in ("true_ARG","estimated_ARG"):
 for lo,hi in BINS:
  z=[r for r in rows if r["model"]==model and r["d0_lo"]==lo]; e=np.array([r["bias"] for r in z]);
  summary.append({"model":model,"d0_lo":lo,"d0_hi":hi,"mean_sites_per_replicate":float(np.mean([r["sites"] for r in z[::10]])),"bias":float(e.mean()),"mae":float(abs(e).mean()),"rmse":float(np.sqrt(np.mean(e**2)))})
(out/"run.json").write_text(json.dumps(summary,indent=2)+"\n")
fig,axes=plt.subplots(1,2,figsize=(11,4.5),sharey=True); labels=[f"{a}" if a==b else f"{a}-{b}" for a,b in BINS]
for ax,model in zip(axes,("true_ARG","estimated_ARG")):
 z=[r for r in summary if r["model"]==model]; ax.axhline(0,color="k",lw=1); ax.plot(labels,[r["bias"] for r in z],"o-"); ax.set(title=model.replace("_"," "),xlabel="Present derived count d0",ylabel="Mean MAP bias (generations)"); ax.tick_params(axis="x",rotation=45)
fig.tight_layout(); fig.savefig(out/"bias_by_allele_frequency.png",dpi=220)
print(json.dumps(summary,indent=2))
