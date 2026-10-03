#!/usr/bin/env python3
"""Plot within-chromosome distances between consecutive VCF records."""
import argparse, json
from array import array
from pathlib import Path
import numpy as np

def main(a):
    positions={}; sites=0
    with a.vcf.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("#"): continue
            fields=line.split("\t",2); chrom=fields[0]; pos=int(fields[1]); sites+=1
            positions.setdefault(chrom,array("q")).append(pos)
    pieces=[]; duplicates=0; by_chrom={}
    for chrom, values in positions.items():
        ordered=np.sort(np.frombuffer(values,dtype=np.int64)); by_chrom[chrom]=len(ordered)
        delta=np.diff(ordered); duplicates += int(np.sum(delta==0)); pieces.append(delta)
    d=np.concatenate(pieces)
    positive=d[d>0]
    summary={"vcf":str(a.vcf.resolve()),"sites":sites,"within_chromosome_intervals":len(d),
             "duplicate_position_intervals":duplicates,
             "chromosome_site_counts":by_chrom,
             "distance_bp":{"mean":float(d.mean()),"median":float(np.median(d)),
                 "q01_q05_q25_q75_q95_q99":[float(x) for x in np.quantile(d,[.01,.05,.25,.75,.95,.99])],
                 "maximum":int(d.max())}}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.summary.write_text(json.dumps(summary,indent=2)+"\n")
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.6))
    axes[0].hist(d[d<=1000],bins=np.arange(0,1010,10),color="#3268a8")
    axes[0].set(xlabel="Distance to preceding SNP (bp)",ylabel="Number of intervals",
                title="Distances up to 1 kb")
    bins=np.logspace(0,np.log10(max(positive.max(),2)),80)
    axes[1].hist(positive,bins=bins,color="#3268a8")
    axes[1].set_xscale("log"); axes[1].set_yscale("log")
    axes[1].set(xlabel="Distance to preceding SNP (bp; log scale)",ylabel="Number of intervals (log scale)",
                title="Full positive-distance distribution")
    fig.suptitle(f"Ancient VCF SNP spacing ({sites:,} sites; chromosome boundaries excluded)")
    fig.tight_layout(); fig.savefig(a.output,dpi=200); plt.close(fig)
    print(json.dumps(summary,indent=2))

if __name__=="__main__":
    p=argparse.ArgumentParser(); p.add_argument("--vcf",type=Path,required=True); p.add_argument("--output",type=Path,required=True); p.add_argument("--summary",type=Path,required=True); main(p.parse_args())
