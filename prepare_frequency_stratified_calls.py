#!/usr/bin/env python3
"""Add fixed observed-modern allele counts to prepared insertion calls."""
import argparse,csv,gzip
from pathlib import Path
import numpy as np

def main():
 p=argparse.ArgumentParser(); p.add_argument("--root",type=Path,default=Path("msprime_variable_ne_error")); p.add_argument("--replicate",type=int,required=True); a=p.parse_args(); label=f"replicate_{a.replicate:03d}"
 src=a.root/"insertion_posterior"/"prepared"/label/"calls"/"chrchr1.npz"; data=np.load(src)
 wanted=set(data["position"].tolist()); by_id={}
 with (a.root/"singer"/label/"position_map.tsv").open() as h:
  for r in csv.DictReader(h,delimiter="\t"):
   p0=int(r["singer_pos"])
   if p0 in wanted: by_id[r["site_id"]]=p0
 counts={}
 with gzip.open(a.root/"simulations"/label/"all_samples.vcf.gz","rt") as h:
  for line in h:
   if line.startswith("#CHROM"):
    names=line.rstrip().split("\t")[9:]; modern=[i for i,n in enumerate(names) if n.startswith("modern_")]; continue
   if line.startswith("#"): continue
   f=line.rstrip().split("\t"); pos=by_id.get(f[2])
   if pos is not None:
    gt=[f[9+i].split(":",1)[0] for i in modern]; counts[pos]=(sum(x=="1" for x in gt),len(gt))
 pos=data["position"]; out=a.root/"insertion_posterior"/"prepared_frequency"/label/"calls"; out.mkdir(parents=True,exist_ok=False)
 np.savez_compressed(out/"chrchr1.npz",**{k:data[k] for k in data.files},
  panel_alt_count=np.array([counts[int(x)][0] for x in pos],dtype=np.uint8),panel_called=np.array([counts[int(x)][1] for x in pos],dtype=np.uint8))
 (out.parent/"draw_manifest.tsv").write_text((src.parent.parent/"draw_manifest.tsv").read_text())
if __name__=="__main__": main()
