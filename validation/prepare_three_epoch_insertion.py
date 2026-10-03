#!/usr/bin/env python3
"""Prepare common calls and a retained-draw manifest for the 3x50 Mb suite."""
import csv, gzip, json
from pathlib import Path
import numpy as np
import tszip
import sys
sys.path[:0] = [str(Path(__file__).resolve().parents[1] / d) for d in ("pipeline", "validation")]
from betabinom_real_data import orientation
from insertion_real_data import chromosome_layout

R=Path("msprime_three_epoch")

def calls_for(chrom):
    rep=int(chrom[3:]); mapped={}
    with (R/"singer"/"input"/f"{chrom}.position_map.tsv").open() as h:
        for x in csv.DictReader(h,delimiter="\t"): mapped[x["site_id"]]=int(x["singer_pos"])
    rows={}; names=None
    with gzip.open(R/"simulations"/f"replicate_{rep:03d}"/"all_samples.vcf.gz","rt") as h:
        for line in h:
            if line.startswith("##"): continue
            f=line.rstrip().split("\t")
            if line.startswith("#CHROM"):
                allnames=f[9:]; keep=[i for i,n in enumerate(allnames) if n.startswith("ancient_")]
                names=[allnames[i] for i in keep]; continue
            pos=mapped.get(f[2])
            if pos is None: continue
            gi=f[8].split(":").index("GT"); gt=[f[9+i].split(":")[gi] for i in keep]
            if any(x not in {"0","1"} for x in gt): raise ValueError(f"bad GT at {chrom}:{pos}")
            rows[pos]=(f[3],f[4],np.asarray(gt,dtype=np.uint8))
    return names,rows

def main():
    out=R/"insertion_three_chrom"/"prepared"; out.mkdir(parents=True,exist_ok=False)
    calls={}; samples=None
    for chrom in ("chr1","chr2","chr3"):
        names,rows=calls_for(chrom)
        if samples is None: samples=names
        elif names != samples: raise ValueError("ancient sample order differs")
        calls[chrom]=rows
    common={c:set(v) for c,v in calls.items()}
    draw_dir=R/"argtest"/"out"/"combined"
    draws=[(i,draw_dir/f"run.combined.{i}.tsz") for i in range(50,100)]
    for draw,path in draws:
        ts=tszip.decompress(path); layout=chromosome_layout(ts)
        eligible={c:set() for c in calls}
        for site in ts.sites():
            if len(site.mutations)!=1: continue
            for chrom,(offset,length) in layout.items():
                if offset <= site.position < offset+length: break
            else: continue
            local=int(site.position-offset); record=calls.get(chrom,{}).get(local)
            if record is None: continue
            mut=site.mutations[0]; tree=ts.at(site.position)
            if tree.parent(mut.node)==-1: continue
            is_alt,_=orientation(site,mut,record[0],record[1])
            if is_alt is not None: eligible[chrom].add(local)
        for chrom in common: common[chrom] &= eligible[chrom]
        print(draw,{c:len(x) for c,x in common.items()},flush=True)
    call_dir=out/"calls"; call_dir.mkdir()
    for chrom in calls:
        pos=np.asarray(sorted(common[chrom]),dtype=np.int64); rec=calls[chrom]
        np.savez_compressed(call_dir/f"{chrom}.npz",position=pos,
          ref=np.asarray([rec[x][0] for x in pos]),alt=np.asarray([rec[x][1] for x in pos]),
          observed_alt=np.stack([rec[x][2] for x in pos]),
          called=np.ones((len(pos),len(samples)),dtype=np.uint8),samples=np.asarray(samples))
    with (out/"draw_manifest.tsv").open("w") as h:
        h.write("draw_id\ttsz_path\n")
        for i,p in draws: h.write(f"{i}\t{p.resolve()}\n")
    (out/"prepare.json").write_text(json.dumps({"draws":[i for i,_ in draws],"samples":samples,
      "common_sites":{c:len(x) for c,x in common.items()}},indent=2)+"\n")

if __name__=="__main__": main()
