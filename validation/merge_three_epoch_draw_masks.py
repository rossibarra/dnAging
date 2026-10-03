#!/usr/bin/env python3
"""Intersect draw masks and write final multichromosome call inputs."""
import json
from pathlib import Path
import numpy as np
from prepare_three_epoch_insertion import R, calls_for

def main():
    out=R/"insertion_three_chrom"/"prepared"; out.mkdir(parents=True,exist_ok=True)
    call_dir=out/"calls"; call_dir.mkdir(exist_ok=False); samples=None; summary={}
    for chrom in ("chr1","chr2","chr3"):
        names,calls=calls_for(chrom)
        if samples is None: samples=names
        elif samples!=names: raise ValueError("ancient sample order differs")
        common=None
        for draw in range(50,100):
            x=np.load(R/"insertion_three_chrom"/"mask_shards"/f"draw_{draw}"/f"eligible_{chrom}.npy")
            current=set(x.tolist()); common=current if common is None else common & current
        pos=np.asarray(sorted(common),dtype=np.int64)
        np.savez_compressed(call_dir/f"{chrom}.npz",position=pos,
          ref=np.asarray([calls[x][0] for x in pos]),alt=np.asarray([calls[x][1] for x in pos]),
          observed_alt=np.stack([calls[x][2] for x in pos]),
          called=np.ones((len(pos),len(samples)),dtype=np.uint8),samples=np.asarray(samples))
        summary[chrom]=len(pos)
    draw_dir=R/"argtest"/"out"/"combined"
    with (out/"draw_manifest.tsv").open("w") as h:
        h.write("draw_id\ttsz_path\n")
        for i in range(50,100): h.write(f"{i}\t{(draw_dir/f'run.combined.{i}.tsz').resolve()}\n")
    (out/"prepare.json").write_text(json.dumps({"draws":list(range(50,100)),"samples":samples,
      "common_sites":summary},indent=2)+"\n")
    print(summary,flush=True)
if __name__=="__main__": main()
