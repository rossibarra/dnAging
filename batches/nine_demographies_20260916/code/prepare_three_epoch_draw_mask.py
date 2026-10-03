#!/usr/bin/env python3
"""Compute insertion-eligible sites for one posterior ARG draw."""
import argparse
from pathlib import Path
import numpy as np
import tszip
from betabinom_real_data import orientation
from insertion_real_data import chromosome_layout
from prepare_three_epoch_insertion import R, calls_for

def main():
    p=argparse.ArgumentParser(); p.add_argument("--draw",type=int,required=True); a=p.parse_args()
    calls={c:calls_for(c)[1] for c in ("chr1","chr2","chr3")}
    ts=tszip.decompress(R/"argtest"/"out"/"combined"/f"run.combined.{a.draw}.tsz")
    layout=chromosome_layout(ts); eligible={c:[] for c in calls}; tree=ts.first()
    ordered=sorted(layout.items(),key=lambda x:x[1][0])
    for site in ts.sites():
        if len(site.mutations)!=1: continue
        chrom=None
        for name,(offset,length) in ordered:
            if offset <= site.position < offset+length: chrom=name; break
        if chrom not in calls: continue
        local=int(site.position-layout[chrom][0]); rec=calls[chrom].get(local)
        if rec is None: continue
        mut=site.mutations[0]; tree.seek(site.position)
        if tree.parent(mut.node)==-1: continue
        is_alt,_=orientation(site,mut,rec[0],rec[1])
        if is_alt is not None: eligible[chrom].append(local)
    out=R/"insertion_three_chrom"/"mask_shards"/f"draw_{a.draw}"; out.mkdir(parents=True,exist_ok=False)
    for chrom,pos in eligible.items(): np.save(out/f"eligible_{chrom}.npy",np.asarray(pos,dtype=np.int64))
    print(a.draw,{c:len(x) for c,x in eligible.items()},flush=True)
if __name__=="__main__": main()
