#!/usr/bin/env python3
"""Prepare one variable-Ne simulation for posterior-ARG insertion inference."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from pathlib import Path

import numpy as np
import tszip

import sys
sys.path[:0] = [str(Path(__file__).resolve().parents[1] / d) for d in ("pipeline", "validation")]
from betabinom_real_data import orientation
from insertion_real_data import chromosome_layout
from insertion_likelihood import PiecewiseConstantNe


DRAW_RE = re.compile(r"\.combined\.(\d+)\.tsz$")


def numbered_draws(path: Path):
    out=[]
    for tree in path.glob("*.tsz"):
        match=DRAW_RE.search(tree.name)
        if match: out.append((int(match.group(1)), tree.resolve()))
    return sorted(out)


def read_calls(vcf: Path, position_map: Path):
    mapped={}
    with position_map.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            mapped[row["site_id"]]=int(row["singer_pos"])
    samples=None; rows={}
    with gzip.open(vcf,"rt") as handle:
        for line in handle:
            if line.startswith("##"): continue
            fields=line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                names=fields[9:]
                keep=[i for i,n in enumerate(names) if n.startswith("ancient_")]
                samples=[names[i] for i in keep]
                continue
            if fields[2] not in mapped: continue
            gt_index=fields[8].split(":").index("GT")
            gt=[fields[9+i].split(":")[gt_index] for i in keep]
            if any(x not in {"0","1"} for x in gt):
                raise ValueError(f"non-haploid/missing ancient GT at site {fields[2]}")
            rows[mapped[fields[2]]]=(fields[3],fields[4],np.asarray(gt,dtype=np.uint8))
    if not samples: raise ValueError(f"no ancient samples in {vcf}")
    return samples,rows


def eligible_positions(tree_path: Path, calls, chrom: str):
    ts=tszip.decompress(tree_path); layout=chromosome_layout(ts)
    if set(layout) != {chrom}:
        raise ValueError(f"expected one chromosome {chrom}; found {sorted(layout)} in {tree_path}")
    offset,length=layout[chrom]; eligible=set(); positions=np.asarray(ts.tables.sites.position)
    for local,(ref,alt,_gt) in calls.items():
        global_pos=offset+local
        site_id=int(np.searchsorted(positions,global_pos))
        if site_id >= len(positions) or positions[site_id] != global_pos: continue
        site=ts.site(site_id)
        if len(site.mutations) != 1: continue
        mut=site.mutations[0]; tree=ts.at(global_pos)
        if tree.parent(mut.node) == -1: continue
        is_alt,_=orientation(site,mut,ref,alt)
        if is_alt is not None: eligible.add(local)
    return eligible


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=Path("msprime_variable_ne_error"))
    p.add_argument("--replicate",type=int,required=True)
    p.add_argument("--keep-last",type=int,default=50)
    p.add_argument("--output-root",type=Path,default=Path("msprime_variable_ne_error/insertion_posterior"))
    args=p.parse_args(); rep=f"replicate_{args.replicate:03d}"
    sim=args.root/"simulations"/rep; arg=args.root/"argtest"/rep
    draws=numbered_draws(arg/"out"/"combined")
    if len(draws) < args.keep_last:
        raise SystemExit(f"{rep}: only {len(draws)} draws; need {args.keep_last}")
    # A completed 0..99 chain retains 50..99. Selecting by numeric rank also
    # fails safely for missing internal draws because continuity is checked.
    selected=draws[-args.keep_last:]
    ids=[x[0] for x in selected]
    if ids != list(range(ids[0],ids[0]+args.keep_last)):
        raise SystemExit(f"{rep}: retained draw IDs are not contiguous: {ids}")
    ne=arg/"ne"/"estimated_ne.tsv"
    PiecewiseConstantNe.from_tsv(ne,"estimated")
    samples,calls=read_calls(
        sim/"all_samples.vcf.gz", args.root/"singer"/rep/"position_map.tsv")
    first=tszip.decompress(selected[0][1]); layout=chromosome_layout(first)
    if len(layout) != 1: raise SystemExit(f"{rep}: expected one chromosome, found {sorted(layout)}")
    chrom=next(iter(layout)); common=set(calls)
    for draw_id,path in selected:
        common &= eligible_positions(path,calls,chrom)
        print(f"{rep} draw {draw_id}: common eligible sites {len(common)}",flush=True)
    positions=np.asarray(sorted(common),dtype=np.int64)
    ref=np.asarray([calls[x][0] for x in positions]); alt=np.asarray([calls[x][1] for x in positions])
    observed=np.stack([calls[x][2] for x in positions])
    called=np.ones_like(observed,dtype=np.uint8)
    out=args.output_root/"prepared"/rep; call_dir=out/"calls"; call_dir.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(call_dir/f"chr{chrom}.npz",position=positions,ref=ref,alt=alt,
                        observed_alt=observed,called=called,samples=np.asarray(samples))
    with (out/"draw_manifest.tsv").open("w") as handle:
        handle.write("draw_id\ttsz_path\n")
        for draw_id,path in selected: handle.write(f"{draw_id}\t{path}\n")
    summary={"replicate":args.replicate,"draw_ids":ids,"burnin_discarded":len(draws)-len(selected),
             "n_draws":len(selected),"chromosome":chrom,"samples":samples,
             "candidate_sites":len(calls),"common_eligible_sites":len(common),"ne_table":str(ne.resolve())}
    (out/"prepare.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2),flush=True)


if __name__ == "__main__": main()
