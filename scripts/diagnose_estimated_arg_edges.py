#!/usr/bin/env python3
"""Compare true mutation times with inferred mutation-bearing edge intervals."""
import csv, gzip, json
from pathlib import Path
import numpy as np
import tszip

R=Path("msprime_variable_ne_error")
for rep in range(1,11):
    label=f"replicate_{rep:03d}"
    truth={}
    with gzip.open(R/"simulations"/label/"mutation_truth.tsv.gz","rt") as h:
        for row in csv.DictReader(h,delimiter="\t"): truth[row["site_id"]]=float(row["mutation_time"])
    mapped={}
    with (R/"singer"/label/"position_map.tsv").open() as h:
        for row in csv.DictReader(h,delimiter="\t"): mapped[int(row["singer_pos"])]=truth[row["site_id"]]
    common=set(np.load(R/"insertion_posterior"/"prepared"/label/"calls"/"chrchr1.npz")["position"].tolist())
    ts=tszip.decompress(R/"argtest"/label/"out"/"combined"/"run.combined.99.tsz")
    contained=[]; ratios=[]; lower_ratios=[]; upper_ratios=[]
    for site in ts.sites():
        pos=int(site.position)
        if pos not in common or pos not in mapped or len(site.mutations)!=1: continue
        mut=site.mutations[0]; tree=ts.at(site.position); parent=tree.parent(mut.node)
        if parent == -1: continue
        lo=tree.time(mut.node); hi=tree.time(parent); t=mapped[pos]
        contained.append(lo <= t <= hi); ratios.append(np.sqrt(max(lo,1e-9)*hi)/t)
        lower_ratios.append(lo/t); upper_ratios.append(hi/t)
    print(json.dumps({"replicate":rep,"sites":len(contained),"true_time_in_inferred_edge":float(np.mean(contained)),
        "median_geomean_edge_time_over_true":float(np.median(ratios)),
        "median_child_time_over_true":float(np.median(lower_ratios)),
        "median_parent_time_over_true":float(np.median(upper_ratios))}))
