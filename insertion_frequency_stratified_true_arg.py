#!/usr/bin/env python3
"""Insertion likelihood split by present derived count on each true modern ARG."""
import argparse, json
from pathlib import Path
import numpy as np
import tskit
from insertion_likelihood import PiecewiseConstantNe, derived_probability_uniform_edge_grid
from msprime_insertion_validation import vcf_call_matrix

def main():
    p=argparse.ArgumentParser(); p.add_argument("--root",type=Path,default=Path("msprime_variable_ne_error")); p.add_argument("--replicate",type=int,required=True); p.add_argument("--output",type=Path,required=True); p.add_argument("--epsilon",type=float,default=.01); a=p.parse_args()
    label=f"replicate_{a.replicate:03d}"; sim=a.root/"simulations"/label
    ts=tskit.load(sim/"known_modern_arg.trees"); names,calls=vcf_call_matrix(sim/"all_samples.vcf.gz",26)
    ne=PiecewiseConstantNe.from_tsv(sim/"constant_ne_epochs.tsv","variable")
    grid=np.arange(0,15000+25,50,dtype=float); ll=np.zeros((27,len(names),len(grid))); sites=np.zeros(27,dtype=int)
    for site in ts.sites():
        if len(site.mutations)!=1 or site.id not in calls: continue
        mut=site.mutations[0]; tree=ts.at(site.position); parent=tree.parent(mut.node)
        if parent == tskit.NULL: continue
        d0,genotype=calls[site.id]
        if tree.num_samples(mut.node) != d0 or not 1 <= d0 <= 25: continue
        q=derived_probability_uniform_edge_grid(tree,mut.node,grid,tree.time(mut.node),tree.time(parent),ne)
        q=np.clip(a.epsilon+(1-2*a.epsilon)*q,1e-300,1-1e-15)
        ll[d0] += genotype[:,None]*np.log(q)+(1-genotype[:,None])*np.log1p(-q); sites[d0]+=1
    a.output.mkdir(parents=True,exist_ok=False); np.save(a.output/"grid.npy",grid); np.save(a.output/"ll_by_d0.npy",ll); np.save(a.output/"sites_by_d0.npy",sites)
    (a.output/"run.json").write_text(json.dumps({"replicate":a.replicate,"samples":names,"sites_by_d0":sites.tolist(),"epsilon":a.epsilon},indent=2)+"\n")
if __name__=="__main__": main()
