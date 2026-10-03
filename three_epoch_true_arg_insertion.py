#!/usr/bin/env python3
"""Insertion inference for three independent known modern ARGs."""
import json
from pathlib import Path
import numpy as np
import tskit
from betabinom_real_data import summarize
from insertion_likelihood import PiecewiseConstantNe, derived_probability_uniform_edge_grid
from msprime_insertion_validation import vcf_call_matrix
from insertion_real_data import block_bootstrap_maps

R=Path("msprime_three_epoch")
def main():
    grid=np.arange(0,10000+25,50,dtype=float); total=None; samples=None; counts={}; blocks=[]
    ne=PiecewiseConstantNe.from_tsv(R/"demography"/"three_epoch.tsv","three_epoch")
    for rep in range(1,4):
        sim=R/"simulations"/f"replicate_{rep:03d}"; ts=tskit.load(sim/"known_modern_arg.trees")
        names,calls=vcf_call_matrix(sim/"all_samples.vcf.gz",100)
        if samples is None: samples=names; total=np.zeros((len(samples),len(grid)))
        elif names != samples: raise ValueError("sample order differs")
        used=0; block_ll=np.zeros((10,len(samples),len(grid)))
        for site in ts.sites():
            if len(site.mutations)!=1 or site.id not in calls: continue
            mut=site.mutations[0]; tree=ts.at(site.position); parent=tree.parent(mut.node)
            if parent==tskit.NULL: continue
            d0,gt=calls[site.id]
            if tree.num_samples(mut.node)!=d0: continue
            q=derived_probability_uniform_edge_grid(tree,mut.node,grid,tree.time(mut.node),tree.time(parent),ne)
            q=np.clip(q,1e-300,1-1e-15)
            contribution=gt[:,None]*np.log(q)+(1-gt[:,None])*np.log1p(-q)
            total += contribution; block_ll[min(int(site.position//5_000_000),9)] += contribution; used+=1
        blocks.append(block_ll)
        counts[f"chr{rep}"]=used; print(rep,used,flush=True)
    out=R/"insertion_three_chrom"/"true_arg_true_ne"; out.mkdir(parents=True,exist_ok=False)
    np.save(out/"grid.npy",grid); np.save(out/"ll_marginal.npy",total)
    draw_blocks=np.asarray(blocks)[None,:,:,:,:]
    np.save(out/"ll_by_block_and_draw.npy",draw_blocks)
    boot=block_bootstrap_maps(draw_blocks,grid,1000,20260914); np.save(out/"block_bootstrap_maps.npy",boot)
    lo,hi=np.quantile(boot,[.025,.975],axis=0)
    with (out/"block_bootstrap_ci.tsv").open("w") as h:
        h.write("sample\tbootstrap_map_ci95_lower\tbootstrap_map_ci95_upper\n")
        for s,a,b in zip(samples,lo,hi): h.write(f"{s}\t{a:.6g}\t{b:.6g}\n")
    (out/"samples.txt").write_text("\n".join(samples)+"\n")
    with (out/"ages_table.tsv").open("w") as h:
        h.write("sample\tmap_T\tmean_T\tmedian_T\tci95_lower_T\tci95_upper_T\n")
        for s,ll in zip(samples,total): h.write(s+"\t"+"\t".join(f"{x:.6g}" for x in summarize(grid,ll))+"\n")
    (out/"run.json").write_text(json.dumps({"epsilon":0,"sites_by_chromosome":counts},indent=2)+"\n")
if __name__=="__main__": main()
