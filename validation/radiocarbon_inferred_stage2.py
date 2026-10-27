#!/usr/bin/env python3
"""Assemble inferred-ARG epsilon likelihoods and run nested block calibration."""
from __future__ import annotations

import argparse
import csv
import json
import re
from argparse import Namespace
from pathlib import Path

import numpy as np

import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/d) for d in ("pipeline","validation")]
from insertion_real_data import logmeanexp, marginalize_independent_chromosomes
import radiocarbon_calibration as calibration
from radiocarbon_oracle_stage1 import write_radiocarbon


EPS_RE=re.compile(r"epsilon_([0-9.]+)$")
DRAW_RE=re.compile(r"draw_(\d+)$")


def load_parts(parts_root, expected_epsilons, expected_draws):
    epsilon_dirs=[]
    for path in parts_root.glob("epsilon_*"):
        match=EPS_RE.match(path.name)
        if match: epsilon_dirs.append((float(match.group(1)),path))
    epsilon_dirs.sort()
    if len(epsilon_dirs)!=expected_epsilons:
        raise ValueError(f"expected {expected_epsilons} epsilon directories; found {len(epsilon_dirs)}")
    epsilon=np.asarray([item[0] for item in epsilon_dirs])
    all_ll=[]; all_blocks=[]; samples=None; grid=None; chromosomes=None
    for epsilon_value,epsilon_dir in epsilon_dirs:
        draws=[]
        for path in epsilon_dir.glob("draw_*"):
            match=DRAW_RE.match(path.name)
            if match: draws.append((int(match.group(1)),path))
        draws.sort()
        if len(draws)!=expected_draws:
            raise ValueError(f"epsilon {epsilon_value:g}: expected {expected_draws} draws; found {len(draws)}")
        draw_ll=[]; draw_blocks=[]; draw_ids=[]
        for draw_id,path in draws:
            names=(path/"samples.txt").read_text().split()
            current_grid=np.load(path/"grid.npy")
            metadata=json.loads((path/"run.json").read_text())
            if samples is None:
                samples=names; grid=current_grid; chromosomes=metadata["chromosomes"]
            if samples!=names or not np.array_equal(grid,current_grid):
                raise ValueError(f"sample order or age grid differs in {path}")
            if chromosomes!=metadata["chromosomes"]:
                raise ValueError(f"chromosomes differ in {path}")
            if not np.isclose(metadata["epsilon"],epsilon_value):
                raise ValueError(f"epsilon metadata mismatch in {path}")
            draw_ids.append(draw_id)
            draw_ll.append(np.load(path/"ll_by_chrom.npy"))
            draw_blocks.append(np.load(path/"ll_by_block.npy"))
        if epsilon_dirs.index((epsilon_value,epsilon_dir))==0:
            reference_draw_ids=draw_ids
        elif draw_ids!=reference_draw_ids:
            raise ValueError("draw IDs differ among epsilon candidates")
        all_ll.append(np.stack(draw_ll)); all_blocks.append(np.stack(draw_blocks))
    return epsilon,grid,samples,chromosomes,reference_draw_ids,np.stack(all_ll),np.stack(all_blocks)


def resampled_likelihood(block_ll,weights):
    """Return epsilon x sample x age after block sums and ARG-draw marginalization."""
    n_epsilon,n_draw,n_chrom,_n_block,n_sample,n_age=block_ll.shape
    result=np.zeros((n_epsilon,n_sample,n_age),dtype=float)
    for e in range(n_epsilon):
        for chrom in range(n_chrom):
            by_draw=np.tensordot(block_ll[e,:,chrom],weights[chrom],axes=([1],[0]))
            result[e]+=logmeanexp(by_draw,axis=0)
    return result


def block_bootstrap(block_ll,epsilon_grid,age_grid,samples,rc,n_folds,n_bootstrap,seed):
    dated=[sample for sample in samples if sample in rc]
    sample_index={sample:i for i,sample in enumerate(samples)}
    fold_map=calibration.automatic_folds(dated,rc,n_folds)
    folds=sorted(set(fold_map.values()))
    n_chrom=block_ll.shape[2]; n_block=block_ll.shape[3]
    rng=np.random.default_rng(seed)
    epsilon_draws=np.empty(n_bootstrap)
    maps=np.empty((n_bootstrap,len(dated)))
    dated_index={sample:i for i,sample in enumerate(dated)}
    for bootstrap in range(n_bootstrap):
        weights=np.vstack([np.bincount(rng.integers(0,n_block,n_block),minlength=n_block)
                           for _ in range(n_chrom)])
        ll=resampled_likelihood(block_ll,weights)
        profiles={sample:calibration.radiocarbon_log_profile(
            age_grid,ll[:,sample_index[sample]],*rc[sample]) for sample in dated}
        all_profile=np.sum([profiles[sample] for sample in dated],axis=0)
        epsilon_draws[bootstrap]=epsilon_grid[int(np.argmax(all_profile))]
        for fold in folds:
            training=[sample for sample in dated if fold_map[sample]!=fold]
            held_out=[sample for sample in dated if fold_map[sample]==fold]
            profile=np.sum([profiles[sample] for sample in training],axis=0)
            best=int(np.argmax(profile))
            for sample in held_out:
                maps[bootstrap,dated_index[sample]]=age_grid[int(np.argmax(
                    ll[best,sample_index[sample]]))]
        if (bootstrap+1)%10==0:
            print(f"block bootstrap {bootstrap+1}/{n_bootstrap}",flush=True)
    return dated,fold_map,epsilon_draws,maps


def merge(args):
    epsilon,grid,samples,chromosomes,draw_ids,draw_ll,block_ll=load_parts(
        args.parts,args.expected_epsilons,args.expected_draws)
    posterior=np.stack([marginalize_independent_chromosomes(x) for x in draw_ll])
    metadata=json.loads((args.run_root/"simulations/replicate_001/metadata.json").read_text())
    true_ages=np.asarray(metadata["true_ancient_ages"],dtype=float)
    if len(true_ages)!=len(samples): raise ValueError("true-age and sample counts differ")
    args.output.mkdir(parents=True,exist_ok=False)
    likelihood=args.output/"likelihood"; likelihood.mkdir()
    np.save(likelihood/"epsilon_grid.npy",epsilon); np.save(likelihood/"grid.npy",grid)
    np.save(likelihood/"ll_marginal_epsilon.npy",posterior)
    (likelihood/"samples.txt").write_text("\n".join(samples)+"\n")
    np.save(args.output/"ll_by_block_epsilon_draw.npy",block_ll)
    radiocarbon=args.output/"synthetic_radiocarbon.tsv"
    write_radiocarbon(radiocarbon,samples,true_ages,grid,args.rc_sd)
    calibration.run(Namespace(
        likelihood_dir=likelihood,radiocarbon=radiocarbon,folds=None,
        n_folds=args.n_folds,unknown_samples=None,age_prior=None,
        bootstrap_replicates=args.individual_bootstrap_replicates,
        bootstrap_seed=args.seed,output=args.output/"calibration"))
    rc=calibration.read_radiocarbon(radiocarbon)
    dated,fold_map,epsilon_boot,maps=block_bootstrap(
        block_ll,epsilon,grid,samples,rc,args.n_folds,args.block_bootstrap_replicates,args.seed)
    np.save(args.output/"block_bootstrap_epsilon.npy",epsilon_boot)
    np.save(args.output/"block_bootstrap_maps.npy",maps)
    true_by_sample={sample:float(age) for sample,age in zip(samples,true_ages)}
    lo,median,hi=np.quantile(maps,[.025,.5,.975],axis=0)
    with (args.output/"block_bootstrap_ages.tsv").open("w",newline="") as handle:
        writer=csv.writer(handle,delimiter="\t",lineterminator="\n")
        writer.writerow(["sample","true_T","bootstrap_median_map_T",
                         "bootstrap_ci95_lower_T","bootstrap_ci95_upper_T","contains_true_95"])
        for sample,a,b,c in zip(dated,lo,median,hi):
            truth=true_by_sample[sample]
            writer.writerow([sample,f"{truth:.12g}",f"{b:.12g}",f"{a:.12g}",f"{c:.12g}",int(a<=truth<=c)])
    result={"stage":"estimated SINGER topology + POLEGON dates, true Ne, error-free calls",
            "draw_ids":draw_ids,"chromosomes":chromosomes,
            "epsilon_grid":epsilon.tolist(),"radiocarbon_sd":args.rc_sd,
            "block_size":5_000_000,"block_bootstrap_replicates":args.block_bootstrap_replicates,
            "block_bootstrap_epsilon_median":float(np.median(epsilon_boot)),
            "block_bootstrap_epsilon_ci95":[float(x) for x in np.quantile(epsilon_boot,[.025,.975])],
            "block_bootstrap_age_coverage":float(np.mean((lo<=true_ages[1:])&(true_ages[1:]<=hi))),
            "folds":fold_map}
    (args.output/"stage2_run.json").write_text(json.dumps(result,indent=2)+"\n")


def parse_args(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parts",type=Path,required=True)
    parser.add_argument("--run-root",type=Path,required=True)
    parser.add_argument("--expected-epsilons",type=int,default=11)
    parser.add_argument("--expected-draws",type=int,default=10)
    parser.add_argument("--rc-sd",type=float,default=200)
    parser.add_argument("--n-folds",type=int,default=5)
    parser.add_argument("--individual-bootstrap-replicates",type=int,default=2000)
    parser.add_argument("--block-bootstrap-replicates",type=int,default=200)
    parser.add_argument("--seed",type=int,default=20261012)
    parser.add_argument("--output",type=Path,required=True)
    return parser.parse_args(argv)


if __name__=="__main__": merge(parse_args())
