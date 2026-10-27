#!/usr/bin/env python3
"""Stage 1 oracle test for radiocarbon-calibrated insertion epsilon.

Each likelihood task applies one reproducible set of symmetric flips to the
error-free ancient calls, then evaluates one candidate epsilon on the known
modern ARG and true demography across three chromosomes.  The merge command
stacks candidate epsilons, creates synthetic radiocarbon densities, and runs
held-out calibration through the production calibration code.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from argparse import Namespace
from pathlib import Path

import numpy as np
import tskit

import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/d) for d in ("pipeline","validation")]
from insertion_likelihood import PiecewiseConstantNe, derived_probability_uniform_edge_grid
from msprime_insertion_validation import vcf_call_matrix
import radiocarbon_calibration as calibration


def simulation_dirs(run_root):
    paths=sorted((run_root/"simulations").glob("replicate_*"))
    if len(paths)!=3:
        raise ValueError(f"expected three chromosome replicates under {run_root}; found {len(paths)}")
    return paths


def injected_calls(sim, n_modern, injected_epsilon, injection_seed, chromosome_index):
    names,calls=vcf_call_matrix(sim/"all_samples.vcf.gz",n_modern)
    site_ids=sorted(calls)
    genotype=np.stack([calls[site][1] for site in site_ids])
    rng=np.random.default_rng(injection_seed+chromosome_index)
    flips=rng.random(genotype.shape)<injected_epsilon
    observed=np.where(flips,1-genotype,genotype).astype(np.int8)
    result={site:(calls[site][0],observed[i]) for i,site in enumerate(site_ids)}
    return names,result,int(flips.sum()),int(flips.size)


def likelihood(args):
    if not 0<=args.injected_epsilon<.5 or not 0<=args.fit_epsilon<.5:
        raise ValueError("epsilons must satisfy 0 <= epsilon < 0.5")
    grid=np.arange(args.age_min,args.age_max+args.age_step/2,args.age_step,dtype=float)
    total=None; samples=None; true_ages=None; counts={}; flipped=0; allele_total=0
    for chromosome_index,sim in enumerate(simulation_dirs(args.run_root),start=1):
        metadata=json.loads((sim/"metadata.json").read_text())
        if true_ages is None: true_ages=np.asarray(metadata["true_ancient_ages"],dtype=float)
        elif not np.array_equal(true_ages,np.asarray(metadata["true_ancient_ages"],dtype=float)):
            raise ValueError("ancient ages differ among chromosome replicates")
        ts=tskit.load(sim/"known_modern_arg.trees")
        names,calls,n_flipped,n_alleles=injected_calls(
            sim,args.n_modern,args.injected_epsilon,args.injection_seed,chromosome_index)
        flipped+=n_flipped; allele_total+=n_alleles
        if samples is None:
            samples=names; total=np.zeros((len(samples),len(grid)),dtype=float)
        elif samples!=names: raise ValueError("ancient sample order differs among chromosomes")
        ne=PiecewiseConstantNe.from_tsv(
            sim/"constant_ne_epochs.tsv",metadata["ne_table_series"])
        used=0; skipped=0
        for site in ts.sites():
            if len(site.mutations)!=1 or site.id not in calls:
                skipped+=1; continue
            mutation=site.mutations[0]; tree=ts.at(site.position)
            d0,observed=calls[site.id]
            if tree.num_samples(mutation.node)!=d0:
                skipped+=1; continue
            parent=tree.parent(mutation.node)
            if parent==tskit.NULL:
                skipped+=1; continue
            p=derived_probability_uniform_edge_grid(
                tree,mutation.node,grid,tree.time(mutation.node),tree.time(parent),ne)
            q=np.clip(args.fit_epsilon+(1-2*args.fit_epsilon)*p,1e-300,1-1e-15)
            total += observed[:,None]*np.log(q)[None,:]
            total += (1-observed)[:,None]*np.log1p(-q)[None,:]
            used+=1
        counts[f"chr{chromosome_index}"]={"used":used,"skipped":skipped}
        print(f"chr{chromosome_index}: {used} sites",flush=True)
    if len(samples)!=len(true_ages): raise ValueError("sample and true-age counts differ")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    temporary=args.output.with_name(args.output.name+f".tmp.{os.environ.get('SLURM_JOB_ID','local')}")
    with temporary.open("wb") as handle:
        np.savez_compressed(handle,fit_epsilon=args.fit_epsilon,Tgrid=grid,
                            samples=np.asarray(samples),true_ages=true_ages,
                            log_likelihood=total,injected_epsilon=args.injected_epsilon,
                            injection_seed=args.injection_seed,flipped=flipped,
                            allele_total=allele_total,
                            counts_json=json.dumps(counts,sort_keys=True))
    os.replace(temporary,args.output)


def write_radiocarbon(path,samples,true_ages,grid,sd):
    if sd<=0: raise ValueError("radiocarbon SD must be positive")
    with path.open("w",newline="") as handle:
        writer=csv.writer(handle,delimiter="\t",lineterminator="\n")
        writer.writerow(["sample","age_generations","density"])
        for sample,true_age in zip(samples,true_ages):
            if true_age==0: continue
            density=np.exp(-.5*((grid-true_age)/sd)**2)
            for age,value in zip(grid,density):
                writer.writerow([sample,f"{age:.12g}",f"{value:.17g}"])


def merge(args):
    paths=sorted(args.parts.glob("epsilon_*.npz"))
    if not paths: raise ValueError(f"no epsilon parts in {args.parts}")
    records=[]
    for path in paths:
        with np.load(path,allow_pickle=False) as part:
            records.append({key:part[key].copy() for key in part.files})
    records.sort(key=lambda part:float(part["fit_epsilon"]))
    epsilon=np.asarray([float(part["fit_epsilon"]) for part in records])
    if np.any(np.diff(epsilon)<=0): raise ValueError("candidate epsilons are duplicated")
    reference=records[0]
    for part in records[1:]:
        for key in ("Tgrid","samples","true_ages"):
            if not np.array_equal(part[key],reference[key]):
                raise ValueError(f"{key} differs among epsilon parts")
        for key in ("injected_epsilon","injection_seed","flipped","allele_total"):
            if part[key]!=reference[key]: raise ValueError(f"{key} differs among epsilon parts")
    args.output.mkdir(parents=True,exist_ok=False)
    likelihood_dir=args.output/"likelihood"; likelihood_dir.mkdir()
    np.save(likelihood_dir/"grid.npy",reference["Tgrid"])
    np.save(likelihood_dir/"epsilon_grid.npy",epsilon)
    np.save(likelihood_dir/"ll_marginal_epsilon.npy",
            np.stack([part["log_likelihood"] for part in records]))
    samples=reference["samples"].astype(str).tolist()
    (likelihood_dir/"samples.txt").write_text("\n".join(samples)+"\n")
    radiocarbon=args.output/"synthetic_radiocarbon.tsv"
    write_radiocarbon(radiocarbon,samples,reference["true_ages"],reference["Tgrid"],args.rc_sd)
    calibration.run(Namespace(
        likelihood_dir=likelihood_dir,radiocarbon=radiocarbon,folds=None,
        n_folds=args.n_folds,unknown_samples=None,age_prior=None,
        bootstrap_replicates=args.bootstrap_replicates,
        bootstrap_seed=args.bootstrap_seed,output=args.output/"calibration"))
    summary={"stage":"oracle radiocarbon epsilon recovery",
             "injected_epsilon":float(reference["injected_epsilon"]),
             "injection_seed":int(reference["injection_seed"]),
             "flipped_alleles":int(reference["flipped"]),
             "alleles":int(reference["allele_total"]),
             "realized_flip_rate":float(reference["flipped"])/int(reference["allele_total"]),
             "candidate_epsilons":epsilon.tolist(),"radiocarbon_sd":args.rc_sd,
             "excluded_from_radiocarbon":[samples[i] for i,t in enumerate(reference["true_ages"]) if t==0]}
    (args.output/"oracle_run.json").write_text(json.dumps(summary,indent=2)+"\n")


def parse_args(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest="command",required=True)
    one=sub.add_parser("likelihood")
    one.add_argument("--run-root",type=Path,required=True)
    one.add_argument("--injected-epsilon",type=float,required=True)
    one.add_argument("--fit-epsilon",type=float,required=True)
    one.add_argument("--injection-seed",type=int,default=2026101101)
    one.add_argument("--n-modern",type=int,default=100)
    one.add_argument("--age-min",type=float,default=0)
    one.add_argument("--age-max",type=float,default=10000)
    one.add_argument("--age-step",type=float,default=50)
    one.add_argument("--output",type=Path,required=True)
    one.set_defaults(func=likelihood)
    combine=sub.add_parser("merge")
    combine.add_argument("--parts",type=Path,required=True)
    combine.add_argument("--rc-sd",type=float,default=200)
    combine.add_argument("--n-folds",type=int,default=5)
    combine.add_argument("--bootstrap-replicates",type=int,default=2000)
    combine.add_argument("--bootstrap-seed",type=int,default=20261011)
    combine.add_argument("--output",type=Path,required=True)
    combine.set_defaults(func=merge)
    return parser.parse_args(argv)


def main():
    args=parse_args(); args.func(args)


if __name__=="__main__": main()
