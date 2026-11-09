#!/usr/bin/env python3
"""Controlled true-ARG recovery experiments for joint scale/epsilon calibration.

The input-time-scale argument emulates an ARG whose node and mutation-edge
times have already been multiplied by that value. Candidate correction scales
are then applied on top, so a compressed input scale of 0.8 should recover a
correction near 1.25. Candidate ages and the supplied Ne history remain in real
generations.
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
from insertion_real_data import parse_epsilon_grid, parse_time_scale_grid
from radiocarbon_oracle_stage1 import simulation_dirs, injected_calls, write_radiocarbon
import joint_scale_epsilon_calibration as calibration


def run(args):
    if not 0<args.input_time_scale or not np.isfinite(args.input_time_scale):
        raise ValueError("input-time-scale must be finite and positive")
    if not 0<=args.injected_epsilon<.5:
        raise ValueError("injected-epsilon must satisfy 0 <= epsilon < 0.5")
    scale_grid=parse_time_scale_grid(args.time_scale_grid)
    epsilon_grid=parse_epsilon_grid(args.epsilon_grid)
    age_grid=np.arange(args.age_min,args.age_max+args.age_step/2,
                       args.age_step,dtype=float)
    total=None; samples=None; true_ages=None; counts={}; flipped=0; allele_total=0
    for chromosome_index,sim in enumerate(simulation_dirs(args.run_root),start=1):
        metadata=json.loads((sim/"metadata.json").read_text())
        current_ages=np.asarray(metadata["true_ancient_ages"],dtype=float)
        if true_ages is None:
            true_ages=current_ages
        elif not np.array_equal(true_ages,current_ages):
            raise ValueError("ancient ages differ among chromosome replicates")
        ts=tskit.load(sim/"known_modern_arg.trees")
        names,calls,n_flipped,n_alleles=injected_calls(
            sim,args.n_modern,args.injected_epsilon,args.injection_seed,
            chromosome_index)
        flipped+=n_flipped; allele_total+=n_alleles
        if samples is None:
            samples=names
            total=np.zeros((len(scale_grid),len(epsilon_grid),len(samples),
                            len(age_grid)),dtype=float)
        elif samples!=names:
            raise ValueError("ancient sample order differs among chromosomes")
        ne=PiecewiseConstantNe.from_tsv(
            sim/"constant_ne_epochs.tsv",metadata["ne_table_series"])
        used=0; skipped=0
        effective_scales=scale_grid*args.input_time_scale
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
            p=np.asarray([derived_probability_uniform_edge_grid(
                tree,mutation.node,age_grid,tree.time(mutation.node),
                tree.time(parent),ne,time_scale=scale)
                for scale in effective_scales])
            eps=epsilon_grid[None,:,None]
            q=np.clip(eps+(1-2*eps)*p[:,None,:],1e-300,1-1e-15)
            observed=np.asarray(observed,dtype=float)
            total += observed[None,None,:,None]*np.log(q)[:,:,None,:]
            total += (1-observed)[None,None,:,None]*np.log1p(-q)[:,:,None,:]
            used+=1
        counts[f"chr{chromosome_index}"]={"used":used,"skipped":skipped}
        print(f"chr{chromosome_index}: {used} sites",flush=True)
    if len(samples)!=len(true_ages):
        raise ValueError("sample and true-age counts differ")

    if args.output.exists():
        raise FileExistsError(args.output)
    temporary=args.output.with_name(
        args.output.name+f".tmp.{os.environ.get('SLURM_JOB_ID','local')}"
    )
    temporary.mkdir(parents=True,exist_ok=False)
    likelihood=temporary/"likelihood"; likelihood.mkdir()
    np.save(likelihood/"grid.npy",age_grid)
    np.save(likelihood/"time_scale_grid.npy",scale_grid)
    np.save(likelihood/"epsilon_grid.npy",epsilon_grid)
    np.save(likelihood/"ll_marginal_scale_epsilon.npy",total)
    (likelihood/"samples.txt").write_text("\n".join(samples)+"\n")
    radiocarbon=temporary/"synthetic_radiocarbon.tsv"
    write_radiocarbon(radiocarbon,samples,true_ages,age_grid,args.rc_sd)
    calibration.run(Namespace(
        likelihood_dir=likelihood,radiocarbon=radiocarbon,folds=None,
        n_folds=args.n_folds,unknown_samples=None,age_prior=None,
        bootstrap_replicates=args.bootstrap_replicates,
        bootstrap_seed=args.bootstrap_seed,output=temporary/"calibration"))
    metadata={
        "stage":"controlled true-ARG joint scale/epsilon recovery",
        "run_root":str(args.run_root.resolve()),
        "scenario":args.scenario,
        "input_time_scale":args.input_time_scale,
        "expected_correction_scale":1/args.input_time_scale,
        "injected_epsilon":args.injected_epsilon,
        "injection_seed":args.injection_seed,
        "flipped_alleles":flipped,"alleles":allele_total,
        "realized_flip_rate":flipped/allele_total,
        "time_scale_grid":scale_grid.tolist(),
        "epsilon_grid":epsilon_grid.tolist(),
        "radiocarbon_sd":args.rc_sd,"counts_by_chromosome":counts,
        "excluded_from_radiocarbon":[s for s,t in zip(samples,true_ages) if t==0],
    }
    (temporary/"recovery_run.json").write_text(json.dumps(metadata,indent=2)+"\n")
    os.replace(temporary,args.output)


def summarize(args):
    rows=[]
    for result in sorted(args.results.glob("*/*")):
        recovery_path=result/"recovery_run.json"
        calibration_path=result/"calibration"/"run.json"
        if not recovery_path.exists() or not calibration_path.exists():
            continue
        recovery=json.loads(recovery_path.read_text())
        fitted=json.loads(calibration_path.read_text())
        joint=fitted["joint_calibration"]["joint_map"]
        metrics=fitted["cross_validation_metrics"]["joint_scale_epsilon"]
        run_name=Path(recovery["run_root"]).name
        model=run_name.rsplit("_r",1)[0]
        expected_scale=float(recovery["expected_correction_scale"])
        expected_epsilon=float(recovery["injected_epsilon"])
        scale_error=float(joint["time_scale"])-expected_scale
        epsilon_error=float(joint["epsilon"])-expected_epsilon
        rows.append({
            "run":run_name,"model":model,"scenario":recovery["scenario"],
            "expected_scale":expected_scale,"estimated_scale":joint["time_scale"],
            "scale_error":scale_error,"expected_epsilon":expected_epsilon,
            "estimated_epsilon":joint["epsilon"],"epsilon_error":epsilon_error,
            "mean_bias":metrics["mean_bias"],"mae":metrics["mae"],
            "rmse":metrics["rmse"],"regression_slope":metrics.get("regression_slope"),
            "scale_recovered":abs(scale_error)<=args.scale_tolerance,
            "epsilon_recovered":abs(epsilon_error)<=args.epsilon_tolerance,
        })
    if len(rows)!=args.expected_results:
        raise ValueError(f"found {len(rows)} complete results; expected {args.expected_results}")
    args.output.mkdir(parents=True,exist_ok=False)
    fields=list(rows[0])
    with (args.output/"recovery_summary.tsv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t",lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    failures=[row for row in rows if not (row["scale_recovered"] and row["epsilon_recovered"])]
    summary={"n_results":len(rows),"n_recovered":len(rows)-len(failures),
             "n_failed":len(failures),"scale_tolerance":args.scale_tolerance,
             "epsilon_tolerance":args.epsilon_tolerance,"failures":failures}
    (args.output/"run.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2),flush=True)
    if args.fail_on_unrecovered and failures:
        raise SystemExit(f"{len(failures)} recovery conditions failed tolerance gates")


def summarize_inferred(args):
    rows=[]
    for result in sorted(path for path in args.results.iterdir() if path.is_dir()):
        calibration_path=result/"calibration"/"run.json"
        if not calibration_path.exists():
            continue
        fitted=json.loads(calibration_path.read_text())
        joint=fitted["joint_calibration"]
        metrics=fitted["cross_validation_metrics"]["joint_scale_epsilon"]
        rows.append({
            "run":result.name,"model":result.name.rsplit("_r",1)[0],
            "estimated_scale":joint["joint_map"]["time_scale"],
            "estimated_epsilon":joint["joint_map"]["epsilon"],
            "scale_boundary_mass":joint["time_scale"]["boundary_mass"],
            "epsilon_boundary_mass":joint["epsilon"]["boundary_mass"],
            "mean_bias":metrics["mean_bias"],"mae":metrics["mae"],
            "rmse":metrics["rmse"],"regression_slope":metrics.get("regression_slope"),
        })
    if len(rows)!=args.expected_results:
        raise ValueError(f"found {len(rows)} inferred results; expected {args.expected_results}")
    args.output.mkdir(parents=True,exist_ok=False)
    with (args.output/"inferred_summary.tsv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter="\t",
                              lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    (args.output/"run.json").write_text(json.dumps(
        {"n_results":len(rows),"results":rows},indent=2)+"\n")


def parse_args(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest="command",required=True)
    one=sub.add_parser("run")
    one.add_argument("--run-root",type=Path,required=True)
    one.add_argument("--scenario",required=True)
    one.add_argument("--input-time-scale",type=float,required=True)
    one.add_argument("--injected-epsilon",type=float,required=True)
    one.add_argument("--injection-seed",type=int,required=True)
    one.add_argument("--time-scale-grid",default="0.7,0.8,0.9,1.0,1.1,1.2,1.25,1.3,1.4")
    one.add_argument("--epsilon-grid",default="0:0.06:0.01")
    one.add_argument("--n-modern",type=int,default=100)
    one.add_argument("--age-min",type=float,default=0)
    one.add_argument("--age-max",type=float,default=10000)
    one.add_argument("--age-step",type=float,default=100)
    one.add_argument("--rc-sd",type=float,default=200)
    one.add_argument("--n-folds",type=int,default=5)
    one.add_argument("--bootstrap-replicates",type=int,default=2000)
    one.add_argument("--bootstrap-seed",type=int,default=20261029)
    one.add_argument("--output",type=Path,required=True)
    one.set_defaults(func=run)
    combine=sub.add_parser("summarize")
    combine.add_argument("--results",type=Path,required=True)
    combine.add_argument("--expected-results",type=int,default=45)
    combine.add_argument("--scale-tolerance",type=float,default=.05)
    combine.add_argument("--epsilon-tolerance",type=float,default=.01)
    combine.add_argument("--fail-on-unrecovered",action="store_true")
    combine.add_argument("--output",type=Path,required=True)
    combine.set_defaults(func=summarize)
    inferred=sub.add_parser("summarize-inferred")
    inferred.add_argument("--results",type=Path,required=True)
    inferred.add_argument("--expected-results",type=int,default=9)
    inferred.add_argument("--output",type=Path,required=True)
    inferred.set_defaults(func=summarize_inferred)
    return parser.parse_args(argv)


def main():
    args=parse_args(); args.func(args)


if __name__=="__main__":
    main()
