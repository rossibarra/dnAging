#!/usr/bin/env python3
"""Diagnose age-dependent ARG time-scale preferences on completed joint grids.

This is a diagnostic, not a new dating model. It integrates each dated sample's
full radiocarbon density, reports individual preferred scales, and compares a
single global scale with separate young/middle/old scales under one shared
epsilon. The likelihood gain is a composite-likelihood diagnostic and is not
assigned a chi-square p-value.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import linregress, spearmanr

import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/d) for d in ("pipeline","validation")]

import joint_scale_epsilon_calibration as joint
import radiocarbon_calibration as base


STRATA=(
    ("young",-np.inf,3000.0),
    ("middle",3000.0,6500.0),
    ("old",6500.0,np.inf),
)


def stratum_name(age):
    for name,lower,upper in STRATA:
        if lower < age <= upper:
            return name
    raise AssertionError(age)


def fixed_epsilon_scale(profile, scale_grid, epsilon_grid, epsilon):
    matches=np.flatnonzero(np.isclose(epsilon_grid,epsilon,rtol=0,atol=1e-12))
    if len(matches)!=1:
        raise ValueError(f"fixed epsilon {epsilon:g} is absent from grid")
    values=np.asarray(profile[:,matches[0]],dtype=float)
    log_weight=values+base.normalized_log_weights(scale_grid)
    log_weight-=base.logsumexp(log_weight)
    mass=np.exp(log_weight)
    return {
        "map":float(scale_grid[int(np.argmax(log_weight))]),
        "mean":float(np.sum(scale_grid*mass)),
    }


def stratified_common_epsilon(profiles, scale_grid, epsilon_grid):
    """Compare raw composite maxima for one scale versus one per age stratum."""
    names=list(profiles)
    arrays=[np.asarray(profiles[name],dtype=float) for name in names]
    global_profile=np.sum(arrays,axis=0)
    global_index=np.unravel_index(np.argmax(global_profile),global_profile.shape)
    global_max=float(global_profile[global_index])
    stratified_by_epsilon=np.asarray([
        sum(float(np.max(profile[:,j])) for profile in arrays)
        for j in range(len(epsilon_grid))
    ])
    epsilon_index=int(np.argmax(stratified_by_epsilon))
    scale_indices=[int(np.argmax(profile[:,epsilon_index])) for profile in arrays]
    return {
        "global_scale":float(scale_grid[global_index[0]]),
        "global_epsilon":float(epsilon_grid[global_index[1]]),
        "stratified_epsilon":float(epsilon_grid[epsilon_index]),
        "stratified_scales":dict(zip(
            names,(float(scale_grid[i]) for i in scale_indices))),
        "composite_loglik_gain":float(stratified_by_epsilon[epsilon_index]-global_max),
    }


def trend_row(label, model, rows):
    age=np.asarray([row["radiocarbon_mean_age"] for row in rows],dtype=float)
    scale=np.asarray([row["scale_posterior_mean"] for row in rows],dtype=float)
    rho,pvalue=spearmanr(age,scale)
    fit=linregress(age/1000.0,scale)
    maps=np.asarray([row["joint_scale_map"] for row in rows],dtype=float)
    return {
        "group":label,"model":model,"n":len(rows),
        "spearman_rho":float(rho),"spearman_pvalue":float(pvalue),
        "linear_scale_change_per_1000_generations":float(fit.slope),
        "linear_r_squared":float(fit.rvalue**2),
        "map_at_scale_boundary_fraction":float(np.mean(
            (maps==rows[0]["scale_grid_min"])|(maps==rows[0]["scale_grid_max"]))),
        "n_unique_joint_scale_maps":int(len(np.unique(maps))),
    }


def write_rows(path, rows):
    with path.open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter="\t",
                              lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def run(args):
    fixed_epsilons=[float(x) for x in args.fixed_epsilons.split(",")]
    per_sample=[]; profile_records=[]
    result_dirs=sorted(path for path in args.results.iterdir() if path.is_dir())
    if len(result_dirs)!=args.expected_results:
        raise ValueError(f"found {len(result_dirs)} result directories; expected {args.expected_results}")
    for result in result_dirs:
        likelihood=result/"merged"
        calibration_meta=json.loads((result/"calibration"/"run.json").read_text())
        age_grid=np.load(likelihood/"grid.npy")
        scale_grid=np.load(likelihood/"time_scale_grid.npy")
        epsilon_grid=np.load(likelihood/"epsilon_grid.npy")
        ll=np.load(likelihood/"ll_marginal_scale_epsilon.npy",mmap_mode="r")
        samples=(likelihood/"samples.txt").read_text().split()
        sample_index={sample:i for i,sample in enumerate(samples)}
        rc=base.read_radiocarbon(Path(calibration_meta["radiocarbon_input"]))
        model=result.name.rsplit("_r",1)[0]
        for sample,(rc_age,rc_log_weight) in rc.items():
            profile=joint.joint_radiocarbon_log_profile(
                age_grid,ll[:,:,sample_index[sample],:],rc_age,rc_log_weight)
            summary=joint.joint_summary(profile,scale_grid,epsilon_grid)
            raw_index=np.unravel_index(np.argmax(profile),profile.shape)
            age=float(np.sum(rc_age*np.exp(rc_log_weight)))
            row={
                "run":result.name,"model":model,"sample":sample,
                "radiocarbon_mean_age":age,"age_stratum":stratum_name(age),
                "joint_scale_map":summary["joint_map"]["time_scale"],
                "joint_epsilon_map":summary["joint_map"]["epsilon"],
                "raw_scale_mle":float(scale_grid[raw_index[0]]),
                "raw_epsilon_mle":float(epsilon_grid[raw_index[1]]),
                "scale_posterior_mean":summary["time_scale"]["mean"],
                "scale_ci95_lower":summary["time_scale"]["ci95"][0],
                "scale_ci95_upper":summary["time_scale"]["ci95"][1],
                "scale_grid_min":float(scale_grid[0]),
                "scale_grid_max":float(scale_grid[-1]),
            }
            for epsilon in fixed_epsilons:
                fixed=fixed_epsilon_scale(profile,scale_grid,epsilon_grid,epsilon)
                suffix=f"{epsilon:g}".replace(".","p")
                row[f"scale_map_epsilon_{suffix}"]=fixed["map"]
                row[f"scale_mean_epsilon_{suffix}"]=fixed["mean"]
            per_sample.append(row)
            profile_records.append({"run":result.name,"model":model,
                                    "stratum":row["age_stratum"],"profile":profile,
                                    "scale_grid":scale_grid,"epsilon_grid":epsilon_grid})

    trend=[]
    for run in sorted({row["run"] for row in per_sample}):
        rows=[row for row in per_sample if row["run"]==run]
        trend.append(trend_row(run,rows[0]["model"],rows))
    for model in sorted({row["model"] for row in per_sample}):
        rows=[row for row in per_sample if row["model"]==model]
        trend.append(trend_row(f"{model}_pooled",model,rows))
    trend.append(trend_row("all_runs","all",per_sample))

    stratum_rows=[]; comparisons=[]
    groups=[(run,run.rsplit("_r",1)[0]) for run in sorted(
        {record["run"] for record in profile_records})]
    groups += [(f"{model}_pooled",model) for model in sorted(
        {record["model"] for record in profile_records})]
    groups += [("all_runs","all")]
    for label,model in groups:
        if label=="all_runs":
            records=profile_records
        elif label.endswith("_pooled"):
            records=[r for r in profile_records if r["model"]==model]
        else:
            records=[r for r in profile_records if r["run"]==label]
        scale_grid=records[0]["scale_grid"]; epsilon_grid=records[0]["epsilon_grid"]
        profiles={}
        for name,_,_ in STRATA:
            chosen=[r["profile"] for r in records if r["stratum"]==name]
            profiles[name]=np.sum(chosen,axis=0)
            summary=joint.joint_summary(profiles[name],scale_grid,epsilon_grid)
            raw_index=np.unravel_index(np.argmax(profiles[name]),profiles[name].shape)
            row={"group":label,"model":model,"age_stratum":name,
                 "n_samples":len(chosen),
                 "joint_scale_map":summary["joint_map"]["time_scale"],
                 "joint_epsilon_map":summary["joint_map"]["epsilon"],
                 "raw_scale_mle":float(scale_grid[raw_index[0]]),
                 "raw_epsilon_mle":float(epsilon_grid[raw_index[1]])}
            for epsilon in fixed_epsilons:
                fixed=fixed_epsilon_scale(profiles[name],scale_grid,epsilon_grid,epsilon)
                suffix=f"{epsilon:g}".replace(".","p")
                row[f"scale_map_epsilon_{suffix}"]=fixed["map"]
            stratum_rows.append(row)
        comparison=stratified_common_epsilon(profiles,scale_grid,epsilon_grid)
        comparisons.append({"group":label,"model":model,
                            "global_scale":comparison["global_scale"],
                            "global_epsilon":comparison["global_epsilon"],
                            "stratified_epsilon":comparison["stratified_epsilon"],
                            **{f"{name}_scale":comparison["stratified_scales"][name]
                               for name,_,_ in STRATA},
                            "composite_loglik_gain":comparison["composite_loglik_gain"]})

    args.output.mkdir(parents=True,exist_ok=False)
    write_rows(args.output/"per_sample_scale_preferences.tsv",per_sample)
    write_rows(args.output/"age_trends.tsv",trend)
    write_rows(args.output/"age_stratum_profiles.tsv",stratum_rows)
    write_rows(args.output/"age_model_comparison.tsv",comparisons)
    metadata={"n_runs":len(result_dirs),"n_dated_profiles":len(per_sample),
              "age_strata":[{
                  "name":n,
                  "lower_exclusive":None if not np.isfinite(lo) else lo,
                  "upper_inclusive":None if not np.isfinite(hi) else hi,
              } for n,lo,hi in STRATA],
              "fixed_epsilons":fixed_epsilons,
              "warning":"composite log-likelihood gains are diagnostic, not chi-square calibrated"}
    (args.output/"run.json").write_text(json.dumps(metadata,indent=2)+"\n")


def parse_args(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results",type=Path,required=True)
    parser.add_argument("--expected-results",type=int,default=9)
    parser.add_argument("--fixed-epsilons",default="0,0.01,0.02")
    parser.add_argument("--output",type=Path,required=True)
    return parser.parse_args(argv)


def main():
    run(parse_args())


if __name__=="__main__":
    main()
