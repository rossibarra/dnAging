#!/usr/bin/env python3
"""Cross-validate insertion epsilon on dated samples and date undated samples.

The likelihood directory must be a merged ``insertion_real_data.py`` result
created with ``run-draw --epsilon-grid``.  Radiocarbon input is a long TSV with
columns ``sample``, ``age_generations``, and ``density``.  Supplying the full
calibrated density is preferred; a one-row point mass is accepted for tests and
special cases but should not replace a calibrated radiocarbon distribution.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def logsumexp(values, axis=None):
    values=np.asarray(values,dtype=float)
    peak=np.max(values,axis=axis,keepdims=True)
    out=np.squeeze(peak,axis=axis)+np.log(np.sum(np.exp(values-peak),axis=axis))
    return out


def quadrature_weights(grid):
    grid=np.asarray(grid,dtype=float)
    if grid.ndim!=1 or len(grid)==0 or np.any(~np.isfinite(grid)):
        raise ValueError("integration grid must be a nonempty finite vector")
    if len(grid)==1:
        return np.ones(1)
    delta=np.diff(grid)
    if np.any(delta<=0):
        raise ValueError("integration grid must be strictly increasing")
    weights=np.empty(len(grid),dtype=float)
    weights[0]=delta[0]/2; weights[-1]=delta[-1]/2
    weights[1:-1]=(delta[:-1]+delta[1:])/2
    return weights


def normalized_log_weights(grid, density=None):
    weights=quadrature_weights(grid)
    if density is not None:
        density=np.asarray(density,dtype=float)
        if density.shape!=weights.shape or np.any(~np.isfinite(density)) or np.any(density<0):
            raise ValueError("density must be finite, nonnegative, and match its grid")
        weights*=density
    total=weights.sum()
    if not total>0:
        raise ValueError("density has zero integrated mass")
    result=np.full(len(weights),-np.inf)
    positive=weights>0
    result[positive]=np.log(weights[positive]/total)
    return result


def read_radiocarbon(path):
    profiles={}
    with path.open(newline="") as handle:
        reader=csv.DictReader(handle,delimiter="\t")
        required={"sample","age_generations","density"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"{path} must contain columns {sorted(required)}")
        for row in reader:
            sample=row["sample"].strip()
            if not sample:
                raise ValueError("radiocarbon TSV contains an empty sample ID")
            profiles.setdefault(sample,[]).append(
                (float(row["age_generations"]),float(row["density"])))
    if not profiles:
        raise ValueError("radiocarbon TSV contains no samples")
    result={}
    for sample,rows in profiles.items():
        rows=sorted(rows)
        age=np.asarray([x[0] for x in rows]); density=np.asarray([x[1] for x in rows])
        if len(np.unique(age))!=len(age):
            raise ValueError(f"radiocarbon ages are duplicated for {sample}")
        log_weight=normalized_log_weights(age,density)
        result[sample]=(age,log_weight)
    return result


def radiocarbon_log_profile(age_grid, ll, rc_age, rc_log_weight):
    """Integrate one sample's (epsilon, age) likelihood over its RC density."""
    if rc_age[0]<age_grid[0] or rc_age[-1]>age_grid[-1]:
        raise ValueError("radiocarbon support falls outside the genetic age grid")
    interpolated=np.vstack([np.interp(rc_age,age_grid,row) for row in ll])
    return logsumexp(interpolated+rc_log_weight[None,:],axis=1)


def epsilon_log_posterior(profile, epsilon_grid):
    """Normalized composite profile; useful diagnostically, not a valid LD-aware CI."""
    log_prior=normalized_log_weights(epsilon_grid)
    joint=np.asarray(profile)+log_prior
    return joint-logsumexp(joint)


def marginal_age_loglik(ll, epsilon_profile, epsilon_grid):
    log_eps=epsilon_log_posterior(epsilon_profile,epsilon_grid)
    return logsumexp(ll+log_eps[:,None],axis=0)


def marginal_age_loglik_weights(ll, log_epsilon_weights):
    """Average an age likelihood over externally estimated epsilon weights."""
    return logsumexp(ll+np.asarray(log_epsilon_weights)[:,None],axis=0)


def bootstrap_epsilon(sample_profiles, epsilon_grid, n_bootstrap, seed):
    """Resample dated individuals and return MLE draws plus discrete log weights."""
    profiles=np.asarray(sample_profiles,dtype=float)
    if profiles.ndim!=2 or len(profiles)==0:
        raise ValueError("epsilon bootstrap requires one or more sample profiles")
    if n_bootstrap<1: raise ValueError("bootstrap-replicates must be positive")
    rng=np.random.default_rng(seed)
    draws=np.empty(n_bootstrap,dtype=float); counts=np.zeros(len(epsilon_grid),dtype=np.int64)
    for b in range(n_bootstrap):
        selected=rng.integers(0,len(profiles),len(profiles))
        best=int(np.argmax(np.sum(profiles[selected],axis=0)))
        draws[b]=epsilon_grid[best]; counts[best]+=1
    log_weight=np.full(len(epsilon_grid),-np.inf)
    positive=counts>0
    log_weight[positive]=np.log(counts[positive]/n_bootstrap)
    return draws,log_weight


def summarize_density(grid, log_likelihood, log_prior=None):
    lp=np.asarray(log_likelihood,dtype=float).copy()
    if log_prior is not None: lp+=log_prior
    peak=np.max(lp); density=np.exp(lp-peak)
    weight=quadrature_weights(grid)
    density/=np.sum(density*weight)
    mass=density*weight; cdf=np.cumsum(mass); cdf/=cdf[-1]
    quantile=lambda p: float(np.interp(p,cdf,grid))
    return {"map_T":float(grid[int(np.argmax(lp))]),
            "mean_T":float(np.sum(grid*mass)),
            "median_T":quantile(.5),"ci95_lower_T":quantile(.025),
            "ci95_upper_T":quantile(.975)}


def weighted_quantile(grid, log_weight, probability):
    mass=np.exp(log_weight); cdf=np.cumsum(mass); cdf/=cdf[-1]
    return float(np.interp(probability,cdf,grid))


def accuracy_metrics(rows, estimate_key="map_T"):
    error=np.asarray([float(row[estimate_key])-float(row["radiocarbon_mean_T"])
                      for row in rows])
    result={"n":len(rows),"mean_bias":float(np.mean(error)),
            "mae":float(np.mean(np.abs(error))),"rmse":float(np.sqrt(np.mean(error**2)))}
    if estimate_key=="map_T":
        result["map_inside_radiocarbon_95_fraction"]=float(np.mean(
            [row["map_in_radiocarbon_95"] for row in rows]))
        result["radiocarbon_mean_inside_genetic_95_fraction"]=float(np.mean(
            [row["radiocarbon_mean_in_genetic_95"] for row in rows]))
        if len(rows)>=2:
            rc=np.asarray([row["radiocarbon_mean_T"] for row in rows],dtype=float)
            estimate=np.asarray([row[estimate_key] for row in rows],dtype=float)
            if np.ptp(rc)>0:
                slope,intercept=np.polyfit(rc,estimate,1)
                result.update(regression_slope=float(slope),regression_intercept=float(intercept))
    return result


def epsilon_summary(profile, epsilon_grid, bootstrap_draws):
    return {"mle":float(epsilon_grid[int(np.argmax(profile))]),
            "sample_bootstrap_mean":float(np.mean(bootstrap_draws)),
            "sample_bootstrap_ci95":[float(x) for x in np.quantile(
                bootstrap_draws,[.025,.975])],
            "bootstrap_boundary_fraction":float(np.mean(
                (bootstrap_draws==epsilon_grid[0])|(bootstrap_draws==epsilon_grid[-1]))),
            "at_boundary":bool(np.argmax(profile) in (0,len(profile)-1))}


def interpolate_epsilon(ll, epsilon_grid, epsilon):
    if not epsilon_grid[0]<=epsilon<=epsilon_grid[-1]:
        return None
    hi=int(np.searchsorted(epsilon_grid,epsilon,side="right"))
    if hi==0: return ll[0]
    if hi==len(epsilon_grid): return ll[-1]
    lo=hi-1
    fraction=(epsilon-epsilon_grid[lo])/(epsilon_grid[hi]-epsilon_grid[lo])
    return (1-fraction)*ll[lo]+fraction*ll[hi]


def read_fold_map(path, dated):
    if path is None: return None
    with path.open(newline="") as handle:
        reader=csv.DictReader(handle,delimiter="\t")
        if not {"sample","fold"}.issubset(reader.fieldnames or []):
            raise ValueError("fold TSV must contain sample and fold columns")
        mapping={row["sample"]:row["fold"] for row in reader}
    missing=sorted(set(dated)-set(mapping))
    if missing: raise ValueError(f"fold TSV omits dated samples: {missing}")
    return {sample:mapping[sample] for sample in dated}


def automatic_folds(dated, rc_profiles, n_folds):
    if n_folds in (0,1) or n_folds>=len(dated):
        return {sample:f"loo_{i+1:03d}" for i,sample in enumerate(dated)}
    if n_folds<2: raise ValueError("n-folds must be 0 (LOO) or at least 2")
    central=[]
    for sample in dated:
        age,log_weight=rc_profiles[sample]
        central.append((float(np.sum(age*np.exp(log_weight))),sample))
    central.sort()
    return {sample:f"fold_{i%n_folds+1:02d}" for i,(_,sample) in enumerate(central)}


def write_profile(path, epsilon_grid, columns):
    with path.open("w",newline="") as handle:
        writer=csv.writer(handle,delimiter="\t",lineterminator="\n")
        writer.writerow(["epsilon",*columns])
        for i,epsilon in enumerate(epsilon_grid):
            writer.writerow([f"{epsilon:.12g}",*[f"{values[i]-np.max(values):.12g}" for values in columns.values()]])


def run(args):
    likelihood_dir=args.likelihood_dir
    age_grid=np.load(likelihood_dir/"grid.npy")
    epsilon_grid=np.load(likelihood_dir/"epsilon_grid.npy")
    ll=np.load(likelihood_dir/"ll_marginal_epsilon.npy",mmap_mode="r")
    samples=(likelihood_dir/"samples.txt").read_text().split()
    expected=(len(epsilon_grid),len(samples),len(age_grid))
    if ll.shape!=expected: raise ValueError(f"epsilon likelihood shape {ll.shape}; expected {expected}")
    if len(set(samples))!=len(samples): raise ValueError("likelihood sample IDs are not unique")
    rc=read_radiocarbon(args.radiocarbon)
    absent=sorted(set(rc)-set(samples))
    if absent: raise ValueError(f"dated samples absent from likelihood: {absent}")
    dated=[sample for sample in samples if sample in rc]
    if len(dated)<2: raise ValueError("at least two dated samples are required for cross-validation")
    index={sample:i for i,sample in enumerate(samples)}
    sample_profiles={sample:radiocarbon_log_profile(
        age_grid,ll[:,index[sample],:],*rc[sample]) for sample in dated}
    fold_map=read_fold_map(args.folds,dated)
    if fold_map is None: fold_map=automatic_folds(dated,rc,args.n_folds)
    folds=sorted(set(fold_map.values()))
    if len(folds)<2: raise ValueError("cross-validation requires at least two folds")

    age_log_prior=None
    if args.age_prior:
        prior=np.loadtxt(args.age_prior,skiprows=1)
        if prior.ndim!=2 or prior.shape[1]<2: raise ValueError("age prior needs age and density columns")
        density=np.interp(age_grid,prior[:,0],prior[:,1],left=0,right=0)
        age_log_prior=normalized_log_weights(age_grid,density)-np.log(quadrature_weights(age_grid))

    all_profile=np.sum([sample_profiles[s] for s in dated],axis=0)
    all_bootstrap,all_epsilon_log_weight=bootstrap_epsilon(
        [sample_profiles[s] for s in dated],epsilon_grid,
        args.bootstrap_replicates,args.bootstrap_seed)
    fold_profiles={}
    cv_rows=[]
    for fold in folds:
        training=[s for s in dated if fold_map[s]!=fold]
        held_out=[s for s in dated if fold_map[s]==fold]
        if not training or not held_out: raise ValueError(f"fold {fold} has empty training or validation set")
        profile=np.sum([sample_profiles[s] for s in training],axis=0)
        fold_profiles[fold]=profile
        _draws,fold_epsilon_log_weight=bootstrap_epsilon(
            [sample_profiles[s] for s in training],epsilon_grid,
            args.bootstrap_replicates,args.bootstrap_seed+len(fold_profiles))
        for sample in held_out:
            sample_ll=ll[:,index[sample],:]
            summary=summarize_density(age_grid,marginal_age_loglik_weights(
                sample_ll,fold_epsilon_log_weight),age_log_prior)
            rc_age,rc_log_weight=rc[sample]
            rc_lo=weighted_quantile(rc_age,rc_log_weight,.025)
            rc_hi=weighted_quantile(rc_age,rc_log_weight,.975)
            row={"sample":sample,"fold":fold,"n_training":len(training),
                 "radiocarbon_mean_T":float(np.sum(rc_age*np.exp(rc_log_weight))),
                 "radiocarbon_ci95_lower_T":rc_lo,"radiocarbon_ci95_upper_T":rc_hi,
                 **summary}
            row["signed_error_T"]=row["map_T"]-row["radiocarbon_mean_T"]
            row["abs_error_T"]=abs(row["signed_error_T"])
            row["map_in_radiocarbon_95"]=rc_lo<=row["map_T"]<=rc_hi
            row["radiocarbon_mean_in_genetic_95"]=(
                row["ci95_lower_T"]<=row["radiocarbon_mean_T"]<=row["ci95_upper_T"])
            for fixed in (.01,.05):
                fixed_ll=interpolate_epsilon(sample_ll,epsilon_grid,fixed)
                row[f"map_T_epsilon_{fixed:g}"]=("" if fixed_ll is None else
                    summarize_density(age_grid,fixed_ll,age_log_prior)["map_T"])
            cv_rows.append(row)

    if args.unknown_samples:
        unknown=args.unknown_samples.read_text().split()
        bad=sorted((set(unknown)-set(samples))|(set(unknown)&set(dated)))
        if bad: raise ValueError(f"unknown-sample list has absent or dated IDs: {bad}")
    else:
        unknown=[sample for sample in samples if sample not in rc]
    unknown_rows=[]
    for sample in unknown:
        summary=summarize_density(age_grid,marginal_age_loglik_weights(
            ll[:,index[sample],:],all_epsilon_log_weight),age_log_prior)
        unknown_rows.append({"sample":sample,**summary})

    args.output.mkdir(parents=True,exist_ok=False)
    write_profile(args.output/"epsilon_profiles.tsv",epsilon_grid,
                  {"all_relative_log_likelihood":all_profile,
                   **{f"{fold}_training_relative_log_likelihood":p
                      for fold,p in fold_profiles.items()}})
    fields=["sample","fold","n_training","radiocarbon_mean_T",
            "radiocarbon_ci95_lower_T","radiocarbon_ci95_upper_T","map_T","mean_T",
            "median_T","ci95_lower_T","ci95_upper_T","signed_error_T","abs_error_T",
            "map_in_radiocarbon_95","radiocarbon_mean_in_genetic_95",
            "map_T_epsilon_0.01","map_T_epsilon_0.05"]
    with (args.output/"cross_validated_ages.tsv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t",lineterminator="\n")
        writer.writeheader(); writer.writerows(cv_rows)
    fields_unknown=["sample","map_T","mean_T","median_T","ci95_lower_T","ci95_upper_T"]
    with (args.output/"unknown_ages.tsv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields_unknown,delimiter="\t",lineterminator="\n")
        writer.writeheader(); writer.writerows(unknown_rows)
    with (args.output/"epsilon_bootstrap.tsv").open("w",newline="") as handle:
        writer=csv.writer(handle,delimiter="\t",lineterminator="\n")
        writer.writerow(["epsilon","sample_bootstrap_probability_mass"])
        probability=np.exp(all_epsilon_log_weight)
        writer.writerows((f"{e:.12g}",f"{p:.12g}") for e,p in zip(epsilon_grid,probability))
    cv_metrics={"profiled_epsilon":accuracy_metrics(cv_rows)}
    for fixed in ("0.01","0.05"):
        key=f"map_T_epsilon_{fixed}"
        if all(row[key]!="" for row in cv_rows):
            cv_metrics[f"fixed_epsilon_{fixed}"]=accuracy_metrics(cv_rows,key)
    metadata={"format_version":1,"likelihood_dir":str(likelihood_dir.resolve()),
              "radiocarbon_input":str(args.radiocarbon.resolve()),
              "n_dated":len(dated),"n_unknown":len(unknown),"folds":fold_map,
              "epsilon":epsilon_summary(all_profile,epsilon_grid,all_bootstrap),
              "epsilon_uncertainty":{"method":"dated-individual bootstrap",
                                     "replicates":args.bootstrap_replicates,
                                     "seed":args.bootstrap_seed},
              "cross_validation_metrics":cv_metrics,
              "epsilon_interpretation":"effective per-allele discordance, not sequencing error",
              "radiocarbon_use":"training densities integrated; held-out densities excluded from genetic posteriors",
              "interval_warning":"genetic intervals propagate dated-individual-bootstrap epsilon uncertainty but are not linkage-calibrated block-bootstrap intervals",
              "age_prior":"uniform" if args.age_prior is None else str(args.age_prior.resolve())}
    (args.output/"run.json").write_text(json.dumps(metadata,indent=2)+"\n")


def parse_args(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--likelihood-dir",type=Path,required=True)
    parser.add_argument("--radiocarbon",type=Path,required=True,
                        help="long TSV: sample, age_generations, density")
    parser.add_argument("--folds",type=Path,help="optional TSV: sample, fold")
    parser.add_argument("--n-folds",type=int,default=0,
                        help="0 uses leave-one-out; otherwise age-stratified K-fold")
    parser.add_argument("--unknown-samples",type=Path,
                        help="optional one-ID-per-line subset; default is every undated sample")
    parser.add_argument("--age-prior",type=Path,
                        help="optional TSV with header then age_generations, density")
    parser.add_argument("--bootstrap-replicates",type=int,default=2000,
                        help="dated-individual bootstrap replicates [2000]")
    parser.add_argument("--bootstrap-seed",type=int,default=20261011)
    parser.add_argument("--output",type=Path,required=True)
    return parser.parse_args(argv)


def main():
    run(parse_args())


if __name__=="__main__": main()
