#!/usr/bin/env python3
"""Cross-validate joint ARG time-scale and epsilon calibration.

The likelihood directory must contain ``ll_marginal_scale_epsilon.npy`` with
logical dimensions time scale x epsilon x sample x age.  Radiocarbon input is
the same long-form calibrated-density TSV used by ``radiocarbon_calibration``.
Held-out dates are used only for scoring, never for their fold's calibration
surface or genetic age posterior.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

import radiocarbon_calibration as base


def joint_radiocarbon_log_profile(age_grid, ll, rc_age, rc_log_weight):
    """Integrate one sample's (scale, epsilon, age) likelihood over RC age."""
    ll=np.asarray(ll,dtype=float)
    if ll.ndim!=3:
        raise ValueError("joint sample likelihood must have scale, epsilon, age axes")
    flat=ll.reshape((-1,ll.shape[-1]))
    profile=base.radiocarbon_log_profile(
        age_grid,flat,rc_age,rc_log_weight)
    return profile.reshape(ll.shape[:2])


def joint_log_posterior(profile, time_scale_grid, epsilon_grid):
    """Normalize a joint profile under uniform-density grid priors."""
    profile=np.asarray(profile,dtype=float)
    expected=(len(time_scale_grid),len(epsilon_grid))
    if profile.shape!=expected or np.any(~np.isfinite(profile)):
        raise ValueError(f"joint profile must be finite with shape {expected}")
    log_prior=(base.normalized_log_weights(time_scale_grid)[:,None]
               +base.normalized_log_weights(epsilon_grid)[None,:])
    posterior=profile+log_prior
    return posterior-base.logsumexp(posterior)


def marginal_age_loglik_weights(ll, log_joint_weights):
    """Average a joint likelihood over externally estimated calibration weights."""
    ll=np.asarray(ll,dtype=float)
    log_joint_weights=np.asarray(log_joint_weights,dtype=float)
    if ll.ndim!=3 or ll.shape[:2]!=log_joint_weights.shape:
        raise ValueError("joint likelihood and calibration weights have incompatible shapes")
    return base.logsumexp(ll+log_joint_weights[:,:,None],axis=(0,1))


def bootstrap_joint(sample_profiles, time_scale_grid, epsilon_grid,
                    n_bootstrap, seed):
    """Resample dated individuals and return joint posterior-mode draws."""
    profiles=np.asarray(sample_profiles,dtype=float)
    expected=(len(time_scale_grid),len(epsilon_grid))
    if profiles.ndim!=3 or len(profiles)==0 or profiles.shape[1:]!=expected:
        raise ValueError("joint bootstrap requires sample x scale x epsilon profiles")
    if n_bootstrap<1:
        raise ValueError("bootstrap-replicates must be positive")
    log_prior=(base.normalized_log_weights(time_scale_grid)[:,None]
               +base.normalized_log_weights(epsilon_grid)[None,:])
    rng=np.random.default_rng(seed)
    draws=np.empty((n_bootstrap,2),dtype=float)
    for b in range(n_bootstrap):
        selected=rng.integers(0,len(profiles),len(profiles))
        index=np.unravel_index(
            np.argmax(np.sum(profiles[selected],axis=0)+log_prior),expected)
        draws[b]=(time_scale_grid[index[0]],epsilon_grid[index[1]])
    return draws


def joint_summary(profile, time_scale_grid, epsilon_grid):
    """Summarize joint and marginal calibration behavior."""
    log_weight=joint_log_posterior(profile,time_scale_grid,epsilon_grid)
    mass=np.exp(log_weight)
    best=np.unravel_index(np.argmax(log_weight),log_weight.shape)
    scale_log=base.logsumexp(log_weight,axis=1)
    epsilon_log=base.logsumexp(log_weight,axis=0)
    scale_mean=float(np.sum(np.asarray(time_scale_grid)*np.exp(scale_log)))
    epsilon_mean=float(np.sum(np.asarray(epsilon_grid)*np.exp(epsilon_log)))
    scale_delta=np.asarray(time_scale_grid)[:,None]-scale_mean
    epsilon_delta=np.asarray(epsilon_grid)[None,:]-epsilon_mean
    scale_var=float(np.sum(mass*scale_delta**2))
    epsilon_var=float(np.sum(mass*epsilon_delta**2))
    correlation=None
    if scale_var>0 and epsilon_var>0:
        correlation=float(np.sum(mass*scale_delta*epsilon_delta)
                          /np.sqrt(scale_var*epsilon_var))
    edge=np.zeros(mass.shape,dtype=bool)
    edge[[0,-1],:]=True; edge[:,[0,-1]]=True
    scale_edges=np.unique([0,len(time_scale_grid)-1])
    epsilon_edges=np.unique([0,len(epsilon_grid)-1])
    return {
        "joint_map":{"time_scale":float(time_scale_grid[best[0]]),
                     "epsilon":float(epsilon_grid[best[1]])},
        "time_scale":{"mean":scale_mean,
                      "ci95":[base.weighted_quantile(time_scale_grid,scale_log,.025),
                              base.weighted_quantile(time_scale_grid,scale_log,.975)],
                      "boundary_mass":float(np.sum(mass[scale_edges,:])),
                      "map_at_boundary":bool(best[0] in (0,len(time_scale_grid)-1))},
        "epsilon":{"mean":epsilon_mean,
                   "ci95":[base.weighted_quantile(epsilon_grid,epsilon_log,.025),
                           base.weighted_quantile(epsilon_grid,epsilon_log,.975)],
                   "boundary_mass":float(np.sum(mass[:,epsilon_edges])),
                   "map_at_boundary":bool(best[1] in (0,len(epsilon_grid)-1))},
        "joint_boundary_mass":float(np.sum(mass[edge])),
        "parameter_correlation":correlation,
    }


def write_joint_profile(path, time_scale_grid, epsilon_grid, columns):
    with path.open("w",newline="") as handle:
        writer=csv.writer(handle,delimiter="\t",lineterminator="\n")
        writer.writerow(["time_scale","epsilon",*columns])
        maxima={name:np.max(values) for name,values in columns.items()}
        for i,time_scale in enumerate(time_scale_grid):
            for j,epsilon in enumerate(epsilon_grid):
                writer.writerow([
                    f"{time_scale:.12g}",f"{epsilon:.12g}",
                    *[f"{values[i,j]-maxima[name]:.12g}"
                      for name,values in columns.items()]])


def write_marginal(path, name, grid, log_weight):
    with path.open("w",newline="") as handle:
        writer=csv.writer(handle,delimiter="\t",lineterminator="\n")
        writer.writerow([name,"posterior_probability_mass"])
        writer.writerows((f"{x:.12g}",f"{p:.12g}")
                         for x,p in zip(grid,np.exp(log_weight)))


def run(args):
    likelihood_dir=args.likelihood_dir
    age_grid=np.load(likelihood_dir/"grid.npy")
    time_scale_grid=np.load(likelihood_dir/"time_scale_grid.npy")
    epsilon_grid=np.load(likelihood_dir/"epsilon_grid.npy")
    ll=np.load(likelihood_dir/"ll_marginal_scale_epsilon.npy",mmap_mode="r")
    samples=(likelihood_dir/"samples.txt").read_text().split()
    expected=(len(time_scale_grid),len(epsilon_grid),len(samples),len(age_grid))
    if ll.shape!=expected:
        raise ValueError(f"joint likelihood shape {ll.shape}; expected {expected}")
    if np.any(np.asarray(time_scale_grid)<=0):
        raise ValueError("time-scale grid values must be positive")
    if np.any((np.asarray(epsilon_grid)<0)|(np.asarray(epsilon_grid)>=.5)):
        raise ValueError("epsilon grid values must satisfy 0 <= epsilon < 0.5")
    if len(set(samples))!=len(samples):
        raise ValueError("likelihood sample IDs are not unique")

    rc_profiles=base.read_radiocarbon(args.radiocarbon)
    absent=sorted(set(rc_profiles)-set(samples))
    if absent:
        raise ValueError(f"dated samples absent from likelihood: {absent}")
    dated=[sample for sample in samples if sample in rc_profiles]
    if len(dated)<2:
        raise ValueError("at least two dated samples are required for cross-validation")
    index={sample:i for i,sample in enumerate(samples)}
    sample_profiles={sample:joint_radiocarbon_log_profile(
        age_grid,ll[:,:,index[sample],:],*rc_profiles[sample]) for sample in dated}
    fold_map=base.read_fold_map(args.folds,dated)
    if fold_map is None:
        fold_map=base.automatic_folds(dated,rc_profiles,args.n_folds)
    folds=sorted(set(fold_map.values()))
    if len(folds)<2:
        raise ValueError("cross-validation requires at least two folds")

    age_log_prior=None
    if args.age_prior:
        prior=np.loadtxt(args.age_prior,skiprows=1)
        if prior.ndim!=2 or prior.shape[1]<2:
            raise ValueError("age prior needs age and density columns")
        density=np.interp(age_grid,prior[:,0],prior[:,1],left=0,right=0)
        age_log_prior=(base.normalized_log_weights(age_grid,density)
                       -np.log(base.quadrature_weights(age_grid)))

    all_profile=np.sum([sample_profiles[s] for s in dated],axis=0)
    all_log_weight=joint_log_posterior(all_profile,time_scale_grid,epsilon_grid)
    all_bootstrap=bootstrap_joint(
        [sample_profiles[s] for s in dated],time_scale_grid,epsilon_grid,
        args.bootstrap_replicates,args.bootstrap_seed)
    fold_profiles={}; fold_summaries={}; cv_rows=[]
    for fold in folds:
        training=[s for s in dated if fold_map[s]!=fold]
        held_out=[s for s in dated if fold_map[s]==fold]
        if not training or not held_out:
            raise ValueError(f"fold {fold} has empty training or validation set")
        profile=np.sum([sample_profiles[s] for s in training],axis=0)
        fold_profiles[fold]=profile
        fold_summaries[fold]=joint_summary(profile,time_scale_grid,epsilon_grid)
        fold_log_weight=joint_log_posterior(profile,time_scale_grid,epsilon_grid)
        for sample in held_out:
            summary=base.summarize_density(
                age_grid,marginal_age_loglik_weights(
                    ll[:,:,index[sample],:],fold_log_weight),age_log_prior)
            rc_age,rc_log_weight=rc_profiles[sample]
            rc_lo=base.weighted_quantile(rc_age,rc_log_weight,.025)
            rc_hi=base.weighted_quantile(rc_age,rc_log_weight,.975)
            row={"sample":sample,"fold":fold,"n_training":len(training),
                 "fold_time_scale_map":fold_summaries[fold]["joint_map"]["time_scale"],
                 "fold_epsilon_map":fold_summaries[fold]["joint_map"]["epsilon"],
                 "radiocarbon_mean_T":float(np.sum(rc_age*np.exp(rc_log_weight))),
                 "radiocarbon_ci95_lower_T":rc_lo,
                 "radiocarbon_ci95_upper_T":rc_hi,**summary}
            row["signed_error_T"]=row["map_T"]-row["radiocarbon_mean_T"]
            row["abs_error_T"]=abs(row["signed_error_T"])
            row["map_in_radiocarbon_95"]=rc_lo<=row["map_T"]<=rc_hi
            row["radiocarbon_mean_in_genetic_95"]=(
                row["ci95_lower_T"]<=row["radiocarbon_mean_T"]<=row["ci95_upper_T"])
            cv_rows.append(row)

    if args.unknown_samples:
        unknown=args.unknown_samples.read_text().split()
        bad=sorted((set(unknown)-set(samples))|(set(unknown)&set(dated)))
        if bad:
            raise ValueError(f"unknown-sample list has absent or dated IDs: {bad}")
    else:
        unknown=[sample for sample in samples if sample not in rc_profiles]
    unknown_rows=[]
    for sample in unknown:
        summary=base.summarize_density(
            age_grid,marginal_age_loglik_weights(
                ll[:,:,index[sample],:],all_log_weight),age_log_prior)
        unknown_rows.append({"sample":sample,**summary})

    args.output.mkdir(parents=True,exist_ok=False)
    write_joint_profile(args.output/"joint_scale_epsilon_profile.tsv",
                        time_scale_grid,epsilon_grid,
                        {"all_relative_log_posterior":all_log_weight,
                         **{f"{fold}_training_relative_log_likelihood":profile
                            for fold,profile in fold_profiles.items()}})
    scale_log=base.logsumexp(all_log_weight,axis=1)
    epsilon_log=base.logsumexp(all_log_weight,axis=0)
    write_marginal(args.output/"scale_marginal.tsv","time_scale",
                   time_scale_grid,scale_log)
    write_marginal(args.output/"epsilon_marginal.tsv","epsilon",
                   epsilon_grid,epsilon_log)
    cv_fields=["sample","fold","n_training","fold_time_scale_map",
               "fold_epsilon_map","radiocarbon_mean_T",
               "radiocarbon_ci95_lower_T","radiocarbon_ci95_upper_T",
               "map_T","mean_T","median_T","ci95_lower_T","ci95_upper_T",
               "signed_error_T","abs_error_T","map_in_radiocarbon_95",
               "radiocarbon_mean_in_genetic_95"]
    with (args.output/"cross_validated_ages.tsv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=cv_fields,delimiter="\t",
                              lineterminator="\n")
        writer.writeheader(); writer.writerows(cv_rows)
    unknown_fields=["sample","map_T","mean_T","median_T",
                    "ci95_lower_T","ci95_upper_T"]
    with (args.output/"unknown_ages.tsv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=unknown_fields,delimiter="\t",
                              lineterminator="\n")
        writer.writeheader(); writer.writerows(unknown_rows)
    with (args.output/"calibration_bootstrap.tsv").open("w",newline="") as handle:
        writer=csv.writer(handle,delimiter="\t",lineterminator="\n")
        writer.writerow(["replicate","time_scale","epsilon"])
        writer.writerows((i,f"{scale:.12g}",f"{epsilon:.12g}")
                         for i,(scale,epsilon) in enumerate(all_bootstrap))

    calibration=joint_summary(all_profile,time_scale_grid,epsilon_grid)
    calibration["dated_sample_bootstrap"]={
        "replicates":args.bootstrap_replicates,"seed":args.bootstrap_seed,
        "time_scale_mean":float(np.mean(all_bootstrap[:,0])),
        "time_scale_ci95":[float(x) for x in np.quantile(all_bootstrap[:,0],[.025,.975])],
        "epsilon_mean":float(np.mean(all_bootstrap[:,1])),
        "epsilon_ci95":[float(x) for x in np.quantile(all_bootstrap[:,1],[.025,.975])],
    }
    metadata={
        "format_version":1,"likelihood_dir":str(likelihood_dir.resolve()),
        "radiocarbon_input":str(args.radiocarbon.resolve()),
        "n_dated":len(dated),"n_unknown":len(unknown),"folds":fold_map,
        "priors":{"time_scale":"uniform density on supplied grid",
                  "epsilon":"uniform density on supplied grid",
                  "age":"uniform" if args.age_prior is None else str(args.age_prior.resolve())},
        "joint_calibration":calibration,"fold_calibration":fold_summaries,
        "cross_validation_metrics":{"joint_scale_epsilon":base.accuracy_metrics(cv_rows)},
        "epsilon_interpretation":"effective per-allele discordance, not sequencing error",
        "time_scale_interpretation":"multiplicative correction applied inside the insertion likelihood",
        "radiocarbon_use":"training densities integrated; held-out densities used only for scoring",
        "interval_warning":"genetic intervals marginalize the joint calibration surface but do not yet include genomic block-bootstrap uncertainty",
    }
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
    parser.add_argument("--bootstrap-replicates",type=int,default=2000)
    parser.add_argument("--bootstrap-seed",type=int,default=20261028)
    parser.add_argument("--output",type=Path,required=True)
    return parser.parse_args(argv)


def main():
    run(parse_args())


if __name__=="__main__":
    main()
