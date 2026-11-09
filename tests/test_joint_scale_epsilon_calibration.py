import csv
import json
from argparse import Namespace

import numpy as np

import joint_scale_epsilon_calibration as joint
import radiocarbon_calibration as base


def test_one_scale_matches_epsilon_only_marginalization():
    age=np.array([0.,10.,20.])
    epsilon=np.array([0.,.1])
    profile=np.array([0.,-2.])
    ll=np.array([[0.,-1.,-4.],[-3.,0.,-3.]])
    epsilon_only=base.marginal_age_loglik(ll,profile,epsilon)
    log_joint=joint.joint_log_posterior(profile[None,:],np.array([1.]),epsilon)
    observed=joint.marginal_age_loglik_weights(ll[None,:,:],log_joint)
    np.testing.assert_allclose(observed,epsilon_only)


def test_joint_cross_validation_excludes_held_out_dates(tmp_path):
    likelihood=tmp_path/"likelihood"; likelihood.mkdir()
    age=np.array([0.,10.]); scale=np.array([1.,2.]); epsilon=np.array([0.,.1])
    samples=["young","old","unknown"]
    ll=np.full((2,2,3,2),-8.)
    # Each dated sample favors a different calibration corner at its RC age.
    ll[0,0,0,0]=4
    ll[1,1,1,1]=4
    ll[:,:,2,1]=2
    np.save(likelihood/"grid.npy",age)
    np.save(likelihood/"time_scale_grid.npy",scale)
    np.save(likelihood/"epsilon_grid.npy",epsilon)
    np.save(likelihood/"ll_marginal_scale_epsilon.npy",ll)
    (likelihood/"samples.txt").write_text("\n".join(samples)+"\n")
    radiocarbon=tmp_path/"radiocarbon.tsv"
    radiocarbon.write_text(
        "sample\tage_generations\tdensity\n"
        "young\t0\t1\n"
        "old\t10\t1\n")
    output=tmp_path/"output"
    args=Namespace(likelihood_dir=likelihood,radiocarbon=radiocarbon,
                   folds=None,n_folds=0,unknown_samples=None,age_prior=None,
                   bootstrap_replicates=50,bootstrap_seed=11,output=output)
    joint.run(args)

    with (output/"cross_validated_ages.tsv").open() as handle:
        rows={row["sample"]:row for row in csv.DictReader(handle,delimiter="\t")}
    # young is calibrated only by old, and old only by young.
    assert float(rows["young"]["fold_time_scale_map"])==2
    assert float(rows["young"]["fold_epsilon_map"])==.1
    assert float(rows["old"]["fold_time_scale_map"])==1
    assert float(rows["old"]["fold_epsilon_map"])==0
    with (output/"unknown_ages.tsv").open() as handle:
        unknown=list(csv.DictReader(handle,delimiter="\t"))
    assert len(unknown)==1 and unknown[0]["sample"]=="unknown"
    assert float(unknown[0]["map_T"])==10
    metadata=json.loads((output/"run.json").read_text())
    assert metadata["n_dated"]==2
    assert metadata["cross_validation_metrics"]["joint_scale_epsilon"]["n"]==2
    assert (output/"joint_scale_epsilon_profile.tsv").exists()
    assert (output/"scale_marginal.tsv").exists()
    assert (output/"epsilon_marginal.tsv").exists()
    assert len((output/"calibration_bootstrap.tsv").read_text().splitlines())==51


def test_joint_bootstrap_is_deterministic():
    profiles=np.arange(3*2*2,dtype=float).reshape(3,2,2)
    scale=np.array([.8,1.2]); epsilon=np.array([0.,.1])
    first=joint.bootstrap_joint(profiles,scale,epsilon,25,17)
    second=joint.bootstrap_joint(profiles,scale,epsilon,25,17)
    np.testing.assert_array_equal(first,second)


def test_joint_summary_interval_contains_concentrated_map_grid_point():
    scale=np.array([.8,.9,1.,1.1])
    epsilon=np.array([0.,.02,.04])
    profile=np.full((len(scale),len(epsilon)),-1000.)
    profile[2,1]=0.
    summary=joint.joint_summary(profile,scale,epsilon)
    assert summary["joint_map"]=={"time_scale":1.,"epsilon":.02}
    assert summary["time_scale"]["ci95"]==[1.,1.]
    assert summary["epsilon"]["ci95"]==[.02,.02]
