import numpy as np

import age_dependent_scale_diagnostic as diagnostic


def test_stratified_scale_test_detects_opposing_age_preferences():
    scale=np.array([.7,1.,1.3])
    epsilon=np.array([0.,.01])
    profiles={
        "young":np.array([[10.,9.],[2.,2.],[0.,0.]]),
        "middle":np.array([[0.,0.],[10.,9.],[0.,0.]]),
        "old":np.array([[0.,0.],[2.,2.],[10.,9.]]),
    }
    result=diagnostic.stratified_common_epsilon(profiles,scale,epsilon)
    assert result["stratified_epsilon"]==0.
    assert result["stratified_scales"]=={"young":.7,"middle":1.,"old":1.3}
    assert result["composite_loglik_gain"]>0


def test_age_strata_include_boundary_ages_once():
    assert diagnostic.stratum_name(3000)=="young"
    assert diagnostic.stratum_name(3000.1)=="middle"
    assert diagnostic.stratum_name(6500)=="middle"
    assert diagnostic.stratum_name(6500.1)=="old"
