import numpy as np

from chromosome_scale_profile import profile_row, subset_profiles


def test_subset_profiles_sums_only_selected_chromosomes(monkeypatch):
    ll = np.zeros((2, 1, 3, 1, 2))
    ll[:, :, 0] = 1
    ll[:, :, 1] = 2
    ll[:, :, 2] = 8
    seen = []

    def fake_profile(age_grid, values, rc_age, rc_log_weight):
        seen.append(values.copy())
        return np.sum(values, axis=-1)

    monkeypatch.setattr(
        "chromosome_scale_profile.joint.joint_radiocarbon_log_profile", fake_profile
    )
    result = subset_profiles(
        ll, np.array([0, 1]), ["dated"],
        {"dated": (np.array([0]), np.array([0]))}, [0, 2]
    )
    assert np.all(seen[0] == 9)
    assert result.shape == (1, 2, 1)


def test_profile_row_reports_requested_subset():
    profiles = np.array([[[0.0], [10.0]], [[0.0], [5.0]]])
    row = profile_row(
        "chr2", ["chr2"], profiles, np.array([0.7, 0.9]), np.array([0.01])
    )
    assert row["chromosomes"] == "chr2"
    assert row["joint_scale_map"] == 0.9
    assert row["joint_epsilon_map"] == 0.01
