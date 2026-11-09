import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import tskit

from singer_polegon_diagnostic import (
    calibration_metrics,
    chunk_inputs,
    ratio_summary,
    tree_sequence_counts,
)


def test_chunk_inputs_selects_raw_nodes_and_checks_polegon(tmp_path):
    chunks = tmp_path / "chr1" / "chunks"
    chunks.mkdir(parents=True)
    for chunk in ("000002", "000001"):
        (chunks / f"{chunk}.yaml").write_text("singer: {}\n")
        for label in ("recombs", "nodes", "muts", "branches", "polegon"):
            (chunks / f"{chunk}_{label}_50.txt").write_text("0\n")
    selected = chunk_inputs(tmp_path / "chr1", "chr1", 50)
    assert [Path(path).name for path in selected["params"]] == ["000001.yaml", "000002.yaml"]
    assert all("_nodes_50.txt" in path and "_polegon_" not in path for path in selected["nodes"])


def dated_tree(parent_time):
    tables = tskit.TableCollection(10)
    child = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    parent = tables.nodes.add_row(time=parent_time)
    tables.edges.add_row(0, 10, parent=parent, child=child)
    tables.sort()
    return tables.tree_sequence()


def test_ratio_summary_and_tree_counts():
    assert ratio_summary(np.array([2, 4, 8]), np.array([1, 2, 4]))["median"] == 2
    result = tree_sequence_counts(dated_tree(10), dated_tree(20))
    assert result["raw_singer"]["edges"] == 1
    assert result["polegon"]["edges"] == 1


def test_calibration_metrics_extracts_joint_cv(tmp_path):
    payload = {
        "joint_calibration": {"joint_map": {"time_scale": 0.7, "epsilon": 0.01}},
        "cross_validation_metrics": {
            "joint_scale_epsilon": {
                "mean_bias": -10, "mae": 20, "rmse": 30, "regression_slope": 0.8
            }
        },
    }
    (tmp_path / "run.json").write_text(json.dumps(payload))
    assert calibration_metrics(tmp_path) == {
        "time_scale": 0.7,
        "epsilon": 0.01,
        "mean_bias": -10,
        "mae": 20,
        "rmse": 30,
        "regression_slope": 0.8,
    }
