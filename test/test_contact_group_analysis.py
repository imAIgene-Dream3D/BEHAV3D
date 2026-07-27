import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/behav3d-mpl-cache")

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg", force=True)

from behav3d.analysis.behavior.track.contact_grouping import (
    compute_track_contact_features,
    merge_track_contact_features_into_obs,
)
from behav3d.analysis.behavior.track.visualization.plots.reports import (
    save_track_contact_group_analysis,
)

# (sample_name, organoid_line, treatment) — one row per sample; organoid_line/treatment are
# sample-level metadata (2x2: {lineX, lineY} x {drugA, drugB}).
_SAMPLES = [
    ("sample_1", "lineX", "drugA"),
    ("sample_2", "lineY", "drugB"),
    ("sample_3", "lineX", "drugB"),
    ("sample_4", "lineY", "drugA"),
]
# (ClusterID, has_contact_bout) patterns, alternated per sample so each sample's *first* track
# (relevant for `save_track_condition_comparison_report`'s per-sample condition dedup) differs
# between samples — otherwise every sample would collapse to the same condition_col level.
_CLUSTER_PATTERN_A = [("A", True), ("B", False), ("A", False), ("B", True)]
_CLUSTER_PATTERN_B = [("A", False), ("B", True), ("A", True), ("B", False)]
_N_TIMEPOINTS = 10
_MIN_BOUT_LENGTH = 5

_TRACKS = [
    (sample_name, organoid_line, treatment, cluster, has_contact)
    for i, (sample_name, organoid_line, treatment) in enumerate(_SAMPLES)
    for cluster, has_contact in (_CLUSTER_PATTERN_A if i % 2 == 0 else _CLUSTER_PATTERN_B)
]


def _build_adata_tracks():
    obs = pd.DataFrame(
        {
            "sample_name": [t[0] for t in _TRACKS],
            "TrackID": list(range(len(_TRACKS))),
            "organoid_line": [t[1] for t in _TRACKS],
            "treatment": [t[2] for t in _TRACKS],
            "ClusterID": [t[3] for t in _TRACKS],
            "position_t_min": 0,
            "position_t_max": _N_TIMEPOINTS - 1,
        }
    )
    obs.index = [str(i) for i in range(len(_TRACKS))]
    return ad.AnnData(X=np.zeros((len(_TRACKS), 1)), obs=obs)


def _build_df_timepoints():
    rows = []
    for track_id, (sample_name, _line, _treatment, _cluster, has_contact) in enumerate(_TRACKS):
        contact_values = np.zeros(_N_TIMEPOINTS, dtype=int)
        if has_contact:
            contact_values[: _MIN_BOUT_LENGTH + 1] = 1
        for t in range(_N_TIMEPOINTS):
            rows.append(
                {
                    "sample_name": sample_name,
                    "TrackID": track_id,
                    "position_t": t,
                    "macro_contact": int(contact_values[t]),
                }
            )
    return pd.DataFrame(rows)


def test_contact_group_analysis_extra_group_cols(tmp_path):
    adata_tracks = _build_adata_tracks()
    df_timepoints = _build_df_timepoints()

    result = save_track_contact_group_analysis(
        adata_tracks,
        df_timepoints,
        tmp_path,
        contact_col="macro_contact",
        min_bout_length=_MIN_BOUT_LENGTH,
        sample_col="sample_name",
        class_col="ClusterID",
        extra_group_cols=["organoid_line"],
        verbose=False,
    )

    assert result["extra_group_cols"] == ["organoid_line"]
    assert Path(result["pdf_path"]).exists()

    # No group_x/group_y — cluster-stack grid pages are paginated per organoid_line value only
    # (one single-panel page per value), no true 2D grid.
    assert result["cluster_stack_grid"]["n_pages"] == 2
    grid_csv = pd.read_csv(result["cluster_stack_grid"]["csv_path"])
    assert "organoid_line" in grid_csv.columns
    assert set(grid_csv["organoid_line"].unique()) == {"lineX", "lineY"}
    assert set(grid_csv["ClusterID"].unique()) == {"A", "B"}

    reverse_csv = pd.read_csv(result["contact_rate_by_cluster"]["csv_path"])
    assert reverse_csv is not None

    diff_csv_path = result["condition_comparison"]["csv_path"]
    assert diff_csv_path is not None
    diff_csv = pd.read_csv(diff_csv_path)
    assert set(diff_csv["group"].unique()) == {"lineX", "lineY"}


def test_contact_group_analysis_group_x_only(tmp_path):
    adata_tracks = _build_adata_tracks()
    df_timepoints = _build_df_timepoints()

    result = save_track_contact_group_analysis(
        adata_tracks,
        df_timepoints,
        tmp_path,
        contact_col="macro_contact",
        min_bout_length=_MIN_BOUT_LENGTH,
        sample_col="sample_name",
        class_col="ClusterID",
        group_x="organoid_line",
        verbose=False,
    )

    # Only one axis set — condition_comparison stays in the old pooled-column layout.
    assert result["condition_comparison"]["csv_path"] is not None
    diff_csv = pd.read_csv(result["condition_comparison"]["csv_path"])
    assert set(diff_csv["group"].unique()) == {"lineX", "lineY"}

    # Cluster-stack grid: single facet axis (organoid_line), true grid with 1 dimension.
    assert result["cluster_stack_grid"]["n_pages"] == 1
    grid_csv = pd.read_csv(result["cluster_stack_grid"]["csv_path"])
    assert set(grid_csv["organoid_line"].unique()) == {"lineX", "lineY"}


def test_contact_group_analysis_group_x_and_group_y(tmp_path):
    adata_tracks = _build_adata_tracks()
    df_timepoints = _build_df_timepoints()

    result = save_track_contact_group_analysis(
        adata_tracks,
        df_timepoints,
        tmp_path,
        contact_col="macro_contact",
        min_bout_length=_MIN_BOUT_LENGTH,
        sample_col="sample_name",
        class_col="ClusterID",
        group_x="organoid_line",
        group_y="treatment",
        verbose=False,
    )

    # Both axes set — cluster-stack grid is a true 2D grid, one page (no extra_group_cols).
    assert result["cluster_stack_grid"]["n_pages"] == 1
    grid_csv = pd.read_csv(result["cluster_stack_grid"]["csv_path"])
    assert set(grid_csv["organoid_line"].unique()) == {"lineX", "lineY"}
    assert set(grid_csv["treatment"].unique()) == {"drugA", "drugB"}

    # condition_comparison switches to the true group_x x group_y 2D grid (binary condition_col).
    diff_csv_path = result["condition_comparison"]["csv_path"]
    assert diff_csv_path is not None
    diff_csv = pd.read_csv(diff_csv_path)
    assert "organoid_line" in diff_csv.columns
    assert "treatment" in diff_csv.columns
    assert "group" not in diff_csv.columns
    assert set(diff_csv["organoid_line"].unique()) == {"lineX", "lineY"}
    assert set(diff_csv["treatment"].unique()) == {"drugA", "drugB"}


def test_contact_group_analysis_without_extra_group_cols(tmp_path):
    adata_tracks = _build_adata_tracks()
    df_timepoints = _build_df_timepoints()

    result = save_track_contact_group_analysis(
        adata_tracks,
        df_timepoints,
        tmp_path,
        contact_col="macro_contact",
        min_bout_length=_MIN_BOUT_LENGTH,
        sample_col="sample_name",
        class_col="ClusterID",
        verbose=False,
    )

    assert result["extra_group_cols"] == []
    assert result["condition_comparison"]["csv_path"] is not None
    diff_csv = pd.read_csv(result["condition_comparison"]["csv_path"])
    assert set(diff_csv["group"].unique()) == {"(all)"}

    # No group_x/group_y/extra_group_cols at all — the cluster-stack grid page is skipped
    # entirely; only the plain per-sample contact-rate page is produced.
    assert result["cluster_stack_grid"] == {"n_pages": 0, "csv_path": None}


def _group_col():
    return "macro_contact_group"


def test_track_missing_from_df_timepoints_raises():
    """A classified track absent from df_timepoints entirely (e.g. stale CSV, ID typo) must raise
    a hard error rather than silently resolving to "no_contact" — that would misrepresent a
    stale/mismatched input as a real negative result. See _assert_tracks_have_window_data."""
    obs = pd.DataFrame(
        {
            "sample_name": ["sample_1", "sample_1"],
            "TrackID": [0, 1],
            "position_t_min": [0, 0],
            "position_t_max": [9, 9],
        }
    )
    obs.index = ["0", "1"]
    adata_tracks = ad.AnnData(X=np.zeros((2, 1)), obs=obs)

    # Only TrackID 0 has any rows in df_timepoints; TrackID 1 is entirely absent.
    df_timepoints = pd.DataFrame(
        {
            "sample_name": ["sample_1"] * 10,
            "TrackID": [0] * 10,
            "position_t": list(range(10)),
            "macro_contact": [0] * 10,
        }
    )

    with pytest.raises(ValueError, match="have no matching timepoints in the filtered track-features CSV"):
        compute_track_contact_features(
            df_timepoints,
            adata_tracks,
            contact_col="macro_contact",
            min_bout_length=5,
            verbose=True,
        )


def test_track_present_in_csv_but_no_timepoints_in_window_raises():
    """A track present in df_timepoints but with zero timepoints falling inside its classified
    [position_t_min, position_t_max] window must also raise. adata_tracks is built directly from
    this CSV, so position_t_min/max are always the min/max of real rows for that track — zero
    overlap can only mean the CSV and h5ad are out of sync, never a legitimate "no contact"
    result. See _assert_tracks_have_window_data."""
    obs = pd.DataFrame(
        {
            "sample_name": ["sample_1", "sample_1"],
            "TrackID": [0, 1],
            "position_t_min": [0, 100],
            "position_t_max": [9, 109],
        }
    )
    obs.index = ["0", "1"]
    adata_tracks = ad.AnnData(X=np.zeros((2, 1)), obs=obs)

    # TrackID 1 is present in df_timepoints, but only at t=0-9 — outside its [100, 109] window.
    df_timepoints = pd.DataFrame(
        {
            "sample_name": ["sample_1"] * 20,
            "TrackID": [0] * 10 + [1] * 10,
            "position_t": list(range(10)) + list(range(10)),
            "macro_contact": [0] * 20,
        }
    )

    with pytest.raises(ValueError, match="have no matching timepoints in the filtered track-features CSV"):
        compute_track_contact_features(
            df_timepoints,
            adata_tracks,
            contact_col="macro_contact",
            min_bout_length=5,
            verbose=True,
        )


def test_split_long_track_windows_get_independent_contact_groups():
    """When split_long_tracks produces multiple obs rows per TrackID (distinguished by
    trajectory_window_id), each window must be evaluated against its own slice of timepoints
    rather than being collapsed into a single row."""
    obs = pd.DataFrame(
        {
            "sample_name": ["sample_1", "sample_1"],
            "TrackID": [0, 0],
            "trajectory_window_id": [0, 1],
            "position_t_min": [0, 10],
            "position_t_max": [9, 19],
        }
    )
    obs.index = ["0", "1"]
    adata_tracks = ad.AnnData(X=np.zeros((2, 1)), obs=obs)

    # Window 0 (t=0-9): a long contact bout. Window 1 (t=10-19): no contact.
    contact_values = [1] * 6 + [0] * 4 + [0] * 10
    df_timepoints = pd.DataFrame(
        {
            "sample_name": ["sample_1"] * 20,
            "TrackID": [0] * 20,
            "position_t": list(range(20)),
            "macro_contact": contact_values,
        }
    )

    features = compute_track_contact_features(
        df_timepoints,
        adata_tracks,
        contact_col="macro_contact",
        min_bout_length=5,
    )
    assert len(features) == 2
    row0 = features.loc[("sample_1", "0", "0")]
    row1 = features.loc[("sample_1", "0", "1")]
    assert row0["macro_contact_group"] == "contact"
    assert row1["macro_contact_group"] == "no_contact"

    merge_track_contact_features_into_obs(
        adata_tracks, features, contact_col="macro_contact", min_bout_length=5
    )
    merged = adata_tracks.obs.set_index("trajectory_window_id")[_group_col()]
    assert merged.loc["0"] == "contact"
    assert merged.loc["1"] == "no_contact"


def test_track_id_dtype_drift_still_matches():
    """TrackID as int in adata_tracks.obs vs. float-string ("5.0") in a re-read CSV — common after
    a pandas round-trip through a numeric column with NaNs elsewhere — must still join correctly
    instead of silently failing to match."""
    obs = pd.DataFrame(
        {
            "sample_name": ["sample_1"],
            "TrackID": [5],
            "position_t_min": [0],
            "position_t_max": [9],
        }
    )
    obs.index = ["0"]
    adata_tracks = ad.AnnData(X=np.zeros((1, 1)), obs=obs)

    df_timepoints = pd.DataFrame(
        {
            "sample_name": ["sample_1"] * 10,
            "TrackID": ["5.0"] * 10,
            "position_t": list(range(10)),
            "macro_contact": [1] * 6 + [0] * 4,
        }
    )

    features = compute_track_contact_features(
        df_timepoints,
        adata_tracks,
        contact_col="macro_contact",
        min_bout_length=5,
    )
    assert features.loc[("sample_1", "5"), "macro_contact_group"] == "contact"
    assert features.loc[("sample_1", "5"), "macro_contact_max_bout_length"] == 6
