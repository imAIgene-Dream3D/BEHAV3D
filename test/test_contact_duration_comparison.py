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
    build_target_class_lookup_from_track_adata,
    touching_column_name,
)
from behav3d.analysis.behavior.track.visualization.plots.contact_duration_report import (
    save_track_contact_duration_comparison,
    _build_comparisons,
)

_N_TIMEPOINTS = 30
_MIN_BOUT_LENGTH = 3

# (sample, class, macro_track_id, [bout lengths for its T-cell tracks]) -- s4 has no "plastic"
# macrophage at all, so paired-by-sample comparisons involving "plastic" must drop s4.
_TOUCH_SPEC = [
    ("s1", "round", 100, [5, 7]),
    ("s2", "round", 100, [6, 8]),
    ("s3", "round", 100, [4, 6]),
    ("s4", "round", 100, [5, 7]),
    ("s1", "elongated", 101, [9, 11]),
    ("s2", "elongated", 101, [10, 12]),
    ("s3", "elongated", 101, [8, 10]),
    ("s4", "elongated", 101, [9, 13]),
    ("s1", "plastic", 102, [13, 15]),
    ("s2", "plastic", 102, [14, 16]),
    ("s3", "plastic", 102, [12, 14]),
]
_SAMPLES = ["s1", "s2", "s3", "s4"]
_MACRO_CLASS_BY_ID = {100: "round", 101: "elongated", 102: "plastic"}


def _build_adata_tracks():
    rows = []
    track_id = 0
    for sample, _cls, _macro_id, bout_lengths in _TOUCH_SPEC:
        for _bout_len in bout_lengths:
            rows.append({
                "sample_name": sample, "TrackID": track_id,
                "position_t_min": 0, "position_t_max": _N_TIMEPOINTS - 1,
            })
            track_id += 1
    obs = pd.DataFrame(rows)
    obs.index = [str(i) for i in range(len(obs))]
    return ad.AnnData(X=np.zeros((len(obs), 1)), obs=obs)


def _build_df_timepoints():
    rows = []
    track_id = 0
    for sample, _cls, macro_id, bout_lengths in _TOUCH_SPEC:
        for bout_len in bout_lengths:
            contact = [1] * bout_len + [0] * (_N_TIMEPOINTS - bout_len)
            touching = [str(macro_id)] * bout_len + [""] * (_N_TIMEPOINTS - bout_len)
            for t in range(_N_TIMEPOINTS):
                rows.append({
                    "sample_name": sample, "TrackID": track_id, "position_t": t,
                    "macro_contact": contact[t], "touching_macros": touching[t],
                })
            track_id += 1
    return pd.DataFrame(rows)


def _build_target_track_adata():
    """One row per (sample, macrophage TrackID) touched anywhere in _TOUCH_SPEC, classified by
    its morphology class -- s4 simply never appears with macro_id=102 ("plastic")."""
    rows = []
    for sample in _SAMPLES:
        macro_ids = sorted({macro_id for s, _c, macro_id, _b in _TOUCH_SPEC if s == sample})
        for macro_id in macro_ids:
            rows.append({
                "sample_name": sample, "TrackID": macro_id, "ClusterID": _MACRO_CLASS_BY_ID[macro_id],
            })
    obs = pd.DataFrame(rows)
    obs.index = [str(i) for i in range(len(obs))]
    return ad.AnnData(X=np.zeros((len(obs), 1)), obs=obs)


def _common_kwargs():
    adata_target = _build_target_track_adata()
    lookup = build_target_class_lookup_from_track_adata(adata_target, class_col="ClusterID")
    return dict(
        contact_col="macro_contact",
        min_bout_length=_MIN_BOUT_LENGTH,
        target_class_lookup=lookup,
        touching_col=touching_column_name("macro"),
        time_varying=False,
        target_cell_type_label="macro",
    )


def test_build_comparisons_pairwise_plus_one_vs_rest():
    comparisons = _build_comparisons(["round", "elongated", "plastic"])
    assert len(comparisons) == 6
    pairwise = {(a, b) for a, classes_a, b, classes_b in comparisons if len(classes_a) == 1 and len(classes_b) == 1}
    assert pairwise == {("round", "elongated"), ("round", "plastic"), ("elongated", "plastic")}
    one_vs_rest = [row for row in comparisons if len(row[1]) > 1 or len(row[3]) > 1]
    assert len(one_vs_rest) == 3
    for label_a, classes_a, label_b, classes_b in one_vs_rest:
        assert classes_a == [label_a]
        assert set(classes_b) == {"round", "elongated", "plastic"} - {label_a}

    # 2 classes -> just the one pairwise entry, no duplicate "vs rest".
    assert len(_build_comparisons(["round", "elongated"])) == 1


def _find_row(csv, label_a, label_b):
    """Look up the comparison row for {label_a, label_b} regardless of which side
    ``_build_comparisons`` (sorted by ``_mixed_label_sort_key``) put in group_a vs group_b."""
    match = csv[
        ((csv["group_a"] == label_a) & (csv["group_b"] == label_b))
        | ((csv["group_a"] == label_b) & (csv["group_b"] == label_a))
    ]
    assert len(match) == 1, f"expected exactly one {label_a} vs {label_b} row, found {len(match)}"
    return match.iloc[0]


def _n_for(row, label):
    return row["n_a"] if row["group_a"] == label else row["n_b"]


def test_welch_mode_pdf_csv_and_group_sizes(tmp_path):
    adata_tracks = _build_adata_tracks()
    df_timepoints = _build_df_timepoints()

    result = save_track_contact_duration_comparison(
        adata_tracks, df_timepoints, tmp_path,
        test_mode="welch", minutes_per_frame=2.0, comparisons_per_page=4, verbose=False,
        **_common_kwargs(),
    )

    assert result["n_comparisons"] == 6
    assert result["n_pages"] == 2  # 6 comparisons, 4 per page -> 2 pages
    assert Path(result["pdf_path"]).exists()
    assert Path(result["csv_path"]).exists()

    csv = pd.read_csv(result["csv_path"])
    assert len(csv) == 6
    assert set(csv["test_mode"]) == {"welch"}
    assert csv["pairing_col"].isna().all()

    round_vs_elong = _find_row(csv, "round", "elongated")
    # round touched by 2 tracks/sample x 4 samples = 8; elongated likewise 8.
    assert _n_for(round_vs_elong, "round") == 8
    assert _n_for(round_vs_elong, "elongated") == 8

    round_vs_plastic = _find_row(csv, "round", "plastic")
    # plastic only touched in s1-s3 -> 2 tracks x 3 samples = 6.
    assert _n_for(round_vs_plastic, "plastic") == 6

    # minutes = timepoints * minutes_per_frame for every row.
    assert csv["mean_a_minutes"].to_numpy() == pytest.approx(csv["mean_a_timepoints"].to_numpy() * 2.0)
    assert csv["mean_b_minutes"].to_numpy() == pytest.approx(csv["mean_b_timepoints"].to_numpy() * 2.0)


def test_paired_mode_drops_incomplete_pairing_units(tmp_path):
    adata_tracks = _build_adata_tracks()
    df_timepoints = _build_df_timepoints()

    result = save_track_contact_duration_comparison(
        adata_tracks, df_timepoints, tmp_path,
        test_mode="paired", pairing_col="sample_name", minutes_per_frame=None, verbose=False,
        **_common_kwargs(),
    )
    assert result["n_comparisons"] == 6
    assert result["pairing_col"] == "sample_name"

    csv = pd.read_csv(result["csv_path"])
    assert set(csv["test_mode"]) == {"paired"}
    assert (csv["pairing_col"] == "sample_name").all()
    # No time metadata given -> minutes columns are all NaN, but nothing errors.
    assert csv["mean_a_minutes"].isna().all()

    round_vs_elong = _find_row(csv, "round", "elongated")
    assert round_vs_elong["n_a"] == 4  # all 4 samples have both round and elongated touches
    assert round_vs_elong["n_b"] == 4

    round_vs_plastic = _find_row(csv, "round", "plastic")
    # s4 has no "plastic" touch -> dropped from the paired comparison.
    assert round_vs_plastic["n_a"] == 3
    assert round_vs_plastic["n_b"] == 3


def test_paired_mode_requires_pairing_col():
    adata_tracks = _build_adata_tracks()
    df_timepoints = _build_df_timepoints()
    with pytest.raises(ValueError, match="pairing_col"):
        save_track_contact_duration_comparison(
            adata_tracks, df_timepoints, Path("unused"),
            test_mode="paired", pairing_col=None, **_common_kwargs(),
        )


def test_invalid_test_mode_raises():
    adata_tracks = _build_adata_tracks()
    df_timepoints = _build_df_timepoints()
    with pytest.raises(ValueError, match="test_mode"):
        save_track_contact_duration_comparison(
            adata_tracks, df_timepoints, Path("unused"),
            test_mode="bogus", **_common_kwargs(),
        )


def test_fewer_than_two_classes_writes_placeholder(tmp_path):
    # Restrict the target-class lookup to just "round" so only one touched class remains.
    adata_target = _build_target_track_adata()
    adata_target = adata_target[adata_target.obs["ClusterID"] == "round"].copy()
    lookup = build_target_class_lookup_from_track_adata(adata_target, class_col="ClusterID")

    adata_tracks = _build_adata_tracks()
    df_timepoints = _build_df_timepoints()

    result = save_track_contact_duration_comparison(
        adata_tracks, df_timepoints, tmp_path,
        contact_col="macro_contact", min_bout_length=_MIN_BOUT_LENGTH,
        target_class_lookup=lookup, touching_col=touching_column_name("macro"),
        time_varying=False, target_cell_type_label="macro",
        test_mode="welch",
    )
    assert result["n_comparisons"] == 0
    assert result["n_pages"] == 1
    assert Path(result["pdf_path"]).exists()
    csv = pd.read_csv(result["csv_path"])
    assert len(csv) == 0
