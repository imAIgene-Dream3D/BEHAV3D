import os

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
    build_target_class_lookup_from_state_adata,
    build_target_class_lookup_from_track_adata,
    compute_track_contact_target_class_features,
    merge_track_target_class_group_into_obs,
    contact_col_target_cell_type,
    touching_column_name,
)

_N_TIMEPOINTS = 20
_MIN_BOUT_LENGTH = 5


def test_contact_col_target_cell_type():
    assert contact_col_target_cell_type("macro_contact") == "macro"
    assert contact_col_target_cell_type("macro_contact_on_distance") == "macro"
    with pytest.raises(ValueError):
        contact_col_target_cell_type("macro_distance")


def test_touching_column_name():
    assert touching_column_name("macro") == "touching_macros"


def _build_reference_tracks():
    """4 T cell tracks: 0 contacts only round macros, 1 only elongated, 2 both (switching
    mid-bout, so the *overall* min-bout-length is met by a run that spans two different target
    classes), 3 never contacts."""
    obs = pd.DataFrame(
        {
            "sample_name": ["s1"] * 4,
            "TrackID": [0, 1, 2, 3],
            "position_t_min": 0,
            "position_t_max": _N_TIMEPOINTS - 1,
        }
    )
    obs.index = [str(i) for i in range(4)]
    return ad.AnnData(X=np.zeros((4, 1)), obs=obs)


def _build_reference_timepoints():
    rows = []

    def add(track_id, contact_flags, touching):
        for t in range(_N_TIMEPOINTS):
            rows.append(
                {
                    "sample_name": "s1",
                    "TrackID": track_id,
                    "position_t": t,
                    "macro_contact": int(contact_flags[t]),
                    "touching_macros": touching[t],
                }
            )

    add(0, [1] * 6 + [0] * (_N_TIMEPOINTS - 6), ["100"] * 6 + [""] * (_N_TIMEPOINTS - 6))
    add(1, [1] * 6 + [0] * (_N_TIMEPOINTS - 6), ["101"] * 6 + [""] * (_N_TIMEPOINTS - 6))
    # 3 timepoints touching macro 100 (round), immediately followed by 3 touching macro 101
    # (elongated) — a single contiguous 6-timepoint "any contact" bout (>= min_bout_length),
    # but neither class alone reaches the threshold on its own.
    add(
        2,
        [1] * 6 + [0] * (_N_TIMEPOINTS - 6),
        ["100"] * 3 + ["101"] * 3 + [""] * (_N_TIMEPOINTS - 6),
    )
    add(3, [0] * _N_TIMEPOINTS, [""] * _N_TIMEPOINTS)
    return pd.DataFrame(rows)


def _build_target_states_adata():
    """Target (macro) per-timepoint states: TrackID 100 = "round" always, 101 = "elongated"."""
    obs = pd.DataFrame(
        {
            "sample_name": ["s1"] * (2 * _N_TIMEPOINTS),
            "TrackID": [100] * _N_TIMEPOINTS + [101] * _N_TIMEPOINTS,
            "position_t": list(range(_N_TIMEPOINTS)) * 2,
            "behavioral_state": ["round"] * _N_TIMEPOINTS + ["elongated"] * _N_TIMEPOINTS,
        }
    )
    return ad.AnnData(X=np.zeros((len(obs), 1)), obs=obs)


def _build_target_track_classification_adata():
    """Target (macro) whole-track classification: TrackID 100 = "round", 101 = "elongated"."""
    obs = pd.DataFrame(
        {
            "sample_name": ["s1", "s1"],
            "TrackID": [100, 101],
            "ClusterID": ["round", "elongated"],
        }
    )
    obs.index = ["0", "1"]
    return ad.AnnData(X=np.zeros((2, 1)), obs=obs)


def _run_overall_contact_grouping(adata_tracks, df_timepoints):
    features = compute_track_contact_features(
        df_timepoints, adata_tracks, contact_col="macro_contact", min_bout_length=_MIN_BOUT_LENGTH,
    )
    merge_track_contact_features_into_obs(
        adata_tracks, features, contact_col="macro_contact", min_bout_length=_MIN_BOUT_LENGTH,
    )


def _bucket_by_track_id(adata_tracks, col):
    return adata_tracks.obs.set_index("TrackID")[col].to_dict()


@pytest.mark.parametrize("time_varying", [True, False])
def test_target_class_group_bucketing(time_varying):
    adata_tracks = _build_reference_tracks()
    df_timepoints = _build_reference_timepoints()
    _run_overall_contact_grouping(adata_tracks, df_timepoints)

    if time_varying:
        adata_states = _build_target_states_adata()
        lookup = build_target_class_lookup_from_state_adata(adata_states, state_col="behavioral_state")
    else:
        adata_target_tracks = _build_target_track_classification_adata()
        lookup = build_target_class_lookup_from_track_adata(adata_target_tracks, class_col="ClusterID")

    long_df, group_df = compute_track_contact_target_class_features(
        df_timepoints, adata_tracks, lookup,
        contact_col="macro_contact", touching_col="touching_macros", time_varying=time_varying,
    )
    merge_track_target_class_group_into_obs(adata_tracks, group_df, contact_col="macro_contact")

    buckets = _bucket_by_track_id(adata_tracks, "macro_contact_target_class_group")
    assert buckets["0"] == "round"
    assert buckets["1"] == "elongated"
    assert buckets["2"] == "Multiple classes"
    assert buckets["3"] == "no_contact"

    # Track 2 touched both classes -> 2 rows in the long-format per-class table; tracks 0/1 -> 1 each.
    long_flat = long_df.reset_index()
    assert set(long_flat.loc[long_flat["TrackID"] == "0", "target_class"]) == {"round"}
    assert set(long_flat.loc[long_flat["TrackID"] == "1", "target_class"]) == {"elongated"}
    assert set(long_flat.loc[long_flat["TrackID"] == "2", "target_class"]) == {"round", "elongated"}
    assert len(long_flat.loc[long_flat["TrackID"] == "3"]) == 0

    # Track 0 touched macro 100 (round) at every one of its 6 contact timepoints, contiguously.
    row0 = long_flat[(long_flat["TrackID"] == "0") & (long_flat["target_class"] == "round")].iloc[0]
    assert row0["macro_contact_class_mean_fraction"] == pytest.approx(6 / _N_TIMEPOINTS)
    assert row0["macro_contact_class_max_bout_length"] == 6

    # Track 2 split its 6-timepoint bout 3/3 between round and elongated.
    row2_round = long_flat[(long_flat["TrackID"] == "2") & (long_flat["target_class"] == "round")].iloc[0]
    row2_elong = long_flat[(long_flat["TrackID"] == "2") & (long_flat["target_class"] == "elongated")].iloc[0]
    assert row2_round["macro_contact_class_max_bout_length"] == 3
    assert row2_elong["macro_contact_class_max_bout_length"] == 3


def test_missing_touching_column_raises_clear_error():
    adata_tracks = _build_reference_tracks()
    df_timepoints = _build_reference_timepoints().drop(columns=["touching_macros"])
    _run_overall_contact_grouping(adata_tracks, df_timepoints)

    adata_states = _build_target_states_adata()
    lookup = build_target_class_lookup_from_state_adata(adata_states, state_col="behavioral_state")

    with pytest.raises(KeyError, match="touching_macros"):
        compute_track_contact_target_class_features(
            df_timepoints, adata_tracks, lookup,
            contact_col="macro_contact", touching_col="touching_macros", time_varying=True,
        )


def test_build_target_class_lookup_from_state_adata_missing_column_raises():
    adata_states = _build_target_states_adata()
    with pytest.raises(KeyError):
        build_target_class_lookup_from_state_adata(adata_states, state_col="does_not_exist")


def test_build_target_class_lookup_from_track_adata_missing_column_raises():
    adata_target_tracks = _build_target_track_classification_adata()
    with pytest.raises(KeyError):
        build_target_class_lookup_from_track_adata(adata_target_tracks, class_col="does_not_exist")
