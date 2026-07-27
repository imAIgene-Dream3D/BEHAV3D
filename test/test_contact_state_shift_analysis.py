import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/behav3d-mpl-cache")

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg", force=True)

from behav3d.analysis.behavior.track.contact_state_shift import (
    _first_qualifying_bout_bounds,
    compute_contact_bout_windows,
    compute_track_state_shift_features,
    summarize_state_shift_track_fractions,
    BEFORE_LABEL,
    AFTER_LABEL,
)
from behav3d.analysis.behavior.track.contact_grouping import compute_track_contact_features
from behav3d.analysis.behavior.track.visualization.plots.contact_state_shift_report import (
    save_track_contact_state_shift_report,
)

_N_TIMEPOINTS = 40
_MIN_BOUT_LENGTH = 5
_STATE_COL = "behavioral_state"

# (TrackID, contact_group, bout_start) — bout occupies [bout_start, bout_start + _MIN_BOUT_LENGTH - 1]
# for "contact" tracks; "no_contact" tracks never have any contact timepoints.
_TRACK_SPECS = [
    (0, "contact", 10),
    (1, "contact", 20),
    (2, "contact", 30),
    (3, "no_contact", None),
    (4, "no_contact", None),
    (5, "no_contact", None),
]


def _state_for_contact_track(t, bout_start):
    bout_end = bout_start + _MIN_BOUT_LENGTH - 1
    if t < bout_start:
        return "static"
    if t <= bout_end:
        return "scanning"
    return "engaging"


def _build_adata_tracks(track_id_dtype=int):
    obs = pd.DataFrame(
        {
            "sample_name": ["sample_1"] * len(_TRACK_SPECS),
            "TrackID": [track_id_dtype(t[0]) for t in _TRACK_SPECS],
            "position_t_min": 0,
            "position_t_max": _N_TIMEPOINTS - 1,
        }
    )
    obs.index = [str(t[0]) for t in _TRACK_SPECS]
    return ad.AnnData(X=np.zeros((len(_TRACK_SPECS), 1)), obs=obs)


def _build_df_timepoints(track_id_dtype=int):
    rows = []
    for track_id, group, bout_start in _TRACK_SPECS:
        contact = np.zeros(_N_TIMEPOINTS, dtype=int)
        if group == "contact":
            bout_end = bout_start + _MIN_BOUT_LENGTH - 1
            contact[bout_start : bout_end + 1] = 1
        for t in range(_N_TIMEPOINTS):
            rows.append({
                "sample_name": "sample_1",
                "TrackID": track_id_dtype(track_id),
                "position_t": t,
                "macro_contact": int(contact[t]),
            })
    return pd.DataFrame(rows)


def _build_adata_states(track_id_dtype=int):
    rows = []
    for track_id, group, bout_start in _TRACK_SPECS:
        for t in range(_N_TIMEPOINTS):
            if group == "contact":
                state = _state_for_contact_track(t, bout_start)
            else:
                state = "static"
            rows.append({
                "sample_name": "sample_1",
                "TrackID": track_id_dtype(track_id),
                "position_t": t,
                _STATE_COL: state,
            })
    obs = pd.DataFrame(rows)
    obs.index = [str(i) for i in range(len(obs))]
    return ad.AnnData(X=np.zeros((len(obs), 1)), obs=obs)


def test_first_qualifying_bout_uses_earliest_run():
    times = list(range(20))
    # Two qualifying runs (length >= 5): [2..6] and [12..17]; must return the first.
    is_contact = [0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0]
    start_t, end_t = _first_qualifying_bout_bounds(times, is_contact, min_bout_length=5)
    assert (start_t, end_t) == (2, 6)


def test_first_qualifying_bout_skips_short_runs():
    times = list(range(10))
    is_contact = [1, 1, 0, 1, 1, 1, 1, 1, 0, 0]  # first run length 2 (too short), second length 5
    start_t, end_t = _first_qualifying_bout_bounds(times, is_contact, min_bout_length=5)
    assert (start_t, end_t) == (3, 7)


def test_first_qualifying_bout_none_when_no_run_qualifies():
    times = list(range(10))
    is_contact = [1, 1, 0, 1, 1, 1, 0, 0, 0, 0]
    start_t, end_t = _first_qualifying_bout_bounds(times, is_contact, min_bout_length=5)
    assert (start_t, end_t) == (None, None)


def test_shared_min_bout_length_matches_existing_contact_grouping():
    """The contact/no_contact split produced here must agree with
    contact_grouping.compute_track_contact_features for the same min_bout_length, since both
    analyses are meant to share one definition of 'a real contact bout'."""
    adata_tracks = _build_adata_tracks()
    df_timepoints = _build_df_timepoints()

    windows = compute_contact_bout_windows(
        df_timepoints, adata_tracks, contact_col="macro_contact", min_bout_length=_MIN_BOUT_LENGTH,
        window_mode="fixed", fixed_window_length=5,
    )
    existing = compute_track_contact_features(
        df_timepoints, adata_tracks, contact_col="macro_contact", min_bout_length=_MIN_BOUT_LENGTH,
    )

    for track_id, group, _bout_start in _TRACK_SPECS:
        key = ("sample_1", str(track_id))
        assert windows.loc[key, "contact_group"] == group
        assert existing.loc[key, "macro_contact_group"] == group


def test_fixed_window_mode_bounds():
    adata_tracks = _build_adata_tracks()
    df_timepoints = _build_df_timepoints()

    windows = compute_contact_bout_windows(
        df_timepoints, adata_tracks, contact_col="macro_contact", min_bout_length=_MIN_BOUT_LENGTH,
        window_mode="fixed", fixed_window_length=5,
    )
    row = windows.loc[("sample_1", "0")]  # bout_start=10, bout_end=14
    assert row["before_start_t"] == 5
    assert row["before_end_t"] == 10
    assert row["after_start_t"] == 14
    assert row["after_end_t"] == 19
    assert row["before_n_timepoints"] == 5
    assert row["after_n_timepoints"] == 5
    assert not row["excluded"]


def test_full_window_mode_bounds():
    adata_tracks = _build_adata_tracks()
    df_timepoints = _build_df_timepoints()

    windows = compute_contact_bout_windows(
        df_timepoints, adata_tracks, contact_col="macro_contact", min_bout_length=_MIN_BOUT_LENGTH,
        window_mode="full",
    )
    row = windows.loc[("sample_1", "0")]  # bout_start=10, bout_end=14, track window [0, 39]
    assert row["before_start_t"] == 0
    assert row["before_end_t"] == 10
    assert row["after_start_t"] == 14
    assert row["after_end_t"] == 39
    assert row["before_n_timepoints"] == 10
    assert row["after_n_timepoints"] == 25


def test_min_window_timepoints_excludes_short_windows():
    adata_tracks = _build_adata_tracks()
    df_timepoints = _build_df_timepoints()

    # bout_start=10 with fixed_window_length=20 -> before window would need [-10, 10), but is
    # clipped to [0, 10) => only 10 timepoints, which is fine; force a short window instead by
    # requiring more timepoints than available before t=10.
    windows = compute_contact_bout_windows(
        df_timepoints, adata_tracks, contact_col="macro_contact", min_bout_length=_MIN_BOUT_LENGTH,
        window_mode="fixed", fixed_window_length=5, min_window_timepoints=6,
    )
    row = windows.loc[("sample_1", "0")]  # before/after windows both have exactly 5 timepoints < 6
    assert row["excluded"]
    assert row["excluded_reason"] == "both_too_short"


def test_null_reference_reproducible_with_seed():
    adata_tracks = _build_adata_tracks()
    df_timepoints = _build_df_timepoints()

    windows_a = compute_contact_bout_windows(
        df_timepoints, adata_tracks, contact_col="macro_contact", min_bout_length=_MIN_BOUT_LENGTH,
        window_mode="fixed", fixed_window_length=5, null_seed=42,
    )
    windows_b = compute_contact_bout_windows(
        df_timepoints, adata_tracks, contact_col="macro_contact", min_bout_length=_MIN_BOUT_LENGTH,
        window_mode="fixed", fixed_window_length=5, null_seed=42,
    )
    windows_c = compute_contact_bout_windows(
        df_timepoints, adata_tracks, contact_col="macro_contact", min_bout_length=_MIN_BOUT_LENGTH,
        window_mode="fixed", fixed_window_length=5, null_seed=7,
    )

    no_contact_keys = [("sample_1", str(t[0])) for t in _TRACK_SPECS if t[1] == "no_contact"]
    bout_starts_a = windows_a.loc[no_contact_keys, "bout_start_t"]
    bout_starts_b = windows_b.loc[no_contact_keys, "bout_start_t"]
    bout_starts_c = windows_c.loc[no_contact_keys, "bout_start_t"]

    assert bout_starts_a.tolist() == bout_starts_b.tolist()
    assert bout_starts_a.tolist() != bout_starts_c.tolist()


def test_state_shift_features_and_report(tmp_path):
    adata_tracks = _build_adata_tracks()
    df_timepoints = _build_df_timepoints()
    adata_states = _build_adata_states()

    result = save_track_contact_state_shift_report(
        adata_tracks, df_timepoints, adata_states, tmp_path,
        contact_col="macro_contact",
        min_bout_length=_MIN_BOUT_LENGTH,
        state_col=_STATE_COL,
        window_mode="fixed",
        fixed_window_length=5,
        min_window_timepoints=3,
        verbose=False,
    )

    assert result["n_contact_tracks"] == 3
    assert result["n_no_contact_tracks"] == 3
    assert result["n_excluded_tracks"] == 0
    assert Path(result["pdf_path"]).exists()
    assert Path(result["track_windows_csv"]).exists()

    diff_csv = pd.read_csv(result["diff_bars_csv"])
    contact_rows = diff_csv[diff_csv["contact_group"] == "contact"].set_index("state")
    # Contact tracks: fully "static" before the bout, fully "engaging" after it.
    assert contact_rows.loc["static", "diff"] == pytest.approx(-1.0)
    assert contact_rows.loc["engaging", "diff"] == pytest.approx(1.0)
    assert contact_rows.loc["scanning", "diff"] == pytest.approx(0.0)

    no_contact_rows = diff_csv[diff_csv["contact_group"] == "no_contact"].set_index("state")
    # No-contact tracks: constant "static" the entire time -> no shift, regardless of where the
    # null reference split point landed.
    assert no_contact_rows.loc["static", "diff"] == pytest.approx(0.0)
    assert no_contact_rows.loc["engaging", "diff"] == pytest.approx(0.0)

    stacked_csv = pd.read_csv(result["stacked_composition_csv"])
    contact_before = stacked_csv[
        (stacked_csv["contact_group"] == "contact") & (stacked_csv["period"] == BEFORE_LABEL)
    ].set_index("state")
    assert contact_before.loc["static", "proportion"] == pytest.approx(1.0)


def test_track_missing_from_csv_raises():
    """A track present in adata_tracks.obs but entirely absent from df_timepoints (stale filtered
    CSV) must raise a hard error instead of being silently dropped from the output."""
    adata_tracks = _build_adata_tracks()
    df_timepoints = _build_df_timepoints()
    # Drop all rows for TrackID 0 from the CSV, simulating a re-run of filtering that removed it.
    df_timepoints = df_timepoints[df_timepoints["TrackID"] != 0]

    with pytest.raises(ValueError, match="have no matching timepoints in the filtered track-features CSV"):
        compute_contact_bout_windows(
            df_timepoints, adata_tracks, contact_col="macro_contact",
            min_bout_length=_MIN_BOUT_LENGTH, window_mode="fixed", fixed_window_length=5,
        )


def test_track_present_in_csv_but_no_timepoints_in_window_raises():
    """A track present in df_timepoints but with zero timepoints inside its classified
    [position_t_min, position_t_max] window must also raise from compute_contact_bout_windows —
    previously this track would just silently vanish from the output (no left-merge-back here,
    unlike compute_track_contact_features). Since adata_tracks is built directly from this CSV,
    zero overlap can only mean the two are out of sync."""
    adata_tracks = _build_adata_tracks()
    df_timepoints = _build_df_timepoints()
    # Shift TrackID 0's rows in the CSV entirely outside its obs window [0, _N_TIMEPOINTS - 1].
    mask = df_timepoints["TrackID"] == 0
    df_timepoints = df_timepoints.copy()
    df_timepoints.loc[mask, "position_t"] = df_timepoints.loc[mask, "position_t"] + 1000

    with pytest.raises(ValueError, match="have no matching timepoints in the filtered track-features CSV"):
        compute_contact_bout_windows(
            df_timepoints, adata_tracks, contact_col="macro_contact",
            min_bout_length=_MIN_BOUT_LENGTH, window_mode="fixed", fixed_window_length=5,
        )


def test_usable_track_missing_from_states_adata_raises():
    """A non-excluded track present in the CSV and adata_tracks, but with zero rows at all in
    adata_states.obs (e.g. the states h5ad wasn't regenerated for this track), must raise a hard
    error instead of silently contributing zero before/after rows."""
    adata_tracks = _build_adata_tracks()
    df_timepoints = _build_df_timepoints()
    adata_states = _build_adata_states()
    # Drop all rows for TrackID 0 from the states h5ad only — CSV/adata_tracks stay consistent.
    adata_states = adata_states[adata_states.obs["TrackID"] != 0].copy()

    with pytest.raises(ValueError, match="no matching timepoints at all in adata_states.obs"):
        compute_track_state_shift_features(
            df_timepoints, adata_tracks, adata_states,
            contact_col="macro_contact", min_bout_length=_MIN_BOUT_LENGTH, state_col=_STATE_COL,
            window_mode="fixed", fixed_window_length=5,
        )


def test_extra_track_in_states_adata_not_in_csv_raises():
    """A track present in adata_states.obs but entirely absent from the filtered CSV (e.g. the
    states h5ad wasn't regenerated after filtering dropped that track) must raise a hard error.
    adata_tracks/df_timepoints are kept mutually consistent (both drop TrackID 0) so this
    specifically exercises the adata_states-vs-CSV check, not the adata_tracks-vs-CSV one."""
    adata_tracks = _build_adata_tracks()
    df_timepoints = _build_df_timepoints()
    adata_states = _build_adata_states()
    adata_tracks = adata_tracks[adata_tracks.obs["TrackID"] != 0].copy()
    df_timepoints = df_timepoints[df_timepoints["TrackID"] != 0]

    with pytest.raises(ValueError, match="missing entirely from the filtered track-features CSV"):
        compute_track_state_shift_features(
            df_timepoints, adata_tracks, adata_states,
            contact_col="macro_contact", min_bout_length=_MIN_BOUT_LENGTH, state_col=_STATE_COL,
            window_mode="fixed", fixed_window_length=5,
        )


def test_dtype_drift_across_all_three_sources():
    """TrackID as int in adata_tracks.obs vs. '5.0'-style string in a re-read CSV vs. int in the
    states adata — the join must still succeed for all three sources."""
    adata_tracks = _build_adata_tracks(track_id_dtype=int)
    df_timepoints = _build_df_timepoints(track_id_dtype=lambda x: f"{float(x)}")
    adata_states = _build_adata_states(track_id_dtype=int)

    features = compute_track_state_shift_features(
        df_timepoints, adata_tracks, adata_states,
        contact_col="macro_contact", min_bout_length=_MIN_BOUT_LENGTH, state_col=_STATE_COL,
        window_mode="fixed", fixed_window_length=5,
    )
    assert features["n_contact_tracks"] == 3
    assert features["n_no_contact_tracks"] == 3
    assert len(features["state_timepoints"]) > 0
