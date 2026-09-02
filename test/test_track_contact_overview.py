import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/behav3d-mpl-cache")

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg", force=True)

from behav3d.analysis.behavior.track.visualization.plots.contact_state_shift_report import (
    save_track_contact_overview_report,
    _compute_contact_bar_segments,
)

_N_TIMEPOINTS = 30
_MIN_BOUT_LENGTH = 5
_STATE_COL = "behavioral_state"

# (sample_name, TrackID, is_contact, bout_start) — bout occupies
# [bout_start, bout_start + _MIN_BOUT_LENGTH - 1] for contact tracks.
_SPECS = [
    ("sample_A", 0, True, 2),
    ("sample_A", 1, True, 12),
    ("sample_A", 2, True, 22),
    ("sample_B", 0, True, 2),
    ("sample_B", 1, True, 12),
    ("sample_B", 2, True, 22),
    ("sample_B", 3, False, None),
]


def _build_adata_tracks(specs):
    obs = pd.DataFrame(
        {
            "sample_name": [s[0] for s in specs],
            "TrackID": [s[1] for s in specs],
            "position_t_min": 0,
            "position_t_max": _N_TIMEPOINTS - 1,
        }
    )
    obs.index = [f"{s[0]}_{s[1]}" for s in specs]
    return ad.AnnData(X=np.zeros((len(specs), 1)), obs=obs)


def _build_df_timepoints(specs, contact_col="macro_contact"):
    rows = []
    for sample_name, track_id, is_contact, bout_start in specs:
        contact = np.zeros(_N_TIMEPOINTS, dtype=int)
        if is_contact:
            bout_end = bout_start + _MIN_BOUT_LENGTH - 1
            contact[bout_start : bout_end + 1] = 1
        for t in range(_N_TIMEPOINTS):
            rows.append(
                {
                    "sample_name": sample_name,
                    "TrackID": track_id,
                    "position_t": t,
                    contact_col: int(contact[t]),
                }
            )
    return pd.DataFrame(rows)


def _state_for_track(t, bout_start):
    if bout_start is None:
        return "static"
    bout_end = bout_start + _MIN_BOUT_LENGTH - 1
    if t < bout_start:
        return "static"
    if t <= bout_end:
        return "scanning"
    return "engaging"


def _build_adata_states(specs):
    rows = []
    for sample_name, track_id, _is_contact, bout_start in specs:
        for t in range(_N_TIMEPOINTS):
            rows.append(
                {
                    "sample_name": sample_name,
                    "TrackID": track_id,
                    "position_t": t,
                    _STATE_COL: _state_for_track(t, bout_start),
                }
            )
    obs = pd.DataFrame(rows)
    obs.index = [str(i) for i in range(len(obs))]
    return ad.AnnData(X=np.zeros((len(obs), 1)), obs=obs)


def test_compute_contact_bar_segments_only_marks_qualifying_bouts_green():
    times = list(range(20))
    # Qualifying bout (length 5, >= min_bout_length): [10..14]. Sub-threshold blip: [2..3] (len 2).
    is_contact = [0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0]
    segments = _compute_contact_bar_segments(times, is_contact, min_bout_length=5)

    green = [s for s in segments if s[2] == "#2ca02c"]
    grey = [s for s in segments if s[2] == "#bbbbbb"]
    assert len(green) == 1
    assert green[0][0] == 10.0
    assert green[0][1] == pytest.approx(5.0)
    # The short blip at [2,3] must render grey, same as true no-contact stretches.
    assert len(grey) == len(segments) - 1
    total_width = sum(s[1] for s in segments)
    assert total_width == pytest.approx(float(times[-1] - times[0] + 1))


def test_track_contact_overview_filters_to_contact_tracks_only(tmp_path):
    adata_tracks = _build_adata_tracks(_SPECS)
    df_timepoints = _build_df_timepoints(_SPECS)
    adata_states = _build_adata_states(_SPECS)

    result = save_track_contact_overview_report(
        adata_tracks, df_timepoints, adata_states, tmp_path,
        contact_col="macro_contact",
        min_bout_length=_MIN_BOUT_LENGTH,
        state_col=_STATE_COL,
        rows_per_page=6,
        verbose=False,
    )

    n_contact_tracks = sum(1 for s in _SPECS if s[2])
    assert result["n_tracks"] == n_contact_tracks
    assert result["n_samples"] == 2
    assert Path(result["pdf_path"]).exists()
    assert Path(result["pdf_path"]).name == "track_contact_overview.pdf"
    assert Path(result["pdf_path"]).parent.name == "macro_contact"


def test_track_contact_overview_no_contact_tracks_writes_placeholder(tmp_path):
    adata_tracks = _build_adata_tracks(_SPECS)
    df_timepoints = _build_df_timepoints(_SPECS)
    adata_states = _build_adata_states(_SPECS)

    # A threshold longer than any bout in the fixture means zero tracks qualify as "contact".
    result = save_track_contact_overview_report(
        adata_tracks, df_timepoints, adata_states, tmp_path,
        contact_col="macro_contact",
        min_bout_length=_N_TIMEPOINTS,
        state_col=_STATE_COL,
        verbose=False,
    )

    assert result["n_tracks"] == 0
    assert result["n_samples"] == 0
    assert Path(result["pdf_path"]).exists()


def test_track_contact_overview_pages_never_mix_samples(tmp_path):
    pypdf = pytest.importorskip("pypdf")

    adata_tracks = _build_adata_tracks(_SPECS)
    df_timepoints = _build_df_timepoints(_SPECS)
    adata_states = _build_adata_states(_SPECS)

    result = save_track_contact_overview_report(
        adata_tracks, df_timepoints, adata_states, tmp_path,
        contact_col="macro_contact",
        min_bout_length=_MIN_BOUT_LENGTH,
        state_col=_STATE_COL,
        rows_per_page=2,
        verbose=False,
    )

    # sample_A has 3 contact tracks, sample_B has 3 contact tracks; with rows_per_page=2 each
    # sample needs its own 2-page split (2 + 1). If pages were naively chunked across the flat
    # list of all 6 contact tracks instead of per-sample, this would collapse to exactly 3 pages
    # (2+2+2) with the last sample_A track sharing a page with the first sample_B track.
    assert result["n_tracks"] == 6
    assert result["n_samples"] == 2

    reader = pypdf.PdfReader(result["pdf_path"])
    assert len(reader.pages) == 4

    page_texts = [page.extract_text() or "" for page in reader.pages]
    expected_samples = ["sample_A", "sample_A", "sample_B", "sample_B"]
    for text, expected_sample in zip(page_texts, expected_samples):
        assert expected_sample in text
        other_sample = "sample_B" if expected_sample == "sample_A" else "sample_A"
        assert other_sample not in text
