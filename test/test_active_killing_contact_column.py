"""
Tests for the ``contact_column`` choice on Active Killing Analysis.

Active Killing's contact-event gate historically only ever looked at the pixel/mask-based
``{target}_contact`` column, silently ignoring the distance-based ``{target}_contact_on_distance``
column that Feature Extraction also produces from the user-configurable Contact Threshold (µm).
See [[test_active_killing_staleness]] for the broader Active Killing test conventions this file
follows.
"""
import os
import time

import pandas as pd
import pytest

from behav3d.features.advanced_timepoint_features import (
    identify_contact_events_global,
    run_active_killing_analysis,
)


def _immune_df():
    """Four timepoints of one immune track, always touching organoid TrackID 5, but with
    the pixel-based and distance-based contact flags diverging at t=2 and t=3 -- isolating
    the effect of ``contact_column`` to exactly those two timepoints."""
    return pd.DataFrame({
        "sample_name": ["s1"] * 4,
        "TrackID": [0] * 4,
        "position_t": [0, 1, 2, 3],
        "organoid_contact": [True, True, False, False],
        "organoid_contact_on_distance": [True, True, True, True],
        "touching_organoids": ["5", "5", "5", "5"],
    })


class TestIdentifyContactEventsGlobalContactColumn:
    def test_default_uses_pixel_contact(self):
        events = identify_contact_events_global(_immune_df(), ["organoid"], min_contact_duration=1)

        assert len(events) == 1
        assert events.iloc[0]["contact_duration"] == 2
        assert events.iloc[0]["contact_end_t"] == 1

    def test_contact_on_distance_respects_distance_threshold(self):
        events = identify_contact_events_global(
            _immune_df(), ["organoid"], min_contact_duration=1, contact_column="contact_on_distance",
        )

        assert len(events) == 1
        assert events.iloc[0]["contact_duration"] == 4
        assert events.iloc[0]["contact_end_t"] == 3

    def test_invalid_contact_column_raises(self):
        with pytest.raises(ValueError, match="contact_column"):
            identify_contact_events_global(_immune_df(), ["organoid"], contact_column="bogus")

    def test_missing_column_returns_empty_with_warning(self, capsys):
        df = _immune_df().drop(columns=["organoid_contact_on_distance"])

        events = identify_contact_events_global(df, ["organoid"], contact_column="contact_on_distance")

        assert events.empty
        out = capsys.readouterr().out
        assert "No 'contact_on_distance' columns found" in out


def _touch_immune_csv(path, mtime=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    _immune_df().to_csv(path, index=False)
    if mtime is not None:
        os.utime(path, (mtime, mtime))


def _touch_organoid_csv(path, mtime=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({
        "sample_name": ["s1"] * 4,
        "TrackID": [5] * 4,
        "position_t": [0, 1, 2, 3],
        "mean_dead_dye": [10.0] * 4,
    }).to_csv(path, index=False)
    if mtime is not None:
        os.utime(path, (mtime, mtime))


class TestRunActiveKillingAnalysisContactColumn:
    """End-to-end: the panel-facing entry point must actually forward contact_column
    through to identify_contact_events_global rather than only accepting it."""

    def test_default_contact_column_is_pixel_based(self, tmp_path, capsys):
        now = time.time()
        _touch_immune_csv(
            tmp_path / "analysis" / "tcell" / "track_features" / "BEHAV3D_tcell_combined_track_features.csv", now,
        )
        _touch_organoid_csv(
            tmp_path / "analysis" / "organoid" / "track_features" / "BEHAV3D_organoid_combined_track_features.csv", now,
        )

        _, _, stats = run_active_killing_analysis(
            metadata=pd.DataFrame(),
            output_dir=tmp_path,
            immune_cell_type="tcell",
            target_cell_types=["organoid"],
            save_results=False,
        )

        assert stats["contact_column"] == "contact"
        assert stats["total_contact_timepoints"] == 2
        out = capsys.readouterr().out
        assert "Contact column: contact" in out

    def test_contact_on_distance_widens_the_contact_window(self, tmp_path):
        now = time.time()
        _touch_immune_csv(
            tmp_path / "analysis" / "tcell" / "track_features" / "BEHAV3D_tcell_combined_track_features.csv", now,
        )
        _touch_organoid_csv(
            tmp_path / "analysis" / "organoid" / "track_features" / "BEHAV3D_organoid_combined_track_features.csv", now,
        )

        _, _, stats = run_active_killing_analysis(
            metadata=pd.DataFrame(),
            output_dir=tmp_path,
            immune_cell_type="tcell",
            target_cell_types=["organoid"],
            contact_column="contact_on_distance",
            save_results=False,
        )

        assert stats["contact_column"] == "contact_on_distance"
        assert stats["total_contact_timepoints"] == 4
