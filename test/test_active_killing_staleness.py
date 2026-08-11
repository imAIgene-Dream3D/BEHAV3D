"""
Tests for staleness detection across Feature Extraction -> Active Killing -> Filtering.

If Feature Extraction is rerun for a cell type, any previously-produced Active Killing
or Filtering output for that cell type was computed from track features that no longer
exist. See [[test_contact_group_analysis]] for the same "stale CSV" hard-error convention
already used for contact_grouping.py.

Active Killing always reads the RAW (unfiltered) combined_track_features.csv, never the
filtered one -- see the comment in run_active_killing_analysis(). Filtering, in turn,
reads Active Killing's advanced-features CSV as its own input when one exists (see
find_advanced_features_csv / filter_tracks' df_input_path), so Active Killing preferring
the filtered CSV would make the two steps each other's upstream dependency: a staleness
check on either side could deadlock, with Filtering refusing to run because Active
Killing looks stale and Active Killing refusing to run because Filtering looks stale.
Always starting Active Killing from raw breaks that cycle; instead, Active Killing warns
the caller (print + stats["filtering_needs_rerun_for"]) that any existing filtered CSV is
now behind and Filtering should be re-run to pick up its output.
"""
import os
import time

import pandas as pd
import pytest

from behav3d.features.advanced_timepoint_features import (
    StaleDataError,
    find_advanced_features_csv,
    run_active_killing_analysis,
)
from behav3d.features.timepoint_features import run_feature_extraction


def _touch(path, mtime):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"TrackID": [0]}).to_csv(path, index=False)
    os.utime(path, (mtime, mtime))


def _touch_track_csv(path, mtime):
    """Like _touch, but with the minimal columns run_active_killing_analysis
    actually reads (sample_name, TrackID, position_t) instead of just TrackID."""
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"sample_name": ["s1"], "TrackID": [0], "position_t": [0]}).to_csv(path, index=False)
    os.utime(path, (mtime, mtime))


def _raw_csv_path(output_dir, cell_type):
    return output_dir / "analysis" / cell_type / "track_features" / f"BEHAV3D_{cell_type}_combined_track_features.csv"


def _filtered_csv_path(output_dir, cell_type):
    return output_dir / "analysis" / cell_type / "track_features" / f"BEHAV3D_{cell_type}_combined_track_features_filtered.csv"


def _advanced_csv_path(output_dir, cell_type, subfolder="organoid"):
    return output_dir / "analysis" / cell_type / "active_killing" / subfolder / f"BEHAV3D_{cell_type}_advanced_track_features.csv"


class TestFindAdvancedFeaturesCsv:
    def test_raises_when_advanced_csv_predates_raw(self, tmp_path):
        """Feature Extraction rerun after Active Killing must not let Filtering silently
        pick up the now-stale advanced-features CSV."""
        now = time.time()
        advanced = _advanced_csv_path(tmp_path, "tcell")
        raw = _raw_csv_path(tmp_path, "tcell")
        _touch(advanced, now - 100)
        _touch(raw, now)

        with pytest.raises(StaleDataError, match="Feature Extraction was rerun"):
            find_advanced_features_csv(tmp_path, "tcell")

    def test_returns_path_when_advanced_csv_is_newer(self, tmp_path):
        now = time.time()
        advanced = _advanced_csv_path(tmp_path, "tcell")
        raw = _raw_csv_path(tmp_path, "tcell")
        _touch(raw, now - 100)
        _touch(advanced, now)

        assert find_advanced_features_csv(tmp_path, "tcell") == advanced

    def test_returns_path_when_no_raw_csv_exists_yet(self, tmp_path):
        advanced = _advanced_csv_path(tmp_path, "tcell")
        _touch(advanced, time.time())

        assert find_advanced_features_csv(tmp_path, "tcell") == advanced

    def test_returns_none_when_no_active_killing_dir(self, tmp_path):
        assert find_advanced_features_csv(tmp_path, "tcell") is None


class TestRunActiveKillingAnalysisAlwaysUsesRaw:
    """Active Killing must never read the filtered CSV -- staleness of the
    filtered CSV relative to raw (in either direction) is irrelevant to it
    and must never raise or change which file gets loaded."""

    def test_immune_raw_used_even_when_filtered_predates_raw(self, tmp_path, capsys):
        now = time.time()
        raw = _raw_csv_path(tmp_path, "tcell")
        filtered = _filtered_csv_path(tmp_path, "tcell")
        _touch_track_csv(filtered, now - 100)
        _touch_track_csv(raw, now)
        _touch_track_csv(_raw_csv_path(tmp_path, "organoid"), now)

        run_active_killing_analysis(
            metadata=pd.DataFrame(),
            output_dir=tmp_path,
            immune_cell_type="tcell",
            target_cell_types=["organoid"],
            save_results=False,
        )

        out = capsys.readouterr().out
        assert f"Loading immune cell tracks from {raw}" in out

    def test_immune_raw_used_even_when_filtered_is_newer(self, tmp_path, capsys):
        now = time.time()
        raw = _raw_csv_path(tmp_path, "tcell")
        filtered = _filtered_csv_path(tmp_path, "tcell")
        _touch_track_csv(raw, now - 100)
        _touch_track_csv(filtered, now)
        _touch_track_csv(_raw_csv_path(tmp_path, "organoid"), now)

        run_active_killing_analysis(
            metadata=pd.DataFrame(),
            output_dir=tmp_path,
            immune_cell_type="tcell",
            target_cell_types=["organoid"],
            save_results=False,
        )

        out = capsys.readouterr().out
        assert f"Loading immune cell tracks from {raw}" in out

    def test_target_raw_used_regardless_of_filtered_mtime(self, tmp_path, capsys):
        now = time.time()
        _touch_track_csv(_raw_csv_path(tmp_path, "tcell"), now)
        target_raw = _raw_csv_path(tmp_path, "organoid")
        target_filtered = _filtered_csv_path(tmp_path, "organoid")
        _touch_track_csv(target_filtered, now - 100)
        _touch_track_csv(target_raw, now)

        run_active_killing_analysis(
            metadata=pd.DataFrame(),
            output_dir=tmp_path,
            immune_cell_type="tcell",
            target_cell_types=["organoid"],
            save_results=False,
        )

        out = capsys.readouterr().out
        assert f"Loading organoid tracks from {target_raw}" in out

    def test_immune_filenotfound_when_only_filtered_exists(self, tmp_path):
        """No raw CSV means Active Killing has nothing to read, even though a
        filtered CSV is sitting right there -- it is never used as a fallback."""
        _touch_track_csv(_filtered_csv_path(tmp_path, "tcell"), time.time())

        with pytest.raises(FileNotFoundError):
            run_active_killing_analysis(
                metadata=pd.DataFrame(),
                output_dir=tmp_path,
                immune_cell_type="tcell",
                target_cell_types=["organoid"],
                save_results=False,
            )

    def test_immune_filenotfound_when_no_tracks_exist(self, tmp_path):
        """Sanity check: absence of input data still raises the pre-existing
        FileNotFoundError, not StaleDataError - staleness only applies when a
        filtered file actually exists but is out of date."""
        with pytest.raises(FileNotFoundError):
            run_active_killing_analysis(
                metadata=pd.DataFrame(),
                output_dir=tmp_path,
                immune_cell_type="tcell",
                target_cell_types=["organoid"],
                save_results=False,
            )


class TestRunActiveKillingAnalysisRerunFilteringWarning:
    """After a run, Active Killing should tell the caller which cell types'
    filtered CSVs are now behind, both via a printed warning and via
    stats["filtering_needs_rerun_for"], so the GUI can surface a popup."""

    def test_warns_and_reports_stats_when_filtered_csv_exists(self, tmp_path, capsys):
        now = time.time()
        _touch_track_csv(_raw_csv_path(tmp_path, "tcell"), now)
        _touch_track_csv(_filtered_csv_path(tmp_path, "tcell"), now)
        _touch_track_csv(_raw_csv_path(tmp_path, "organoid"), now)

        # save_results=True: the warning is about a *freshly written* advanced
        # CSV outrunning the existing filtered CSV, so it only fires when this
        # run actually produced new output.
        _, _, stats = run_active_killing_analysis(
            metadata=pd.DataFrame(),
            output_dir=tmp_path,
            immune_cell_type="tcell",
            target_cell_types=["organoid"],
            save_results=True,
        )

        assert stats["filtering_needs_rerun_for"] == ["tcell"]
        out = capsys.readouterr().out
        assert "Re-run Filtering" in out
        assert "tcell" in out

    def test_no_warning_when_no_filtered_csv_exists(self, tmp_path, capsys):
        now = time.time()
        _touch_track_csv(_raw_csv_path(tmp_path, "tcell"), now)
        _touch_track_csv(_raw_csv_path(tmp_path, "organoid"), now)

        _, _, stats = run_active_killing_analysis(
            metadata=pd.DataFrame(),
            output_dir=tmp_path,
            immune_cell_type="tcell",
            target_cell_types=["organoid"],
            save_results=True,
        )

        assert stats["filtering_needs_rerun_for"] == []
        out = capsys.readouterr().out
        assert "Re-run Filtering" not in out

    def test_immune_filenotfound_when_no_tracks_exist(self, tmp_path):
        """Sanity check: absence of input data still raises the pre-existing
        FileNotFoundError, not StaleDataError - staleness only applies when a
        filtered file actually exists but is out of date."""
        with pytest.raises(FileNotFoundError):
            run_active_killing_analysis(
                metadata=pd.DataFrame(),
                output_dir=tmp_path,
                immune_cell_type="tcell",
                target_cell_types=["organoid"],
                save_results=False,
            )


class TestRunFeatureExtractionProactiveWarning:
    def test_warns_when_filtered_and_active_killing_outputs_already_exist(self, tmp_path, capsys):
        now = time.time()
        _touch(_filtered_csv_path(tmp_path, "tcell"), now - 100)
        _touch(_advanced_csv_path(tmp_path, "tcell"), now - 100)

        run_feature_extraction(metadata=pd.DataFrame(columns=["sample_name"]), output_dir=tmp_path, cell_type="tcell")

        out = capsys.readouterr().out
        assert "WARNING: Feature Extraction for 'tcell' was just rerun" in out
        assert "Filtering" in out
        assert "Active Killing" in out

    def test_no_warning_when_no_downstream_output_exists(self, tmp_path, capsys):
        run_feature_extraction(metadata=pd.DataFrame(columns=["sample_name"]), output_dir=tmp_path, cell_type="tcell")

        out = capsys.readouterr().out
        assert "WARNING: Feature Extraction" not in out
