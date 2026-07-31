"""
Tests for staleness detection across Feature Extraction -> Active Killing -> Filtering.

If Feature Extraction is rerun for a cell type, any previously-produced Active Killing
or Filtering output for that cell type was computed from track features that no longer
exist. See [[test_contact_group_analysis]] for the same "stale CSV" hard-error convention
already used for contact_grouping.py.
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


class TestRunActiveKillingAnalysisStaleness:
    def test_raises_when_immune_filtered_csv_predates_raw(self, tmp_path):
        now = time.time()
        raw = _raw_csv_path(tmp_path, "tcell")
        filtered = _filtered_csv_path(tmp_path, "tcell")
        _touch(filtered, now - 100)
        _touch(raw, now)

        with pytest.raises(StaleDataError, match="Re-run Filtering"):
            run_active_killing_analysis(
                metadata=pd.DataFrame(),
                output_dir=tmp_path,
                immune_cell_type="tcell",
                target_cell_types=["organoid"],
                save_results=False,
            )

    def test_raises_when_target_filtered_csv_predates_raw(self, tmp_path):
        now = time.time()
        # Immune tracks resolve cleanly (raw only, no filtered) so the loop reaches the target check.
        _touch(_raw_csv_path(tmp_path, "tcell"), now - 200)

        target_raw = _raw_csv_path(tmp_path, "organoid")
        target_filtered = _filtered_csv_path(tmp_path, "organoid")
        _touch(target_filtered, now - 100)
        _touch(target_raw, now)

        with pytest.raises(StaleDataError, match="Re-run Filtering"):
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
