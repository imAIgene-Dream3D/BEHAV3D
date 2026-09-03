"""
Regression tests for the AUC (trapezoidal-integration) calculations that were
touched by the ``np.trapz`` -> ``np.trapezoid`` migration (NumPy 2.4 removed
``np.trapz`` entirely, see [[state_composition]]). ``np.trapezoid`` is meant to
be a drop-in replacement, so these tests check actual numeric correctness
against hand-computed trapezoidal-rule values -- not just "it doesn't crash".
"""
import numpy as np
import pandas as pd
import pytest

from behav3d.analysis.invasiveness_analysis import _collapse_series
from behav3d.analysis.behavior.state.visualization.plots.state_composition import (
    _compute_relative_auc_table,
)


class TestCollapseSeriesAuc:
    """``_collapse_series(..., summary_stat="auc")`` -- per-sample AUC/span used
    as one of the invasiveness-analysis summary stats."""

    def test_basic_trapezoid(self):
        times = np.array([0, 1, 2, 3], dtype=float)
        values = np.array([0, 2, 2, 0], dtype=float)
        # trapezoid areas: (0+2)/2*1 + (2+2)/2*1 + (2+0)/2*1 = 1 + 2 + 1 = 4
        # span = 3 - 0 = 3 -> result = 4/3
        result = _collapse_series(times, values, "auc")
        assert result == pytest.approx(4.0 / 3.0)

    def test_unsorted_times_and_nan_masking(self):
        times = np.array([2, 0, 3, 1], dtype=float)
        values = np.array([2, 0, 0, np.nan], dtype=float)
        # sorted by time: t=[0,1,2,3], y=[0,nan,2,0]; NaN dropped -> t=[0,2,3], y=[0,2,0]
        # trapezoid: (0+2)/2*2 + (2+0)/2*1 = 2 + 1 = 3; span = 3 - 0 = 3 -> result = 1.0
        result = _collapse_series(times, values, "auc")
        assert result == pytest.approx(1.0)

    def test_single_valid_point_after_masking_returns_that_value(self):
        times = np.array([0, 1], dtype=float)
        values = np.array([5, np.nan], dtype=float)
        # only one non-NaN point remains -> short-circuits before any trapezoid call
        result = _collapse_series(times, values, "auc")
        assert result == pytest.approx(5.0)

    def test_zero_span_falls_back_to_mean(self):
        times = np.array([1, 1], dtype=float)
        values = np.array([3, 7], dtype=float)
        # span == 0 -> short-circuits to nanmean instead of dividing by zero
        result = _collapse_series(times, values, "auc")
        assert result == pytest.approx(5.0)


class TestComputeRelativeAucTable:
    """``_compute_relative_auc_table`` -- per-sample/per-state trapezoidal AUC
    of relative-composition curves, used by the state-composition report."""

    @staticmethod
    def _sample_matrix():
        return pd.DataFrame(
            {"stateA": [0.0, 2.0, 2.0, 0.0], "stateB": [1.0, 1.0, 1.0, 1.0]},
            index=pd.Index([0, 1, 2, 3], name="position_t"),
        )

    @staticmethod
    def _row(table, sample, state):
        match = table[(table["sample_name"] == sample) & (table["state_id"] == state)]
        assert len(match) == 1, f"expected exactly one row for ({sample}, {state})"
        return match.iloc[0]

    def test_auc_values_missing_state_and_single_timepoint_guard(self):
        mat = self._sample_matrix()
        single_point = pd.DataFrame(
            {"stateA": [5.0], "stateB": [2.0]}, index=pd.Index([0], name="position_t")
        )

        table = _compute_relative_auc_table(
            {"sample1": mat, "sample2": single_point},
            state_order=["stateA", "stateB", "stateC"],
            include_pooled_summary=False,
        )

        # stateA: trapezoid area over [0,2,2,0] at x=[0,1,2,3] -> 1+2+1 = 4
        assert self._row(table, "sample1", "stateA")["auc_relative"] == pytest.approx(4.0)
        # stateB: constant height-1 rectangle, width 3 -> 3.0
        assert self._row(table, "sample1", "stateB")["auc_relative"] == pytest.approx(3.0)
        # stateC absent from mat.columns -> treated as an all-zero curve -> 0.0
        assert self._row(table, "sample1", "stateC")["auc_relative"] == pytest.approx(0.0)
        # sample2 has only 1 timepoint -> the n_timepoints >= 2 guard forces 0.0
        assert self._row(table, "sample2", "stateA")["auc_relative"] == pytest.approx(0.0)
        assert self._row(table, "sample2", "stateA")["n_timepoints"] == 1

    def test_pooled_summary_row_uses_relative_pooled_matrix(self):
        mat = self._sample_matrix()
        table = _compute_relative_auc_table(
            {"sample1": mat},
            state_order=["stateA"],
            include_pooled_summary=True,
            relative_pooled=mat,
            pooled_sample_name="__all__",
        )
        pooled = self._row(table, "__all__", "stateA")
        assert pooled["auc_relative"] == pytest.approx(4.0)
