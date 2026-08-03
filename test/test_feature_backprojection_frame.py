"""
Tests for `backproject_feature_at_timepoint`, the shared single-frame
feature-to-label mapping primitive used by both the quick napari-plugin
and notebook feature-backprojection previews (current-timepoint-only,
no full-experiment zarr written).
"""
import numpy as np
import pandas as pd
import pytest

from behav3d.analysis.backprojection import backproject_feature_at_timepoint


def test_basic_mapping_2d():
    labels = np.array([[0, 1, 1], [2, 2, 0]])
    df = pd.DataFrame({"TrackID": [1, 2], "elongation": [1.5, 2.5]})
    mapped, ids = backproject_feature_at_timepoint(labels, df, "elongation")
    assert ids.tolist() == [1, 2]
    np.testing.assert_allclose(mapped, [[0, 1.5, 1.5], [2.5, 2.5, 0]])


def test_3d_frame_shape_preserved():
    labels = np.zeros((2, 3, 3), dtype=int)
    labels[0, 0, 0] = 1
    labels[1, 1, 1] = 2
    df = pd.DataFrame({"TrackID": [1, 2], "nr_dead_mask_pixels": [0, 5]})
    mapped, ids = backproject_feature_at_timepoint(labels, df, "nr_dead_mask_pixels")
    assert mapped.shape == labels.shape
    assert mapped[0, 0, 0] == 0  # legitimate zero value, not background-confused
    assert mapped[1, 1, 1] == 5
    assert set(ids.tolist()) == {1, 2}


def test_missing_track_reported_via_ids():
    labels = np.array([[1, 3]])  # label 3 has no feature row
    df = pd.DataFrame({"TrackID": [1], "elongation": [2.0]})
    mapped, ids = backproject_feature_at_timepoint(labels, df, "elongation")
    assert 3 not in ids.tolist()
    assert 1 in ids.tolist()


def test_original_trackid_preferred():
    labels = np.array([[10, 20]])  # segmentation label ids
    df = pd.DataFrame({
        "TrackID": [1, 2],             # post-split renumbered ids
        "original_TrackID": [10, 20],  # pixel-matching ids
        "elongation": [1.0, 2.0],
    })
    mapped, ids = backproject_feature_at_timepoint(labels, df, "elongation")
    np.testing.assert_allclose(mapped, [[1.0, 2.0]])
    assert set(ids.tolist()) == {10, 20}


def test_non_numeric_column_raises():
    labels = np.array([[1]])
    df = pd.DataFrame({"TrackID": [1], "category": ["alive"]})
    with pytest.raises(ValueError):
        backproject_feature_at_timepoint(labels, df, "category")


def test_missing_feature_column_raises():
    labels = np.array([[1]])
    df = pd.DataFrame({"TrackID": [1]})
    with pytest.raises(ValueError):
        backproject_feature_at_timepoint(labels, df, "nope")


def test_missing_track_column_raises():
    labels = np.array([[1]])
    df = pd.DataFrame({"elongation": [1.0]})
    with pytest.raises(ValueError):
        backproject_feature_at_timepoint(labels, df, "elongation")


def test_empty_dataframe_returns_background():
    labels = np.array([[1, 2], [0, 1]])
    df = pd.DataFrame({"TrackID": pd.Series([], dtype=int), "elongation": pd.Series([], dtype=float)})
    mapped, ids = backproject_feature_at_timepoint(labels, df, "elongation")
    assert ids.size == 0
    np.testing.assert_array_equal(mapped, np.zeros_like(labels))
