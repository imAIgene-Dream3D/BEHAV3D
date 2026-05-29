from pathlib import Path
import shutil
import uuid

import numpy as np

from behav3d.io.formats.zarr import save_as_zarr
from behav3d.io.images import load_zarr
from behav3d.preprocessing.segmentation.multicolor_segment_processing import (
    _erode_candidate,
    _normalize_segmentation_labels_for_view,
    calculate_multicolor_overlap,
    process_multicolor_segments,
)


def _make_case_dir(name):
    root = Path(__file__).resolve().parent / ".tmp_multicolor_segment_processing"
    root.mkdir(exist_ok=True)
    case_dir = root / f"{name}_{uuid.uuid4().hex}"
    case_dir.mkdir()
    return case_dir


def _cleanup_case_dir(path):
    if path.exists():
        shutil.rmtree(path)


def _write_zarr(path, arr):
    save_as_zarr(arr.astype(np.uint16), path)
    return Path(path)


def _make_empty(shape=(1, 5, 8, 8)):
    return np.zeros(shape, dtype=np.uint16)


def test_full_containment_removes_loser():
    case_dir = _make_case_dir("full_containment")
    try:
        img1 = _make_empty()
        img2 = _make_empty()

        img1[:, 1:4, 1:6, 1:6] = 1
        img2[:, 1:4, 2:4, 2:4] = 7

        in1 = _write_zarr(case_dir / "img1_segments.zarr", img1)
        in2 = _write_zarr(case_dir / "img2_segments.zarr", img2)

        out1, out2 = process_multicolor_segments([in1, in2], erosion_pixels=1)

        overlap_stats = calculate_multicolor_overlap([out1, out2])
        assert np.array_equal(np.asarray(load_zarr(out1)), img1)
        assert np.count_nonzero(np.asarray(load_zarr(out2))) == 0
        assert overlap_stats["overlap_any_count"] == 0
    finally:
        _cleanup_case_dir(case_dir)


def test_partial_overlap_single_survivor_gets_union():
    case_dir = _make_case_dir("single_survivor")
    try:
        img1 = _make_empty()
        img2 = _make_empty()

        img1[:, 1:4, 1:6, 1:6] = 11
        img2[:, 1:4, 3:7, 3:6] = 22

        in1 = _write_zarr(case_dir / "img1_segments.zarr", img1)
        in2 = _write_zarr(case_dir / "img2_segments.zarr", img2)

        out1, out2 = process_multicolor_segments([in1, in2], erosion_pixels=1)

        res1 = np.asarray(load_zarr(out1))
        res2 = np.asarray(load_zarr(out2))
        expected_union = (img1 > 0) | (img2 > 0)

        assert np.array_equal(res1 > 0, expected_union)
        assert np.count_nonzero(res2) == 0
    finally:
        _cleanup_case_dir(case_dir)


def test_partial_overlap_multiple_survivors_repartition_without_overlap():
    case_dir = _make_case_dir("multiple_survivors")
    try:
        img1 = _make_empty()
        img2 = _make_empty()

        img1[:, 1:4, 1:5, 1:5] = 3
        img2[:, 1:4, 3:7, 3:7] = 4

        in1 = _write_zarr(case_dir / "img1_segments.zarr", img1)
        in2 = _write_zarr(case_dir / "img2_segments.zarr", img2)

        out1, out2 = process_multicolor_segments([in1, in2], erosion_pixels=0)

        res1 = np.asarray(load_zarr(out1))
        res2 = np.asarray(load_zarr(out2))
        processed_union = (res1 > 0) | (res2 > 0)
        original_union = (img1 > 0) | (img2 > 0)

        assert not np.any((res1 > 0) & (res2 > 0))
        assert np.array_equal(processed_union, original_union)
    finally:
        _cleanup_case_dir(case_dir)


def test_seed_erosion_does_not_apply_min_size_filter():
    case_dir = _make_case_dir("seed_stage_no_size_filter")
    try:
        img1 = _make_empty(shape=(1, 8, 10, 10))
        img2 = _make_empty(shape=(1, 8, 10, 10))

        img1[:, 1:7, 1:8, 1:8] = 5
        img2[:, 1:7, 4:9, 4:9] = 6

        shared = (img1[0] > 0) & (img2[0] > 0)
        seed1 = _erode_candidate((img1[0] > 0) & ~shared, 1)
        seed2 = _erode_candidate((img2[0] > 0) & ~shared, 1)

        seed1_size = int(np.count_nonzero(seed1))
        seed2_size = int(np.count_nonzero(seed2))
        assert seed1_size > seed2_size > 0

        min_size = seed2_size + 1
        assert int(np.count_nonzero(img2[0])) > min_size

        in1 = _write_zarr(case_dir / "img1_segments.zarr", img1)
        in2 = _write_zarr(case_dir / "img2_segments.zarr", img2)

        out1, out2 = process_multicolor_segments(
            [in1, in2],
            erosion_pixels=1,
            min_size=min_size,
        )

        res1 = np.asarray(load_zarr(out1))
        res2 = np.asarray(load_zarr(out2))
        processed_union = (res1 > 0) | (res2 > 0)
        original_union = (img1 > 0) | (img2 > 0)

        assert np.count_nonzero(res1) > 0
        assert np.count_nonzero(res2) > 0
        assert not np.any((res1 > 0) & (res2 > 0))
        assert np.array_equal(processed_union, original_union)
    finally:
        _cleanup_case_dir(case_dir)


def test_final_size_filter_removes_small_resolved_segment():
    case_dir = _make_case_dir("final_size_filter")
    try:
        img1 = _make_empty(shape=(1, 6, 8, 8))
        img2 = _make_empty(shape=(1, 6, 8, 8))

        img1[:, 1:5, 1:6, 1:6] = 10
        img2[:, 1:5, 4:6, 4:6] = 20

        min_size = int(np.count_nonzero(img2[0])) + 1

        in1 = _write_zarr(case_dir / "img1_segments.zarr", img1)
        in2 = _write_zarr(case_dir / "img2_segments.zarr", img2)

        out1, out2 = process_multicolor_segments(
            [in1, in2],
            erosion_pixels=0,
            min_size=min_size,
        )

        res1 = np.asarray(load_zarr(out1))
        res2 = np.asarray(load_zarr(out2))

        assert np.count_nonzero(res1) > 0
        assert np.count_nonzero(res2) == 0
    finally:
        _cleanup_case_dir(case_dir)


def test_calculate_overlap_returns_pairwise_and_global_counts():
    case_dir = _make_case_dir("overlap_stats")
    try:
        a = _make_empty(shape=(2, 2, 4, 4))
        b = _make_empty(shape=(2, 2, 4, 4))
        c = _make_empty(shape=(2, 2, 4, 4))
        d = _make_empty(shape=(2, 2, 4, 4))

        a[:, :, 0:2, 0:2] = 1
        b[:, :, 1:3, 1:3] = 2
        c[:, :, 1:2, 1:2] = 3
        d[:, :, 3:4, 3:4] = 4

        paths = [
            _write_zarr(case_dir / "a.zarr", a),
            _write_zarr(case_dir / "b.zarr", b),
            _write_zarr(case_dir / "c.zarr", c),
            _write_zarr(case_dir / "d.zarr", d),
        ]

        stats = calculate_multicolor_overlap(paths)

        assert stats["pairwise_overlap_matrix"].shape == (4, 4)
        assert stats["pairwise_overlap_by_timepoint"].shape == (2, 4, 4)
        assert stats["pairwise_overlap_matrix"][0, 1] == 4
        assert stats["pairwise_overlap_matrix"][0, 2] == 4
        assert stats["pairwise_overlap_matrix"][1, 2] == 4
        assert stats["pairwise_overlap_matrix"][0, 3] == 0
        assert stats["overlap_all_count"] == 0
        assert stats["overlap_any_count"] == int(np.count_nonzero(stats["overlap_mask"]))
        named = {item["name"]: item for item in stats["pairwise_named_stats"]}
        assert named["a-b"]["overlapping_count"] == 4
        assert named["a-b"]["non_overlapping_count"] == 24
    finally:
        _cleanup_case_dir(case_dir)


def test_normalize_labels_for_view_accepts_singleton_channel_axis():
    raw_shape = (2, 3, 4, 5, 6)
    labels = np.zeros((2, 1, 4, 5, 6), dtype=np.uint8)

    normalized = _normalize_segmentation_labels_for_view(labels, raw_shape, layer_name="demo")

    assert normalized.shape == (2, 4, 5, 6)


def test_process_supports_sibling_and_in_place_modes():
    case_dir = _make_case_dir("output_modes")
    try:
        img1 = _make_empty()
        img2 = _make_empty()
        img1[:, 1:4, 1:6, 1:6] = 1
        img2[:, 1:4, 2:4, 2:4] = 2

        sibling_in1 = _write_zarr(case_dir / "sibling1.zarr", img1)
        sibling_in2 = _write_zarr(case_dir / "sibling2.zarr", img2)
        sibling_outs = process_multicolor_segments([sibling_in1, sibling_in2], erosion_pixels=1)
        assert sibling_outs[0].name == "sibling1_multicolor_processed.zarr"
        assert sibling_outs[1].exists()
        assert sibling_in1.exists()

        inplace_in1 = _write_zarr(case_dir / "inplace1.zarr", img1)
        inplace_in2 = _write_zarr(case_dir / "inplace2.zarr", img2)
        inplace_outs = process_multicolor_segments(
            [inplace_in1, inplace_in2],
            overwrite=True,
            erosion_pixels=1,
        )
        assert inplace_outs == [inplace_in1, inplace_in2]
        assert not (case_dir / "inplace1_multicolor_processed.zarr").exists()
    finally:
        _cleanup_case_dir(case_dir)
