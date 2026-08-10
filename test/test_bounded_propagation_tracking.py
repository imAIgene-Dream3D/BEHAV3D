from pathlib import Path
import shutil
import uuid

import numpy as np

from behav3d.io.formats.zarr import save_as_zarr
from behav3d.io.images import load_zarr
from behav3d.preprocessing.tracking.bounded_propagation_tracking import (
    propagate_tracks_bounded,
    _prune_markers_to_best_region,
    _label_leftover_regions,
)


def _make_case_dir(name):
    root = Path(__file__).resolve().parent / ".tmp_bounded_propagation_tracking"
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


def _run(case_dir, t_seg, **kwargs):
    """Write ``t_seg`` (T,Z,Y,X), run propagate_tracks_bounded, return frame 1."""
    in_path = _write_zarr(case_dir / "segments.zarr", t_seg)
    tracked_img = case_dir / "tracked.zarr"
    tracked_csv = case_dir / "tracked.csv"
    kwargs.setdefault("dilation_nr_pixels", 1)
    kwargs.setdefault("segment_size_min", 1)
    kwargs.setdefault("min_overlap_fraction", 0.0)
    propagate_tracks_bounded(
        segments_path=in_path,
        tracked_img_outpath=tracked_img,
        tracked_csv_outpath=tracked_csv,
        **kwargs,
    )
    result = np.asarray(load_zarr(tracked_img))
    return result[1]


# ---------------------------------------------------------------------------
# End-to-end tests via propagate_tracks_bounded
# ---------------------------------------------------------------------------

def test_motivating_scenario_split_marker_leftover_becomes_one_new_track():
    """comp_A -> A, comp_B -> B, and the shared/merged region -> ONE new
    track (not stolen by A/B, not split into two new tracks)."""
    case_dir = _make_case_dir("motivating")
    try:
        t_seg = np.zeros((2, 3, 10, 45), dtype=np.uint16)

        # t=0: previous tracks A (1) and B (2), each reaching partway toward
        # the middle so the merged region at t=1 overlaps both.
        t_seg[0, :, 2:8, 0:16] = 1
        t_seg[0, :, 2:8, 24:40] = 2

        # t=1: three disconnected raw-mask regions separated by gap>=4
        # (safe against the dilation_nr_pixels=1 gap-closing pass).
        t_seg[1, :, 2:8, 0:10] = 9    # comp_A: only overlaps A (overlap 180 vs 36)
        t_seg[1, :, 2:8, 14:26] = 9   # comp_AB: overlaps both A and B, but weakly
        t_seg[1, :, 2:8, 30:40] = 9   # comp_B: only overlaps B (overlap 180 vs 36)

        frame1 = _run(case_dir, t_seg)

        assert np.all(frame1[:, 2:8, 0:10] == 1)
        assert np.all(frame1[:, 2:8, 30:40] == 2)

        ab_labels = set(np.unique(frame1[:, 2:8, 14:26]).tolist())
        assert ab_labels == {3}, f"expected exactly one new track (3), got {ab_labels}"

        assert set(np.unique(frame1).tolist()) == {0, 1, 2, 3}
    finally:
        _cleanup_case_dir(case_dir)


def test_single_track_continuation_unaffected():
    """One previous track, one overlapping current region: whole region
    inherits the track's label unchanged, no new track created."""
    case_dir = _make_case_dir("continuation")
    try:
        t_seg = np.zeros((2, 3, 10, 20), dtype=np.uint16)
        t_seg[0, :, 2:8, 5:15] = 1
        t_seg[1, :, 2:8, 6:16] = 5  # slight shift, still one connected blob

        frame1 = _run(case_dir, t_seg)

        mask1 = t_seg[1] != 0
        assert np.all(frame1[mask1] == 1)
        assert np.array_equal(frame1 > 0, mask1)
        assert set(np.unique(frame1).tolist()) == {0, 1}
    finally:
        _cleanup_case_dir(case_dir)


def test_shared_best_region_still_splits_via_watershed():
    """Both tracks' old footprints sit inside ONE current-frame blob (their
    only candidate region), so nothing gets pruned and watershed still
    splits it into the two tracks as before."""
    case_dir = _make_case_dir("shared_split")
    try:
        t_seg = np.zeros((2, 3, 10, 30), dtype=np.uint16)
        t_seg[0, :, 2:8, 0:10] = 1
        t_seg[0, :, 2:8, 15:25] = 2
        t_seg[1, :, 2:8, 0:25] = 9  # one connected blob spanning both old footprints

        frame1 = _run(case_dir, t_seg)

        labels = set(np.unique(frame1).tolist())
        assert labels == {0, 1, 2}, f"expected split preserved with no leftover, got {labels}"
        assert np.any(frame1 == 1) and np.any(frame1 == 2)
    finally:
        _cleanup_case_dir(case_dir)


def test_min_overlap_fraction_gates_whole_region_not_fragment():
    """Track's overlap fraction against the whole current region (not a
    post-split fragment) determines whether it keeps the region."""
    case_dir = _make_case_dir("min_overlap")
    try:
        t_seg = np.zeros((2, 3, 10, 20), dtype=np.uint16)
        t_seg[0, :, 2:8, 0:3] = 1     # footprint size 3*6*3=54
        t_seg[1, :, 2:8, 0:20] = 9    # one region, size 3*6*20=360; overlap frac = 54/360 = 0.15

        # 0.0 (default, passes): track 1 keeps the whole region.
        frame1_default = _run(case_dir, t_seg, min_overlap_fraction=0.0)
        assert set(np.unique(frame1_default).tolist()) == {0, 1}
    finally:
        _cleanup_case_dir(case_dir)

    case_dir = _make_case_dir("min_overlap_gated")
    try:
        t_seg = np.zeros((2, 3, 10, 20), dtype=np.uint16)
        t_seg[0, :, 2:8, 0:3] = 1
        t_seg[1, :, 2:8, 0:20] = 9

        # 0.2 (above the 0.15 actual fraction, fails): track 1 is pruned to
        # nothing everywhere; the whole region becomes one leftover -> a
        # brand-new track (next_label = max(seg_prev at t=0) + 1 = 2).
        frame1_gated = _run(case_dir, t_seg, min_overlap_fraction=0.2)
        assert 1 not in np.unique(frame1_gated)
        assert set(np.unique(frame1_gated).tolist()) == {0, 2}
    finally:
        _cleanup_case_dir(case_dir)


def test_2d_flat_leftover_is_removed_even_though_size_filter_would_pass():
    """A leftover component that's only 1 voxel thick along Z is removed by
    segment_2d_filter, even though its voxel count alone clears
    segment_size_min."""
    case_dir = _make_case_dir("flat_leftover")
    try:
        t_seg = np.zeros((2, 3, 10, 20), dtype=np.uint16)
        # No previous tracks at all -> everything at t=1 is unclaimed and
        # becomes a leftover. The blob spans only Z=1 (a single slice) but
        # is large in Y/X (6*20=120 voxels), comfortably above
        # segment_size_min=1.
        t_seg[1, 1:2, 2:8, 0:20] = 9

        frame1 = _run(case_dir, t_seg, segment_size_min=1)

        assert set(np.unique(frame1).tolist()) == {0}, "flat (2D) leftover should be filtered out"
    finally:
        _cleanup_case_dir(case_dir)


# ---------------------------------------------------------------------------
# Focused unit tests on the pruning / leftover-labelling helpers directly
# ---------------------------------------------------------------------------

def test_prune_markers_to_best_region_picks_argmax_per_track():
    # Two regions, track 1 overlaps region 1 more than region 2.
    regions = np.zeros((1, 10, 10), dtype=np.int64)
    regions[:, :, 0:4] = 1
    regions[:, :, 6:10] = 2

    seg_prev = np.zeros((1, 10, 10), dtype=np.uint16)
    seg_prev[:, :, 0:3] = 1  # mostly inside region 1 (3 px) vs none in region 2
    seg_prev[:, :, 6:7] = 1  # a sliver also touching region 2 (1 px)

    pruned = _prune_markers_to_best_region(regions, seg_prev, min_overlap_fraction=0.0)

    # Kept: track 1's original marker pixels that fall within its chosen
    # best region (region 1) -> exactly its X[0:3] slice, unchanged.
    assert np.array_equal(pruned == 1, (seg_prev == 1) & (regions == 1))
    assert not np.any(pruned[:, :, 6:7] == 1)  # sliver in region 2 pruned away


def test_prune_markers_to_best_region_gates_on_min_overlap_fraction():
    regions = np.ones((1, 10, 10), dtype=np.int64)  # single region, size 100
    seg_prev = np.zeros((1, 10, 10), dtype=np.uint16)
    seg_prev[:, :, 0:5] = 1  # overlap 50 -> fraction 0.5

    kept = _prune_markers_to_best_region(regions, seg_prev, min_overlap_fraction=0.4)
    assert np.any(kept == 1)

    dropped = _prune_markers_to_best_region(regions, seg_prev, min_overlap_fraction=0.6)
    assert not np.any(dropped == 1)


def test_label_leftover_regions_gives_one_id_per_disconnected_patch():
    mask = np.zeros((1, 10, 10), dtype=bool)
    mask[:, :, 0:3] = True
    mask[:, :, 6:9] = True  # disconnected from the first patch

    new_seg = np.zeros_like(mask, dtype=np.uint16)  # nothing claimed by watershed
    out, next_label = _label_leftover_regions(mask, new_seg, next_label=5)

    labels = set(np.unique(out).tolist()) - {0}
    assert labels == {5, 6}
    assert next_label == 7
