"""Correctness check for multi-worker ``create_group_tracked_segments``.

Builds a tiny synthetic "group" (two member cell types, two samples) and
asserts that running the build with ``n_workers=1`` (each sample's
timepoints processed serially) and ``n_workers=2`` (each sample's
timepoints split across a ``ThreadPoolExecutor``, writing to fixed zarr
indices out of order) produce byte-identical tracked images and tracks
CSVs — i.e. per-timepoint parallel dispatch is a correctness-preserving
refactor of the original sequential-append loop, not just "doesn't crash".
"""
from pathlib import Path
import shutil
import uuid

import numpy as np
import pandas as pd

from behav3d.io.formats.zarr import save_as_zarr
from behav3d.io.images import load_image
from behav3d.analysis.grouping import (
    create_group_tracked_segments,
    group_tracked_csv_path,
    merged_track_features_csv,
)

SAMPLES = ["sample1", "sample2"]
MEMBERS = ["tcell", "nkcell"]  # both immune-category -> "im_" prefix columns
T = 3
SHAPE = (T, 1, 4, 4)  # T, Z, Y, X

# Disjoint regions per member so combining them never overwrites a pixel
# written by another member.
REGIONS = {
    "tcell": (slice(0, 2), slice(0, 2)),
    "nkcell": (slice(2, 4), slice(2, 4)),
}
# Each member has two tracks (label ids 1 and 2) mapped to distinct group
# TrackIDs so the merged output never collides.
GROUP_TRACKID = {
    "tcell": {1: 1, 2: 2},
    "nkcell": {1: 3, 2: 4},
}


def _case_dir(name):
    root = Path(__file__).resolve().parent / ".tmp_group_tracked_segments_parallel"
    root.mkdir(exist_ok=True)
    case_dir = root / f"{name}_{uuid.uuid4().hex}"
    case_dir.mkdir()
    return case_dir


def _cleanup(path):
    if path.exists():
        shutil.rmtree(path)


def _make_member_image(cell_type, sample_idx):
    """Two disjoint labels (1, 2) at a fixed spatial region, constant over T.

    ``sample_idx`` offsets the region within its half so the two samples'
    images differ (not strictly required, but avoids an accidentally
    trivial test where every sample is pixel-identical).
    """
    arr = np.zeros(SHAPE, dtype=np.uint16)
    ry, rx = REGIONS[cell_type]
    y0, x0 = ry.start, rx.start
    # label 1 in the left/top half of the region, label 2 in the other half
    if sample_idx == 0:
        arr[:, 0, y0, x0] = 1
        arr[:, 0, y0, x0 + 1] = 2
    else:
        arr[:, 0, y0 + 1, x0] = 1
        arr[:, 0, y0 + 1, x0 + 1] = 2
    return arr


def _build_inputs(case_dir):
    """Write member tracked zarrs + metadata + a merged CSV. Returns
    (metadata_df, group_id)."""
    group_id = "immune_merged"
    member_paths = {s: {} for s in SAMPLES}

    for sidx, sample_name in enumerate(SAMPLES):
        sample_dir = case_dir / "images" / sample_name
        sample_dir.mkdir(parents=True)
        for cell_type in MEMBERS:
            arr = _make_member_image(cell_type, sidx)
            path = sample_dir / f"{sample_name}_{cell_type}_tracked.zarr"
            save_as_zarr(arr, path)
            member_paths[sample_name][cell_type] = path

    metadata_rows = []
    for sample_name in SAMPLES:
        row = {
            "sample_name": sample_name,
            "pixel_distance_xy": 1,
            "pixel_distance_z": 1,
        }
        for cell_type in MEMBERS:
            row[f"im_{cell_type}_tracks_image_path"] = str(member_paths[sample_name][cell_type])
        metadata_rows.append(row)
    metadata = pd.DataFrame(metadata_rows)

    merged_rows = []
    for sample_name in SAMPLES:
        for cell_type in MEMBERS:
            for origin_track_id, group_track_id in GROUP_TRACKID[cell_type].items():
                for t in range(T):
                    merged_rows.append({
                        "TrackID": group_track_id,
                        "origin_cell_type": cell_type,
                        "origin_TrackID": origin_track_id,
                        "sample_name": sample_name,
                        "position_t": t,
                    })
    merged_df = pd.DataFrame(merged_rows)

    out_csv = merged_track_features_csv(case_dir, group_id)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    merged_df.to_csv(out_csv, index=False)

    return metadata, group_id


def test_parallel_matches_serial():
    case_dir = _case_dir("group_tracked_segments")
    try:
        metadata, group_id = _build_inputs(case_dir)

        out_serial = case_dir / "out_serial"
        out_parallel = case_dir / "out_parallel"
        for out_dir in (out_serial, out_parallel):
            merged_src = merged_track_features_csv(case_dir, group_id)
            merged_dst = merged_track_features_csv(out_dir, group_id)
            merged_dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(merged_src, merged_dst)

        written_serial = create_group_tracked_segments(
            out_serial, group_id, metadata, n_workers=1
        )
        written_parallel = create_group_tracked_segments(
            out_parallel, group_id, metadata, n_workers=2
        )

        assert set(written_serial.keys()) == set(SAMPLES)
        assert set(written_parallel.keys()) == set(SAMPLES)

        for sample_name in SAMPLES:
            img_serial = np.asarray(load_image(written_serial[sample_name]["image"]))
            img_parallel = np.asarray(load_image(written_parallel[sample_name]["image"]))
            assert np.array_equal(img_serial, img_parallel), sample_name

            csv_serial = pd.read_csv(written_serial[sample_name]["csv"])
            csv_parallel = pd.read_csv(written_parallel[sample_name]["csv"])
            pd.testing.assert_frame_equal(
                csv_serial.sort_values("SegmentID").reset_index(drop=True),
                csv_parallel.sort_values("SegmentID").reset_index(drop=True),
            )

            # Sanity: the group TrackIDs actually made it into the combined image.
            expected_ids = {tid for m in GROUP_TRACKID.values() for tid in m.values()}
            assert set(np.unique(img_serial)) - {0} == expected_ids

            # And the tracks CSV path matches the deterministic convention.
            assert written_serial[sample_name]["csv"] == group_tracked_csv_path(
                out_serial, sample_name, group_id
            )
    finally:
        _cleanup(case_dir)
