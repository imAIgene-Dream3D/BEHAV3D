"""Reporter Propagation: static-object tracking from the single best detection.

For near-static objects whose segmentation is unreliable or intermittent over
time (present at some timepoints, missing at others, "flickering"), this
method does not attempt frame-to-frame linking at all. Instead it:

1. Pools every segmented instance from every timepoint of the video.
2. Groups together any instances that spatially overlap each other, no
   matter how far apart in time they occur (so any amount of flickering is
   tolerated with no gap/frame parameter needed — grouping is purely
   spatial).
3. Within each group, picks the single LARGEST instance as that object's
   canonical mask.
4. Stamps that one mask onto every timepoint of the video, unchanged.

This is intentionally not a variant of ``propagation_tracking.py`` (which
does frame-to-frame watershed propagation) — it shares no code path with it
other than the generic per-sample I/O plumbing.
"""

from pathlib import Path
import shutil

import numpy as np
import pandas as pd
from tqdm import tqdm

from behav3d.io.images import load_image, append_to_zarr
from behav3d.preprocessing.tracking import convert_tracked_image_to_csv
from behav3d.preprocessing.tracking.propagation_tracking import (
    _resolve_segments_column,
    _resolve_tracking_output_columns,
    _ensure_metadata_output_columns,
)


def _mask_bbox(mask):
    """Return (mins, maxs) per axis for a boolean mask's bounding box."""
    coords = np.nonzero(mask)
    mins = tuple(int(c.min()) for c in coords)
    maxs = tuple(int(c.max()) for c in coords)
    return mins, maxs


def _bboxes_overlap(bbox_a, bbox_b):
    mins_a, maxs_a = bbox_a
    mins_b, maxs_b = bbox_b
    return all(
        a_min <= b_max and b_min <= a_max
        for a_min, a_max, b_min, b_max in zip(mins_a, maxs_a, mins_b, maxs_b)
    )


def _group_static_instances(seg, segment_size_min=100, min_overlap_fraction=0.1):
    """Pool segment instances across all timepoints and group by spatial overlap.

    Returns a list of dicts, one per group: {"best_mask": bool array,
    "best_size": int}, each representing one static object's canonical mask
    (the largest instance ever detected within that spatially-overlapping
    group of instances).

    Merge decisions are tested against each group's CURRENT BEST (largest
    so far) mask, not an ever-growing union of every instance ever merged
    into the group. Testing against a growing union is a classic "chain
    drift" bug: small overlaps accumulate across many timepoints until the
    union has drifted far from where the group started, and instances end
    up merged together mainly because they each touched *some* edge of the
    sprawling union — not because they actually resemble the group's real
    object. That silently produces groups whose single "best" mask covers
    only a small fraction of everything that got merged into them (most
    instances end up with zero overlap with the final propagated mask).
    Anchoring every comparison to the current best mask keeps groups tied to
    the most complete, concrete representative observed so far.
    """
    groups = []  # each: {"best_mask": bool arr, "best_size": int, "bbox": (mins,maxs)}

    T = seg.shape[0]
    for t in tqdm(range(T), total=T, desc="Pooling instances"):
        t_seg = np.asarray(seg[t])
        labels = np.unique(t_seg)
        labels = labels[labels != 0]

        for label in labels:
            inst_mask = t_seg == label
            inst_size = int(inst_mask.sum())
            if inst_size < segment_size_min:
                continue
            inst_bbox = _mask_bbox(inst_mask)

            matched = []
            for gi, grp in enumerate(groups):
                if not _bboxes_overlap(inst_bbox, grp["bbox"]):
                    continue
                intersection = int(np.logical_and(inst_mask, grp["best_mask"]).sum())
                if intersection == 0:
                    continue
                smaller_size = min(inst_size, grp["best_size"])
                if smaller_size == 0:
                    continue
                if (intersection / smaller_size) >= min_overlap_fraction:
                    matched.append(gi)

            if not matched:
                groups.append({
                    "best_mask": inst_mask.copy(),
                    "best_size": inst_size,
                    "bbox": inst_bbox,
                })
                continue

            primary_gi = matched[0]
            primary = groups[primary_gi]
            # Merge any additional matched groups into the primary one
            # (this instance bridges them into a single object identity).
            for gi in matched[1:]:
                other = groups[gi]
                if other["best_size"] > primary["best_size"]:
                    primary["best_mask"] = other["best_mask"]
                    primary["best_size"] = other["best_size"]
                    primary["bbox"] = other["bbox"]
            for gi in sorted(matched[1:], reverse=True):
                del groups[gi]

            if inst_size > primary["best_size"]:
                primary["best_mask"] = inst_mask.copy()
                primary["best_size"] = inst_size
                primary["bbox"] = inst_bbox

    return groups


def reporter_propagate_tracks(
    segments_path,
    tracked_img_outpath,
    tracked_csv_outpath,
    element_size_x=1,
    element_size_y=1,
    element_size_z=1,
    segment_size_min=100,
    min_overlap_fraction=0.1,
    **kwargs,
):
    seg = load_image(segments_path)
    T = seg.shape[0]

    if tracked_img_outpath.exists():
        shutil.rmtree(tracked_img_outpath)

    groups = _group_static_instances(
        seg, segment_size_min=segment_size_min, min_overlap_fraction=min_overlap_fraction,
    )

    frame_shape = np.asarray(seg[0]).shape
    if not groups:
        # No qualifying objects anywhere in the video — write an empty
        # canonical (all-background) frame T times and bail with an empty CSV.
        canonical = np.zeros(frame_shape, dtype=np.uint16)
        for _ in range(T):
            append_to_zarr(np.expand_dims(canonical, axis=0), tracked_img_outpath)
        df_tracks = convert_tracked_image_to_csv(
            img_path=tracked_img_outpath,
            outpath=tracked_csv_outpath,
            element_size_x=element_size_x,
            element_size_y=element_size_y,
            element_size_z=element_size_z,
        )
        return df_tracks

    canonical = np.zeros(frame_shape, dtype=np.uint16)
    for track_id, grp in enumerate(groups, start=1):
        canonical[grp["best_mask"]] = track_id

    for _ in tqdm(range(T), total=T, desc="Stamping canonical mask"):
        append_to_zarr(np.expand_dims(canonical, axis=0), tracked_img_outpath)

    df_tracks = convert_tracked_image_to_csv(
        img_path=tracked_img_outpath,
        outpath=tracked_csv_outpath,
        element_size_x=element_size_x,
        element_size_y=element_size_y,
        element_size_z=element_size_z,
    )
    return df_tracks


def run_reporter_propagation_tracking(
    metadata,
    output_dir,
    cell_type,
    overwrite=False,
    segment_size_min=100,
    min_overlap_fraction=0.1,
    progress_cb=None,
    **kwargs,
):
    """Run Reporter Propagation tracking for near-static, unreliably-segmented objects.

    Parameters
    ----------
    metadata : pd.DataFrame
        DataFrame containing sample information
    output_dir : str or Path
        Root output directory
    cell_type : str
        Name of the cell type to track
    overwrite : bool
        Whether to overwrite existing tracking results
    segment_size_min : int
        Minimum segment size in voxels; smaller instances are ignored
        entirely (treated as noise).
    min_overlap_fraction : float
        Minimum spatial overlap (as a fraction of the smaller segment's
        size) required for two segments to be grouped as the same object.
        0 = any shared pixel counts.
    progress_cb : callable, optional
        Called as ``progress_cb(current, total, label)`` from inside the
        per-sample loop so a GUI can drive a progress bar. ``None`` is a
        no-op (default).
    **kwargs : dict
    """
    total_samples = len(metadata)
    for i, (idx, sample) in enumerate(metadata.iterrows()):
        sample_name = sample['sample_name']
        if progress_cb is not None:
            try:
                progress_cb(i, total_samples, f"{cell_type} / {sample_name}")
            except Exception:
                pass
        print(f"Tracking sample: {sample_name}")

        element_size_x = sample["pixel_distance_xy"]
        element_size_y = sample["pixel_distance_xy"]
        element_size_z = sample["pixel_distance_z"]

        tracked_img_outdir = Path(output_dir, "images", sample_name)
        tracked_csv_outdir = Path(output_dir, "trackdata", sample_name, cell_type)

        segments_col, segments_path = _resolve_segments_column(sample, cell_type)

        if pd.isna(segments_path) or segments_path is None:
            print(f"Warning: No segmentation data found for {sample_name}, {cell_type}. Skipping tracking.")
            continue

        segments_path = Path(segments_path)
        if not segments_path.exists():
            print(f"Warning: Segmentation file not found: {segments_path}. Skipping tracking.")
            continue

        tracked_img_outpath = Path(tracked_img_outdir, f"{sample_name}_{cell_type}_tracked.zarr")
        tracked_csv_outpath = Path(tracked_csv_outdir, f"{sample_name}_{cell_type}_tracks.csv")

        if not tracked_img_outdir.exists():
            tracked_img_outdir.mkdir(parents=True)
        if not tracked_csv_outdir.exists():
            tracked_csv_outdir.mkdir(parents=True)

        if (
            (not tracked_csv_outpath.exists() or not tracked_img_outpath.exists())
            or overwrite
        ):
            reporter_propagate_tracks(
                segments_path=segments_path,
                tracked_img_outpath=tracked_img_outpath,
                tracked_csv_outpath=tracked_csv_outpath,
                element_size_x=element_size_x,
                element_size_y=element_size_y,
                element_size_z=element_size_z,
                segment_size_min=segment_size_min,
                min_overlap_fraction=min_overlap_fraction,
            )
        else:
            print("Tracking already exists... Provide overwrite=True to overwrite... Loading existing tracking data")

        img_col, csv_col = _resolve_tracking_output_columns(cell_type, segments_col)
        _ensure_metadata_output_columns(metadata, img_col, csv_col)
        metadata.at[idx, img_col] = str(tracked_img_outpath)
        metadata.at[idx, csv_col] = str(tracked_csv_outpath)

    if progress_cb is not None:
        try:
            progress_cb(total_samples, total_samples, f"{cell_type} done")
        except Exception:
            pass
    return metadata
