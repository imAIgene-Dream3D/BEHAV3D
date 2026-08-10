"""Bounded Propagation: watershed propagation that never lets a
single track ID span more than one physically disconnected region.

Like ``propagation_tracking.py``, each timepoint's mask is watershed-seeded
by the previous timepoint's tracked labels — but the decision about which
track gets to seed which piece of geometry is made BEFORE watershed runs,
not cleaned up after the fact. The raw current-frame mask is first
connected-component-labelled into independent "regions". Every
previous-frame track then picks, by simple pixel-overlap argmax, the single
region its old footprint overlaps most (several different tracks are
allowed to pick the same region — that's what preserves genuine merge/split
behaviour). Every track's marker pixels lying outside its chosen region are
pruned to zero before the existing two-pass watershed runs. As a direct,
structural consequence: a region only one track chose is flooded to that
track's label untouched; a region several tracks chose is still split
between them by watershed exactly as before; and a region no track chose
(all its markers were pruned) is left completely unlabeled by watershed. A
final pass connected-component-labels that unlabeled leftover foreground and
hands each disconnected leftover patch its own brand-new track id. A track
ID can therefore never span disconnected geometry — and a genuinely single
leftover blob is never accidentally splintered into more than one new track
just because several old tracks each grazed it.
"""

from pathlib import Path
import shutil

import numpy as np
import pandas as pd
from skimage.measure import label as bp_label
from tqdm import tqdm

from behav3d.io.images import load_image, append_to_zarr
from behav3d.preprocessing.segmentation import segment_size_filter, segment_2d_filter
from behav3d.preprocessing.tracking import convert_tracked_image_to_csv
from behav3d.preprocessing.tracking.propagation_tracking import (
    _resolve_segments_column,
    _resolve_tracking_output_columns,
    _ensure_metadata_output_columns,
    _watershed_from_markers,
)


def _prune_markers_to_best_region(regions, seg_prev_tp, min_overlap_fraction=0.0):
    """Restrict every previous-frame track's marker pixels to the single
    current-frame connected region it overlaps most, before watershed runs.

    ``regions`` is the raw current-frame mask connected-component-labelled
    independently of any track (``skimage.measure.label`` on the mask,
    default connectivity). For every nonzero previous-frame track label T in
    ``seg_prev_tp``, this finds the region T's old footprint overlaps most
    by exact pixel count, and zeroes every one of T's marker pixels that
    falls outside that one region. Downstream watershed can then never let
    T's seed leak into a region some other track's old footprint explains
    better.

    A track only keeps a home region if that region qualifies: the
    track/region pixel-overlap count divided by the REGION's total pixel
    count (the whole raw current-frame region, before any watershed
    splitting) must be >= ``min_overlap_fraction``. A track that fails this
    for every region it touches is pruned to zero everywhere; any region it
    was the sole candidate for is then picked up by the leftover-labelling
    pass instead of being assigned to it.

    Multiple different tracks MAY end up sharing the same best region — this
    is intentional, not a bug: watershed still splits that region between
    them exactly as before.

    Ties (a track with exactly equal overlap against two regions) are
    broken deterministically toward the lower region id.

    Returns a new array, same shape/dtype as ``seg_prev_tp``, with every
    non-home-region pixel zeroed. ``seg_prev_tp`` itself is not mutated.
    """
    pruned = np.zeros_like(seg_prev_tp)

    both = (regions != 0) & (seg_prev_tp != 0)
    if not np.any(both):
        return pruned

    # Region sizes (denominator for min_overlap_fraction).
    reg_ids, reg_counts = np.unique(regions[regions != 0], return_counts=True)
    region_sizes = {int(r): int(c) for r, c in zip(reg_ids, reg_counts)}

    # Exact pixel-overlap count for every (track, region) pair that shares
    # at least one pixel — one vectorized pass.
    pairs = np.stack([seg_prev_tp[both], regions[both]], axis=1)
    uniq_pairs, counts = np.unique(pairs, axis=0, return_counts=True)
    # uniq_pairs is lexicographically sorted ascending by (track, region),
    # so within a track, candidate regions are visited in ascending region
    # id order -> keeping only on strict improvement reproduces
    # np.argmax's first-occurrence tie-break (lowest region id wins).

    best_region = {}  # track -> (region_id, overlap)
    for (t, r), overlap in zip(uniq_pairs, counts):
        t, r, overlap = int(t), int(r), int(overlap)
        if overlap / region_sizes[r] < min_overlap_fraction:
            continue  # doesn't fill enough of the region to qualify
        current = best_region.get(t)
        if current is None or overlap > current[1]:
            best_region[t] = (r, overlap)

    for t, (r, _overlap) in best_region.items():
        keep = (seg_prev_tp == t) & (regions == r)
        pruned[keep] = t

    return pruned


def _label_leftover_regions(mask, new_seg, next_label):
    """Give every disconnected patch of watershed-unclaimed foreground its
    own brand-new track id.

    After watershed floods every region that some track(s) chose as their
    best home, any ``mask`` area belonging to a region NO track chose (all
    its markers were pruned in ``_prune_markers_to_best_region``) is left
    labelled 0 even though it's foreground. This connected-component-labels
    that leftover foreground and assigns each disconnected leftover
    component a fresh, incrementing ``next_label`` — so a genuinely single
    leftover blob (e.g. a merge region no old track claimed as its best)
    becomes exactly ONE new track instead of splintering into as many
    pieces as tracks that grazed it.

    Returns ``(new_seg, next_label)``; does not mutate ``new_seg`` in place.
    """
    leftover = mask & (new_seg == 0)
    out = new_seg.copy()
    if not np.any(leftover):
        return out, next_label

    leftover_components = bp_label(leftover)
    for cid in np.unique(leftover_components):
        if cid == 0:
            continue
        out[leftover_components == cid] = next_label
        next_label += 1
    return out, next_label


def _run_single_timepoint_propagation_bounded(
    t_seg,
    seg_prev_tp,
    next_label,
    dilation_nr_pixels=2,
    segment_size_min=20,
    min_overlap_fraction=0.0,
    ):
    mask = t_seg != 0
    seg_prev_tp[mask==0]=0

    regions = bp_label(mask)
    pruned_markers = _prune_markers_to_best_region(
        regions, seg_prev_tp, min_overlap_fraction=min_overlap_fraction,
    )

    # Gap-closing dilation is always disabled here (regardless of the
    # dilation_nr_pixels argument): it would let watershed fuse an
    # unclaimed leftover fragment into a neighboring track's label across
    # a real background gap, producing one track ID that spans two
    # disconnected pixel blobs.
    new_seg = _watershed_from_markers(mask, pruned_markers, dilation_nr_pixels=0)
    new_seg, next_label = _label_leftover_regions(mask, new_seg, next_label)

    new_seg = segment_size_filter(new_seg, size_min=segment_size_min)
    new_seg = segment_2d_filter(new_seg)
    return new_seg, next_label


def propagate_tracks_bounded(
    segments_path,
    tracked_img_outpath,
    tracked_csv_outpath,
    element_size_x=1,
    element_size_y=1,
    element_size_z=1,
    dilation_nr_pixels=2,
    segment_size_min=20,
    min_overlap_fraction=0.0,
    **kwargs
    ):

    seg = load_image(segments_path)

    if tracked_img_outpath.exists():
        shutil.rmtree(tracked_img_outpath)

    seg_prev_tp = np.asarray(seg[0])
    next_label = int(seg_prev_tp.max()) + 1
    for t, t_seg in tqdm(enumerate(seg), total=seg.shape[0]):
        t_seg = np.asarray(t_seg)

        if t==0:
            t_tracked_seg = t_seg
        else:
            t_tracked_seg, next_label = _run_single_timepoint_propagation_bounded(
                t_seg,
                seg_prev_tp,
                next_label,
                dilation_nr_pixels=dilation_nr_pixels,
                segment_size_min=segment_size_min,
                min_overlap_fraction=min_overlap_fraction,
            )
        seg_prev_tp = t_tracked_seg.copy()
        t_tracked_seg = np.expand_dims(t_tracked_seg, axis=0)
        append_to_zarr(t_tracked_seg, tracked_img_outpath)

    df_tracks = convert_tracked_image_to_csv(
            img_path=tracked_img_outpath,
            outpath=tracked_csv_outpath,
            element_size_x=element_size_x,
            element_size_y=element_size_y,
            element_size_z=element_size_z
        )
    return df_tracks


def run_bounded_propagation_tracking(
    metadata,
    output_dir,
    cell_type,
    overwrite=False,
    dilation_nr_pixels=2,
    segment_size_min=20,
    min_overlap_fraction=0.0,
    progress_cb=None,
    **kwargs
    ):
    """Run Bounded Propagation tracking on any cell type.

    Behaves like propagation tracking, except a track ID can never span more
    than one physically disconnected region. Before watershed runs, the raw
    current-frame mask is split into its own connected regions and every
    previous-frame track picks whichever region its old footprint overlaps
    most; a track's marker is discarded everywhere outside that one region.
    Watershed then floods each region from whichever track(s) chose it — a
    region no track picked is left over and gets a brand-new track id.

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
    dilation_nr_pixels : int
        Number of pixels to dilate the mask for watershed propagation
    segment_size_min : int
        Minimum segment size in voxels (smaller segments are filtered out)
        after watershed. Note this is NOT the same job as a segmentation-stage
        size filter: the raw input mask has typically already been filtered
        for whole-object size before it reaches tracking. This filter instead
        cleans up small/degenerate fragments that watershed splitting itself
        can introduce (e.g. a merge region divided between two tracks, or a
        small leftover sliver spun off into a new track), so it's usually
        set much lower than a segmentation-stage minimum. Segments are also
        always additionally filtered for being 2D/flat (via
        ``segment_2d_filter``), regardless of this value.
    min_overlap_fraction : float
        Minimum fraction of a current-frame connected region's area that a
        previous-frame track's old footprint must fill for that track to be
        allowed to claim the region as its home (denominator: the whole raw
        region, before any watershed splitting). ``0.0`` (default) means any
        positive overlap qualifies. Note: a track sharing a region with
        others will generally read a smaller fraction here than it would
        against just its own post-split share, so tuned non-default values
        may need lowering.
    progress_cb : callable, optional
        Called as ``progress_cb(current, total, label)`` from inside the
        per-sample loop so a GUI can drive a progress bar. ``None`` is
        a no-op (default).
    **kwargs : dict
    """
    total_samples = len(metadata)
    for i, (idx, sample) in enumerate(metadata.iterrows()):
        sample_name=sample['sample_name']
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
            (
                not tracked_csv_outpath.exists() or
                not tracked_img_outpath.exists()
            ) or overwrite
            ):
            propagate_tracks_bounded(
                segments_path=segments_path,
                tracked_img_outpath=tracked_img_outpath,
                tracked_csv_outpath=tracked_csv_outpath,
                element_size_x=element_size_x,
                element_size_y=element_size_y,
                element_size_z=element_size_z,
                dilation_nr_pixels=dilation_nr_pixels,
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
