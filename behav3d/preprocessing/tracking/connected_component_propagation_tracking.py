"""Connected-Component Propagation: watershed propagation that never lets a
single track ID span more than one physically disconnected region.

Like ``propagation_tracking.py``, each timepoint's mask is watershed-seeded
by the previous timepoint's tracked labels. But afterwards, every physically
disconnected fragment of the result — regardless of which label watershed
happened to paint it with — is re-identified against ALL previous-frame
tracks' footprints via a stable matching (tracks "propose" to fragments in
order of pixel overlap; a fragment keeps its best suitor and rejects it for
someone with more overlap, bumping the loser back to its next-best
candidate). A fragment nobody has (sufficient) overlap with becomes a
brand-new track. A track ID can therefore never span disconnected geometry.
"""

from pathlib import Path
import shutil

import numpy as np
import pandas as pd
from skimage.measure import label as cc_label
from tqdm import tqdm

from behav3d.io.images import load_image, append_to_zarr
from behav3d.preprocessing.segmentation import segment_size_filter
from behav3d.preprocessing.tracking import convert_tracked_image_to_csv
from behav3d.preprocessing.tracking.propagation_tracking import (
    _resolve_segments_column,
    _resolve_tracking_output_columns,
    _ensure_metadata_output_columns,
    _watershed_from_markers,
)


def _resolve_fragment_identities(new_seg, seg_prev_tp, next_label, min_overlap_fraction=0.0):
    """Reassign every connected fragment to whichever previous-frame track
    has the largest overlap with it (stable matching across ALL tracks, not
    just the label watershed originally painted the fragment with), using
    exact pixel-overlap counts.

    A track only qualifies as a candidate for a fragment if its
    previous-frame footprint fills at least ``min_overlap_fraction`` of that
    fragment. A fragment with no qualifying candidate becomes a brand-new
    track. Returns ``(new_seg, next_label)``.
    """
    if not new_seg.any():
        return new_seg, next_label

    # 1. Physically-connected fragments, per-label connectivity (never let
    #    a neighbouring different label's pixels bridge two fragments of
    #    the same label together), each given a unique running fragment id.
    frag_ids = np.zeros_like(new_seg, dtype=np.int64)
    frag_masks = {}
    next_fid = 1
    for lbl in np.unique(new_seg):
        if lbl == 0:
            continue
        components = cc_label(new_seg == lbl)
        for cid in np.unique(components):
            if cid == 0:
                continue
            comp_mask = components == cid
            frag_ids[comp_mask] = next_fid
            frag_masks[next_fid] = comp_mask
            next_fid += 1
    frag_sizes = {fid: int(mask.sum()) for fid, mask in frag_masks.items()}

    # 2. Exact pixel-overlap count for every (fragment, previous track) pair
    #    that shares at least one pixel — one vectorized pass.
    both = (frag_ids != 0) & (seg_prev_tp != 0)
    track_candidates = {}
    if np.any(both):
        pairs = np.stack([frag_ids[both], seg_prev_tp[both]], axis=1)
        uniq_pairs, counts = np.unique(pairs, axis=0, return_counts=True)
        for (fid, t), overlap in zip(uniq_pairs, counts):
            fid, t, overlap = int(fid), int(t), int(overlap)
            if overlap / frag_sizes[fid] < min_overlap_fraction:
                continue  # doesn't fill enough of the fragment to qualify
            track_candidates.setdefault(t, []).append((fid, overlap))
    for t in track_candidates:
        track_candidates[t].sort(key=lambda x: -x[1])

    # 3. Gale-Shapley deferred acceptance: tracks propose to fragments in
    #    order of overlap size; a fragment always holds its best suitor.
    next_choice = {t: 0 for t in track_candidates}
    held_by = {}  # fragment_id -> (track, overlap)
    free = list(track_candidates.keys())
    while free:
        t = free.pop()
        choices = track_candidates[t]
        idx = next_choice[t]
        if idx >= len(choices):
            continue  # exhausted candidates -> track ends here
        fid, overlap = choices[idx]
        next_choice[t] = idx + 1
        current = held_by.get(fid)
        if current is None or overlap > current[1]:
            if current is not None:
                free.append(current[0])
            held_by[fid] = (t, overlap)
        else:
            free.append(t)

    # 4. Apply the assignment; unclaimed fragments get fresh IDs.
    out = np.zeros_like(new_seg)
    for fid, mask in frag_masks.items():
        holder = held_by.get(fid)
        if holder is None:
            out[mask] = next_label
            next_label += 1
        else:
            out[mask] = holder[0]
    return out, next_label


def _run_single_timepoint_propagation_connected(
    t_seg,
    seg_prev_tp,
    next_label,
    dilation_nr_pixels=2,
    segment_size_min=100,
    min_overlap_fraction=0.0,
    ):
    mask = t_seg != 0
    seg_prev_tp[mask==0]=0
    new_seg = _watershed_from_markers(mask, seg_prev_tp, dilation_nr_pixels)
    new_seg, next_label = _resolve_fragment_identities(
        new_seg, seg_prev_tp, next_label, min_overlap_fraction=min_overlap_fraction,
    )
    new_seg = segment_size_filter(new_seg, size_min=segment_size_min)
    return new_seg, next_label


def propagate_tracks_connected(
    segments_path,
    tracked_img_outpath,
    tracked_csv_outpath,
    element_size_x=1,
    element_size_y=1,
    element_size_z=1,
    dilation_nr_pixels=2,
    segment_size_min=100,
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
            t_tracked_seg, next_label = _run_single_timepoint_propagation_connected(
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


def run_connected_component_propagation_tracking(
    metadata,
    output_dir,
    cell_type,
    overwrite=False,
    dilation_nr_pixels=2,
    segment_size_min=100,
    min_overlap_fraction=0.0,
    progress_cb=None,
    **kwargs
    ):
    """Run Connected-Component Propagation tracking on any cell type.

    Behaves like propagation tracking, except a track ID can never span more
    than one physically disconnected region: every disconnected fragment is
    resolved against ALL previous-frame tracks' footprints (not just the
    label watershed produced it with), so a splitting object hands its
    leftover pieces either to whichever other track fits best, or to a new
    track if nobody's overlap is decisive.

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
    min_overlap_fraction : float
        Minimum fraction of a fragment's area that a previous-frame track's
        footprint must fill for that track to be allowed to claim it.
        ``0.0`` (default) means any positive overlap qualifies.
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
            propagate_tracks_connected(
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
