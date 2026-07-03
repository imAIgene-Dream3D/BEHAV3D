"""
ConvPaint batch inference queue (UNIFIED MULTI-CLASS MODE).

Loads the trained unified ConvPaint model + label-map sidecar, runs it on
every sample × timepoint exactly **once**, and writes per-cell-type
``segments`` / ``mask`` zarrs in the same directory layout that the rest of
the BEHAV3D pipeline expects.

Output layout (unchanged from before)::

    images/{sample_name}/{sample_name}_{cell_type}_segments.zarr
    images/{sample_name}/{sample_name}_{cell_type}_mask.zarr
    images/{sample_name}/{sample_name}_mask_dead.zarr

Two strategies are supported:

* ``ConvPaint Mask + EDT/Watershed`` (default) — uses ``segment`` to obtain
  multi-class labels per voxel, then ``mask = labels == k`` per cell type
  and EDT/watershed.
* ``ConvPaint Probability + Watershed`` — uses ``predict_probas`` to obtain a
  ``(n_classes, Z, Y, X)`` volume per frame, slices ``probas[k-1]`` for each
  cell type, and runs the seeded watershed.

Death stays a separate binary ConvPaint model (``ConvPaintModel_Dead.pkl``)
because death is a state, not a cell-type identity.
"""

import sys
import time
import shutil
import warnings
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import zarr
from tqdm import tqdm

from behav3d.io.images import _ensure_zarr, load_image
from behav3d.core.metadata import (
    detect_organoid_types_from_metadata,
    detect_immune_cell_types_from_metadata,
    detect_other_cell_types_from_metadata,
    has_dead_channel,
)
from behav3d.preprocessing.segmentation.convpaint_postprocess import (
    segment_per_z,
    predict_probas_per_z,
    mask_to_instances,
    probability_to_instances,
)
from behav3d.preprocessing.segmentation.convpaint_label_map import (
    BACKGROUND_LABEL,
    cell_type_order,
    death_input_channels_path,
    load_input_channels,
    load_label_map,
    unified_input_channels_path,
    unified_label_map_path,
    unified_model_path,
)


STRATEGY_EDT = "ConvPaint Mask + EDT/Watershed"
STRATEGY_PEAK_EDT = "ConvPaint Mask + Peak EDT/Watershed"
STRATEGY_PROB = "ConvPaint Probability + Watershed"
DEFAULT_STRATEGY = STRATEGY_EDT


def _normalize_strategy(value):
    val = str(value or "").strip()
    if val == STRATEGY_PROB:
        return STRATEGY_PROB
    if val == STRATEGY_PEAK_EDT:
        return STRATEGY_PEAK_EDT
    if val and val != STRATEGY_EDT:
        print(
            f"  \u26a0\ufe0f Unknown convpaint_strategy '{val}' "
            f"\u2014 falling back to '{STRATEGY_EDT}'."
        )
    return STRATEGY_EDT


def _resolve_effective_strategy(cell_type, global_strategy, per_ct_strategies):
    """Return the effective ConvPaint strategy for *cell_type*."""
    if per_ct_strategies:
        explicit = per_ct_strategies.get(cell_type)
        if explicit:
            return _normalize_strategy(explicit)
    return _normalize_strategy(global_strategy)


def _normalize_class_weights(value):
    val = str(value or "").strip()
    if val.lower() in ("", "none", "null"):
        return None
    if val.lower() == "balanced":
        return "Balanced"
    if val.lower() == "sqrtbalanced":
        return "SqrtBalanced"
    return val


def _default_device():
    """Return ``'cuda:0'`` when an NVIDIA GPU is present, otherwise ``'auto'``."""
    try:
        import torch
        if torch.cuda.is_available():
            return "cuda:0"
    except Exception:
        pass
    return "auto"


def _set_torch_device(device_str):
    """Set the active PyTorch CUDA device before running ConvPaint."""
    if not device_str or device_str in ("auto", "cpu"):
        return
    try:
        import torch
        if device_str.startswith("cuda:"):
            idx = int(device_str.split(":")[1])
            torch.cuda.set_device(idx)
            print(
                f"ConvPaint: Selected GPU device {idx} ({torch.cuda.get_device_name(idx)})"
            )
    except Exception as e:
        print(f"\u26a0\ufe0f Could not set CUDA device '{device_str}': {e}")


# ---------------------------------------------------------------------------
# Loading helpers
# ---------------------------------------------------------------------------

def _load_unified_model(pixelclass_dir, provided_path=None):
    """Load the trained unified ConvPaint model.

    ``provided_path`` overrides the default location when supplied. Returns
    ``None`` if no model can be loaded.
    """
    from napari_convpaint import ConvpaintModel

    if provided_path:
        p = Path(provided_path)
        if p.exists() and p.suffix.lower() == ".pkl":
            print(f"  \u2705 Using provided unified ConvPaint model \u2192 {p.name}")
            return ConvpaintModel(model_path=str(p))

    p = unified_model_path(pixelclass_dir)
    if p.exists():
        print(f"  \u2705 Found unified ConvPaint model \u2192 {p.name}")
        return ConvpaintModel(model_path=str(p))

    print(f"  \u274c No unified ConvPaint model found in {pixelclass_dir}")
    return None


def _load_unified_label_map(pixelclass_dir, provided_path=None):
    """Load the unified label-map sidecar (next to the .pkl by default)."""
    if provided_path:
        p = Path(provided_path)
        if p.exists():
            return load_label_map(p)

    p = unified_label_map_path(pixelclass_dir)
    if p.exists():
        return load_label_map(p)
    raise FileNotFoundError(
        f"Could not find unified label-map sidecar at {p}. "
        "Train the unified ConvPaint classifier before running batch inference."
    )


def _load_death_model(pixelclass_dir, provided_path=None):
    """Load the binary death ConvPaint model (separate from the unified one)."""
    from napari_convpaint import ConvpaintModel

    if provided_path:
        p = Path(provided_path)
        if p.exists() and p.suffix.lower() == ".pkl":
            print(f"  \u2705 Using provided death ConvPaint model \u2192 {p.name}")
            return ConvpaintModel(model_path=str(p))

    pixelclass_dir = Path(pixelclass_dir)
    candidates = [
        "ConvPaintModel_Dead.pkl",
        "ConvPaintModel_dead.pkl",
    ]
    for fname in candidates:
        p = pixelclass_dir / fname
        if p.exists():
            print(f"  \u2705 Found death ConvPaint model \u2192 {p.name}")
            return ConvpaintModel(model_path=str(p))

    print(f"  \u26a0\ufe0f No death ConvPaint model found in {pixelclass_dir}")
    return None


def _open_zarr_output(path, dtype, shape, overwrite):
    path = Path(path)
    if overwrite and path.exists():
        shutil.rmtree(path)
    if path.exists():
        return zarr.open(str(path), mode="r+")
    return zarr.open(
        str(path), mode="w", shape=shape, dtype=dtype,
        chunks=(1,) + shape[1:],
    )


def _resolve_unified_clf_path(*legacy_dicts, unified_path=None):
    """Pick the unified model path, prioritising explicit user input.

    The legacy per-cell-type ``clf_*_paths`` dicts are accepted only to keep
    older notebooks/configs from breaking; if any of them point to existing
    .pkl files we surface a one-time deprecation warning.
    """
    if unified_path:
        return str(unified_path)

    for d in legacy_dicts:
        if not d:
            continue
        for v in d.values():
            if v and Path(v).exists() and str(v).lower().endswith(".pkl"):
                warnings.warn(
                    "Per-cell-type ConvPaint classifier paths are no longer "
                    "supported in unified mode; supply a single "
                    "'clf_unified_path' instead. Ignoring legacy entry.",
                    DeprecationWarning,
                    stacklevel=3,
                )
                break
    return None


def _slice_frame_channels(frame, channel_indices, *, label):
    """Slice ``(C, Z, Y, X)`` frames by sidecar channels with safe fallback."""
    if not channel_indices:
        return frame
    valid = []
    n_channels = int(frame.shape[0])
    for idx in channel_indices:
        try:
            idx = int(idx)
        except (TypeError, ValueError):
            continue
        if 0 <= idx < n_channels and idx not in valid:
            valid.append(idx)
    if not valid:
        print(
            f"  \u26a0\ufe0f {label}: input-channel sidecar has no valid indices "
            f"for frame with {n_channels} channels; using all channels."
        )
        return frame
    return frame[valid]


# ---------------------------------------------------------------------------
# Main queue
# ---------------------------------------------------------------------------

def run_convpaint_segmentation(
    output_dir,
    metadata,
    metadata_csv_path,
    timepoint_range=None,
    clf_unified_path=None,
    clf_organoid_paths=None,  # legacy, ignored except for deprecation warning
    clf_immune_paths=None,    # legacy
    clf_other_paths=None,     # legacy
    clf_death_path=None,
    convpaint_config=None,
    organoid_types=None,
    immune_types=None,
    other_types=None,
    only_segment=False,
    overwrite_existing=False,
    n_workers=1,
    convpaint_strategy=DEFAULT_STRATEGY,
    per_ct_convpaint_strategies=None,
    only_cell_types=None,
    **_kwargs,
):
    """Unified-mode ConvPaint inference queue.

    Runs the unified multi-class classifier exactly once per timepoint and
    writes per-cell-type ``segments`` / ``mask`` zarrs alongside the
    optional binary ``mask_dead.zarr``. Drop-in replacement for the legacy
    per-cell-type queue: the on-disk output layout is unchanged so
    downstream tracking/analysis works without modification.

    ``per_ct_convpaint_strategies`` (plugin only): optional
    ``{cell_type: strategy}`` mapping. When set, each cell type uses its own
    strategy regardless of the global one. Mirrors APOC's
    ``per_ct_strategies`` parameter. Cell types missing from the mapping
    fall back to ``convpaint_strategy``.
    """
    convpaint_strategy = _normalize_strategy(convpaint_strategy)
    per_ct_norm = {
        str(ct): _normalize_strategy(strat)
        for ct, strat in (per_ct_convpaint_strategies or {}).items()
    }
    if per_ct_norm:
        summary = ", ".join(f"{k}:{v}" for k, v in per_ct_norm.items())
        print(
            f"Running ConvPaint segmentation pipeline (global strategy: {convpaint_strategy};"
            f" per-cell-type: {summary})"
        )
    else:
        print(f"Running ConvPaint segmentation pipeline (strategy: {convpaint_strategy})")

    output_dir = Path(output_dir)
    pixelclass_dir = output_dir / "images" / "PixelClassification"

    organoid_types = detect_organoid_types_from_metadata(metadata)
    immune_types = detect_immune_cell_types_from_metadata(metadata)
    other_types = detect_other_cell_types_from_metadata(metadata)
    all_cell_types = organoid_types + immune_types + other_types
    _has_death = has_dead_channel(metadata)

    # Resolve the unified model path (with a friendly fallback that warns on
    # legacy per-cell-type clf path dicts).
    resolved_unified = _resolve_unified_clf_path(
        clf_organoid_paths, clf_immune_paths, clf_other_paths,
        unified_path=clf_unified_path,
    )

    # Load unified model + label map up-front.
    unified_model = None
    label_map = None
    if not only_segment:
        unified_model = _load_unified_model(pixelclass_dir, resolved_unified)
    try:
        label_map = _load_unified_label_map(pixelclass_dir)
    except FileNotFoundError as exc:
        if not only_segment:
            print(f"  \u274c {exc}")
            return metadata
        else:
            label_map = None  # only_segment mode tolerates a missing map

    # Determine which cell types to process: must be in metadata AND in label map.
    if label_map is not None:
        active_cell_types = [
            ct for ct in all_cell_types
            if ct in label_map.get("celltype_to_label", {})
        ]
        skipped = [ct for ct in all_cell_types if ct not in active_cell_types]
        if skipped:
            print(
                f"  \u26a0\ufe0f Skipping (not in label map): {skipped}"
            )
    else:
        # only_segment mode without a label map: process everything that has a
        # cached mask on disk (per-sample check happens later).
        active_cell_types = list(all_cell_types)

    # Optional plugin-supplied filter: restrict to a subset (e.g. Run [CT]).
    if only_cell_types:
        only_set = {str(c) for c in only_cell_types}
        active_cell_types = [ct for ct in active_cell_types if ct in only_set]

    # Death model is independent.
    model_death = None
    if _has_death and not only_segment:
        model_death = _load_death_model(pixelclass_dir, clf_death_path)

    unified_input_channels = load_input_channels(
        unified_input_channels_path(pixelclass_dir)
    )
    death_input_channels = load_input_channels(
        death_input_channels_path(pixelclass_dir)
    )
    if unified_input_channels:
        print(f"  Unified ConvPaint input channels: {unified_input_channels}")
    if death_input_channels:
        print(f"  Death ConvPaint input channels: {death_input_channels}")

    if not only_segment and unified_model is None and model_death is None:
        print("  \u26a0\ufe0f No ConvPaint models found. Nothing to do.")
        return metadata

    if not active_cell_types and not model_death:
        print("  \u26a0\ufe0f No active cell types to process.")
        return metadata

    cfg = convpaint_config or {}

    # Allow config to override tile_image / use_dask on loaded models.
    use_dask = bool(cfg.get("convpaint_use_dask", False))
    tile_image = bool(cfg.get("convpaint_tile_image", False))
    if use_dask and not tile_image:
        print(
            "  \u2139\ufe0f 'use_dask' enabled but 'tile_image' is not \u2014 "
            "Dask only parallelizes tiled predictions, so it will be a no-op."
        )
    for model in (unified_model, model_death):
        if model is None:
            continue
        if tile_image:
            try:
                model.set_params(tile_image=True)
            except Exception:
                pass
        if use_dask:
            try:
                model.set_params(use_dask=True)
            except Exception as exc:
                print(
                    f"  \u26a0\ufe0f Could not set use_dask on a ConvPaint model: {exc}"
                )

    sample_names = metadata['sample_name'].unique()

    for sample_name in sample_names:
        t0 = time.time()
        img_outdir = output_dir / "images" / sample_name
        done_dir = img_outdir / ".done_markers"

        sample_row = metadata[metadata['sample_name'] == sample_name].iloc[0]
        raw_image_path = sample_row.get('raw_image_path')
        if not raw_image_path or not Path(raw_image_path).exists():
            print(f"  \u26a0\ufe0f Raw image not found for {sample_name}")
            continue

        raw_image_path = Path(raw_image_path)
        
        axis_order = sample_row.get('dimension_order', "TCZYX")
        if not isinstance(axis_order, str) or not axis_order:
            axis_order = "TCZYX"

        img = load_image(raw_image_path, axis_order=axis_order)
        n_timepoints = img.shape[0]

        if timepoint_range is not None:
            if isinstance(timepoint_range, (range, list)):
                t_range = list(timepoint_range)
            else:
                s, e = timepoint_range
                t_range = list(range(s, e + 1))
        else:
            t_range = list(range(n_timepoints))

        # Check existing segmentations for skip logic (cell-type level)
        active_cts_for_sample = list(active_cell_types)
        has_death_for_sample = _has_death

        if not overwrite_existing:
            # Check cell types
            remaining_cts = []
            for ct in active_cts_for_sample:
                seg_path = img_outdir / f"{sample_name}_{ct}_segments.zarr"
                if seg_path.exists():
                    all_done = all((done_dir / f"{ct}_t{t:04d}.done").exists() for t in t_range)
                    if all_done:
                        continue  # Skip this cell type for this sample
                remaining_cts.append(ct)
            active_cts_for_sample = remaining_cts

            # Check death mask
            if has_death_for_sample and model_death:
                death_path = img_outdir / f"{sample_name}_mask_dead.zarr"
                if death_path.exists():
                    all_done = all((done_dir / f"dead_t{t:04d}.done").exists() for t in t_range)
                    if all_done:
                        has_death_for_sample = False

            if not active_cts_for_sample and (not has_death_for_sample or not model_death):
                print(f"  \u23ED\uFE0F Skipping {sample_name} (all selected outputs already exist)")
                continue

        _ensure_zarr(raw_image_path, label=f"Raw image for '{sample_name}'")

        img_outdir.mkdir(parents=True, exist_ok=True)
        done_dir.mkdir(parents=True, exist_ok=True)
        spatial_shape = img.shape[2:]
        out_shape = (n_timepoints,) + spatial_shape

        zarr_segs = {}
        zarr_masks = {}
        for ct in active_cts_for_sample:
            zarr_segs[ct] = _open_zarr_output(
                img_outdir / f"{sample_name}_{ct}_segments.zarr",
                "uint16", out_shape, overwrite_existing if not only_segment else True,
            )
            if only_segment:
                mask_path = img_outdir / f"{sample_name}_{ct}_mask.zarr"
                if not mask_path.exists():
                    raise FileNotFoundError(
                        f"Cannot 'Only Resegment' \u2014 mask not found: {mask_path}"
                    )
                zarr_masks[ct] = zarr.open(str(mask_path), mode="r")
            else:
                zarr_masks[ct] = _open_zarr_output(
                    img_outdir / f"{sample_name}_{ct}_mask.zarr",
                    "uint16", out_shape, overwrite_existing,
                )

        zarr_death = None
        if has_death_for_sample and model_death and not only_segment:
            zarr_death = _open_zarr_output(
                img_outdir / f"{sample_name}_mask_dead.zarr",
                "uint16", out_shape, overwrite_existing,
            )

        def _load_tp(img_obj, t_idx):
            return np.asarray(img_obj[t_idx])

        pbar = tqdm(
            t_range,
            desc=f"  \u23F1\uFE0F {sample_name}",
            leave=True, unit="tp",
            file=sys.stdout, dynamic_ncols=True,
        )

        with ThreadPoolExecutor(max_workers=1) as executor:
            future = executor.submit(_load_tp, img, t_range[0])

            saved_device = cfg.get("convpaint_device", _default_device())
            if saved_device and saved_device not in ("auto", "cpu"):
                _set_torch_device(saved_device)
                fe_device = "gpu"
            elif saved_device == "cpu":
                fe_device = "cpu"
            else:
                fe_device = "auto"

            for i, t in enumerate(pbar):
                t_img = future.result()  # (C, Z, Y, X)

                if i + 1 < len(t_range):
                    future = executor.submit(_load_tp, img, t_range[i + 1])

                # Filter active cell types for this timepoint
                active_cts_for_timepoint = []
                for ct in active_cts_for_sample:
                    if not overwrite_existing and (done_dir / f"{ct}_t{t:04d}.done").exists():
                        continue
                    active_cts_for_timepoint.append(ct)

                # Bucket cell types by their effective strategy so each
                # strategy can run in its own (single) inference call.
                strat_buckets = {}
                for ct in active_cts_for_timepoint:
                    eff = _resolve_effective_strategy(
                        ct, convpaint_strategy, per_ct_norm,
                    )
                    strat_buckets.setdefault(eff, []).append(ct)

                if only_segment:
                    # Re-run only the post-processing using cached masks.
                    for eff_strategy, cts in strat_buckets.items():
                        _resegment_timepoint(
                            t, t_img, cts, cfg, label_map,
                            zarr_masks, zarr_segs, eff_strategy,
                        )
                else:
                    if active_cts_for_timepoint:
                        unified_frame = _slice_frame_channels(
                            t_img, unified_input_channels, label="Unified ConvPaint"
                        )
                        for eff_strategy, cts in strat_buckets.items():
                            if eff_strategy == STRATEGY_PROB:
                                _process_probability_timepoint(
                                    unified_model, unified_frame, t, cts,
                                    label_map, cfg, zarr_masks, zarr_segs,
                                    fe_device=fe_device, diag=(i == 0),
                                )
                            else:
                                _process_edt_timepoint(
                                    unified_model, unified_frame, t, cts,
                                    label_map, cfg, zarr_masks, zarr_segs,
                                    fe_device=fe_device,
                                    diag=(i == 0),
                                    marker_strategy=(
                                        "peak" if eff_strategy == STRATEGY_PEAK_EDT
                                        else "threshold"
                                    ),
                                )
                
                # Mark timepoint as done for these cell types
                for ct in active_cts_for_timepoint:
                    (done_dir / f"{ct}_t{t:04d}.done").touch()

                if zarr_death is not None and model_death is not None:
                    if not overwrite_existing and (done_dir / f"dead_t{t:04d}.done").exists():
                        pass
                    else:
                        _t0_death = time.time()
                        death_frame = _slice_frame_channels(
                            t_img, death_input_channels, label="Death ConvPaint"
                        )
                        death_seg = segment_per_z(
                            model_death, death_frame, fe_use_device=fe_device,
                        )
                        death_mask = (death_seg >= 2).astype(np.uint16)
                        zarr_death[t] = death_mask
                        (done_dir / f"dead_t{t:04d}.done").touch()
                        if i == 0:
                            print(
                                f"    \u23F1 death segment_per_z + write: "
                                f"{time.time() - _t0_death:.2f}s"
                            )

        # Write paths back into the metadata table.
        row_idx = metadata.index[
            metadata['sample_name'] == sample_name
        ].tolist()[0]

        for ct in active_cts_for_sample:
            if ct in organoid_types:
                prefix = 'or'
            elif ct in immune_types:
                prefix = 'im'
            elif ct in other_types:
                prefix = 'ot'
            else:
                continue

            path_col = f'{prefix}_{ct}_segments_image_path'
            seg_path = img_outdir / f"{sample_name}_{ct}_segments.zarr"
            if seg_path.exists():
                if path_col not in metadata.columns:
                    metadata[path_col] = pd.NA
                if metadata[path_col].dtype != 'object':
                    metadata[path_col] = metadata[path_col].astype('object')
                metadata.at[row_idx, path_col] = str(seg_path)

        if zarr_death is not None:
            death_path = img_outdir / f"{sample_name}_mask_dead.zarr"
            if death_path.exists():
                if 'dead_mask_path' not in metadata.columns:
                    metadata['dead_mask_path'] = pd.NA
                if metadata['dead_mask_path'].dtype != 'object':
                    metadata['dead_mask_path'] = metadata['dead_mask_path'].astype('object')
                metadata.at[row_idx, 'dead_mask_path'] = str(death_path)

        elapsed = time.time() - t0
        print(f"  \u2705 {sample_name} done in {elapsed:.1f}s")

    print("\n\u2705 ConvPaint queue finished successfully. "
          "Updated *_segments_image_path columns in metadata.")

    new_md = metadata.copy()
    new_md['preprocess_pixelclass_done'] = True
    new_md['preprocess_pixelclass_engine'] = "ConvPaint"
    return new_md


# ---------------------------------------------------------------------------
# Strategy implementations (one inference call per timepoint).
# ---------------------------------------------------------------------------

def _process_edt_timepoint(
    unified_model, frame, t, active_cell_types, label_map, cfg,
    zarr_masks, zarr_segs, fe_device=None, diag=False,
    marker_strategy="threshold",
):
    """One ``segment`` call → derive per-cell-type binary masks → EDT/watershed."""
    _t0 = time.time()
    multi_labels = segment_per_z(unified_model, frame, fe_use_device=fe_device)
    multi_labels = np.asarray(multi_labels)
    if diag:
        print(
            f"    \u23F1 unified segment_per_z: {time.time() - _t0:.2f}s  "
            f"shape={multi_labels.shape}"
        )

    for ct in active_cell_types:
        _t0 = time.time()
        k = int(label_map["celltype_to_label"][ct])
        mask_out = (multi_labels == k).astype(np.uint16)
        zarr_masks[ct][t] = mask_out

        fill_holes = bool(cfg.get(f"{ct}_fill_holes", True))
        opening_nr_pixels = int(cfg.get(f"{ct}_opening_nr_pixels", 0))
        edt_thr = float(cfg.get(f"{ct}_edt_threshold", 1.0))
        segment_size_min = int(cfg.get(f"{ct}_segment_size_min", 10))
        _raw_peak_dist = cfg.get(f"{ct}_peak_min_distance", 0)
        peak_min_distance = None if (not _raw_peak_dist or float(_raw_peak_dist) <= 0) else float(_raw_peak_dist)
        peak_min_ratio = float(cfg.get(f"{ct}_peak_min_ratio", 0.35))

        instances = mask_to_instances(
            mask_out,
            edt_thr=edt_thr,
            opening_nr_pixels=opening_nr_pixels,
            fill_holes=fill_holes,
            segment_size_min=segment_size_min,
            marker_strategy=marker_strategy,
            peak_min_distance=peak_min_distance,
            peak_min_ratio=peak_min_ratio,
        )
        zarr_segs[ct][t] = instances.astype(np.uint16)
        if diag:
            print(f"    \u23F1 {ct} mask->instances + writes: {time.time() - _t0:.2f}s")


def _process_probability_timepoint(
    unified_model, frame, t, active_cell_types, label_map, cfg,
    zarr_masks, zarr_segs, fe_device=None, diag=False,
):
    """One ``predict_probas`` call → slice per cell type → seeded watershed."""
    from behav3d.preprocessing.segmentation.segmentation_utils import postprocess_mask

    _t0 = time.time()
    probas = predict_probas_per_z(unified_model, frame, fe_use_device=fe_device)
    probas = np.asarray(probas)
    if diag:
        print(
            f"    \u23F1 unified predict_probas: {time.time() - _t0:.2f}s  "
            f"shape={probas.shape}"
        )

    for ct in active_cell_types:
        k = int(label_map["celltype_to_label"][ct])
        prob_axis = k - 1  # class indices [1..N+1] -> axis [0..N]
        if not (0 <= prob_axis < probas.shape[0]):
            # Should not happen with a correctly trained model + label map.
            print(
                f"    \u26a0\ufe0f Probability slice missing for {ct} "
                f"(class {k}, axis {prob_axis} out of {probas.shape[0]}). "
                "Skipping."
            )
            continue

        prob_map = probas[prob_axis].astype(np.float32)
        opening_nr_pixels = int(cfg.get(f"{ct}_opening_nr_pixels", 0))
        mask_thr = float(cfg.get(f"{ct}_prob_mask_threshold", 0.5))
        seed_thr = float(cfg.get(f"{ct}_prob_seed_threshold", 0.8))
        segment_size_min = int(cfg.get(f"{ct}_segment_size_min", 10))

        if diag:
            print(
                f"    \u23F1 {ct} prob slice range="
                f"[{prob_map.min():.3f}, {prob_map.max():.3f}]"
            )

        mask_out = postprocess_mask(
            prob_map > mask_thr,
            fill_holes=False,
            opening_nr_pixels=opening_nr_pixels,
        ).astype(np.uint16)
        zarr_masks[ct][t] = mask_out

        instances = probability_to_instances(
            prob_map,
            mask_thr=mask_thr,
            seed_thr=seed_thr,
            opening_nr_pixels=opening_nr_pixels,
            segment_size_min=segment_size_min,
        )
        zarr_segs[ct][t] = instances.astype(np.uint16)


def _resegment_timepoint(
    t, _frame, active_cell_types, cfg, label_map,
    zarr_masks, zarr_segs, strategy,
):
    """Resegment from cached masks (no feature extraction).

    For probability strategy we don't have the original probability map cached
    on disk, so we fall back to EDT-style reseeding using the cached binary
    mask. This matches the previous behaviour for the legacy queue.
    """
    for ct in active_cell_types:
        mask_out = np.asarray(zarr_masks[ct][t])
        opening_nr_pixels = int(cfg.get(f"{ct}_opening_nr_pixels", 0))
        segment_size_min = int(cfg.get(f"{ct}_segment_size_min", 10))

        _raw_peak_dist = cfg.get(f"{ct}_peak_min_distance", 0)
        peak_min_distance = None if (not _raw_peak_dist or float(_raw_peak_dist) <= 0) else float(_raw_peak_dist)
        peak_min_ratio = float(cfg.get(f"{ct}_peak_min_ratio", 0.35))
        if strategy == STRATEGY_PROB:
            seed_thr = float(cfg.get(f"{ct}_prob_seed_threshold", 0.8))
            instances = mask_to_instances(
                mask_out,
                edt_thr=seed_thr,
                opening_nr_pixels=opening_nr_pixels,
                fill_holes=False,
                segment_size_min=segment_size_min,
            )
        else:
            edt_thr = float(cfg.get(f"{ct}_edt_threshold", 1.0))
            fill_holes = bool(cfg.get(f"{ct}_fill_holes", True))
            instances = mask_to_instances(
                mask_out,
                edt_thr=edt_thr,
                opening_nr_pixels=opening_nr_pixels,
                fill_holes=fill_holes,
                segment_size_min=segment_size_min,
                marker_strategy="peak" if strategy == STRATEGY_PEAK_EDT else "threshold",
                peak_min_distance=peak_min_distance,
                peak_min_ratio=peak_min_ratio,
            )
        zarr_segs[ct][t] = instances.astype(np.uint16)
