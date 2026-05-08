"""
ConvPaint batch inference queue.

Standalone module that loads trained ConvpaintModel .pkl files,
runs them on every sample × timepoint, and writes segments / masks to .zarr
in the same directory layout that the rest of the BEHAV3D pipeline expects.

Output layout (identical to APOC):
    images/{sample_name}/{sample_name}_{cell_type}_segments.zarr
    images/{sample_name}/{sample_name}_{cell_type}_mask.zarr
    images/{sample_name}/{sample_name}_mask_dead.zarr
"""

import os
import sys
import time
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor

from behav3d.io.images import _ensure_zarr, load_image, save_as_zarr
from behav3d.core.metadata import (
    detect_organoid_types_from_metadata,
    detect_immune_cell_types_from_metadata,
    detect_other_cell_types_from_metadata,
    has_dead_channel,
)


def _set_torch_device(device_str):
    """Set the active PyTorch CUDA device before running ConvPaint."""
    if not device_str or device_str in ("auto", "cpu"):
        return
    try:
        import torch
        if device_str.startswith("cuda:"):
            idx = int(device_str.split(":")[1])
            torch.cuda.set_device(idx)
            print(f"ConvPaint: Selected GPU device {idx} ({torch.cuda.get_device_name(idx)})")
    except Exception as e:
        print(f"\u26a0\ufe0f Could not set CUDA device '{device_str}': {e}")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_convpaint_model(pixelclass_dir, cell_type, provided_path=None):
    """Load a trained ConvpaintModel .pkl file."""
    from napari_convpaint import ConvpaintModel

    pixelclass_dir = Path(pixelclass_dir)

    # 1. Manually provided path
    if provided_path:
        p = Path(provided_path)
        if p.exists() and p.suffix.lower() == '.pkl':
            print(f"  ✅ {cell_type}: using provided ConvPaint model → {p.name}")
            return ConvpaintModel(model_path=str(p))

    # 2. Well-known filenames
    filenames_to_try = [
        f'ConvPaintModel_{cell_type}.pkl',
        f'ConvPaintModel_{cell_type.capitalize()}.pkl',
        f'ConvPaintModel_{cell_type.lower()}.pkl',
    ]
    for fname in filenames_to_try:
        p = pixelclass_dir / fname
        if p.exists():
            print(f"  ✅ {cell_type}: found ConvPaint model → {p.name}")
            return ConvpaintModel(model_path=str(p))

    # 3. Scan directory for .pkl containing cell type name
    if pixelclass_dir.is_dir():
        ct_lower = cell_type.lower()
        for pkl_file in pixelclass_dir.glob('*.pkl'):
            if ct_lower in pkl_file.stem.lower() and 'convpaint' in pkl_file.stem.lower():
                print(f"  ✅ {cell_type}: found ConvPaint model (scan) → {pkl_file.name}")
                return ConvpaintModel(model_path=str(pkl_file))

    print(f"  ❌ {cell_type}: no ConvPaint .pkl model found in {pixelclass_dir}")
    return None


def _open_zarr_output(path, dtype, shape, overwrite):
    """Open a zarr array for output, creating/recreating as needed."""
    path = Path(path)
    if overwrite and path.exists():
        shutil.rmtree(path)
    if path.exists():
        return zarr.open(str(path), mode="r+")
    else:
        return zarr.open(
            str(path), mode="w", shape=shape, dtype=dtype,
            chunks=(1,) + shape[1:],
        )



# ---------------------------------------------------------------------------
# Main queue
# ---------------------------------------------------------------------------

def run_convpaint_segmentation(
    output_dir,
    metadata,
    metadata_csv_path,
    timepoint_range=None,
    clf_organoid_paths=None,
    clf_immune_paths=None,
    clf_other_paths=None,
    clf_death_path=None,
    convpaint_config=None,
    organoid_types=None,
    immune_types=None,
    other_types=None,
    only_segment=False,
    overwrite_existing=False,
    n_workers=1,
    convpaint_strategy="ConvPaint (Direct Segmentation)",
    **_kwargs,
):
    """
    Fully independent ConvPaint inference queue.

    Writes output .zarr files to the same directory layout as the APOC/CPU
    pipeline so downstream tracking works seamlessly.
    """
    print("Running ConvPaint segmentation pipeline...")

    output_dir = Path(output_dir)
    pixelclass_dir = output_dir / "images" / "PixelClassification"

    # --- Discover cell types ---
    organoid_types = organoid_types or detect_organoid_types_from_metadata(metadata)
    immune_types = immune_types or detect_immune_cell_types_from_metadata(metadata)
    other_types = other_types or detect_other_cell_types_from_metadata(metadata)
    all_cell_types = organoid_types + immune_types + other_types
    _has_death = has_dead_channel(metadata)

    # --- Load models ---
    models = {}
    clf_organoid_paths = clf_organoid_paths or {}
    for ct in organoid_types:
        m = _load_convpaint_model(pixelclass_dir, ct, clf_organoid_paths.get(ct))
        if m:
            models[ct] = m

    clf_immune_paths = clf_immune_paths or {}
    for ct in immune_types:
        m = _load_convpaint_model(pixelclass_dir, ct, clf_immune_paths.get(ct))
        if m:
            models[ct] = m

    clf_other_paths = clf_other_paths or {}
    for ct in other_types:
        m = _load_convpaint_model(pixelclass_dir, ct, clf_other_paths.get(ct))
        if m:
            models[ct] = m

    active_cell_types = [ct for ct in all_cell_types if ct in models]

    model_death = None
    if _has_death:
        model_death = _load_convpaint_model(pixelclass_dir, "dead", clf_death_path)

    print(f"  Loaded ConvPaint models for: {list(models.keys())}"
          + (f" + Death" if model_death else ""))

    if not active_cell_types and not model_death:
        print("  ⚠️ No ConvPaint models found. Nothing to do.")
        return metadata

    # --- Read config for per-cell-type params ---
    cfg = convpaint_config or {}

    # --- Process each sample ---
    sample_names = metadata['sample_name'].unique()

    for sample_name in sample_names:
        t0 = time.time()
        img_outdir = output_dir / "images" / sample_name

        # Skip if all outputs exist
        if not overwrite_existing:
            all_done = True
            for ct in active_cell_types:
                if not (img_outdir / f"{sample_name}_{ct}_segments.zarr").exists():
                    all_done = False
                    break
            if all_done and _has_death and model_death:
                if not (img_outdir / f"{sample_name}_mask_dead.zarr").exists():
                    all_done = False
            if all_done and active_cell_types:
                print(f"  ⏭️ Skipping {sample_name} (all outputs already exist)")
                continue

        # Load raw image
        sample_row = metadata[metadata['sample_name'] == sample_name].iloc[0]
        raw_image_path = sample_row.get('raw_image_path')
        if not raw_image_path or not Path(raw_image_path).exists():
            print(f"  ⚠️ Raw image not found for {sample_name}")
            continue

        raw_image_path = Path(raw_image_path)
        _ensure_zarr(raw_image_path, label=f"Raw image for '{sample_name}'")

        axis_order = sample_row.get('dimension_order', "TCZYX")
        if not isinstance(axis_order, str) or not axis_order:
            axis_order = "TCZYX"

        img = load_image(raw_image_path, axis_order=axis_order)  # (T, C, Z, Y, X)
        n_timepoints = img.shape[0]

        img_outdir.mkdir(parents=True, exist_ok=True)
        spatial_shape = img.shape[2:]  # (Z, Y, X)
        out_shape = (n_timepoints,) + spatial_shape

        # Create output zarr arrays
        zarr_segs = {}
        zarr_masks = {}
        for ct in active_cell_types:
            zarr_segs[ct] = _open_zarr_output(
                img_outdir / f"{sample_name}_{ct}_segments.zarr",
                "uint16", out_shape, overwrite_existing if not only_segment else True,
            )
            if only_segment:
                mask_path = img_outdir / f"{sample_name}_{ct}_mask.zarr"
                if not mask_path.exists():
                    raise FileNotFoundError(
                        f"Cannot 'Only Resegment' — mask not found: {mask_path}"
                    )
                zarr_masks[ct] = zarr.open(str(mask_path), mode="r")
            else:
                zarr_masks[ct] = _open_zarr_output(
                    img_outdir / f"{sample_name}_{ct}_mask.zarr",
                    "uint16", out_shape, overwrite_existing,
                )

        zarr_death = None
        if _has_death and model_death and not only_segment:
            zarr_death = _open_zarr_output(
                img_outdir / f"{sample_name}_mask_dead.zarr",
                "uint16", out_shape, overwrite_existing,
            )

        # Timepoint range
        if timepoint_range is not None:
            if isinstance(timepoint_range, (range, list)):
                t_range = list(timepoint_range)
            else:
                s, e = timepoint_range
                t_range = list(range(s, e + 1))
        else:
            t_range = list(range(n_timepoints))

        def _load_tp(img_obj, t_idx):
            return np.asarray(img_obj[t_idx])

        pbar = tqdm(
            t_range,
            desc=f"  ⏱️ {sample_name}",
            leave=True, unit="tp",
            file=sys.stdout, dynamic_ncols=True,
        )

        with ThreadPoolExecutor(max_workers=1) as executor:
            future = executor.submit(_load_tp, img, t_range[0])

            # Apply GPU device selection once (global setting)
            saved_device = cfg.get("convpaint_device", "auto")
            if saved_device and saved_device not in ("auto", "cpu"):
                _set_torch_device(saved_device)
                fe_device = "gpu"
            elif saved_device == "cpu":
                fe_device = "cpu"
            else:
                fe_device = None

            for i, t in enumerate(pbar):
                t_img = future.result()  # (C, Z, Y, X)

                if i + 1 < len(t_range):
                    future = executor.submit(_load_tp, img, t_range[i + 1])

                for ct in active_cell_types:
                    cp_model = models[ct]

                    # Use all channels (no filtering)
                    frame = t_img

                    # Strategy routing
                    if convpaint_strategy == "ConvPaint Probability Map + Watershed":
                        _process_probability_strategy(
                            cp_model, frame, t, ct, cfg,
                            zarr_masks, zarr_segs, only_segment,
                            fe_device=fe_device,
                            diag=(i == 0),
                        )
                    elif convpaint_strategy == "ConvPaint Mask + EDT/Watershed Resegmentation":
                        _process_edt_strategy(
                            cp_model, frame, t, ct, cfg,
                            zarr_masks, zarr_segs, only_segment,
                            fe_device=fe_device,
                        )
                    else:
                        # Direct segmentation
                        _process_direct_strategy(
                            cp_model, frame, t, ct,
                            zarr_masks, zarr_segs, only_segment,
                            fe_device=fe_device,
                        )

                # Death channel
                if zarr_death is not None and model_death is not None:
                    death_frame = t_img  # use all channels
                    death_seg = np.asarray(
                        model_death.segment(death_frame)
                    )
                    death_mask = (death_seg > 0).astype(np.uint16)
                    zarr_death[t] = death_mask

        # Update metadata
        row_idx = metadata.index[
            metadata['sample_name'] == sample_name
        ].tolist()[0]

        for ct in active_cell_types:
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
        print(f"  ✅ {sample_name} done in {elapsed:.1f}s")

    print("\n✅ ConvPaint queue finished successfully.")

    new_md = metadata.copy()
    new_md['preprocess_pixelclass_done'] = True
    new_md['preprocess_pixelclass_engine'] = "ConvPaint"
    return new_md


# ---------------------------------------------------------------------------
# Strategy implementations
# ---------------------------------------------------------------------------

def _process_direct_strategy(
    cp_model, frame, t, ct, zarr_masks, zarr_segs, only_segment, fe_device=None,
):
    """Direct segmentation: ConvPaint class labels become instance labels."""
    if not only_segment:
        seg_out = np.asarray(cp_model.segment(frame, fe_use_device=fe_device)).astype(np.uint16)
        # ConvPaint labels: 1=background, 2=foreground (matching our annotation)
        # Convert to mask: foreground = class 2
        mask_out = (seg_out >= 2).astype(np.uint16)
        zarr_masks[ct][t] = mask_out
    else:
        mask_out = np.asarray(zarr_masks[ct][t])
        seg_out = mask_out

    # For direct strategy, use connected components as instances
    from skimage.measure import label as sk_label
    instances = sk_label(mask_out.astype(bool)).astype(np.uint16)
    zarr_segs[ct][t] = instances


def _process_edt_strategy(
    cp_model, frame, t, ct, cfg, zarr_masks, zarr_segs, only_segment, fe_device=None,
):
    """Mask + EDT/Watershed resegmentation."""
    from behav3d.preprocessing.segmentation.segmentation_utils import (
        postprocess_mask, segment_mask,
    )

    if not only_segment:
        seg_out = np.asarray(cp_model.segment(frame, fe_use_device=fe_device)).astype(np.uint16)
        mask_out = (seg_out >= 2).astype(np.uint16)
        zarr_masks[ct][t] = mask_out
    else:
        mask_out = np.asarray(zarr_masks[ct][t])

    fill_holes = cfg.get(f"{ct}_fill_holes", True)
    opening_nr_pixels = int(cfg.get(f"{ct}_opening_nr_pixels", 0))
    proc_mask = postprocess_mask(
        mask_out, fill_holes=fill_holes, opening_nr_pixels=opening_nr_pixels,
    )

    edt_thr = float(cfg.get(f"{ct}_edt_threshold", 1.0))
    segment_size_min = int(cfg.get(f"{ct}_segment_size_min", 10))
    seg_refined = segment_mask(
        proc_mask, edt_thr=edt_thr, edt_thr_refined=None,
        segment_size_min=segment_size_min, use_dims=3, n_workers=1,
    )
    zarr_segs[ct][t] = np.asarray(seg_refined).astype(np.uint16)


def _process_probability_strategy(
    cp_model, frame, t, ct, cfg,
    zarr_masks, zarr_segs, only_segment, fe_device=None, diag=False,
):
    """Probability Map + Watershed strategy."""
    from scipy import ndimage as ndi
    from skimage.measure import label as sk_label
    from skimage.segmentation import watershed
    from behav3d.preprocessing.segmentation import (
        segment_size_filter, segment_2d_filter,
    )
    from behav3d.preprocessing.segmentation.segmentation_utils import (
        postprocess_mask, segment_mask,
    )

    opening_nr_pixels = int(cfg.get(f"{ct}_opening_nr_pixels", 0))

    if not only_segment:
        _t0 = time.time()
        # predict_probas returns (n_classes, Z, Y, X) — we want class index 1
        # (foreground probability, since class 0=background in ConvPaint output)
        probas = cp_model.predict_probas(frame, fe_use_device=fe_device)
        # probas shape: (n_classes, Z, Y, X) for single image
        probas = np.asarray(probas)
        # The foreground class is the last class (index -1, i.e. class 2 in binary)
        # ConvPaint class 1=BG, class 2=FG, so index 1 in probas (0-indexed)
        if probas.ndim == 3:
            # Single-class edge case — treat as foreground probability
            prob_map = probas.astype(np.float32)
        else:
            # Multi-class: take the foreground class probability
            prob_map = probas[-1].astype(np.float32)

        if diag:
            print(
                f"    ⏱ {ct} predict (prob): {time.time()-_t0:.2f}s  "
                f"range=[{prob_map.min():.3f}, {prob_map.max():.3f}]"
            )

        _t0 = time.time()
        mask_thr = float(cfg.get(f"{ct}_prob_mask_threshold", 0.5))
        mask_out = postprocess_mask(
            prob_map > mask_thr,
            fill_holes=False,
            opening_nr_pixels=opening_nr_pixels,
        ).astype(np.uint16)
        zarr_masks[ct][t] = mask_out
        if diag:
            print(
                f"    ⏱ {ct} mask threshold + write: {time.time()-_t0:.2f}s"
            )
    else:
        mask_out = np.asarray(zarr_masks[ct][t])
        prob_map = None

    seed_thr = float(cfg.get(f"{ct}_prob_seed_threshold", 0.8))
    segment_size_min = int(cfg.get(f"{ct}_segment_size_min", 10))

    if prob_map is not None:
        _t0 = time.time()
        mask_bool = mask_out.astype(bool)
        cc_labels = sk_label(mask_bool)
        seed_mask = (prob_map > seed_thr) & mask_bool
        segments = np.zeros_like(cc_labels, dtype=np.uint16)
        next_id = 0

        obj_slices = ndi.find_objects(cc_labels)
        for comp_idx, slc in enumerate(obj_slices, start=1):
            if slc is None:
                continue
            comp_mask = cc_labels[slc] == comp_idx
            sub_seeds = sk_label(seed_mask[slc] & comp_mask)
            n_seeds = int(sub_seeds.max())

            if n_seeds <= 1:
                next_id += 1
                segments[slc][comp_mask] = next_id
            else:
                sub_result = watershed(
                    -prob_map[slc], markers=sub_seeds, mask=comp_mask,
                )
                for s in range(1, int(sub_result.max()) + 1):
                    next_id += 1
                    segments[slc][sub_result == s] = next_id

        if diag:
            print(
                f"    ⏱ {ct} per-component watershed: "
                f"{time.time()-_t0:.2f}s"
            )
    else:
        segments = segment_mask(
            mask_out.astype(bool),
            edt_thr=seed_thr,
            segment_size_min=segment_size_min,
            use_dims=3, n_workers=1,
        )

    _t0 = time.time()
    segments = segment_size_filter(segments, size_min=segment_size_min)
    segments = segment_2d_filter(segments)
    if diag:
        print(f"    ⏱ {ct} size+2d filter: {time.time()-_t0:.2f}s")

    zarr_segs[ct][t] = np.asarray(segments).astype(np.uint16)
