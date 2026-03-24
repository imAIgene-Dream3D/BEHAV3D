"""
Independent APOC GPU-Accelerated inference queue.

This module is completely standalone - it does NOT import anything from
napari_pixelclassifier.py. It loads trained .cl ObjectSegmenter classifiers,
runs them on every sample timepoint, and writes segments / masks to .zarr in
the same directory layout that the rest of the BEHAV3D pipeline (tracking etc.)
expects.
"""

import os
import sys
import time
import shutil
from pathlib import Path
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor

# Defeat buggy PyOpenCL compiler caching that causes TypeErrors on some systems
os.environ['PYOPENCL_NO_CACHE'] = '1'
# Suppress compiler output unless there is an actual error
os.environ['PYOPENCL_COMPILER_OUTPUT'] = '0'

import numpy as np
import zarr
import apoc
import pyclesperanto_prototype as cle

from behav3d.io.images import load_image, save_as_zarr
from behav3d.core.metadata import (
    detect_organoid_types_from_metadata,
    detect_immune_cell_types_from_metadata,
    detect_other_cell_types_from_metadata,
    has_dead_channel,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_classifier(pixelclass_dir, cell_type, provided_path=None):
    """Load an APOC ObjectSegmenter .cl file. Case-insensitive lookup."""
    pixelclass_dir = Path(pixelclass_dir)

    # 1. Use manually-provided path if it exists
    if provided_path:
        provided_path = Path(provided_path)
        if provided_path.exists() and provided_path.suffix.lower() == '.cl':
            print(f"  ✅ {cell_type}: using provided classifier → {provided_path.name}")
            return apoc.ObjectSegmenter(opencl_filename=str(provided_path))

    # 2. Try well-known name patterns
    filenames_to_try = [
        f'PixelClassifier_{cell_type}.cl',
        f'PixelClassifier_{cell_type.capitalize()}.cl',
        f'PixelClassifier_{cell_type.lower()}.cl',
        f'PixelClassifier_{cell_type.upper()}.cl',
    ]
    for fname in filenames_to_try:
        p = pixelclass_dir / fname
        if p.exists():
            print(f"  ✅ {cell_type}: found classifier → {p.name}")
            return apoc.ObjectSegmenter(opencl_filename=str(p))

    # 3. Fallback: scan directory for ANY .cl file whose name contains the cell type
    if pixelclass_dir.is_dir():
        ct_lower = cell_type.lower()
        for cl_file in pixelclass_dir.glob('*.cl'):
            if ct_lower in cl_file.stem.lower():
                print(f"  ✅ {cell_type}: found classifier (scan) → {cl_file.name}")
                return apoc.ObjectSegmenter(opencl_filename=str(cl_file))

    # 4. Not found – emit diagnostic warnings
    if pixelclass_dir.is_dir():
        for jf in pixelclass_dir.glob('*.joblib'):
            if cell_type.lower() in jf.stem.lower():
                print(f"  ⚠️ Found {jf.name} but APOC engine needs .cl files. "
                      f"Please (re-)train with APOC.")
                break
    print(f"  ❌ {cell_type}: no .cl classifier found in {pixelclass_dir}")
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
            chunks=(1,) + shape[1:],  # one chunk per timepoint
        )


# ---------------------------------------------------------------------------
# Main queue
# ---------------------------------------------------------------------------

def run_apoc_segmentation(
    output_dir,
    metadata,
    metadata_csv_path,
    timepoint_range=None,
    clf_organoid_paths=None,
    clf_immune_paths=None,
    clf_other_paths=None,
    clf_death_path=None,
    apoc_config=None,
    organoid_types=None,
    immune_types=None,
    other_types=None,
    only_segment=False,
    overwrite_existing=False,
    n_workers=1,
    gpu_device=None,
    **_kwargs,  # absorb unused CPU-specific params (EDT, opening, fill_holes etc.)
):
    """
    Fully independent APOC ObjectSegmenter inference queue.

    Writes output .zarr files to the same directory layout the CPU pipeline
    uses so that the downstream tracking step can seamlessly consume them:
        images/{sample_name}/{sample_name}_{cell_type}_segments.zarr
        images/{sample_name}/{sample_name}_{cell_type}_mask.zarr
        images/{sample_name}/{sample_name}_mask_dead.zarr
    """
    if gpu_device:
        print(f"Selecting GPU device: {gpu_device}")
        cle.select_device(gpu_device)
        
    print("Running fully independent APOC ObjectSegmenter pipeline...")

    output_dir = Path(output_dir)
    pixelclass_dir = output_dir / "images" / "PixelClassification"

    # --- Discover cell types from metadata (same logic as the rest of BEHAV3D) ---
    organoid_types = detect_organoid_types_from_metadata(metadata)
    immune_types = detect_immune_cell_types_from_metadata(metadata)
    other_types = detect_other_cell_types_from_metadata(metadata)
    all_cell_types = organoid_types + immune_types + other_types

    has_death = has_dead_channel(metadata)

    # --- Load classifiers ---
    classifiers = {}
    clf_organoid_paths = clf_organoid_paths or {}
    for ct in organoid_types:
        clf = _load_classifier(pixelclass_dir, ct, clf_organoid_paths.get(ct))
        if clf: classifiers[ct] = clf

    clf_immune_paths = clf_immune_paths or {}
    for ct in immune_types:
        clf = _load_classifier(pixelclass_dir, ct, clf_immune_paths.get(ct))
        if clf: classifiers[ct] = clf

    clf_other_paths = clf_other_paths or {}
    for ct in other_types:
        clf = _load_classifier(pixelclass_dir, ct, clf_other_paths.get(ct))
        if clf: classifiers[ct] = clf

    active_cell_types = [ct for ct in all_cell_types if ct in classifiers]

    clf_death = None
    if has_death:
        death_default = pixelclass_dir / 'PixelClassifier_Death.cl'
        if clf_death_path:
            p = Path(clf_death_path)
            if p.exists() and p.suffix.lower() == '.cl':
                # Copy to canonical location only if it's a different file
                if not death_default.exists() or not p.samefile(death_default):
                    pixelclass_dir.mkdir(parents=True, exist_ok=True)
                    shutil.copy(p, death_default)
        # Also try the directory-scan fallback for Death
        if not death_default.exists() and pixelclass_dir.is_dir():
            for cl_file in pixelclass_dir.glob('*.cl'):
                if 'death' in cl_file.stem.lower() or 'dead' in cl_file.stem.lower():
                    death_default = cl_file
                    break
        if death_default.exists():
            clf_death = apoc.ObjectSegmenter(opencl_filename=str(death_default))

    print(f"  Loaded classifiers for: {list(classifiers.keys())}"
          + (f" + Death" if clf_death else ""))

    # --- Prepare channel indices and post-processing params ---
    def get_indices(ct):
        if not apoc_config: return [0]
        chan_names = apoc_config.get(f"apoc_{ct}_channels", [])
        if not chan_names: return [0]
        indices = []
        for name in chan_names:
            try: indices.append(int(name.split(" ")[-1]))
            except: pass
        return indices if indices else [0]

    clf_channels = {ct: get_indices(ct) for ct in all_cell_types}
    death_channels = get_indices("dead") if clf_death else []

    # --- Process each sample ---
    sample_names = metadata['sample_name'].unique()
    
    for sample_name in sample_names:
        t0 = time.time()
 
        # Output directory
        img_outdir = output_dir / "images" / sample_name
 
        # SKIP LOGIC: If all active outputs already exist and overwrite=False, skip
        if not overwrite_existing:
            all_done = True
            for ct in active_cell_types:
                if not (img_outdir / f"{sample_name}_{ct}_segments.zarr").exists():
                    all_done = False
                    break
            if all_done and has_death and clf_death:
                if not (img_outdir / f"{sample_name}_mask_dead.zarr").exists():
                    all_done = False

            if all_done and active_cell_types:  # guard: don't skip if nothing was checked
                print(f"  ⏭️ Skipping {sample_name} (all outputs already exist)")
                continue

        # Use raw_image_path from metadata (source of truth)
        sample_row = metadata[metadata['sample_name'] == sample_name].iloc[0]
        raw_image_path = sample_row.get('raw_image_path')
        if not raw_image_path:
            print(f"  ⚠️ No raw_image_path found in metadata for {sample_name}")
            continue
        raw_image_path = Path(raw_image_path)
        if not raw_image_path.exists():
            print(f"  ⚠️ Raw image not found for {sample_name}: {raw_image_path}")
            continue

        # Get dimension order from metadata for this sample
        axis_order = sample_row.get('dimension_order', "TCZYX")
        if not isinstance(axis_order, str) or not axis_order:
            axis_order = "TCZYX"

        img = load_image(raw_image_path, axis_order=axis_order)             # lazy load via metadata path
        n_timepoints = img.shape[0]

        # Ensure output directory exists
        img_outdir.mkdir(parents=True, exist_ok=True)

        # Spatial shape for output arrays
        spatial_shape = img.shape[2:]           # (Z, Y, X)  — img is (T, C, Z, Y, X)
        out_shape = (n_timepoints,) + spatial_shape

        # Create output zarr arrays
        zarr_segs = {}
        zarr_masks = {}
        for ct in active_cell_types:
            zarr_segs[ct] = _open_zarr_output(
                img_outdir / f"{sample_name}_{ct}_segments.zarr",
                "uint16", out_shape, overwrite_existing,
            )
            zarr_masks[ct] = _open_zarr_output(
                img_outdir / f"{sample_name}_{ct}_mask.zarr",
                "uint16", out_shape, overwrite_existing,
            )

        zarr_death = None
        if has_death and clf_death:
            zarr_death = _open_zarr_output(
                img_outdir / f"{sample_name}_mask_dead.zarr",
                "uint16", out_shape, overwrite_existing,
            )

        # Timepoint progress via tqdm (per-sample)
        if timepoint_range is not None:
            if isinstance(timepoint_range, (range, list)):
                t_range = list(timepoint_range)
            else:
                # Handle old tuple (start, end)
                s, e = timepoint_range
                t_range = list(range(s, e + 1))
        else:
            t_range = list(range(n_timepoints))

        def _load_tp(img_obj, t_idx):
            return np.asarray(img_obj[t_idx])

        pbar = tqdm(t_range, desc=f"  ⏱️ {sample_name}", leave=True, unit="tp", file=sys.stdout, dynamic_ncols=True)
        
        # Prefetch pipeline via ThreadPoolExecutor
        # Overlaps disk I/O (loading next TP) with GPU work (predicting current TP)
        with ThreadPoolExecutor(max_workers=1) as executor:
            # Prefetch the first timepoint
            future = executor.submit(_load_tp, img, t_range[0])
            
            for i, t in enumerate(pbar):
                # 1. Get current timepoint (blocks if disk I/O not finished)
                t_img = future.result() 
                
                # 2. Prefetch the NEXT timepoint immediately
                if i + 1 < len(t_range):
                    future = executor.submit(_load_tp, img, t_range[i + 1])

                # 3. Process current timepoint (GPU-bound)
                for ct in active_cell_types:
                    indices = clf_channels[ct]
                    imgs = [t_img[i_ch] for i_ch in indices]
                    imgs_to_pass = imgs[0] if len(imgs) == 1 else imgs
                    
                    # Predict gives labels (segments) directly
                    seg = np.asarray(classifiers[ct].predict(image=imgs_to_pass)).astype(np.uint16)
                    mask = (seg > 0).astype(np.uint16)
                    
                    zarr_segs[ct][t] = seg
                    zarr_masks[ct][t] = mask

                if zarr_death is not None:
                    indices = death_channels
                    imgs = [t_img[i_ch] for i_ch in indices]
                    imgs_to_pass = imgs[0] if len(imgs) == 1 else imgs
                    death_mask = (np.asarray(clf_death.predict(image=imgs_to_pass)) > 0).astype(np.uint16)
                    zarr_death[t] = death_mask

        elapsed = time.time() - t0
        print(f"  ✅ {sample_name} done in {elapsed:.1f}s")

    print("\n✅ APOC queue finished successfully.")

    # Tag metadata
    new_md = metadata.copy()
    new_md['preprocess_pixelclass_done'] = True
    new_md['preprocess_pixelclass_engine'] = "APOC"
    return new_md
