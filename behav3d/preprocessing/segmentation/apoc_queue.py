"""
Independent APOC GPU-Accelerated inference queue.

This module is completely standalone - it does NOT import anything from
napari_pixelclassifier.py. It loads trained .cl ObjectSegmenter classifiers,
runs them on every sample timepoint, and writes segments / masks to .zarr in
the same directory layout that the rest of the BEHAV3D pipeline (tracking etc.)
expects.
"""

import time
import shutil
from pathlib import Path

import numpy as np
import zarr
import apoc

from behav3d.io.images import load_image, save_as_zarr


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_classifier(pixelclass_dir, cell_type, provided_path=None):
    """Load an APOC ObjectSegmenter .cl file."""
    default_path = pixelclass_dir / f'PixelClassifier_{cell_type.capitalize()}.cl'
    if provided_path:
        provided_path = Path(provided_path)
        if not provided_path.samefile(default_path):
            shutil.copy(provided_path, default_path)
    clf_path = default_path
    if clf_path.exists():
        return apoc.ObjectSegmenter(opencl_filename=str(clf_path))
    else:
        raise FileNotFoundError(f"APOC classifier not found for {cell_type}: {clf_path}")


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
    only_segment=False,
    overwrite_existing=False,
    n_workers=1,
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
    print("Running fully independent APOC ObjectSegmenter pipeline...")

    output_dir = Path(output_dir)
    pixelclass_dir = output_dir / "images" / "PixelClassification"

    # --- Discover cell types from metadata ---
    organoid_types = (
        metadata['organoid_type'].dropna().unique().tolist()
        if 'organoid_type' in metadata.columns else []
    )
    immune_types = (
        metadata['immune_type'].dropna().unique().tolist()
        if 'immune_type' in metadata.columns else []
    )
    other_types = (
        metadata['other_channel_type'].dropna().unique().tolist()
        if 'other_channel_type' in metadata.columns else []
    )
    all_cell_types = organoid_types + immune_types + other_types

    has_death = (
        'dead_channel' in metadata.columns
        and not metadata['dead_channel'].isnull().all()
    )

    # --- Load classifiers ---
    classifiers = {}
    clf_organoid_paths = clf_organoid_paths or {}
    for ct in organoid_types:
        classifiers[ct] = _load_classifier(pixelclass_dir, ct, clf_organoid_paths.get(ct))

    clf_immune_paths = clf_immune_paths or {}
    for ct in immune_types:
        classifiers[ct] = _load_classifier(pixelclass_dir, ct, clf_immune_paths.get(ct))

    clf_other_paths = clf_other_paths or {}
    for ct in other_types:
        classifiers[ct] = _load_classifier(pixelclass_dir, ct, clf_other_paths.get(ct))

    clf_death = None
    if has_death:
        death_default = pixelclass_dir / 'PixelClassifier_Death.cl'
        if clf_death_path:
            p = Path(clf_death_path)
            if not p.samefile(death_default):
                shutil.copy(p, death_default)
        if death_default.exists():
            clf_death = apoc.ObjectSegmenter(opencl_filename=str(death_default))

    print(f"  Loaded classifiers for: {list(classifiers.keys())}"
          + (f" + Death" if clf_death else ""))

    # --- Process each sample ---
    for sample_name in metadata['sample_name'].unique():
        t0 = time.time()
        print(f"\n▶ Processing sample: {sample_name}")

        # Find raw zarr
        raw_zarr = output_dir / "images" / sample_name / f"{sample_name}.zarr"
        if not raw_zarr.exists():
            print(f"  ⚠️ Raw zarr not found at {raw_zarr}, skipping.")
            continue

        img = load_image(raw_zarr)             # lazy zarr / dask array
        img = np.asarray(img)                   # materialise
        n_timepoints = img.shape[0]

        # Output directory (same as CPU pipeline)
        img_outdir = output_dir / "images" / sample_name
        img_outdir.mkdir(parents=True, exist_ok=True)

        # Spatial shape for output arrays
        spatial_shape = img.shape[2:]           # (Z, Y, X)  — img is (T, C, Z, Y, X)
        out_shape = (n_timepoints,) + spatial_shape

        # Create output zarr arrays
        zarr_segs = {}
        zarr_masks = {}
        for ct in all_cell_types:
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

        # Timepoint range
        if timepoint_range is not None:
            t_range = list(timepoint_range)
        else:
            t_range = list(range(n_timepoints))

        for t_i, t in enumerate(t_range):
            t_img = img[t]                      # (C, Z, Y, X)

            for ct in all_cell_types:
                seg = np.asarray(classifiers[ct].predict(image=t_img)).astype(np.uint16)
                mask = (seg > 0).astype(np.uint16)
                zarr_segs[ct][t] = seg
                zarr_masks[ct][t] = mask

            if zarr_death is not None:
                death_mask = (
                    np.asarray(clf_death.predict(image=t_img)) > 0
                ).astype(np.uint16)
                zarr_death[t] = death_mask

            if (t_i + 1) % 5 == 0 or t_i == len(t_range) - 1:
                print(f"    Timepoint {t_i + 1}/{len(t_range)}")

        elapsed = time.time() - t0
        print(f"  ✅ {sample_name} done in {elapsed:.1f}s")

    print("\n✅ APOC queue finished successfully.")

    # Tag metadata
    new_md = metadata.copy()
    new_md['preprocess_pixelclass_done'] = True
    new_md['preprocess_pixelclass_engine'] = "APOC"
    return new_md
