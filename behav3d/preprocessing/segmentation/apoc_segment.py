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
import pandas as pd
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

def _get_probability_classifier(clf, positive_class_id=2):
    """Create a probability-output classifier from an existing ObjectSegmenter.

    The ObjectSegmenter .cl file contains OpenCL code that outputs integer class
    labels. This function patches the footer to output the probability of the
    foreground class instead (vote fraction s{class-1} / sum_s).

    The patched file is cached per classifier file identity so it is recreated
    when the source .cl changes after retraining.
    """
    import re, tempfile

    if not hasattr(_get_probability_classifier, '_cache'):
        _get_probability_classifier._cache = {}

    src_path = clf.opencl_file
    src_stat = os.stat(src_path)
    cache_key = (
        src_path,
        int(src_stat.st_mtime_ns),
        int(src_stat.st_size),
        int(positive_class_id),
    )
    if cache_key in _get_probability_classifier._cache:
        return _get_probability_classifier._cache[cache_key]

    with open(src_path, 'r') as f:
        content = f.read()

    # The label-mode footer looks like:
    #   float max_s=s0;
    #   int cls=1;
    #   if (max_s < s1) { ... cls=2; }
    #   ...
    #   WRITE_IMAGE (out, POS_out_INSTANCE(x,y,z,0), cls);
    # }
    #
    # Replace with probability-mode footer:
    #   float sum_s = s0 + s1 [+ s2 ...];
    #   WRITE_IMAGE (out, POS_out_INSTANCE(x,y,z,0), s{idx} / sum_s);
    # }

    # Determine number of classes from the header (counts "float s0=0", "float s1=0", etc.)
    class_vars = re.findall(r' float (s\d+)=0;', content)
    n_classes = len(class_vars)
    prob_idx = positive_class_id - 1  # APOC classes are 1-indexed, vars are 0-indexed

    # Build the new footer
    sum_expr = ' + '.join(class_vars)
    new_footer = f' float sum_s={sum_expr};\n'
    new_footer += f' WRITE_IMAGE (out, POS_out_INSTANCE(x,y,z,0), s{prob_idx} / sum_s);\n}}\n'

    # Replace: everything from " float max_s=s0;" to the end of the kernel
    patched = re.sub(
        r' float max_s=s0;.*',
        new_footer,
        content,
        flags=re.DOTALL,
    )

    # Write to a temp .cl file
    tmp_fd, tmp_path = tempfile.mkstemp(suffix='.cl', prefix='apoc_prob_')
    with os.fdopen(tmp_fd, 'w') as f:
        f.write(patched)

    # Load as a PixelClassifier with probability output
    prob_clf = apoc.PixelClassifier(
        opencl_filename=tmp_path,
        overwrite_classname='ObjectSegmenter',  # match the class name in the header
    )
    prob_clf.output_probability_of_class = positive_class_id

    _get_probability_classifier._cache[cache_key] = prob_clf
    return prob_clf

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
    apoc_strategy="APOC (Direct Instance Segmentation)",
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
                "uint16", out_shape, overwrite_existing if not only_segment else True,
            )
            if only_segment:
                # Open mask in read-only mode, we need it to do resegmentation
                mask_path = img_outdir / f"{sample_name}_{ct}_mask.zarr"
                if not mask_path.exists():
                    raise FileNotFoundError(f"Cannot 'Only Resegment' because mask does not exist: {mask_path}")
                zarr_masks[ct] = zarr.open(str(mask_path), mode="r")
            else:
                zarr_masks[ct] = _open_zarr_output(
                    img_outdir / f"{sample_name}_{ct}_mask.zarr",
                    "uint16", out_shape, overwrite_existing,
                )

        zarr_death = None
        if has_death and clf_death and not only_segment:
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

                # 3. Process current timepoint (GPU-bound or CPU-bound if only_segment)
                for ct in active_cell_types:
                    cfg = apoc_config or {}
                    indices = clf_channels[ct]
                    imgs = [t_img[i_ch] for i_ch in indices]
                    imgs_to_pass = imgs[0] if len(imgs) == 1 else imgs

                    if apoc_strategy == "APOC Probability Map + Watershed":
                        # ── Probability Map strategy ──────────────────────
                        from skimage.measure import label as sk_label
                        from skimage.segmentation import watershed
                        from behav3d.preprocessing.segmentation import segment_size_filter, segment_2d_filter
                        from behav3d.preprocessing.segmentation.segmentation_utils import postprocess_mask, segment_mask
                        _diag = (i == 0)  # print timings on first timepoint only

                        opening_nr_pixels = int(cfg.get(f"{ct}_opening_nr_pixels", 0))
                        if not only_segment:
                            _t0 = time.time()
                            prob_clf = _get_probability_classifier(classifiers[ct])
                            prob_map = np.asarray(
                                prob_clf.predict(image=imgs_to_pass)
                            ).astype(np.float32)
                            if _diag: print(f"    ⏱ {ct} predict (prob): {time.time()-_t0:.2f}s  range=[{prob_map.min():.3f}, {prob_map.max():.3f}]")

                            _t0 = time.time()
                            mask_thr = float(cfg.get(f"{ct}_prob_mask_threshold", 0.5))
                            mask_out = postprocess_mask(
                                prob_map > mask_thr,
                                fill_holes=False,
                                opening_nr_pixels=opening_nr_pixels,
                            ).astype(np.uint16)
                            zarr_masks[ct][t] = mask_out
                            if _diag: print(f"    ⏱ {ct} mask threshold + opening + write: {time.time()-_t0:.2f}s")
                        else:
                            mask_out = np.asarray(zarr_masks[ct][t])
                            prob_map = None

                        seed_thr = float(cfg.get(f"{ct}_prob_seed_threshold", 0.8))
                        segment_size_min = int(cfg.get(f"{ct}_segment_size_min", 10))

                        if prob_map is not None:
                            _t0 = time.time()
                            from scipy import ndimage as ndi
                            mask_bool = mask_out.astype(bool)
                            cc_labels = sk_label(mask_bool)
                            seed_mask = (prob_map > seed_thr) & mask_bool
                            segments = np.zeros_like(cc_labels)
                            next_id = 0
                            obj_slices = ndi.find_objects(cc_labels)
                            n_comps = len(obj_slices)
                            n_split = 0

                            for comp_idx, slc in enumerate(obj_slices, start=1):
                                if slc is None:
                                    continue
                                comp_mask = cc_labels[slc] == comp_idx
                                sub_seeds = sk_label(seed_mask[slc] & comp_mask)
                                n_seeds = sub_seeds.max()

                                if n_seeds <= 1:
                                    next_id += 1
                                    segments[slc][comp_mask] = next_id
                                else:
                                    n_split += 1
                                    sub_result = watershed(
                                        -prob_map[slc], markers=sub_seeds, mask=comp_mask
                                    )
                                    for s in range(1, sub_result.max() + 1):
                                        next_id += 1
                                        segments[slc][sub_result == s] = next_id
                            if _diag: print(f"    ⏱ {ct} per-component watershed: {time.time()-_t0:.2f}s ({n_comps} components, {n_split} split)")
                        else:
                            segments = segment_mask(mask_out.astype(bool), edt_thr=seed_thr, segment_size_min=segment_size_min, use_dims=3, n_workers=1)

                        _t0 = time.time()
                        segments = segment_size_filter(segments, size_min=segment_size_min)
                        if _diag: print(f"    ⏱ {ct} size_filter: {time.time()-_t0:.2f}s")
                        _t0 = time.time()
                        segments = segment_2d_filter(segments)
                        if _diag: print(f"    ⏱ {ct} 2d_filter: {time.time()-_t0:.2f}s")
                        zarr_segs[ct][t] = np.asarray(segments).astype(np.uint16)

                    elif apoc_strategy == "APOC Mask + EDT/Watershed Resegmentation":
                        # ── EDT/Watershed strategy ────────────────────────
                        from behav3d.preprocessing.segmentation.segmentation_utils import (
                            postprocess_mask, segment_mask
                        )
                        if not only_segment:
                            seg_out = np.asarray(classifiers[ct].predict(image=imgs_to_pass)).astype(np.uint16)
                            mask_out = (seg_out > 0).astype(np.uint16)
                            zarr_masks[ct][t] = mask_out
                        else:
                            mask_out = np.asarray(zarr_masks[ct][t])

                        fill_holes = cfg.get(f"{ct}_fill_holes", True)
                        opening_nr_pixels = cfg.get(f"{ct}_opening_nr_pixels", 0)
                        proc_mask = postprocess_mask(mask_out, fill_holes=fill_holes, opening_nr_pixels=opening_nr_pixels)

                        edt_thr = cfg.get(f"{ct}_edt_threshold", 1.0)
                        segment_size_min = cfg.get(f"{ct}_segment_size_min", 10)
                        seg_refined = segment_mask(proc_mask, edt_thr=edt_thr, edt_thr_refined=None, segment_size_min=segment_size_min, use_dims=3, n_workers=1)
                        zarr_segs[ct][t] = np.asarray(seg_refined).astype(np.uint16)

                    else:
                        # ── Default: Direct APOC ──────────────────────────
                        if not only_segment:
                            seg_out = np.asarray(classifiers[ct].predict(image=imgs_to_pass)).astype(np.uint16)
                            mask_out = (seg_out > 0).astype(np.uint16)
                            zarr_masks[ct][t] = mask_out
                        else:
                            mask_out = np.asarray(zarr_masks[ct][t])
                            seg_out = mask_out
                            if t == 0:
                                print(f"⚠️ Warning: 'Only resegment' selected but strategy is Direct APOC. Cannot tweak instance rules without EDT method. Resaving {ct} instances.")
                        zarr_segs[ct][t] = seg_out

                if zarr_death is not None:
                    indices = death_channels
                    imgs = [t_img[i_ch] for i_ch in indices]
                    imgs_to_pass = imgs[0] if len(imgs) == 1 else imgs
                    death_mask = (np.asarray(clf_death.predict(image=imgs_to_pass)) > 0).astype(np.uint16)
                    zarr_death[t] = death_mask

        # Update metadata with output paths
        row_idx = metadata.index[metadata['sample_name'] == sample_name].tolist()[0]
        for ct in active_cell_types:
            if ct in organoid_types: prefix = 'or'
            elif ct in immune_types: prefix = 'im'
            elif ct in other_types: prefix = 'ot'
            else: continue
            
            path_col = f'{prefix}_{ct}_segments_image_path'
            seg_path = img_outdir / f"{sample_name}_{ct}_segments.zarr"
            if seg_path.exists():
                # Ensure column exists and has object dtype to prevent warnings
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

    print("\n✅ APOC queue finished successfully.")

    # Tag metadata
    new_md = metadata.copy()
    new_md['preprocess_pixelclass_done'] = True
    new_md['preprocess_pixelclass_engine'] = "APOC"
    return new_md
