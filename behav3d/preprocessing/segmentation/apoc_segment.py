"""
Independent APOC GPU-Accelerated inference queue.

This module is completely standalone - it does NOT import anything from
napari_pixelclassifier.py. It loads trained .cl ObjectSegmenter classifiers,
runs them on every sample timepoint, and writes segments / masks to .zarr in
the same directory layout that the rest of the BEHAV3D pipeline (tracking etc.)
expects.
"""

import os
import re
import sys
import time
import shutil
from pathlib import Path
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor

# Configure PyOpenCL (cache off, build chatter silenced) before
# apoc/pyclesperanto_prototype are imported below.
from behav3d.core.opencl_env import configure_pyopencl
configure_pyopencl()

import numpy as np
import pandas as pd
import zarr
import apoc
import pyclesperanto_prototype as cle

from behav3d.io.images import _ensure_zarr, load_image, save_as_zarr
from behav3d.core.metadata import (
    detect_organoid_types_from_metadata,
    detect_immune_cell_types_from_metadata,
    detect_other_cell_types_from_metadata,
    has_dead_channel,
)
from behav3d.preprocessing.segmentation.segment_journal import (
    file_fingerprint,
    journal_path,
    mark_done,
    new_journal,
    open_output_zarr,
    params_fingerprint,
    plan_output,
    preflight_conflicts,
    raise_for_conflicts,
    write_journal,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_classifier_header_value(opencl_path, key, default=None):
    """Read a single APOC metadata value from the header of a trained .cl file."""
    path = Path(opencl_path)
    if not path.exists():
        return default

    prefix = f"{key} = "
    with path.open() as f:
        line = ""
        count = 0
        while line != "*/" and count < 80:
            line = f.readline()
            if not line:
                break
            count += 1
            line = line.strip()
            if line.startswith(prefix):
                return line[len(prefix):].strip()
    return default


def _parse_channel_indices(value):
    """Parse a comma-separated list of channel indices."""
    if value is None:
        return []
    indices = []
    for token in str(value).split(","):
        token = token.strip()
        if not token:
            continue
        try:
            idx = int(token)
        except (TypeError, ValueError):
            continue
        if idx >= 0 and idx not in indices:
            indices.append(idx)
    return indices


def _parse_channel_names(value):
    """Parse a pipe-separated list of channel names."""
    if value is None:
        return []
    names = []
    for token in str(value).split("|"):
        name = token.strip()
        if name and name not in names:
            names.append(name)
    return names


def _channel_index_from_name(name):
    """Extract the integer index from a canonical `Channel N` layer name."""
    match = re.fullmatch(r"Channel\s+(\d+)", str(name or "").strip())
    if not match:
        return None
    try:
        return int(match.group(1))
    except (TypeError, ValueError):
        return None


def _channel_indices_from_names(names):
    """Convert canonical channel layer names to integer indices."""
    indices = []
    for name in names or []:
        idx = _channel_index_from_name(name)
        if idx is not None and idx not in indices:
            indices.append(idx)
    return indices


def _channel_indices_from_config(apoc_config, cell_type):
    """Read legacy channel selections from saved BEHAV3D config."""
    if not apoc_config:
        return []
    chan_names = apoc_config.get(f"apoc_{cell_type}_channels", [])
    if isinstance(chan_names, str):
        chan_names = [chan_names]
    return _channel_indices_from_names(chan_names)


def _resolve_classifier_channel_indices(opencl_path, cell_type, apoc_config=None):
    """Resolve classifier input channels, preferring embedded .cl metadata."""
    indices = _parse_channel_indices(
        _read_classifier_header_value(opencl_path, "input_channel_indices")
    )
    if indices:
        return indices

    channel_names = _parse_channel_names(
        _read_classifier_header_value(opencl_path, "input_channel_names")
    )
    indices = _channel_indices_from_names(channel_names)
    if indices:
        return indices

    indices = _channel_indices_from_config(apoc_config, cell_type)
    if indices:
        return indices

    return [0]


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


def _cfg_val(cfg, ct, key, default):
    """Look up a per-cell-type instance-segmentation postprocessing param.

    Two callers feed ``apoc_config`` here with different key conventions: the
    napari GUI (``behav3d/napari/_segmentation.py``) writes ``apoc_{ct}_{key}``,
    while the notebook widget (``behav3d/widgets/segmentation.py``) writes the
    unprefixed ``{ct}_{key}`` into its own ``pixel_classifier`` config. Prefer
    the prefixed key and fall back to the unprefixed one so both stay correct.

    A key that is *present but* ``None`` counts as missing: the GUI tab's
    ``get_config()`` emits every post-processing key on every tab, filling in
    ``None`` for the widgets its own strategy never created (an EDT/Watershed
    tab still reports ``peak_min_ratio: None``). Those nulls must not shadow
    the default, or the ``float()``/``int()`` casts downstream crash.
    """
    val = cfg.get(f"apoc_{ct}_{key}")
    if val is None:
        val = cfg.get(f"{ct}_{key}")
    return default if val is None else val


# ---------------------------------------------------------------------------
# Resume fingerprints
# ---------------------------------------------------------------------------
#
# Two fingerprints per cell type, because "Only Resegment" reuses the cached mask
# and rewrites only the instance labels. Changing a post-processing threshold must
# invalidate the segments without throwing away a mask that is still perfectly good.


def _mask_fingerprint(cfg, ct, clf, channels):
    """Hash everything that determines the *mask* for one cell type."""
    return params_fingerprint({
        "engine": "apoc",
        "stage": "mask",
        "classifier": file_fingerprint(getattr(clf, "opencl_file", None)),
        "channels": list(channels or []),
        "prob_mask_threshold": _cfg_val(cfg, ct, "prob_mask_threshold", 0.5),
        "opening_nr_pixels": _cfg_val(cfg, ct, "opening_nr_pixels", 0),
        "fill_holes": _cfg_val(cfg, ct, "fill_holes", True),
    })


def _segments_fingerprint(cfg, ct, clf, channels, strategy):
    """Hash the mask inputs *plus* everything that turns a mask into instances."""
    return params_fingerprint({
        "engine": "apoc",
        "stage": "segments",
        "mask": _mask_fingerprint(cfg, ct, clf, channels),
        "strategy": str(strategy),
        "prob_seed_threshold": _cfg_val(cfg, ct, "prob_seed_threshold", 0.8),
        "edt_threshold": _cfg_val(cfg, ct, "edt_threshold", 1.0),
        "segment_size_min": _cfg_val(cfg, ct, "segment_size_min", 10),
        "peak_min_distance": _cfg_val(cfg, ct, "peak_min_distance", 0),
        "peak_min_ratio": _cfg_val(cfg, ct, "peak_min_ratio", 0.35),
    })


def _death_fingerprint(clf_death, channels):
    """Hash the death mask inputs. The death mask has no instance stage."""
    return params_fingerprint({
        "engine": "apoc",
        "stage": "death",
        "classifier": file_fingerprint(getattr(clf_death, "opencl_file", None)),
        "channels": list(channels or []),
    })


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
    per_ct_strategies=None,  # dict {cell_type: strategy_string} for Advanced mode; None = use global apoc_strategy
    only_cell_types=None,    # list of cell-type strings to process; None = process all
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

    # --- Discover cell types from metadata ---
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

    # Optional cell-type filter (used by per-cell-type Run button in the GUI)
    if only_cell_types is not None:
        active_cell_types = [ct for ct in active_cell_types if ct in only_cell_types]

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
    def get_indices(ct, clf=None):
        opencl_path = getattr(clf, "opencl_file", None)
        if not opencl_path:
            return _channel_indices_from_config(apoc_config, ct) or [0]
        return _resolve_classifier_channel_indices(opencl_path, ct, apoc_config=apoc_config)

    clf_channels = {ct: get_indices(ct, classifiers.get(ct)) for ct in all_cell_types}
    death_channels = get_indices("dead", clf_death) if clf_death else []

    # --- Process each sample ---
    sample_names = metadata['sample_name'].unique()

    # ------------------------------------------------------------------
    # Fingerprint preflight
    # ------------------------------------------------------------------
    # Every mismatch is knowable before any frame is written, so find them all now
    # and fail once rather than aborting the batch half-way through.
    if not overwrite_existing:
        conflicts = []
        for sample_name in sample_names:
            row = metadata[metadata['sample_name'] == sample_name].iloc[0]
            raw_path = row.get('raw_image_path')
            if not raw_path or not Path(raw_path).exists():
                continue
            axis = row.get('dimension_order', "TCZYX")
            if not isinstance(axis, str) or not axis:
                axis = "TCZYX"
            try:
                probe = load_image(Path(raw_path), axis_order=axis)
                shape = (probe.shape[0],) + tuple(probe.shape[2:])
            except Exception:
                continue  # unreadable input is the main loop's problem to report
            img_dir = output_dir / "images" / sample_name
            entries = []
            for ct in active_cell_types:
                channels = clf_channels[ct] or [0]
                strategy = (per_ct_strategies or {}).get(ct, apoc_strategy)
                entries.append((
                    img_dir / f"{sample_name}_{ct}_segments.zarr", shape, "uint16",
                    _segments_fingerprint(apoc_config or {}, ct, classifiers[ct], channels, strategy),
                    f"{ct} segments for {sample_name}",
                ))
                if not only_segment:
                    entries.append((
                        img_dir / f"{sample_name}_{ct}_mask.zarr", shape, "uint16",
                        _mask_fingerprint(apoc_config or {}, ct, classifiers[ct], channels),
                        f"{ct} mask for {sample_name}",
                    ))
            if has_death and clf_death and not only_segment:
                entries.append((
                    img_dir / f"{sample_name}_mask_dead.zarr", shape, "uint16",
                    _death_fingerprint(clf_death, death_channels),
                    f"dead mask for {sample_name}",
                ))
            conflicts.extend(preflight_conflicts(
                entries,
                overwrite_existing=overwrite_existing,
                timepoint_range=timepoint_range,
            ))
        raise_for_conflicts(conflicts, "APOC")

    for sample_name in sample_names:
        t0 = time.time()
 
        # Output directory
        img_outdir = output_dir / "images" / sample_name

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

        # Spatial shape for output arrays — needed before the skip decisions, which
        # compare the shape on disk against what this run would write.
        spatial_shape = img.shape[2:]           # (Z, Y, X)  — img is (T, C, Z, Y, X)
        out_shape = (n_timepoints,) + spatial_shape

        seg_path_of = lambda ct: img_outdir / f"{sample_name}_{ct}_segments.zarr"
        mask_path_of = lambda ct: img_outdir / f"{sample_name}_{ct}_mask.zarr"
        death_path = img_outdir / f"{sample_name}_mask_dead.zarr"

        # ------------------------------------------------------------------
        # Resume / skip / overwrite decision, per output array
        # ------------------------------------------------------------------
        active_cts_for_sample = list(active_cell_types)
        has_death_for_sample = has_death

        seg_fps, mask_fps = {}, {}
        seg_plans, mask_plans = {}, {}
        remaining_cts = []
        for ct in active_cts_for_sample:
            channels = clf_channels[ct] or [0]
            strategy = (per_ct_strategies or {}).get(ct, apoc_strategy)
            mask_fps[ct] = _mask_fingerprint(apoc_config or {}, ct, classifiers[ct], channels)
            seg_fps[ct] = _segments_fingerprint(
                apoc_config or {}, ct, classifiers[ct], channels, strategy
            )
            seg_plans[ct] = plan_output(
                seg_path_of(ct), out_shape, "uint16", seg_fps[ct],
                f"{ct} segments for {sample_name}",
                overwrite_existing=overwrite_existing,
                timepoint_range=timepoint_range,
                requested_timepoints=t_range,
            )
            # Under only_segment the mask is an input, never rewritten, so it gets
            # no plan of its own — it is read as-is and left untouched.
            mask_plans[ct] = None if only_segment else plan_output(
                mask_path_of(ct), out_shape, "uint16", mask_fps[ct],
                f"{ct} mask for {sample_name}",
                overwrite_existing=overwrite_existing,
                timepoint_range=timepoint_range,
                requested_timepoints=t_range,
            )
            mask_complete = only_segment or mask_plans[ct].complete
            if seg_plans[ct].complete and mask_complete:
                continue  # nothing left to do for this cell type
            remaining_cts.append(ct)
        active_cts_for_sample = remaining_cts

        death_plan = None
        if has_death_for_sample and clf_death and not only_segment:
            death_fp = _death_fingerprint(clf_death, death_channels)
            death_plan = plan_output(
                death_path, out_shape, "uint16", death_fp,
                f"dead mask for {sample_name}",
                overwrite_existing=overwrite_existing,
                timepoint_range=timepoint_range,
                requested_timepoints=t_range,
            )
            if death_plan.complete:
                has_death_for_sample = False
                death_plan = None
        else:
            has_death_for_sample = False

        if not active_cts_for_sample and not has_death_for_sample:
            print(f"  ⏭️ Skipping {sample_name} (all selected outputs already complete)")
            continue

        _ensure_zarr(raw_image_path, label=f"Raw image for '{sample_name}'")

        img_outdir.mkdir(parents=True, exist_ok=True)

        # Create output zarr arrays and their journals
        zarr_segs = {}
        zarr_masks = {}
        seg_journals = {}
        mask_journals = {}
        for ct in active_cts_for_sample:
            zarr_segs[ct] = open_output_zarr(
                seg_path_of(ct), "uint16", out_shape, seg_plans[ct].recreate,
            )
            seg_journals[ct] = new_journal(
                seg_fps[ct], out_shape, "uint16", done=seg_plans[ct].done
            )
            write_journal(journal_path(seg_path_of(ct)), seg_journals[ct])

            if only_segment:
                # Open mask in read-only mode, we need it to do resegmentation
                mask_path = mask_path_of(ct)
                if not mask_path.exists():
                    raise FileNotFoundError(f"Cannot 'Only Resegment' because mask does not exist: {mask_path}")
                zarr_masks[ct] = zarr.open(str(mask_path), mode="r")
            else:
                zarr_masks[ct] = open_output_zarr(
                    mask_path_of(ct), "uint16", out_shape, mask_plans[ct].recreate,
                )
                mask_journals[ct] = new_journal(
                    mask_fps[ct], out_shape, "uint16", done=mask_plans[ct].done
                )
                write_journal(journal_path(mask_path_of(ct)), mask_journals[ct])

        zarr_death = None
        death_journal = None
        if has_death_for_sample:
            zarr_death = open_output_zarr(
                death_path, "uint16", out_shape, death_plan.recreate,
            )
            death_journal = new_journal(
                death_fp, out_shape, "uint16", done=death_plan.done
            )
            write_journal(journal_path(death_path), death_journal)

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
                for ct in active_cts_for_sample:
                    # The journals are authoritative: under overwrite they were
                    # started empty, so this gate covers both resume and redo.
                    need_seg = t not in seg_journals[ct]["done"]
                    need_mask = (not only_segment) and t not in mask_journals[ct]["done"]
                    if not need_seg and not need_mask:
                        continue
                    # Instances are derived from the probability map, which is never
                    # persisted, so segments cannot be produced without re-running the
                    # mask prediction. When either is missing, both get rewritten.

                    cfg = apoc_config or {}
                    opencl_path = getattr(classifiers[ct], "opencl_file", None)
                    n_clf_ch_raw = _read_classifier_header_value(opencl_path, "n_image_channels")
                    if n_clf_ch_raw is not None:
                        try:
                            n_clf_ch = int(n_clf_ch_raw)
                            if n_clf_ch != t_img.shape[0]:
                                raise ValueError(
                                    f"Classifier for '{ct}' was trained on an image with {n_clf_ch} channel(s), "
                                    f"but the current image has {t_img.shape[0]} channel(s). "
                                    f"Retrain the classifier on an image with the same number of channels."
                                )
                        except (TypeError, ValueError) as _e:
                            if "channel" in str(_e).lower():
                                raise
                    indices = clf_channels[ct] or [0]
                    imgs = [t_img[i_ch] for i_ch in indices]
                    imgs_to_pass = imgs[0] if len(imgs) == 1 else imgs

                    # Resolve effective strategy: per-cell-type override takes precedence when provided
                    effective_strategy = (per_ct_strategies or {}).get(ct, apoc_strategy)

                    if effective_strategy == "APOC Probability Map + Watershed":
                        # ── Probability Map strategy ──────────────────────
                        from skimage.measure import label as sk_label
                        from skimage.segmentation import watershed
                        from behav3d.preprocessing.segmentation import segment_size_filter, segment_2d_filter
                        from behav3d.preprocessing.segmentation.segmentation_utils import postprocess_mask, segment_mask
                        _diag = (i == 0)  # print timings on first timepoint only

                        opening_nr_pixels = int(_cfg_val(cfg, ct, "opening_nr_pixels", 0))
                        if not only_segment:
                            _t0 = time.time()
                            prob_clf = _get_probability_classifier(classifiers[ct])
                            prob_map = np.asarray(
                                prob_clf.predict(image=imgs_to_pass)
                            ).astype(np.float32)
                            if _diag: print(f"    ⏱ {ct} predict (prob): {time.time()-_t0:.2f}s  range=[{prob_map.min():.3f}, {prob_map.max():.3f}]")

                            _t0 = time.time()
                            mask_thr = float(_cfg_val(cfg, ct, "prob_mask_threshold", 0.5))
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

                        seed_thr = float(_cfg_val(cfg, ct, "prob_seed_threshold", 0.8))
                        segment_size_min = int(_cfg_val(cfg, ct, "segment_size_min", 10))
                        if _diag: print(f"    ⚙ {ct} params: mask_thr={mask_thr if not only_segment else 'n/a'}, seed_thr={seed_thr}, min_size={segment_size_min}, opening_px={opening_nr_pixels}")

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

                    elif effective_strategy in {
                        "APOC Mask + EDT/Watershed Resegmentation",
                        "APOC Mask + Peak EDT/Watershed Resegmentation",
                    }:
                        # ── EDT/Watershed strategy ────────────────────────
                        from behav3d.preprocessing.segmentation.segmentation_utils import (
                            postprocess_mask, segment_mask
                        )
                        _diag = (i == 0)  # print timings/params on first timepoint only
                        if not only_segment:
                            seg_out = np.asarray(classifiers[ct].predict(image=imgs_to_pass)).astype(np.uint16)
                            mask_out = (seg_out > 0).astype(np.uint16)
                            zarr_masks[ct][t] = mask_out
                        else:
                            mask_out = np.asarray(zarr_masks[ct][t])

                        fill_holes = _cfg_val(cfg, ct, "fill_holes", True)
                        opening_nr_pixels = _cfg_val(cfg, ct, "opening_nr_pixels", 0)
                        proc_mask = postprocess_mask(mask_out, fill_holes=fill_holes, opening_nr_pixels=opening_nr_pixels)

                        edt_thr = _cfg_val(cfg, ct, "edt_threshold", 1.0)
                        segment_size_min = _cfg_val(cfg, ct, "segment_size_min", 10)
                        _raw_peak_dist = _cfg_val(cfg, ct, "peak_min_distance", 0)
                        peak_min_distance = None if (not _raw_peak_dist or float(_raw_peak_dist) <= 0) else float(_raw_peak_dist)
                        peak_min_ratio = float(_cfg_val(cfg, ct, "peak_min_ratio", 0.35))
                        if _diag: print(f"    ⚙ {ct} params: edt_thr={edt_thr}, min_size={segment_size_min}, fill_holes={fill_holes}, opening_px={opening_nr_pixels}, peak_min_distance={peak_min_distance}, peak_min_ratio={peak_min_ratio}")
                        seg_refined = segment_mask(
                            proc_mask,
                            edt_thr=edt_thr,
                            edt_thr_refined=None,
                            segment_size_min=segment_size_min,
                            use_dims=3,
                            n_workers=1,
                            marker_strategy=(
                                "peak"
                                if apoc_strategy == "APOC Mask + Peak EDT/Watershed Resegmentation"
                                else "threshold"
                            ),
                            peak_min_distance=peak_min_distance,
                            peak_min_ratio=peak_min_ratio,
                        )
                        zarr_segs[ct][t] = np.asarray(seg_refined).astype(np.uint16)

                    else:
                        # ── Default: Direct APOC (or unknown strategy) ────
                        if i == 0:
                            print(f"    ⚙ {ct}: Direct APOC strategy — no postprocessing parameters applied")
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

                    # Mark this timepoint done, per array, before moving on: the
                    # failure this guards against gives no chance to flush later.
                    if not only_segment:
                        mark_done(mask_journals[ct], journal_path(mask_path_of(ct)), t)
                    mark_done(seg_journals[ct], journal_path(seg_path_of(ct)), t)

                if zarr_death is not None:
                    if t in death_journal["done"]:
                        pass  # Skip timepoint for death mask if already computed
                    else:
                        n_death_ch_raw = _read_classifier_header_value(clf_death.opencl_file, "n_image_channels")
                        if n_death_ch_raw is not None:
                            try:
                                n_death_ch = int(n_death_ch_raw)
                                if n_death_ch != t_img.shape[0]:
                                    raise ValueError(
                                        f"Death classifier was trained on an image with {n_death_ch} channel(s), "
                                        f"but the current image has {t_img.shape[0]} channel(s). "
                                        f"Retrain the classifier on an image with the same number of channels."
                                    )
                            except (TypeError, ValueError) as _e:
                                if "channel" in str(_e).lower():
                                    raise
                        indices = death_channels or [0]
                        imgs = [t_img[i_ch] for i_ch in indices]
                        imgs_to_pass = imgs[0] if len(imgs) == 1 else imgs
                        death_mask = (np.asarray(clf_death.predict(image=imgs_to_pass)) > 0).astype(np.uint16)
                        zarr_death[t] = death_mask
                        mark_done(death_journal, journal_path(death_path), t)

        # Update metadata with output paths
        row_idx = metadata.index[metadata['sample_name'] == sample_name].tolist()[0]
        for ct in active_cts_for_sample:
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
