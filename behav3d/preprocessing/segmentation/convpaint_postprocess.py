"""
Shared ConvPaint post-processing helpers.

Both the live training widget (:mod:`behav3d.preprocessing.segmentation.convpaint_train`)
and the batch queue (:mod:`behav3d.preprocessing.segmentation.convpaint_segment`) use
these helpers to convert pixel-classifier outputs into instance segmentations,
keeping a single source of truth for the segmentation strategies.

`ConvpaintModel` only exposes pixel-level outputs (`segment` returns class labels,
`predict_probas` returns per-class probabilities), so instances must always be
derived in post-processing.  Two strategies are supported:

- ``mask_to_instances``  → EDT/Watershed resegmentation of the binary mask.
  Supports threshold markers and peak EDT markers.
- ``probability_to_instances`` → seeded watershed of the raw probability map.

The "Direct" strategy that previously used connected components on the binary
mask was removed because it cannot separate touching cells and was therefore
misleading to users.
"""

from __future__ import annotations

import io
import contextlib
import numpy as np

_SUPPRESSED_CONVPAINT_MESSAGES = (
    "FE model is designed for imagenet normalization",
    "not declared as 'rgb' (parameter channel_mode)",
)


# ---------------------------------------------------------------------------
# Memory-safe per-Z prediction wrappers (used by both the widget and the batch
# queue).  ConvPaint feature extractors expand each pixel into hundreds of
# features, so processing a full 3-D volume at once would allocate
# ``(n_features, Z, Y, X)`` and easily exceed available RAM on typical
# microscopy data.  Iterating over Z keeps peak memory to one slice.
# ---------------------------------------------------------------------------


def _call_deduplicated(fn, *args, **kwargs):
    """Call ``fn(*args, **kwargs)`` while suppressing duplicate print lines.

    ConvPaint emits informational messages (e.g. the imagenet-normalization
    notice) once per call, which means once per Z-slice when we iterate.
    This wrapper captures stdout for each call and only forwards lines that
    have not been seen before within the current loop, so each message
    appears at most once regardless of Z depth.

    Returns ``(result, seen)`` where ``seen`` is the set of already-printed
    lines; pass it in subsequent calls to accumulate across slices.
    """
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        result = fn(*args, **kwargs)
    return result, buf.getvalue()


def _forward_unique_messages(output, seen):
    """Print unique ConvPaint lines, skipping known noisy notices."""
    for line in output.splitlines():
        if not line:
            continue
        if any(msg in line for msg in _SUPPRESSED_CONVPAINT_MESSAGES):
            continue
        if line in seen:
            continue
        seen.add(line)
        print(line)


def segment_per_z(model, frame, fe_use_device=None):
    """Segment a ``(C, Z, Y, X)`` frame slice-by-slice and stack over Z.

    Returns a ``(Z, Y, X)`` array of class labels (1 = background, 2 = foreground
    for the binary annotations used in BEHAV3D).
    """
    if frame.ndim == 4 and frame.shape[1] > 1:
        z_slices = []
        seen: set[str] = set()
        for z in range(frame.shape[1]):
            seg_z, output = _call_deduplicated(
                model.segment, frame[:, z, :, :], fe_use_device=fe_use_device
            )
            z_slices.append(np.asarray(seg_z))
            _forward_unique_messages(output, seen)
        return np.stack(z_slices, axis=0)
    seg, output = _call_deduplicated(
        model.segment, frame, fe_use_device=fe_use_device
    )
    _forward_unique_messages(output, set())
    return np.asarray(seg)


def predict_probas_per_z(model, frame, fe_use_device=None):
    """Predict probability maps for a ``(C, Z, Y, X)`` frame slice-by-slice.

    Returns a ``(n_classes, Z, Y, X)`` array.
    """
    if frame.ndim == 4 and frame.shape[1] > 1:
        slices = []
        seen: set[str] = set()
        for z in range(frame.shape[1]):
            p, output = _call_deduplicated(
                model.predict_probas, frame[:, z, :, :], fe_use_device=fe_use_device
            )
            slices.append(np.asarray(p))
            _forward_unique_messages(output, seen)
        return np.stack(slices, axis=1)
    probas, output = _call_deduplicated(
        model.predict_probas, frame, fe_use_device=fe_use_device
    )
    _forward_unique_messages(output, set())
    return np.asarray(probas)


# ---------------------------------------------------------------------------
# Instance-segmentation helpers.
# Both helpers accept 3-D ``(Z, Y, X)`` or 4-D ``(T, Z, Y, X)`` arrays and
# always return ``uint16`` instance labels with size + 2D filtering applied.
# ---------------------------------------------------------------------------

def _finalize_segments(segments, segment_size_min):
    from behav3d.preprocessing.segmentation import (
        segment_size_filter,
        segment_2d_filter,
    )

    segments = np.asarray(segments).copy()
    segments = segment_size_filter(segments, size_min=segment_size_min)
    if segments.ndim == 3:
        segments = segment_2d_filter(segments)
    return np.asarray(segments).astype(np.uint16)


def _mask_volume_to_instances(mask, edt_thr, opening_nr_pixels, fill_holes,
                              segment_size_min, marker_strategy="threshold",
                              peak_min_distance=None, peak_min_ratio=0.35):
    from behav3d.preprocessing.segmentation.segmentation_utils import (
        postprocess_mask,
        segment_mask,
    )

    proc_mask = postprocess_mask(
        np.asarray(mask).astype(bool),
        fill_holes=fill_holes,
        opening_nr_pixels=opening_nr_pixels,
    )
    instances = segment_mask(
        proc_mask,
        edt_thr=edt_thr,
        edt_thr_refined=None,
        segment_size_min=segment_size_min,
        use_dims=3,
        n_workers=1,
        marker_strategy=marker_strategy,
        peak_min_distance=peak_min_distance,
        peak_min_ratio=peak_min_ratio,
    )
    return _finalize_segments(instances, segment_size_min)


def mask_to_instances(mask, *, edt_thr, opening_nr_pixels, fill_holes,
                      segment_size_min, marker_strategy="threshold",
                      peak_min_distance=None, peak_min_ratio=0.35):
    """EDT/watershed instance segmentation of a binary mask.

    Parameters mirror the per-cell-type spinners in the napari widget /
    notebook controls.  ``mask`` may be 3-D ``(Z, Y, X)`` or 4-D
    ``(T, Z, Y, X)``.

    ``peak_min_distance`` and ``peak_min_ratio`` are only used when
    ``marker_strategy="peak"``.  Pass ``peak_min_distance=None`` (or ``0``) to
    use the automatic default of ``2 × edt_thr``.
    """
    mask = np.asarray(mask)
    if mask.ndim == 4:
        return np.stack(
            [
                _mask_volume_to_instances(
                    mask[t], edt_thr, opening_nr_pixels, fill_holes,
                    segment_size_min, marker_strategy,
                    peak_min_distance, peak_min_ratio,
                )
                for t in range(mask.shape[0])
            ],
            axis=0,
        )
    return _mask_volume_to_instances(
        mask, edt_thr, opening_nr_pixels, fill_holes, segment_size_min,
        marker_strategy, peak_min_distance, peak_min_ratio,
    )


def _probability_volume_to_instances(prob_map, mask_thr, seed_thr,
                                     opening_nr_pixels, segment_size_min):
    from scipy import ndimage as ndi
    from skimage.measure import label as sk_label
    from skimage.segmentation import watershed
    from behav3d.preprocessing.segmentation.segmentation_utils import (
        postprocess_mask,
    )

    prob_map = np.asarray(prob_map).astype(np.float32)
    mask_out = postprocess_mask(
        prob_map > mask_thr,
        fill_holes=False,
        opening_nr_pixels=opening_nr_pixels,
    ).astype(bool)

    cc_labels = sk_label(mask_out)
    seed_mask = (prob_map > seed_thr) & mask_out
    segments = np.zeros_like(cc_labels, dtype=np.uint16)
    next_id = 0

    for comp_idx, slc in enumerate(ndi.find_objects(cc_labels), start=1):
        if slc is None:
            continue
        comp_mask = cc_labels[slc] == comp_idx
        sub_seeds = sk_label(seed_mask[slc] & comp_mask)
        n_seeds = int(sub_seeds.max())
        if n_seeds <= 1:
            next_id += 1
            segments[slc][comp_mask] = next_id
            continue
        sub_result = watershed(
            -prob_map[slc], markers=sub_seeds, mask=comp_mask,
        )
        for seed_id in range(1, int(sub_result.max()) + 1):
            next_id += 1
            segments[slc][sub_result == seed_id] = next_id

    return _finalize_segments(segments, segment_size_min)


def probability_to_instances(prob_map, *, mask_thr, seed_thr,
                             opening_nr_pixels, segment_size_min):
    """Seeded watershed instance segmentation of a probability volume.

    ``prob_map`` is the foreground-class probability and may be 3-D
    ``(Z, Y, X)`` or 4-D ``(T, Z, Y, X)``.
    """
    prob_map = np.asarray(prob_map)
    if prob_map.ndim == 4:
        return np.stack(
            [
                _probability_volume_to_instances(
                    prob_map[t], mask_thr, seed_thr,
                    opening_nr_pixels, segment_size_min,
                )
                for t in range(prob_map.shape[0])
            ],
            axis=0,
        )
    return _probability_volume_to_instances(
        prob_map, mask_thr, seed_thr, opening_nr_pixels, segment_size_min,
    )


# ---------------------------------------------------------------------------
# Helpers for extracting the foreground-class probability from
# ``ConvpaintModel.predict_probas`` output.
# ---------------------------------------------------------------------------

def foreground_probability(probas):
    """Return the foreground-class probability volume from a probas array.

    ConvPaint encodes class 1 = background, class 2 = foreground (matching
    the BEHAV3D annotation convention).  ``predict_probas`` returns either
    ``(n_classes, Z, Y, X)`` or ``(Z, Y, X)`` (single-class edge case).
    """
    probas = np.asarray(probas)
    if probas.ndim == 3:
        return probas.astype(np.float32)
    return probas[-1].astype(np.float32)
