"""Pure backend primitives for editing 4D tracked-segment volumes.

Each operation works on numpy arrays and returns the modified frames as a
``dict[t, np.ndarray]`` plus an ``OpResult`` describing what changed.  The
:class:`~behav3d.editing.edit_buffer.EditBuffer` wraps these into committed
dirty frames.  The primitives can also be called directly from
scripts/tests with a plain ``EditBuffer``-like object.

Coordinate convention: every volume is ``T x Z x Y x X`` with integer label
IDs and ``0`` reserved for background.  All primitives are restricted to one
or more user-selected labels — they never touch unrelated pixels.

The ``split_label`` propagation reuses
``_run_single_timepoint_propagation`` from
``behav3d.preprocessing.tracking.propagation_tracking`` so that splits behave
identically to propagation tracking (watershed of the label mask seeded by
the previous frame's sub-labels).

Operations that process each timepoint independently (``merge_labels``,
``delete_label``, ``delete_labels``, ``fill_label``) accept an optional
``n_workers`` argument and dispatch frames via :class:`ThreadPoolExecutor`
(numpy/skimage ops release the GIL so true parallel CPU utilisation is
achieved).  When ``n_workers`` is ``None`` the function uses one worker per
frame up to ``os.cpu_count()``.
"""
from __future__ import annotations

import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from functools import partial
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from skimage.filters import threshold_otsu
from skimage.morphology import binary_erosion
from skimage.segmentation import expand_labels, watershed

from behav3d.preprocessing.tracking.propagation_tracking import (
    _run_single_timepoint_propagation,
)


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------
@dataclass
class OpResult:
    """Description of a primitive's effect on the tracked-segments volume.

    Attributes
    ----------
    name:
        Human-readable name of the operation (e.g. ``"split"``).
    new_frames:
        Mapping ``t -> np.ndarray`` containing only the timepoints that
        actually changed.  Frames that were inspected but unchanged are
        deliberately excluded so the buffer's dirty set stays minimal.
    affected_labels:
        Set of label IDs that exist in the result and were touched (existing
        IDs reassigned, new IDs created, etc.).
    new_labels:
        Newly-introduced label IDs (subset of ``affected_labels``).
    summary:
        Short, single-line message intended for the editor's log box.
    """
    name: str
    new_frames: Dict[int, np.ndarray] = field(default_factory=dict)
    affected_labels: List[int] = field(default_factory=list)
    new_labels: List[int] = field(default_factory=list)
    summary: str = ""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _normalize_t_range(buf, t_range) -> Tuple[int, int]:
    """Clamp ``t_range`` to ``[0, T-1]`` and return ``(t0, t1)`` inclusive."""
    T = int(buf.shape[0])
    if t_range is None:
        return 0, T - 1
    t0, t1 = int(t_range[0]), int(t_range[1])
    if t1 < t0:
        t0, t1 = t1, t0
    return max(0, t0), min(T - 1, t1)


def _next_label_id(buf, taken: Iterable[int] = ()) -> int:
    """Return the next available label ID (max + 1 across the whole volume).

    Uses :meth:`~behav3d.editing.edit_buffer.EditBuffer.max_label` when
    available (O(1) after the first call) instead of a full unique-label
    scan, making repeated splits on large datasets much faster.
    """
    taken_max = max((int(x) for x in taken), default=0)
    if hasattr(buf, "max_label"):
        return max(buf.max_label(), taken_max) + 1
    # Fallback for buf objects that don't expose max_label (tests etc.).
    used = set(int(x) for x in taken)
    used.add(0)
    used.update(int(x) for x in buf.unique_labels())
    return max(used) + 1


def lifetime_of(buf, label_id: int) -> Tuple[Optional[int], Optional[int]]:
    """Return ``(t_first, t_last)`` inclusive frames where ``label_id`` exists.

    Delegates to :meth:`~behav3d.editing.edit_buffer.EditBuffer.get_lifetime`
    when available so repeated calls for the same label (e.g. every click)
    are served from a cache instead of scanning all T frames each time.

    Returns ``(None, None)`` if the label is not present anywhere.
    """
    if hasattr(buf, "get_lifetime"):
        return buf.get_lifetime(label_id)
    # Fallback for buf objects that don't expose get_lifetime (tests etc.).
    label_id = int(label_id)
    T = int(buf.shape[0])
    first = None
    last = None
    for t in range(T):
        frame = buf.peek(t)
        if (frame == label_id).any():
            if first is None:
                first = t
            last = t
    return first, last


def _binary_erosion_gpu(mask: np.ndarray, footprint: np.ndarray) -> np.ndarray:
    """GPU-accelerated binary erosion via PyTorch ``conv3d``.

    Erosion is equivalent to a 3D convolution of the float mask with the
    footprint kernel: a voxel is kept eroded (True) only if every voxel in
    its footprint neighbourhood was also True, i.e. ``conv == sum(footprint)``.

    Parameters
    ----------
    mask:
        Boolean array of shape ``(Z, Y, X)``.
    footprint:
        Boolean structuring element, e.g. from :func:`_ellipsoid_footprint`.

    Returns
    -------
    np.ndarray
        Boolean array, same shape as ``mask``.
    """
    import torch
    import torch.nn.functional as F

    device = torch.device("cuda")
    # (1, 1, Z, Y, X) float32 tensors on GPU.
    x = torch.from_numpy(mask.astype(np.float32)).to(device).unsqueeze(0).unsqueeze(0)
    k = torch.from_numpy(footprint.astype(np.float32)).to(device).unsqueeze(0).unsqueeze(0)
    pad = tuple(s // 2 for s in footprint.shape)  # (pZ, pY, pX)
    conv = F.conv3d(x, k, padding=pad).squeeze(0).squeeze(0)
    target = float(footprint.sum())
    result = (conv >= target - 1e-3).cpu().numpy()
    return result


# Resolved once at import time: use GPU path when CUDA is available and
# PyTorch is installed; fall back to skimage CPU path otherwise.
def _binary_erosion_impl(mask: np.ndarray, footprint: np.ndarray) -> np.ndarray:
    """Apply binary erosion to ``mask`` using ``footprint``.

    Routes to :func:`_binary_erosion_gpu` when CUDA is available (RTX/Tesla
    class GPU gives a 5–10× speed-up per frame for large volumes) and falls
    back to :func:`skimage.morphology.binary_erosion` on CPU otherwise.
    """
    return _EROSION_FN(mask, footprint)


def _skimage_erosion(mask: np.ndarray, footprint: np.ndarray) -> np.ndarray:
    return binary_erosion(mask, footprint=footprint)


def _pick_erosion_backend():
    """Return the best available erosion function at import time."""
    try:
        import torch
        if torch.cuda.is_available():
            # Warm the CUDA context with a tiny dry-run so the first real call
            # is not delayed by driver initialisation.
            import torch.nn.functional as F
            _dummy = torch.zeros(1, 1, 1, 1, 1, device="cuda")
            F.conv3d(_dummy, _dummy)
            del _dummy
            return _binary_erosion_gpu
    except Exception:
        pass
    return _skimage_erosion


_EROSION_FN = _pick_erosion_backend()


def _ellipsoid_footprint(r_xy: int, r_z: int) -> np.ndarray:
    """Discrete (2*r_z+1, 2*r_xy+1, 2*r_xy+1) ellipsoid mask for 3D ops."""
    r_xy = max(0, int(r_xy))
    r_z = max(0, int(r_z))
    if r_xy == 0 and r_z == 0:
        return np.ones((1, 1, 1), dtype=bool)
    if r_z == 0:
        # 2D disk replicated in z
        rr = np.arange(-r_xy, r_xy + 1)
        yy, xx = np.meshgrid(rr, rr, indexing="ij")
        disk = (yy * yy + xx * xx) <= r_xy * r_xy
        return disk[None, :, :]
    if r_xy == 0:
        zz = np.arange(-r_z, r_z + 1)
        return (np.abs(zz)[:, None, None] <= r_z)
    zz = np.arange(-r_z, r_z + 1)[:, None, None]
    yy = np.arange(-r_xy, r_xy + 1)[None, :, None]
    xx = np.arange(-r_xy, r_xy + 1)[None, None, :]
    norm = (zz * zz) / max(1, r_z * r_z) + (yy * yy + xx * xx) / max(1, r_xy * r_xy)
    return norm <= 1.0


# ---------------------------------------------------------------------------
# Primitives
# ---------------------------------------------------------------------------
def split_label(
    buf,
    label_id: int,
    seeds_zyx: Sequence[Tuple[int, int, int]],
    t_ref: int,
    t_range: Optional[Tuple[int, int]] = None,
    keep_first_id: bool = True,
    new_ids: Optional[Sequence[int]] = None,
    progress_cb: Optional[Callable[[int, int], None]] = None,
) -> OpResult:
    """Split a single label into N sub-labels seeded by user clicks.

    Workflow at the *reference frame* ``t_ref``:

    1. Restrict to the binary mask ``vol[t_ref] == label_id``.
    2. Drop a marker for each seed (``(z, y, x)``); each seed must lie
       inside the mask.
    3. Run a 3D watershed seeded by those markers.
    4. The sub-label that contains the *first* seed keeps the original
       ``label_id`` (when ``keep_first_id=True``); the other(s) get new IDs
       from ``new_ids`` if supplied or are auto-allocated as
       ``max_label_id + k``.

    Then, propagating frame-by-frame across ``t_range`` (defaulting to the
    full lifetime of ``label_id``), the sub-labels from the previous frame
    are used as markers for ``_run_single_timepoint_propagation`` restricted
    to the next frame's ``label_id`` mask, so each sub-label keeps its
    identity over time.
    """
    label_id = int(label_id)
    if len(seeds_zyx) < 2:
        raise ValueError("split_label needs at least 2 seeds")

    T = int(buf.shape[0])
    if not 0 <= int(t_ref) < T:
        raise IndexError(f"t_ref={t_ref} out of bounds [0, {T - 1}]")

    t0, t1 = _normalize_t_range(buf, t_range)
    if not (t0 <= t_ref <= t1):
        raise ValueError(
            f"t_ref={t_ref} must lie within t_range=[{t0}, {t1}]"
        )

    # Allocate IDs for each seed; the first seed keeps label_id when asked.
    n_seeds = len(seeds_zyx)
    if new_ids is not None:
        new_ids = [int(x) for x in new_ids]
        if len(new_ids) != (n_seeds - 1 if keep_first_id else n_seeds):
            raise ValueError(
                "new_ids length must match the number of children that need a new ID"
            )
    else:
        # Auto-allocate sequentially starting at max+1.
        next_id = _next_label_id(buf, taken=[label_id])
        if keep_first_id:
            new_ids = [next_id + i for i in range(n_seeds - 1)]
        else:
            new_ids = [next_id + i for i in range(n_seeds)]
    if keep_first_id:
        seed_ids = [label_id] + new_ids
    else:
        seed_ids = list(new_ids)

    # ---- Reference frame split ------------------------------------------
    ref = buf.peek(t_ref).copy()
    mask = ref == label_id
    if not mask.any():
        raise ValueError(
            f"label {label_id} is empty at t_ref={t_ref}; nothing to split"
        )

    seeds_img = np.zeros_like(ref)
    for sid, (zz, yy, xx) in zip(seed_ids, seeds_zyx):
        zz = int(np.clip(zz, 0, ref.shape[0] - 1))
        yy = int(np.clip(yy, 0, ref.shape[1] - 1))
        xx = int(np.clip(xx, 0, ref.shape[2] - 1))
        if not mask[zz, yy, xx]:
            raise ValueError(
                f"seed ({zz}, {yy}, {xx}) is not inside label {label_id} at t_ref={t_ref}"
            )
        seeds_img[zz, yy, xx] = sid

    sub = watershed(np.zeros_like(ref, dtype=np.float32), markers=seeds_img, mask=mask)
    if not sub.any():
        raise RuntimeError("watershed produced an empty result; check seed positions")

    new_frames: Dict[int, np.ndarray] = {}
    out_ref = ref.copy()
    out_ref[mask] = sub[mask]
    new_frames[t_ref] = out_ref

    # Progress tracking: reference frame counts as step 1.
    _total_frames = t1 - t0 + 1
    _done = [1]
    if progress_cb:
        progress_cb(1, _total_frames)

    # ---- Propagation forward and backward in time -----------------------
    def _propagate(direction: int, start: int, end: int):
        if direction == 0:
            return
        prev = out_ref.copy()
        step = direction
        t = start
        while (step > 0 and t <= end) or (step < 0 and t >= end):
            frame = buf.peek(t).copy()
            mask_t = frame == label_id
            if not mask_t.any():
                # Track has gaps — stop propagating in this direction.
                break
            # Restrict prev sub-labels to current mask (so they only seed
            # pixels still inside the parent).
            prev_seg = np.where(np.isin(prev, seed_ids), prev, 0).astype(prev.dtype)
            prev_inside = np.where(mask_t, prev_seg, 0)
            # Some labels from prev may not survive into this frame (cell
            # disappears) — in that case, only the surviving children get
            # propagated.  Watershed needs at least one marker.
            if not (prev_inside > 0).any():
                # Fallback: use the previous frame's sub-label centroids
                # projected back into the new mask.
                break
            t_seg_label = mask_t.astype(prev.dtype)
            new_sub = _run_single_timepoint_propagation(
                t_seg_label, prev_inside.copy(),
                dilation_nr_pixels=2,
                segment_size_min=0,
            )
            out = frame.copy()
            out[mask_t] = new_sub[mask_t]
            new_frames[t] = out
            prev = out
            _done[0] += 1
            if progress_cb:
                progress_cb(_done[0], _total_frames)
            t += step

    _propagate(direction=+1, start=t_ref + 1, end=t1)
    _propagate(direction=-1, start=t_ref - 1, end=t0)

    affected = sorted(set(seed_ids))
    return OpResult(
        name="split",
        new_frames=new_frames,
        affected_labels=affected,
        new_labels=[i for i in seed_ids if i != label_id],
        summary=(
            f"Split label {label_id} at t={t_ref} into {affected}; "
            f"propagated across {len(new_frames)} frame(s)"
        ),
    )


def create_label(
    buf,
    seeds_zyx: Sequence[Tuple[int, int, int]],
    t_ref: int,
    t_range: Optional[Tuple[int, int]] = None,
    new_ids: Optional[Sequence[int]] = None,
    image_ref: Optional[np.ndarray] = None,
    image_stack=None,
    progress_cb: Optional[Callable[[int, int], None]] = None,
) -> OpResult:
    """Create N new labels in unlabelled background pixels, seeded by user clicks.

    Workflow at the *reference frame* ``t_ref``:

    1. Restrict to the binary mask ``vol[t_ref] == 0`` (background only).
       When ``image_ref`` is provided the mask is further restricted to voxels
       whose intensity is at or above the automatic Otsu threshold of
       ``image_ref``, so the new label cannot bleed into truly empty space.
    2. Drop one marker per seed ``(z, y, x)``; each seed must lie inside the
       (possibly tightened) background mask.
    3. Run a 3D watershed.  When ``image_ref`` is provided the watershed image
       is ``image_ref.max() - image_ref`` (negated), so bright voxels are
       low-cost basins and boundaries form naturally at dim regions between
       organoids.  Without ``image_ref`` a flat (zeros) image is used, giving
       a Voronoi-style partition of background pixels.
    4. Each seed receives a new, auto-allocated label ID.

    Then, propagating frame-by-frame across ``t_range`` (defaulting to the
    full time axis when omitted), the new sub-labels from the previous frame
    are used as markers for an intensity-guided watershed (when
    ``image_stack`` is provided) or a binary watershed restricted to the
    background pixels of the next frame, so the new labels follow the
    unlabelled region over time.

    Parameters
    ----------
    image_ref:
        Optional 3-D ``(Z, Y, X)`` intensity array for the reference frame
        (e.g. the selected raw fluorescence channel).  When supplied, the
        watershed is guided by the image and the eligible mask is restricted
        to bright (above-Otsu) voxels.
    image_stack:
        Optional 4-D array or dask array of shape ``(T, Z, Y, X)`` for the
        same fluorescence channel across all timepoints.  When supplied,
        each propagation frame uses the same Otsu + negated-image watershed
        strategy as ``t_ref``, preventing the new label from flooding the
        entire background in frames where no image guidance would otherwise
        apply.
    """
    n_seeds = len(seeds_zyx)
    if n_seeds < 1:
        raise ValueError("create_label needs at least 1 seed")

    T = int(buf.shape[0])
    t_ref = int(t_ref)
    if not 0 <= t_ref < T:
        raise IndexError(f"t_ref={t_ref} out of bounds [0, {T - 1}]")

    t0, t1 = _normalize_t_range(buf, t_range)
    if not (t0 <= t_ref <= t1):
        raise ValueError(
            f"t_ref={t_ref} must lie within t_range=[{t0}, {t1}]"
        )

    # Allocate a fresh ID for every seed.
    if new_ids is not None:
        new_ids = [int(x) for x in new_ids]
        if len(new_ids) != n_seeds:
            raise ValueError("new_ids length must match the number of seeds")
    else:
        next_id = _next_label_id(buf)
        new_ids = [next_id + i for i in range(n_seeds)]

    # ---- Reference frame watershed on background mask -------------------
    ref = buf.peek(t_ref).copy()
    bg_mask = ref == 0

    # When an image channel is provided, tighten the mask to only bright
    # (organoid-like) voxels using an automatic Otsu threshold, and use the
    # negated image as the watershed cost so that boundaries form at dim
    # regions between organoids rather than Voronoi-partitioning all background.
    if image_ref is not None:
        img = np.asarray(image_ref, dtype=np.float32)
        otsu_thresh = float(threshold_otsu(img))
        bg_mask = bg_mask & (img >= otsu_thresh)
        ws_image = img.max() - img
    else:
        ws_image = np.zeros_like(ref, dtype=np.float32)

    if not bg_mask.any():
        raise ValueError(
            f"No eligible background pixels at t_ref={t_ref}. "
            "If using a channel, try a different one or check that the seeds "
            "are placed on bright, unlabelled voxels."
        )

    seeds_img = np.zeros_like(ref)
    for nid, (zz, yy, xx) in zip(new_ids, seeds_zyx):
        zz = int(np.clip(zz, 0, ref.shape[0] - 1))
        yy = int(np.clip(yy, 0, ref.shape[1] - 1))
        xx = int(np.clip(xx, 0, ref.shape[2] - 1))
        if not bg_mask[zz, yy, xx]:
            raise ValueError(
                f"seed ({zz}, {yy}, {xx}) is not inside the eligible background "
                f"at t_ref={t_ref}; place seeds only on unlabelled voxels "
                "(and on bright voxels when a channel is selected)"
            )
        seeds_img[zz, yy, xx] = nid

    sub = watershed(ws_image, markers=seeds_img, mask=bg_mask)
    if not sub.any():
        raise RuntimeError("watershed produced an empty result; check seed positions")

    new_frames: Dict[int, np.ndarray] = {}
    out_ref = ref.copy()
    out_ref[bg_mask] = sub[bg_mask]
    new_frames[t_ref] = out_ref

    _total_frames = t1 - t0 + 1
    _done = [1]
    if progress_cb:
        progress_cb(1, _total_frames)

    # ---- Propagation forward and backward in time -----------------------
    def _propagate(direction: int, start: int, end: int):
        if direction == 0:
            return
        prev = out_ref.copy()
        step = direction
        t = start
        while (step > 0 and t <= end) or (step < 0 and t >= end):
            frame = buf.peek(t).copy()
            bg_mask_t = frame == 0
            if not bg_mask_t.any():
                # No background left — stop propagating.
                break

            # When a full image stack is available, apply the same
            # Otsu + negated-image strategy as the reference frame so that
            # watershed boundaries form at dim valleys between objects.
            # This prevents a single marker from flooding the entire
            # background in frames where no intensity guidance is used.
            if image_stack is not None:
                img_t = np.asarray(image_stack[t]).astype(np.float32)
                otsu_thresh = float(threshold_otsu(img_t))
                eligible = bg_mask_t & (img_t >= otsu_thresh)
                ws_image_t = img_t.max() - img_t
            else:
                eligible = bg_mask_t
                ws_image_t = np.zeros(frame.shape, dtype=np.float32)

            if not eligible.any():
                break

            # Keep only the new-label pixels from prev that fall inside the
            # eligible mask (they are the propagation markers).
            prev_seg = np.where(np.isin(prev, new_ids), prev, 0).astype(prev.dtype)
            prev_inside = np.where(eligible, prev_seg, 0)
            if not (prev_inside > 0).any():
                # New labels fully disappeared from the eligible region — stop.
                break

            # Intensity-guided watershed mirrors the reference-frame logic.
            new_sub = watershed(ws_image_t, markers=prev_inside, mask=eligible)
            out = frame.copy()
            out[eligible] = new_sub[eligible]
            new_frames[t] = out
            prev = out
            _done[0] += 1
            if progress_cb:
                progress_cb(_done[0], _total_frames)
            t += step

    _propagate(direction=+1, start=t_ref + 1, end=t1)
    _propagate(direction=-1, start=t_ref - 1, end=t0)

    return OpResult(
        name="create",
        new_frames=new_frames,
        affected_labels=list(new_ids),
        new_labels=list(new_ids),
        summary=(
            f"Created label(s) {new_ids} at t={t_ref}; "
            f"propagated across {len(new_frames)} frame(s)"
        ),
    )


def merge_labels(
    buf,
    label_ids: Sequence[int],
    target_id: Optional[int] = None,
    t_range: Optional[Tuple[int, int]] = None,
    n_workers: Optional[int] = None,
    progress_cb: Optional[Callable[[int, int], None]] = None,
) -> OpResult:
    """Merge several labels into one.

    All pixels equal to any of ``label_ids`` (across ``t_range``, default =
    union of their lifetimes) are rewritten to ``target_id``.  ``target_id``
    defaults to ``min(label_ids)`` so the lowest existing TrackID is kept.

    Frames are processed in parallel via :class:`ThreadPoolExecutor` when
    ``n_workers`` > 1 (numpy ``isin`` releases the GIL).
    """
    ids = [int(x) for x in label_ids]
    if len(ids) < 2:
        raise ValueError("merge_labels needs at least 2 labels")
    target_id = int(target_id) if target_id is not None else min(ids)
    if target_id not in ids:
        ids.append(target_id)

    if t_range is None:
        starts: List[int] = []
        ends: List[int] = []
        for lid in ids:
            f, l = lifetime_of(buf, lid)
            if f is not None:
                starts.append(f)
                ends.append(l)  # type: ignore[arg-type]
        if not starts:
            raise ValueError(f"none of the labels {ids} are present in the volume")
        t_range = (min(starts), max(ends))

    t0, t1 = _normalize_t_range(buf, t_range)
    others = [i for i in ids if i != target_id]

    # Prefetch sequentially — buf.peek() is not thread-safe.
    frames_to_process: List[Tuple[int, np.ndarray]] = []
    for t in range(t0, t1 + 1):
        frame = buf.peek(t)
        if np.isin(frame, others).any():
            frames_to_process.append((t, frame.copy()))

    def _merge_one(t_frame: Tuple[int, np.ndarray]) -> Tuple[int, Optional[np.ndarray]]:
        t, frame = t_frame
        mask = np.isin(frame, others)
        if not mask.any():
            return t, None
        out = frame.copy()
        out[mask] = target_id
        return t, out

    _nw = min(max(1, int(n_workers)) if n_workers else (os.cpu_count() or 1),
              max(1, len(frames_to_process)))
    _total = len(frames_to_process)
    new_frames: Dict[int, np.ndarray] = {}
    with ThreadPoolExecutor(max_workers=_nw) as ex:
        for _i, (t, out) in enumerate(ex.map(_merge_one, frames_to_process)):
            if progress_cb:
                progress_cb(_i + 1, _total)
            if out is not None:
                new_frames[t] = out

    return OpResult(
        name="merge",
        new_frames=new_frames,
        affected_labels=sorted({target_id, *others}),
        new_labels=[],
        summary=(
            f"Merged labels {others} into {target_id} across "
            f"{len(new_frames)} frame(s)"
        ),
    )


def erode_label(
    buf,
    label_id: int,
    r_xy: int = 1,
    r_z: int = 1,
    t_range: Optional[Tuple[int, int]] = None,
    n_workers: Optional[int] = None,
    progress_cb: Optional[Callable[[int, int], None]] = None,
) -> OpResult:
    """Erode a single label by a 3D ellipsoid footprint, per timepoint.

    Erosion is restricted to the label's mask: pixels of *other* labels are
    never touched.  Eroded-away pixels become background (``0``).

    Frames are processed in parallel using a :class:`ThreadPoolExecutor`
    (``skimage.morphology.binary_erosion`` releases the GIL, so all available
    CPU cores are utilised).
    """
    label_id = int(label_id)
    if t_range is None:
        f, l = lifetime_of(buf, label_id)
        if f is None:
            raise ValueError(f"label {label_id} not present anywhere in the volume")
        t_range = (f, l)
    t0, t1 = _normalize_t_range(buf, t_range)

    footprint = _ellipsoid_footprint(r_xy=r_xy, r_z=r_z)

    # Prefetch frames sequentially — buf.peek() is not thread-safe.
    candidate_frames: List[Tuple[int, np.ndarray]] = []
    for t in range(t0, t1 + 1):
        frame = buf.peek(t)
        if (frame == label_id).any():
            candidate_frames.append((t, frame.copy()))

    if not candidate_frames:
        return OpResult(
            name="erode",
            new_frames={},
            affected_labels=[label_id],
            new_labels=[],
            summary=f"Eroded label {label_id}: no frames contained the label",
        )

    def _erode_one(t_frame: Tuple[int, np.ndarray]):
        t, frame = t_frame
        mask = frame == label_id
        eroded = _binary_erosion_impl(mask, footprint)
        removed = mask & ~eroded
        if not removed.any():
            return t, None
        out = frame.copy()
        out[removed] = 0
        return t, out

    new_frames: Dict[int, np.ndarray] = {}
    _nw = min(max(1, int(n_workers)) if n_workers else (os.cpu_count() or 1),
              max(1, len(candidate_frames)))
    _total = len(candidate_frames)
    with ThreadPoolExecutor(max_workers=_nw) as ex:
        for _i, (t, out) in enumerate(ex.map(_erode_one, candidate_frames)):
            if progress_cb:
                progress_cb(_i + 1, _total)
            if out is not None:
                new_frames[t] = out

    return OpResult(
        name="erode",
        new_frames=new_frames,
        affected_labels=[label_id],
        new_labels=[],
        summary=(
            f"Eroded label {label_id} (XY={r_xy}px, Z={r_z}px) across "
            f"{len(new_frames)} frame(s)"
        ),
    )


def _dilate_one_frame(
    t_frame: Tuple[int, np.ndarray],
    label_id: int,
    r_xy: int,
    r_z: int,
) -> Tuple[int, Optional[np.ndarray]]:
    """Dilate ``label_id`` within a single pre-fetched frame.

    Returns ``(t, out_frame)`` if the frame changed, ``(t, None)`` otherwise.
    Kept at module level so it is picklable for :class:`ThreadPoolExecutor`.
    """
    t, frame = t_frame
    out = frame.copy()
    if r_xy > 0:
        for z in range(out.shape[0]):
            expanded = expand_labels(out[z], distance=r_xy)
            # expand_labels assigns every background pixel to its nearest label,
            # which naturally stops at other labels' borders.  We only accept
            # pixels that label_id claimed from background — leaving all other
            # labels untouched.
            new_pixels = (out[z] == 0) & (expanded == label_id)
            out[z][new_pixels] = label_id
    if r_z > 0:
        for _ in range(r_z):
            grown = out.copy()
            for z in range(out.shape[0]):
                bg = grown[z] == 0
                above = grown[z - 1] if z > 0 else None
                below = grown[z + 1] if z < out.shape[0] - 1 else None
                take = np.zeros_like(grown[z], dtype=bool)
                if above is not None:
                    take |= bg & (above == label_id)
                if below is not None:
                    take |= bg & (below == label_id)
                if take.any():
                    out[z][take] = label_id
    if np.array_equal(out, frame):
        return t, None
    return t, out


def dilate_label(
    buf,
    label_id: int,
    r_xy: int = 1,
    r_z: int = 1,
    t_range: Optional[Tuple[int, int]] = None,
    n_workers: Optional[int] = None,
    progress_cb: Optional[Callable[[int, int], None]] = None,
) -> OpResult:
    """Dilate / expand a single label by ``r_xy`` (px) in plane and ``r_z`` (slices).

    Uses :func:`skimage.segmentation.expand_labels` so the expansion never
    crosses into another label — neighbour cells are protected.
    Anisotropic expansion is implemented by running ``expand_labels`` once
    per Z slice for ``r_xy`` and then propagating ownership across
    ``+/- r_z`` slices for ``r_z``.

    Frames are processed in parallel using a :class:`ThreadPoolExecutor`
    (``skimage.segmentation.expand_labels`` releases the GIL, so all available
    CPU cores are utilised).
    """
    label_id = int(label_id)
    if t_range is None:
        f, l = lifetime_of(buf, label_id)
        if f is None:
            raise ValueError(f"label {label_id} not present anywhere in the volume")
        t_range = (f, l)
    t0, t1 = _normalize_t_range(buf, t_range)
    r_xy = max(0, int(r_xy))
    r_z = max(0, int(r_z))
    if r_xy == 0 and r_z == 0:
        return OpResult(name="dilate", summary="no-op (radii are 0)")

    # Prefetch frames sequentially — buf.peek() is not thread-safe.
    candidate_frames: List[Tuple[int, np.ndarray]] = []
    for t in range(t0, t1 + 1):
        frame = buf.peek(t)
        if (frame == label_id).any():
            candidate_frames.append((t, frame.copy()))

    if not candidate_frames:
        return OpResult(
            name="dilate",
            new_frames={},
            affected_labels=[label_id],
            new_labels=[],
            summary=f"Dilated label {label_id}: no frames contained the label",
        )

    _worker = partial(_dilate_one_frame, label_id=label_id, r_xy=r_xy, r_z=r_z)

    new_frames: Dict[int, np.ndarray] = {}
    _nw = min(max(1, int(n_workers)) if n_workers else (os.cpu_count() or 1),
              max(1, len(candidate_frames)))
    _total = len(candidate_frames)
    with ThreadPoolExecutor(max_workers=_nw) as ex:
        for _i, (t, out) in enumerate(ex.map(_worker, candidate_frames)):
            if progress_cb:
                progress_cb(_i + 1, _total)
            if out is not None:
                new_frames[t] = out

    return OpResult(
        name="dilate",
        new_frames=new_frames,
        affected_labels=[label_id],
        new_labels=[],
        summary=(
            f"Dilated label {label_id} (XY={r_xy}px, Z={r_z}px) across "
            f"{len(new_frames)} frame(s)"
        ),
    )


def delete_label(
    buf,
    label_id: int,
    t_range: Optional[Tuple[int, int]] = None,
    n_workers: Optional[int] = None,
    progress_cb: Optional[Callable[[int, int], None]] = None,
) -> OpResult:
    """Set every pixel of ``label_id`` to background (``0``) in ``t_range``.

    Frames are processed in parallel via :class:`ThreadPoolExecutor` when
    ``n_workers`` > 1.
    """
    label_id = int(label_id)
    if t_range is None:
        f, l = lifetime_of(buf, label_id)
        if f is None:
            raise ValueError(f"label {label_id} not present anywhere in the volume")
        t_range = (f, l)
    t0, t1 = _normalize_t_range(buf, t_range)

    # Prefetch sequentially — buf.peek() is not thread-safe.
    frames_to_process: List[Tuple[int, np.ndarray]] = []
    for t in range(t0, t1 + 1):
        frame = buf.peek(t)
        if (frame == label_id).any():
            frames_to_process.append((t, frame.copy()))

    def _del_one(t_frame: Tuple[int, np.ndarray]) -> Tuple[int, np.ndarray]:
        t, frame = t_frame
        out = frame.copy()
        out[frame == label_id] = 0
        return t, out

    _nw = min(max(1, int(n_workers)) if n_workers else (os.cpu_count() or 1),
              max(1, len(frames_to_process)))
    _total = len(frames_to_process)
    new_frames: Dict[int, np.ndarray] = {}
    with ThreadPoolExecutor(max_workers=_nw) as ex:
        for _i, (t, out) in enumerate(ex.map(_del_one, frames_to_process)):
            if progress_cb:
                progress_cb(_i + 1, _total)
            new_frames[t] = out

    return OpResult(
        name="delete",
        new_frames=new_frames,
        affected_labels=[label_id],
        new_labels=[],
        summary=(
            f"Deleted label {label_id} across {len(new_frames)} frame(s)"
        ),
    )


def delete_labels(
    buf,
    label_ids: Sequence[int],
    t_ranges: Optional[Dict[int, Tuple[int, int]]] = None,
    n_workers: Optional[int] = None,
    progress_cb: Optional[Callable[[int, int], None]] = None,
) -> OpResult:
    """Set every pixel of ``label_ids`` to background (``0``) in their respective ranges.

    Frames are processed in parallel via :class:`ThreadPoolExecutor` when
    ``n_workers`` > 1.
    """
    # Resolve per-label time ranges.
    resolved_ranges: Dict[int, Tuple[int, int]] = {}
    starts = []
    ends = []
    for lid in label_ids:
        rng = t_ranges.get(lid) if t_ranges else None
        if rng is None:
            f, l = lifetime_of(buf, lid)
            if f is None:
                continue
            rng = (f, l)
        rng = _normalize_t_range(buf, rng)
        resolved_ranges[lid] = rng
        starts.append(rng[0])
        ends.append(rng[1])

    if not resolved_ranges:
        return OpResult(
            name="delete",
            new_frames={},
            affected_labels=list(label_ids),
            new_labels=[],
            summary=f"Deleted labels {list(label_ids)}: none present in the volume"
        )

    t0_global, t1_global = min(starts), max(ends)

    # Prefetch sequentially — buf.peek() is not thread-safe.
    frames_to_process: List[Tuple[int, np.ndarray]] = []
    for t in range(t0_global, t1_global + 1):
        frame = buf.peek(t)
        active = [lid for lid, (lt0, lt1) in resolved_ranges.items() if lt0 <= t <= lt1]
        if any((frame == lid).any() for lid in active):
            frames_to_process.append((t, frame.copy()))

    def _del_multi(t_frame: Tuple[int, np.ndarray]) -> Tuple[int, Optional[np.ndarray]]:
        t, frame = t_frame
        out = None
        for lid, (lt0, lt1) in resolved_ranges.items():
            if not (lt0 <= t <= lt1):
                continue
            mask = (out if out is not None else frame) == lid
            if mask.any():
                if out is None:
                    out = frame.copy()
                out[mask] = 0
        return t, out

    _nw = min(max(1, int(n_workers)) if n_workers else (os.cpu_count() or 1),
              max(1, len(frames_to_process)))
    _total = len(frames_to_process)
    new_frames: Dict[int, np.ndarray] = {}
    with ThreadPoolExecutor(max_workers=_nw) as ex:
        for _i, (t, out) in enumerate(ex.map(_del_multi, frames_to_process)):
            if progress_cb:
                progress_cb(_i + 1, _total)
            if out is not None:
                new_frames[t] = out

    return OpResult(
        name="delete",
        new_frames=new_frames,
        affected_labels=list(label_ids),
        new_labels=[],
        summary=f"Deleted label(s) {list(label_ids)} across {len(new_frames)} frame(s)"
    )


def fill_label(
    buf,
    label_id: int,
    t_range: Optional[Tuple[int, int]] = None,
    n_workers: Optional[int] = None,
    progress_cb: Optional[Callable[[int, int], None]] = None,
) -> OpResult:
    """Fill internal holes inside a label's 3D mask at each timepoint.

    Uses :func:`scipy.ndimage.binary_fill_holes` to close cavities in the
    binary mask of ``label_id``.  Only voxels that are currently **background**
    (``== 0``) inside the filled region are assigned to ``label_id`` — other
    labels' pixels are never overwritten.

    This is useful for hollow sphere-like organoids that are segmented only on
    the outer shell: after filling the interior becomes solid.

    Frames are processed in parallel via :class:`ThreadPoolExecutor`.

    Parameters
    ----------
    buf:
        An :class:`~behav3d.editing.edit_buffer.EditBuffer` instance.
    label_id:
        The label whose holes are to be filled.
    t_range:
        Inclusive ``(t_first, t_last)`` range.  Defaults to the label's
        full lifetime.
    n_workers:
        Number of worker threads.  Defaults to ``os.cpu_count()``.
    progress_cb:
        Optional ``(current, total)`` progress callback.
    """
    try:
        from scipy.ndimage import binary_fill_holes as _fill_holes
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "fill_label requires scipy.  Install it with: pip install scipy"
        ) from exc

    label_id = int(label_id)
    if t_range is None:
        f, l = lifetime_of(buf, label_id)
        if f is None:
            raise ValueError(f"label {label_id} not present anywhere in the volume")
        t_range = (f, l)
    t0, t1 = _normalize_t_range(buf, t_range)

    # Prefetch sequentially — buf.peek() is not thread-safe.
    frames_to_process: List[Tuple[int, np.ndarray]] = []
    for t in range(t0, t1 + 1):
        frame = buf.peek(t)
        if (frame == label_id).any():
            frames_to_process.append((t, frame.copy()))

    if not frames_to_process:
        return OpResult(
            name="fill",
            new_frames={},
            affected_labels=[label_id],
            new_labels=[],
            summary=f"Fill label {label_id}: no frames contained the label",
        )

    def _fill_one(t_frame: Tuple[int, np.ndarray]) -> Tuple[int, Optional[np.ndarray]]:
        t, frame = t_frame
        mask = frame == label_id
        filled = _fill_holes(mask)
        # New pixels = filled but not in original mask AND background
        new_pixels = filled & ~mask & (frame == 0)
        if not new_pixels.any():
            return t, None
        out = frame.copy()
        out[new_pixels] = label_id
        return t, out

    _nw = min(max(1, int(n_workers)) if n_workers else (os.cpu_count() or 1),
              max(1, len(frames_to_process)))
    _total = len(frames_to_process)
    new_frames: Dict[int, np.ndarray] = {}
    with ThreadPoolExecutor(max_workers=_nw) as ex:
        for _i, (t, out) in enumerate(ex.map(_fill_one, frames_to_process)):
            if progress_cb:
                progress_cb(_i + 1, _total)
            if out is not None:
                new_frames[t] = out

    return OpResult(
        name="fill",
        new_frames=new_frames,
        affected_labels=[label_id],
        new_labels=[],
        summary=(
            f"Filled label {label_id}: plugged holes in {len(new_frames)} frame(s)"
        ),
    )
