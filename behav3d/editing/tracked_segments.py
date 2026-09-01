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
``delete_label``, ``delete_labels``, ``fill_holes``) accept an optional
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
from scipy import ndimage
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


def _seed_components_by_nearest(
    mask_t: np.ndarray,
    prev: np.ndarray,
    seed_ids: Sequence[int],
) -> np.ndarray:
    """Build propagation markers so every mask blob keeps a sub-label.

    The default seeding for split propagation keeps a previous sub-label only
    where its pixels still overlap the current frame's mask.  For small,
    spatially-separate cells that move between frames, a blob can end up with
    *no* overlapping marker; a plain watershed then leaves that whole
    (disconnected) blob as background, silently erasing the label.

    This helper starts from the same overlap markers but additionally rescues
    any sub-label that lost *all* overlap this frame, so it is not erased.  The
    rescue is a strict **one-to-one** match: each lost sub-label is matched to a
    distinct orphan blob (one that received no overlap marker) by nearest
    centroid distance.  Doing it one-to-one is essential — matching each blob to
    its nearest label independently can assign *both* blobs to the same label
    and drop the other, which would delete one of the split IDs across time.

    Blobs already carrying an overlap marker are left untouched, so the
    watershed still splits genuinely-connected regions along a ridge exactly as
    before.  Any leftover orphan blob (a stray same-label fragment, or more
    blobs than labels) is attached to its nearest sub-label so its pixels are
    not dropped.  Only existing mask pixels are ever labelled — a cell that
    truly left the frame has no blob and so is still (correctly) not propagated.
    """
    from scipy.ndimage import label as _cc_label, center_of_mass

    prev_seg = np.where(np.isin(prev, seed_ids), prev, 0).astype(prev.dtype)
    markers = np.where(mask_t, prev_seg, 0).astype(prev.dtype)

    # Previous-frame centroid of each sub-label, computed on the whole prev
    # frame (not just the overlap) so a moved cell still has a reference point.
    centroids: Dict[int, np.ndarray] = {}
    for sid in seed_ids:
        m = prev == sid
        if m.any():
            centroids[int(sid)] = np.array(center_of_mass(m))
    if not centroids:
        return markers

    cc, n_cc = _cc_label(mask_t)
    if n_cc == 0:
        return markers

    comp_centroids: Dict[int, np.ndarray] = {}
    comp_has_marker: Dict[int, bool] = {}
    for comp in range(1, n_cc + 1):
        comp_mask = cc == comp
        comp_has_marker[comp] = bool((markers[comp_mask] > 0).any())
        comp_centroids[comp] = np.array(center_of_mass(comp_mask))

    def _seed_component(comp: int, sid: int) -> None:
        comp_mask = cc == comp
        cz, cy, cx = (int(round(v)) for v in comp_centroids[comp])
        if not comp_mask[cz, cy, cx]:
            # Centroid can fall outside a concave/annular blob; use any voxel.
            zs, ys, xs = np.nonzero(comp_mask)
            cz, cy, cx = int(zs[0]), int(ys[0]), int(xs[0])
        markers[cz, cy, cx] = sid

    # Sub-labels already represented by an overlap marker somewhere this frame.
    present = {int(v) for v in np.unique(markers) if v != 0}
    # Sub-labels that lost all overlap this frame — at risk of being erased.
    missing = [sid for sid in centroids if sid not in present]
    # Components with no overlap marker — available to host a rescued label.
    orphans = [c for c in range(1, n_cc + 1) if not comp_has_marker[c]]

    # 1) Rescue each missing sub-label into a DISTINCT orphan blob, matched
    #    one-to-one by ascending centroid distance (greedy global-nearest), so
    #    two lost labels can never collapse onto the same blob.
    pairs = [
        (float(np.linalg.norm(centroids[sid] - comp_centroids[comp])), sid, comp)
        for sid in missing
        for comp in orphans
    ]
    pairs.sort(key=lambda p: p[0])
    used_sid: set = set()
    used_comp: set = set()
    for _d, sid, comp in pairs:
        if sid in used_sid or comp in used_comp:
            continue
        _seed_component(comp, sid)
        used_sid.add(sid)
        used_comp.add(comp)

    # 2) Any orphan blob still unseeded (stray same-label fragment, or more
    #    blobs than labels) → attach to its nearest sub-label so it is kept.
    for comp in orphans:
        if comp in used_comp:
            continue
        nearest = min(
            centroids,
            key=lambda sid: float(np.linalg.norm(centroids[sid] - comp_centroids[comp])),
        )
        _seed_component(comp, nearest)

    return markers


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
            # Seed the current mask from the previous frame's sub-labels.
            # Overlapping pixels seed directly; any blob that would otherwise
            # be left marker-less (e.g. a small, separate cell that moved) is
            # matched to the nearest previous sub-label so it is not erased.
            prev_inside = _seed_components_by_nearest(mask_t, prev, seed_ids)
            # No sub-label survives anywhere in this frame's mask — stop.
            if not (prev_inside > 0).any():
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

    def _merge_one(t_frame: Tuple[int, np.ndarray]) -> Tuple[int, Optional[np.ndarray]]:
        t, frame = t_frame
        mask = np.isin(frame, others)
        if not mask.any():
            return t, None
        out = frame.copy()
        out[mask] = target_id
        return t, out

    _nw = min(max(1, int(n_workers)) if n_workers else (os.cpu_count() or 1), max(1, t1 - t0 + 1))
    batch_size = _nw * 2
    _total = t1 - t0 + 1
    _processed = 0
    new_frames: Dict[int, np.ndarray] = {}

    with ThreadPoolExecutor(max_workers=_nw) as ex:
        batch = []
        for t in range(t0, t1 + 1):
            frame = buf.peek(t)
            if np.isin(frame, others).any():
                batch.append((t, frame.copy()))
            else:
                _processed += 1
                if progress_cb:
                    progress_cb(_processed, _total)
            
            if len(batch) >= batch_size or t == t1:
                if batch:
                    for _i, (t_out, out) in enumerate(ex.map(_merge_one, batch)):
                        _processed += 1
                        if progress_cb:
                            progress_cb(_processed, _total)
                        if out is not None:
                            new_frames[t_out] = out
                    batch = []

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


def swap_labels(
    buf,
    label_ids: Sequence[int],
    t_range: Optional[Tuple[int, int]] = None,
    n_workers: Optional[int] = None,
    progress_cb: Optional[Callable[[int, int], None]] = None,
) -> OpResult:
    """Swap two label IDs across the chosen time range.

    Fixes identity swaps introduced by the tracker: when two cells lie close
    together the tracker can exchange their TrackIDs from some frame onward.
    Selecting the two labels and the range where they are wrong re-exchanges
    them so each track recovers its correct identity.

    Frames are rewritten when at least one of the two labels exists in that
    timepoint.  If both are present their IDs are swapped; if only one is
    present (e.g. the other has a tracking gap) that label is rewritten to the
    complementary ID so the relabel stays consistent across the gap.
    Timepoints where neither label exists are left unchanged.

    .. note::
        Swapping over the **full union** of both lifetimes is a pure rename
        (``A`` becomes ``B`` and vice-versa everywhere) and does not fix a
        partial swap.  To repair a real ID swap, pass the ``t_range`` that
        starts at the frame where the swap began.

    Frames are processed in parallel via :class:`ThreadPoolExecutor` (numpy
    ``==`` releases the GIL), mirroring :func:`merge_labels`.
    """
    ids = [int(x) for x in label_ids]
    if len(ids) != 2:
        raise ValueError("swap_labels needs exactly 2 labels")
    if ids[0] == ids[1]:
        raise ValueError("swap_labels needs two distinct labels")

    a, b = ids
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

    def _swap_one(t_frame: Tuple[int, np.ndarray]) -> Tuple[int, Optional[np.ndarray]]:
        t, frame = t_frame
        mask_a = frame == a
        mask_b = frame == b
        has_a = mask_a.any()
        has_b = mask_b.any()
        if not has_a and not has_b:
            return t, None
        out = frame.copy()
        if has_a and has_b:
            out[mask_a] = b
            out[mask_b] = a
        elif has_a:
            out[mask_a] = b
        else:
            out[mask_b] = a
        return t, out

    _nw = min(max(1, int(n_workers)) if n_workers else (os.cpu_count() or 1), max(1, t1 - t0 + 1))
    batch_size = _nw * 2
    _total = t1 - t0 + 1
    _processed = 0
    new_frames: Dict[int, np.ndarray] = {}

    with ThreadPoolExecutor(max_workers=_nw) as ex:
        batch = []
        for t in range(t0, t1 + 1):
            frame = buf.peek(t)
            if (frame == a).any() or (frame == b).any():
                batch.append((t, frame.copy()))
            else:
                _processed += 1
                if progress_cb:
                    progress_cb(_processed, _total)

            if len(batch) >= batch_size or t == t1:
                if batch:
                    for _i, (t_out, out) in enumerate(ex.map(_swap_one, batch)):
                        _processed += 1
                        if progress_cb:
                            progress_cb(_processed, _total)
                        if out is not None:
                            new_frames[t_out] = out
                    batch = []

    return OpResult(
        name="swap",
        new_frames=new_frames,
        affected_labels=sorted({a, b}),
        new_labels=[],
        summary=(
            f"Swapped labels {a} and {b} across {len(new_frames)} frame(s)"
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

    _nw = min(max(1, int(n_workers)) if n_workers else (os.cpu_count() or 1), max(1, t1 - t0 + 1))
    batch_size = _nw * 2
    _total = t1 - t0 + 1
    _processed = 0
    _found_any = False
    new_frames: Dict[int, np.ndarray] = {}

    with ThreadPoolExecutor(max_workers=_nw) as ex:
        batch = []
        for t in range(t0, t1 + 1):
            frame = buf.peek(t)
            if (frame == label_id).any():
                _found_any = True
                batch.append((t, frame.copy()))
            else:
                _processed += 1
                if progress_cb:
                    progress_cb(_processed, _total)
            
            if len(batch) >= batch_size or t == t1:
                if batch:
                    for _i, (t_out, out) in enumerate(ex.map(_erode_one, batch)):
                        _processed += 1
                        if progress_cb:
                            progress_cb(_processed, _total)
                        if out is not None:
                            new_frames[t_out] = out
                    batch = []

    if not _found_any:
        return OpResult(
            name="erode",
            new_frames={},
            affected_labels=[label_id],
            new_labels=[],
            summary=f"Eroded label {label_id}: no frames contained the label",
        )

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

    _worker = partial(_dilate_one_frame, label_id=label_id, r_xy=r_xy, r_z=r_z)

    _nw = min(max(1, int(n_workers)) if n_workers else (os.cpu_count() or 1), max(1, t1 - t0 + 1))
    batch_size = _nw * 2
    _total = t1 - t0 + 1
    _processed = 0
    _found_any = False
    new_frames: Dict[int, np.ndarray] = {}

    with ThreadPoolExecutor(max_workers=_nw) as ex:
        batch = []
        for t in range(t0, t1 + 1):
            frame = buf.peek(t)
            if (frame == label_id).any():
                _found_any = True
                batch.append((t, frame.copy()))
            else:
                _processed += 1
                if progress_cb:
                    progress_cb(_processed, _total)
            
            if len(batch) >= batch_size or t == t1:
                if batch:
                    for _i, (t_out, out) in enumerate(ex.map(_worker, batch)):
                        _processed += 1
                        if progress_cb:
                            progress_cb(_processed, _total)
                        if out is not None:
                            new_frames[t_out] = out
                    batch = []

    if not _found_any:
        return OpResult(
            name="dilate",
            new_frames={},
            affected_labels=[label_id],
            new_labels=[],
            summary=f"Dilated label {label_id}: no frames contained the label",
        )

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

    def _del_one(t_frame: Tuple[int, np.ndarray]) -> Tuple[int, np.ndarray]:
        t, frame = t_frame
        out = frame.copy()
        out[frame == label_id] = 0
        return t, out

    _nw = min(max(1, int(n_workers)) if n_workers else (os.cpu_count() or 1), max(1, t1 - t0 + 1))
    batch_size = _nw * 2
    _total = t1 - t0 + 1
    _processed = 0
    new_frames: Dict[int, np.ndarray] = {}

    with ThreadPoolExecutor(max_workers=_nw) as ex:
        batch = []
        for t in range(t0, t1 + 1):
            frame = buf.peek(t)
            if (frame == label_id).any():
                batch.append((t, frame.copy()))
            else:
                _processed += 1
                if progress_cb:
                    progress_cb(_processed, _total)
            
            if len(batch) >= batch_size or t == t1:
                if batch:
                    for _i, (t_out, out) in enumerate(ex.map(_del_one, batch)):
                        _processed += 1
                        if progress_cb:
                            progress_cb(_processed, _total)
                        new_frames[t_out] = out
                    batch = []

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

    _nw = min(max(1, int(n_workers)) if n_workers else (os.cpu_count() or 1), max(1, t1_global - t0_global + 1))
    batch_size = _nw * 2
    _total = t1_global - t0_global + 1
    _processed = 0
    new_frames: Dict[int, np.ndarray] = {}

    with ThreadPoolExecutor(max_workers=_nw) as ex:
        batch = []
        for t in range(t0_global, t1_global + 1):
            frame = buf.peek(t)
            active = [lid for lid, (lt0, lt1) in resolved_ranges.items() if lt0 <= t <= lt1]
            if any((frame == lid).any() for lid in active):
                batch.append((t, frame.copy()))
            else:
                _processed += 1
                if progress_cb:
                    progress_cb(_processed, _total)
            
            if len(batch) >= batch_size or t == t1_global:
                if batch:
                    for _i, (t_out, out) in enumerate(ex.map(_del_multi, batch)):
                        _processed += 1
                        if progress_cb:
                            progress_cb(_processed, _total)
                        if out is not None:
                            new_frames[t_out] = out
                    batch = []

    return OpResult(
        name="delete",
        new_frames=new_frames,
        affected_labels=list(label_ids),
        new_labels=[],
        summary=f"Deleted label(s) {list(label_ids)} across {len(new_frames)} frame(s)"
    )


def fill_holes(
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
            "fill_holes requires scipy.  Install it with: pip install scipy"
        ) from exc

    label_id = int(label_id)
    if t_range is None:
        f, l = lifetime_of(buf, label_id)
        if f is None or l is None:
            raise ValueError(f"label {label_id} not present anywhere in the volume")
        t_range = (f, l)
    t0, t1 = _normalize_t_range(buf, t_range)

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

    _nw = min(max(1, int(n_workers)) if n_workers else (os.cpu_count() or 1), max(1, t1 - t0 + 1))
    batch_size = _nw * 2
    _total = t1 - t0 + 1
    _processed = 0
    _found_any = False
    new_frames: Dict[int, np.ndarray] = {}

    with ThreadPoolExecutor(max_workers=_nw) as ex:
        batch = []
        for t in range(t0, t1 + 1):
            frame = buf.peek(t)
            if (frame == label_id).any():
                _found_any = True
                batch.append((t, frame.copy()))
            else:
                _processed += 1
                if progress_cb:
                    progress_cb(_processed, _total)
            
            if len(batch) >= batch_size or t == t1:
                if batch:
                    for _i, (t_out, out) in enumerate(ex.map(_fill_one, batch)):
                        _processed += 1
                        if progress_cb:
                            progress_cb(_processed, _total)
                        if out is not None:
                            new_frames[t_out] = out
                    batch = []

    if not _found_any:
        return OpResult(
            name="fill_holes",
            new_frames={},
            affected_labels=[label_id],
            new_labels=[],
            summary=f"Fill holes in label {label_id}: no frames contained the label",
        )

    return OpResult(
        name="fill_holes",
        new_frames=new_frames,
        affected_labels=[label_id],
        new_labels=[],
        summary=(
            f"Filled label {label_id}: plugged holes in {len(new_frames)} frame(s)"
        ),
    )


# ---------------------------------------------------------------------------
# Fill temporal gaps (interpolation across missing timepoints)
# ---------------------------------------------------------------------------
def _bbox_slices(mask: np.ndarray):
    """Return a tuple of slices for the tight bounding box of ``mask``.

    ``mask`` must contain at least one True voxel.
    """
    coords = np.array(np.nonzero(mask))
    mins = coords.min(axis=1)
    maxs = coords.max(axis=1) + 1
    return mins, maxs


def _interpolate_mask(
    mask_before: np.ndarray,
    mask_after: np.ndarray,
    num_steps: int,
) -> List[np.ndarray]:
    """Interpolate between two binary masks via signed-distance blending.

    Both masks are aligned to a shared, linearly-interpolated centroid and then
    blended in signed-distance space, so the result is a real shape morph rather
    than a hard switch between the two endpoints.

    Performance: all geometric work (centroid shift, distance transform, blend)
    is restricted to a padded bounding box around the union of the two masks.
    Cells occupy a tiny fraction of a ``Z x Y x X`` volume, so this is typically
    one to two orders of magnitude faster than operating on the full frame, and
    the ``> 0`` threshold near the object is unchanged as long as the padding
    exceeds the centroid displacement plus the number of interpolation steps.

    Returns a list of ``num_steps`` float32 masks (values 0.0 / 1.0), each the
    same shape as the inputs.
    """
    if num_steps <= 0:
        return []

    full_shape = mask_before.shape
    union = mask_before | mask_after
    if not union.any():
        return [np.zeros(full_shape, dtype=np.float32) for _ in range(num_steps)]

    # Exact-equality fast path — morphing a shape into itself is a no-op.
    if np.array_equal(mask_before, mask_after):
        return [mask_before.astype(np.float32) for _ in range(num_steps)]

    com_before = (
        np.array(ndimage.center_of_mass(mask_before)) if mask_before.any() else None
    )
    com_after = (
        np.array(ndimage.center_of_mass(mask_after)) if mask_after.any() else None
    )

    # ---- Restrict to a padded bounding box around the union --------------
    mins, maxs = _bbox_slices(union)
    if com_before is not None and com_after is not None:
        shift_pad = int(np.ceil(np.abs(com_after - com_before).max())) + 2
    else:
        shift_pad = 2
    pad = shift_pad + max(1, int(num_steps))
    lo = np.maximum(mins - pad, 0)
    hi = np.minimum(maxs + pad, full_shape)
    sl = tuple(slice(int(a), int(b)) for a, b in zip(lo, hi))

    mb = mask_before[sl]
    ma = mask_after[sl]

    def signed_distance(mask: np.ndarray) -> np.ndarray:
        inside = ndimage.distance_transform_edt(mask)
        outside = ndimage.distance_transform_edt(~mask)
        return inside - outside

    n = num_steps + 1
    interpolated: List[np.ndarray] = []
    for i in range(1, n):
        alpha = i / n
        if com_before is not None and com_after is not None:
            target_com = (1 - alpha) * com_before + alpha * com_after
            shifted_before = (
                ndimage.shift(mb.astype(np.float32), target_com - com_before, order=0)
                >= 0.5
            )
            shifted_after = (
                ndimage.shift(ma.astype(np.float32), target_com - com_after, order=0)
                >= 0.5
            )
        else:
            shifted_before, shifted_after = mb, ma

        sdf_interp = (1 - alpha) * signed_distance(shifted_before) + alpha * signed_distance(
            shifted_after
        )
        out = np.zeros(full_shape, dtype=np.float32)
        out[sl] = (sdf_interp > 0).astype(np.float32)
        interpolated.append(out)

    return interpolated


def fill_gaps(
    buf,
    label_id: int,
    t_range: Optional[Tuple[int, int]] = None,
    max_gap_size: int = 5,
    n_workers: Optional[int] = None,
    progress_cb: Optional[Callable[[int, int], None]] = None,
) -> OpResult:
    """Fill temporal gaps in a label's lifetime via mask interpolation.

    If a label disappears for one or more consecutive frames and then
    reappears, the missing frames are reconstructed by interpolating the binary
    mask between the last frame before the gap and the first frame after it (see
    :func:`_interpolate_mask`).

    Only background voxels (``== 0``) are ever written — other labels are never
    overwritten.  A gap whose interpolated mask would overlap another label is
    reported as *conflicted* and left untouched (it usually needs a Split).

    Parameters
    ----------
    buf:
        An :class:`~behav3d.editing.edit_buffer.EditBuffer` instance.
    label_id:
        The label whose gaps are to be filled.
    t_range:
        Inclusive ``(t_first, t_last)`` range.  Defaults to the label's full
        lifetime.
    max_gap_size:
        Maximum number of missing timepoints to interpolate.  A gap with
        ``<= max_gap_size`` missing frames is filled; a larger gap is skipped
        and listed in the summary (a long gap usually indicates a tracking /
        merge error, not a true short dropout).
    n_workers:
        Number of worker threads used to process gaps in parallel.  Defaults to
        ``os.cpu_count()``.
    progress_cb:
        Optional ``(current, total)`` progress callback.
    """
    label_id = int(label_id)
    max_gap_size = max(1, int(max_gap_size))

    if t_range is None:
        f, l = lifetime_of(buf, label_id)
        if f is None or l is None:
            raise ValueError(f"label {label_id} not present anywhere in the volume")
        t_range = (f, l)
    t0, t1 = _normalize_t_range(buf, t_range)

    # ---- Pass 1: presence scan (each frame read once, not retained) -----
    _scan_total = t1 - t0 + 1
    present: Dict[int, bool] = {}
    for _i, t in enumerate(range(t0, t1 + 1), start=1):
        present[t] = bool((buf.peek(t) == label_id).any())
        if progress_cb:
            progress_cb(_i, _scan_total + 1)

    # ---- Continuous segments where the label exists --------------------
    segments: List[Tuple[int, int]] = []
    seg_start: Optional[int] = None
    for t in range(t0, t1 + 1):
        if present[t] and seg_start is None:
            seg_start = t
        elif not present[t] and seg_start is not None:
            segments.append((seg_start, t - 1))
            seg_start = None
    if seg_start is not None:
        segments.append((seg_start, t1))

    if len(segments) < 2:
        if progress_cb:
            progress_cb(_scan_total + 1, _scan_total + 1)
        return OpResult(
            name="fill_gaps",
            new_frames={},
            affected_labels=[label_id],
            new_labels=[],
            summary=f"Fill gaps in label {label_id}: no gaps found",
        )

    # ---- Classify gaps by size ----------------------------------------
    gaps: List[Tuple[int, int]] = []
    big_gaps: List[Tuple[int, int]] = []
    for i in range(len(segments) - 1):
        gap_start = segments[i][1] + 1
        gap_end = segments[i + 1][0] - 1
        if gap_start <= gap_end:
            if (gap_end - gap_start + 1) > max_gap_size:
                big_gaps.append((gap_start, gap_end))
            else:
                gaps.append((gap_start, gap_end))

    # ---- Read only the frames referenced by fillable gaps -------------
    needed: set = set()
    for gap_start, gap_end in gaps:
        needed.update(range(gap_start - 1, gap_end + 2))
    frames: Dict[int, np.ndarray] = {t: buf.peek(t) for t in sorted(needed)}

    conflicted_gaps: List[Tuple[int, int]] = []

    def _process_gap(gap: Tuple[int, int]):
        gap_start, gap_end = gap
        mask_before = frames[gap_start - 1] == label_id
        mask_after = frames[gap_end + 1] == label_id
        num = gap_end - gap_start + 1
        interpolated = _interpolate_mask(mask_before, mask_after, num)
        out_frames: Dict[int, np.ndarray] = {}
        for idx, t in enumerate(range(gap_start, gap_end + 1)):
            binary = interpolated[idx] >= 0.5
            frame_t = frames[t]
            if ((frame_t != 0) & (frame_t != label_id) & binary).any():
                return gap, "conflict", {}
            eligible = (frame_t == 0) & binary
            if eligible.any():
                out = frame_t.copy()
                out[eligible] = label_id
                out_frames[t] = out
        return gap, "ok", out_frames

    new_frames: Dict[int, np.ndarray] = {}
    if gaps:
        _nw = min(
            max(1, int(n_workers)) if n_workers else (os.cpu_count() or 1),
            len(gaps),
        )
        _done = 0
        with ThreadPoolExecutor(max_workers=_nw) as ex:
            for gap, status, out_frames in ex.map(_process_gap, gaps):
                _done += 1
                if progress_cb:
                    progress_cb(
                        _scan_total + int(_done / max(1, len(gaps))),
                        _scan_total + 1,
                    )
                if status == "conflict":
                    conflicted_gaps.append(gap)
                else:
                    new_frames.update(out_frames)
    if progress_cb:
        progress_cb(_scan_total + 1, _scan_total + 1)

    filled_gaps = [
        g
        for g in gaps
        if g not in conflicted_gaps
        and any(t in new_frames for t in range(g[0], g[1] + 1))
    ]

    if not new_frames:
        summary = f"Fill gaps in label {label_id}: "
        if conflicted_gaps:
            summary += (
                f"all gaps conflict with other labels (require split): "
                f"{conflicted_gaps}"
            )
        else:
            summary += "no gaps to fill"
        if big_gaps:
            summary += (
                f"; gaps larger than {max_gap_size} timepoint(s) skipped: {big_gaps}"
            )
        return OpResult(
            name="fill_gaps",
            new_frames={},
            affected_labels=[label_id],
            new_labels=[],
            summary=summary,
        )

    summary = (
        f"Filled {len(filled_gaps)} gap(s) in label {label_id} "
        f"across {len(new_frames)} frame(s) via interpolation"
    )
    if big_gaps:
        summary += f"; gaps larger than {max_gap_size} timepoint(s) skipped: {big_gaps}"
    if conflicted_gaps:
        summary += f"; conflicts with other labels (require split): {conflicted_gaps}"

    return OpResult(
        name="fill_gaps",
        new_frames=new_frames,
        affected_labels=[label_id],
        new_labels=[],
        summary=summary,
    )
