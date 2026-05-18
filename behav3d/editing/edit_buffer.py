"""In-memory editing buffer for ``T x Z x Y x X`` tracked-segment volumes.

The buffer wraps a zarr-backed dask array (the on-disk
``*_tracked.zarr``) and:

* Loads frames lazily into a small per-frame numpy cache.
* Stages every primitive's output as a *commit* on a bounded undo/redo
  stack so the editor's Undo button can revert any layer-changing
  action before Save.
* On :meth:`save` writes only the dirty frames back to the original
  zarr (via :func:`behav3d.io.formats.zarr.write_zarr_parallel`) and
  regenerates the matching tracks CSV via
  :func:`behav3d.preprocessing.tracking.convert_tracked_image_to_csv`.

The on-disk file is touched **only** by :meth:`save`; until then every
edit lives in memory and can be discarded with :meth:`discard`.
"""
from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import dask.array as da
import numpy as np

from behav3d.io.formats.zarr import load_zarr, write_zarr_parallel
from behav3d.preprocessing.tracking import convert_tracked_image_to_csv
from behav3d.editing.tracked_segments import OpResult


# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------
@dataclass
class _Commit:
    """One reversible operation: pre/post snapshots of affected frames only."""
    name: str
    summary: str
    # Per-frame numpy snapshots — small because operations only touch one
    # label at a time; the buffer caps the maximum number of commits kept.
    pre: Dict[int, np.ndarray] = field(default_factory=dict)
    post: Dict[int, np.ndarray] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# EditBuffer
# ---------------------------------------------------------------------------
class EditBuffer:
    """Lazy, undoable view of a tracked-segment 4D zarr.

    Parameters
    ----------
    zarr_path:
        Path to the original ``*_tracked.zarr`` (writable).
    csv_path:
        Path to the matching ``*_tracks.csv`` (regenerated on save).
    pixel_size_xy, pixel_size_z:
        Used by :func:`convert_tracked_image_to_csv` when regenerating
        the CSV; usually pulled from the metadata row.
    max_history:
        Maximum number of commits kept on the undo stack.  Older
        commits are dropped (they cannot be undone, but they remain
        applied).
    cache_size:
        Maximum number of un-edited frames kept hot in memory.  Edited
        frames are always kept regardless.
    """

    def __init__(
        self,
        zarr_path: Path,
        csv_path: Optional[Path] = None,
        pixel_size_xy: float = 1.0,
        pixel_size_z: float = 1.0,
        max_history: int = 30,
        cache_size: int = 4,
        max_undoable_frames: int = 50,
    ) -> None:
        self.zarr_path = Path(zarr_path)
        self.csv_path = Path(csv_path) if csv_path is not None else None
        self.pixel_size_xy = float(pixel_size_xy)
        self.pixel_size_z = float(pixel_size_z)
        self.max_history = int(max_history)
        self.cache_size = int(cache_size)
        self.max_undoable_frames = int(max_undoable_frames)

        self._darr: da.Array = load_zarr(self.zarr_path)
        if self._darr.ndim != 4:
            raise ValueError(
                f"Tracked segments zarr must be 4D (T,Z,Y,X); got shape {self._darr.shape}"
            )
        self.shape: Tuple[int, ...] = tuple(int(x) for x in self._darr.shape)
        self.dtype = np.dtype(self._darr.dtype)

        # Mutated frames (kept in memory until save).
        self._dirty: Dict[int, np.ndarray] = {}
        # Hot cache of frames freshly read from disk.
        self._cache: "OrderedDict[int, np.ndarray]" = OrderedDict()
        # Undo/redo
        self._undo: List[_Commit] = []
        self._redo: List[_Commit] = []
        # Listeners that get notified whenever frames are mutated.
        self._listeners: List[Callable[[List[int]], None]] = []

        # --- Performance caches ------------------------------------------
        # Maximum non-zero label ID seen across the whole volume.  None means
        # "not yet computed"; populated lazily by max_label().  Maintained
        # incrementally by apply() so that subsequent splits never rescan the
        # full volume.  Reset to None by undo/redo (conservative) and by
        # save/discard (full reload).
        self._max_label: Optional[int] = None
        # Per-label lifetime cache: label_id → (t_first, t_last).  Populated
        # lazily by lifetime_of(); invalidated for affected labels whenever
        # frames change.
        self._lifetime_cache: Dict[int, Tuple[int, int]] = {}

    # ------------------------------------------------------------------
    # Listener / observer API
    # ------------------------------------------------------------------
    def add_frames_changed_listener(self, fn: Callable[[List[int]], None]) -> None:
        self._listeners.append(fn)

    def _emit(self, frames: List[int]) -> None:
        for fn in list(self._listeners):
            try:
                fn(list(frames))
            except Exception:
                pass

    # ------------------------------------------------------------------
    # Frame access
    # ------------------------------------------------------------------
    def peek(self, t: int) -> np.ndarray:
        """Return the current state of frame ``t`` (dirty or on-disk)."""
        t = int(t)
        if t in self._dirty:
            return self._dirty[t]
        if t in self._cache:
            self._cache.move_to_end(t)
            return self._cache[t]
        frame = np.asarray(self._darr[t])
        self._cache[t] = frame
        if len(self._cache) > self.cache_size:
            # Drop oldest unmutated frame.
            self._cache.popitem(last=False)
        return frame

    def unique_labels(self) -> np.ndarray:
        """Return the set of unique non-zero labels across the whole volume.

        Cached frames + dirty frames are read from memory; all other
        frames stream through dask.
        """
        seen: set = set()
        T = int(self.shape[0])
        for t in range(T):
            frame = self.peek(t)
            ids = np.unique(frame)
            seen.update(int(x) for x in ids if int(x) != 0)
        return np.array(sorted(seen), dtype=self.dtype)

    def max_label(self) -> int:
        """Return the maximum non-zero label ID in the volume.

        The result is cached after the first (full) scan and maintained
        incrementally by :meth:`apply` so that repeated calls — e.g. once
        per split — are O(1) instead of O(T × Z × Y × X).

        The cache is conservatively reset to ``None`` by :meth:`undo` and
        :meth:`redo` (those operations may lower the true maximum) and
        cleared entirely by :meth:`save` / :meth:`discard`.
        """
        if self._max_label is not None:
            return self._max_label
        # First call: full scan.
        T = int(self.shape[0])
        cur_max = 0
        for t in range(T):
            frame = self.peek(t)
            if frame.size:
                cur_max = max(cur_max, int(frame.max()))
        self._max_label = cur_max
        return self._max_label

    # ------------------------------------------------------------------
    # Lifetime cache helpers
    # ------------------------------------------------------------------
    def get_lifetime(self, label_id: int) -> Tuple[Optional[int], Optional[int]]:
        """Return ``(t_first, t_last)`` for ``label_id``, using the cache.

        On a cache miss the full O(T) scan is performed and the result is
        stored.  On subsequent calls (including from different primitives in
        the same editing session) the lookup is O(1).
        """
        label_id = int(label_id)
        if label_id in self._lifetime_cache:
            return self._lifetime_cache[label_id]
        T = int(self.shape[0])
        first: Optional[int] = None
        last: Optional[int] = None
        for t in range(T):
            frame = self.peek(t)
            if (frame == label_id).any():
                if first is None:
                    first = t
                last = t
        if first is not None:
            self._lifetime_cache[label_id] = (first, last)  # type: ignore[assignment]
        return first, last

    def invalidate_lifetime(self, label_id: int) -> None:
        """Remove ``label_id`` from the lifetime cache.

        Called by :meth:`apply` for every label touched by a commit so that
        the next :meth:`get_lifetime` call performs a fresh scan.
        """
        self._lifetime_cache.pop(int(label_id), None)

    # ------------------------------------------------------------------
    # Apply / undo / redo
    # ------------------------------------------------------------------
    def apply(self, op: OpResult) -> None:
        """Apply an :class:`OpResult` from a primitive as a new commit.

        Pre-snapshots are taken from the current state (which may be
        already-edited dirty frames); post-snapshots come from the op.
        Empty ops (no ``new_frames``) are ignored so the undo stack
        stays meaningful.

        When the operation modifies more than ``max_undoable_frames``
        timepoints, pre-snapshots are skipped to avoid an out-of-memory
        crash (N × frame_size bytes for pre + same for post would be
        unaffordable for long tracks).  The frames are still written to
        ``_dirty`` so the edit is applied, but the commit is **not**
        pushed onto the undo stack.  Callers that need to inform the user
        can check the returned bool: ``True`` = undoable, ``False`` = applied
        but not undoable.
        """
        if not op.new_frames:
            return
        large = len(op.new_frames) > self.max_undoable_frames
        commit = _Commit(name=op.name, summary=op.summary)
        for t, post_frame in op.new_frames.items():
            t = int(t)
            # Use copy=False: worker frames are freshly computed and never
            # mutated, so sharing the reference is safe and halves peak RAM.
            commit.post[t] = np.asarray(post_frame).astype(self.dtype, copy=False)
            self._dirty[t] = commit.post[t]
            if not large:
                # Store pre-snapshot only for operations that are small enough
                # to keep in memory for undo.
                commit.pre[t] = self.peek(t).copy()
        if not large:
            self._undo.append(commit)
            # New action invalidates redo history.
            self._redo.clear()
            # Bound history.
            while len(self._undo) > self.max_history:
                self._undo.pop(0)

        # --- Update performance caches -----------------------------------
        # Advance _max_label without rescanning if the new frames contain
        # a higher ID.
        if self._max_label is not None:
            for frame in commit.post.values():
                if frame.size:
                    self._max_label = max(self._max_label, int(frame.max()))
        # Invalidate lifetime cache for every label touched by this op so
        # the next lookup triggers a fresh (but lazy) scan.
        for label_id in op.affected_labels:
            self.invalidate_lifetime(label_id)

        self._emit(list(commit.post.keys()))

    def can_undo(self) -> bool:
        return bool(self._undo)

    def can_redo(self) -> bool:
        return bool(self._redo)

    def undo(self) -> Optional[str]:
        """Revert the latest commit; return its summary or ``None``."""
        if not self._undo:
            return None
        commit = self._undo.pop()
        for t, pre_frame in commit.pre.items():
            t = int(t)
            self._dirty[t] = pre_frame.copy()
            # If the pre-frame matches what's on disk, we can drop it
            # from the dirty set, but reading from disk to compare would
            # cost a full frame read; safer to leave it dirty until save
            # — write_zarr_parallel handles a no-op write fine.
        self._redo.append(commit)
        # Undo may lower the max label or restore lifetimes — reset
        # conservatively so the next access triggers a fresh scan.
        self._max_label = None
        # Invalidate lifetimes for every label that appears in the pre-frames
        # (the state we're reverting to) and the post-frames (what we're
        # undoing away from).
        for frames_dict in (commit.pre, commit.post):
            for frame in frames_dict.values():
                for lbl in np.unique(frame):
                    self.invalidate_lifetime(int(lbl))
        self._emit(list(commit.pre.keys()))
        return f"Undone: {commit.summary}"

    def redo(self) -> Optional[str]:
        if not self._redo:
            return None
        commit = self._redo.pop()
        for t, post_frame in commit.post.items():
            self._dirty[int(t)] = post_frame.copy()
        self._undo.append(commit)
        while len(self._undo) > self.max_history:
            self._undo.pop(0)
        # Re-advance max_label if possible; otherwise reset.
        if self._max_label is not None:
            for frame in commit.post.values():
                if frame.size:
                    self._max_label = max(self._max_label, int(frame.max()))
        for frame in commit.post.values():
            for lbl in np.unique(frame):
                self.invalidate_lifetime(int(lbl))
        self._emit(list(commit.post.keys()))
        return f"Redone: {commit.summary}"

    # ------------------------------------------------------------------
    # State queries
    # ------------------------------------------------------------------
    def is_dirty(self) -> bool:
        return bool(self._dirty)

    def dirty_count(self) -> int:
        return len(self._dirty)

    def history(self) -> List[Tuple[str, str]]:
        return [(c.name, c.summary) for c in self._undo]

    # ------------------------------------------------------------------
    # Save / discard
    # ------------------------------------------------------------------
    def save(self, regenerate_csv: bool = True) -> Tuple[int, Optional[Path]]:
        """Persist dirty frames to the original zarr and, optionally,
        rewrite the tracks CSV.

        Returns a tuple ``(n_frames_written, csv_path_or_None)``.
        """
        if not self._dirty:
            return 0, None
        n = 0
        for t in sorted(self._dirty):
            data = self._dirty[t]
            write_zarr_parallel(self.zarr_path, index=t, data=data)
            n += 1
        # The dirty data is now on disk: clear in-memory mutations and
        # refresh the underlying dask view so subsequent peek()s read
        # the new bytes.
        self._dirty.clear()
        self._cache.clear()
        self._undo.clear()
        self._redo.clear()
        self._max_label = None
        self._lifetime_cache.clear()
        self._darr = load_zarr(self.zarr_path)

        out_csv: Optional[Path] = None
        if regenerate_csv and self.csv_path is not None:
            self.csv_path.parent.mkdir(parents=True, exist_ok=True)
            convert_tracked_image_to_csv(
                img_path=self.zarr_path,
                outpath=self.csv_path,
                element_size_x=self.pixel_size_xy,
                element_size_y=self.pixel_size_xy,
                element_size_z=self.pixel_size_z,
            )
            out_csv = self.csv_path
        return n, out_csv

    def discard(self) -> None:
        """Drop every staged edit; the next :meth:`peek` re-reads disk."""
        self._dirty.clear()
        self._cache.clear()
        self._undo.clear()
        self._redo.clear()
        self._max_label = None
        self._lifetime_cache.clear()
        self._darr = load_zarr(self.zarr_path)
        self._emit(list(range(int(self.shape[0]))))
