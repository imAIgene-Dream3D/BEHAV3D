"""In-memory editing buffer for ``T x Z x Y x X`` tracked-segment volumes.

The buffer wraps a zarr-backed dask array (the on-disk
``*_tracked.zarr``) and:

* Loads frames lazily into a small per-frame numpy cache.
* Stages every primitive's output as a dirty frame in memory.
* On :meth:`save` writes only the dirty frames back to the original
  zarr (via :func:`behav3d.io.formats.zarr.write_zarr_parallel`) and
  regenerates the matching tracks CSV via
  :func:`behav3d.preprocessing.tracking.convert_tracked_image_to_csv`.

.. note::
    Undo/redo has been intentionally removed to minimise RAM usage.
    The :meth:`discard` method reverts **all** unsaved edits back to
    the on-disk state, which is the recommended recovery path.

The on-disk file is touched **only** by :meth:`save`; until then every
edit lives in memory and can be discarded with :meth:`discard`.
"""
from __future__ import annotations

from collections import OrderedDict
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import dask.array as da
import numpy as np

from behav3d.io.formats.zarr import load_zarr, write_zarr_parallel
from behav3d.preprocessing.tracking import convert_tracked_image_to_csv
from behav3d.editing.tracked_segments import OpResult


# ---------------------------------------------------------------------------
# EditBuffer
# ---------------------------------------------------------------------------
class EditBuffer:
    """Lazy, committable view of a tracked-segment 4D zarr.

    Parameters
    ----------
    zarr_path:
        Path to the original ``*_tracked.zarr`` (writable).
    csv_path:
        Path to the matching ``*_tracks.csv`` (regenerated on save).
    pixel_size_xy, pixel_size_z:
        Used by :func:`convert_tracked_image_to_csv` when regenerating
        the CSV; usually pulled from the metadata row.
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
        cache_size: int = 4,
        # Kept for backward-compat — ignored (undo/redo removed).
        max_history: int = 30,
        max_undoable_frames: int = 50,
    ) -> None:
        self.zarr_path = Path(zarr_path)
        self.csv_path = Path(csv_path) if csv_path is not None else None
        self.pixel_size_xy = float(pixel_size_xy)
        self.pixel_size_z = float(pixel_size_z)
        self.cache_size = int(cache_size)

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
        # Listeners that get notified whenever frames are mutated.
        self._listeners: List[Callable[[List[int]], None]] = []

        # --- Performance caches ------------------------------------------
        # Maximum non-zero label ID seen across the whole volume.  None means
        # "not yet computed"; populated lazily by max_label().  Maintained
        # incrementally by commit() so that subsequent splits never rescan the
        # full volume.  Reset to None by save/discard (full reload).
        self._max_label: Optional[int] = None
        # Per-label lifetime cache: label_id → (t_first, t_last).  Populated
        # lazily by get_lifetime() or eagerly by build_lifetime_index();
        # invalidated for affected labels whenever frames change.
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
        incrementally by :meth:`commit` so that repeated calls — e.g. once
        per split — are O(1) instead of O(T × Z × Y × X).

        The cache is cleared by :meth:`save` / :meth:`discard`.
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

        Call :meth:`build_lifetime_index` on editor start to pre-populate
        the cache for all labels in a single pass, avoiding O(T) per-label
        scans on first use.
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

    def build_lifetime_index(
        self,
        progress_cb: Optional[Callable[[int, int], None]] = None,
    ) -> None:
        """Scan the entire volume once and populate the lifetime cache.

        After this call, :meth:`get_lifetime` is O(1) for every label
        present in the volume.  Also pre-computes :meth:`max_label`.

        This is called in a background thread when the editor opens so
        the user does not block on first label click.

        Parameters
        ----------
        progress_cb:
            Optional ``(current_frame, total_frames)`` callback invoked
            after each timepoint is processed.
        is_interrupted:
            Optional callback returning True if the scan should abort early.
        """
        T = int(self.shape[0])
        # first_seen[label] → t_first, last_seen[label] → t_last
        first_seen: Dict[int, int] = {}
        last_seen: Dict[int, int] = {}
        cur_max = 0

        for t in range(T):
            if is_interrupted is not None and is_interrupted():
                return
            frame = self.peek(t)
            if frame.size:
                cur_max = max(cur_max, int(frame.max()))
            for lbl in np.unique(frame):
                lbl = int(lbl)
                if lbl == 0:
                    continue
                if lbl not in first_seen:
                    first_seen[lbl] = t
                last_seen[lbl] = t
            if progress_cb is not None:
                progress_cb(t + 1, T)

        # Bulk-write to cache (don't overwrite entries that were already
        # invalidated / re-populated by edits that raced with the scan).
        for lbl, t_first in first_seen.items():
            if lbl not in self._lifetime_cache:
                self._lifetime_cache[lbl] = (t_first, last_seen[lbl])
        if self._max_label is None:
            self._max_label = cur_max

    def invalidate_lifetime(self, label_id: int) -> None:
        """Remove ``label_id`` from the lifetime cache.

        Called by :meth:`commit` for every label touched by an operation
        so that the next :meth:`get_lifetime` call performs a fresh scan.
        """
        self._lifetime_cache.pop(int(label_id), None)

    # ------------------------------------------------------------------
    # Commit (replaces the old apply/undo/redo)
    # ------------------------------------------------------------------
    def commit(self, op: OpResult) -> None:
        """Apply an :class:`OpResult` from a primitive to the dirty set.

        Unlike the old ``apply()`` method, **no undo snapshot is stored**,
        which keeps RAM usage proportional only to the number of dirty
        frames (not the number of operations performed).  To recover from
        an unwanted edit, use :meth:`discard` to revert all unsaved changes.

        Empty ops (``op.new_frames`` is empty) are ignored.
        """
        if not op.new_frames:
            return

        changed: List[int] = []
        for t, post_frame in op.new_frames.items():
            t = int(t)
            self._dirty[t] = np.asarray(post_frame).astype(self.dtype, copy=False)
            changed.append(t)

        # --- Update performance caches -----------------------------------
        if self._max_label is not None:
            for frame in op.new_frames.values():
                if frame.size:
                    self._max_label = max(self._max_label, int(frame.max()))
        for label_id in op.affected_labels:
            self.invalidate_lifetime(label_id)

        self._emit(changed)

    # backward-compat alias so callers using the old name still work
    def apply(self, op: OpResult) -> None:  # noqa: D401
        """Alias for :meth:`commit` (backward-compatibility)."""
        self.commit(op)

    # ------------------------------------------------------------------
    # State queries
    # ------------------------------------------------------------------
    def is_dirty(self) -> bool:
        return bool(self._dirty)

    def dirty_count(self) -> int:
        return len(self._dirty)

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

        # Open the zarr array ONCE so all per-frame writes share the same
        # store handle.  Closing the store at the end flushes any internal
        # write buffers (important on Windows where per-call open/close can
        # leave chunks unlocked before the OS has flushed them).
        import zarr as _zarr
        try:
            arr = _zarr.open_array(str(self.zarr_path), mode="r+")
        except Exception:
            # Fallback: write frame-by-frame via the helper (pre-existing path).
            for t in sorted(self._dirty):
                write_zarr_parallel(self.zarr_path, index=t, data=self._dirty[t])
            arr = None

        n = 0
        if arr is not None:
            for t in sorted(self._dirty):
                arr[int(t)] = np.asarray(self._dirty[t])
                n += 1
            try:
                arr.store.close()
            except Exception:
                pass
        else:
            n = len(self._dirty)

        # The dirty data is now on disk: clear in-memory mutations and
        # refresh the underlying dask view so subsequent peek()s read
        # the new bytes.
        self._dirty.clear()
        self._cache.clear()
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
        self._max_label = None
        self._lifetime_cache.clear()
        self._darr = load_zarr(self.zarr_path)
        self._emit(list(range(int(self.shape[0]))))
