"""Shared on-disk progress journal for segmentation backends.

A power loss or a kill mid-run leaves a partially filled label zarr. The array is
pre-allocated and written one timepoint at a time, so the frames already in it are
perfectly good - the only thing missing is a record of *which* ones, and a promise
that they were produced with the settings now being requested. That is all this
journal is.

The design started life inside
:mod:`behav3d.preprocessing.segmentation.cellpose_sam_prediction` and was lifted here
so APOC and ConvPaint can share it instead of each inventing their own bookkeeping.
The engines that use it agree on three rules:

* One journal sidecar per output zarr array, named ``<array>.progress.json``.
* A journal vouches for frames **only** under the fingerprint it records. Once the
  settings differ it says nothing about what is on disk, and resuming into it would
  mix two different segmentations in one array.
* Losing a journal costs a recompute. Trusting a damaged one corrupts an output.
  Every failure path here therefore degrades towards "no journal".
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
from pathlib import Path
from typing import NamedTuple, Optional

#: Version tag for the on-disk progress journal, so a future format change can be
#: detected instead of silently misreading an old file as "nothing done yet".
JOURNAL_VERSION = 1

#: Above this size a model file is fingerprinted by (size, mtime) instead of by
#: content. Hashing a large ConvPaint pickle on every sample would cost more than the
#: check is worth; APOC's .cl classifiers are a few KB of text and always get the
#: stronger content hash.
_HASH_SIZE_LIMIT = 4 * 1024 * 1024


def journal_path(out_zarr) -> Path:
    """Return the journal path beside *out_zarr* (``..._segments.progress.json``)."""
    return Path(out_zarr).with_suffix(".progress.json")


def params_fingerprint(payload: dict) -> str:
    """Hash everything that changes what a segmented frame looks like.

    Resuming is only sound when the frames already on disk would have come out the
    same way today. Anything that feeds inference or post-processing therefore
    belongs in *payload*; run bookkeeping (timepoint range, device, thread count,
    tiling) deliberately does not, because none of it changes the resulting labels.
    """
    blob = json.dumps(payload, sort_keys=True, default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:16]


def file_fingerprint(path) -> Optional[str]:
    """Identify a classifier/model file for the purposes of a fingerprint.

    Returns ``None`` when the file is missing, which folds into the payload as a
    distinct value - a run made without a classifier must not match one made with it.
    """
    p = Path(path) if path else None
    if p is None or not p.exists():
        return None
    try:
        stat = p.stat()
        if stat.st_size <= _HASH_SIZE_LIMIT:
            return hashlib.sha256(p.read_bytes()).hexdigest()[:16]
        return f"{stat.st_size}:{stat.st_mtime_ns}"
    except Exception:
        # An unreadable classifier is a problem, but not this function's problem to
        # report; returning None just means "cannot vouch for a resume".
        return None


def read_journal(path) -> Optional[dict]:
    """Return the journal at *path*, or ``None`` if absent/unreadable/stale."""
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except Exception:
        # A truncated journal is the expected casualty of a power cut. Treating it
        # as "no journal" costs a recompute; trusting it could corrupt the output.
        return None
    if not isinstance(data, dict) or data.get("version") != JOURNAL_VERSION:
        return None
    return data


def write_journal(path, data: dict) -> None:
    """Write *data* to *path* atomically.

    Written on every frame of a run whose whole purpose is surviving abrupt power
    loss, so a half-written journal is a real possibility rather than a theoretical
    one: fsync the temp file, then rename over the target.
    """
    path = Path(path)
    tmp = Path(str(path) + ".tmp")
    try:
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(data, fh)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    except Exception:
        # Losing the journal degrades resume to a full recompute; it must never
        # take the segmentation itself down with it.
        try:
            tmp.unlink(missing_ok=True)
        except Exception:
            pass


def new_journal(fingerprint: str, shape, dtype: str = "uint16",
                done=None, extra: Optional[dict] = None) -> dict:
    """Build a fresh journal dict carrying *fingerprint* and *shape*.

    *shape* is the full array shape; ``shape[0]`` is what :func:`journal_state` later
    reports completeness against.
    """
    journal = {
        "version": JOURNAL_VERSION,
        "fingerprint": fingerprint,
        "shape": [int(s) for s in shape],
        "dtype": str(dtype),
        "done": sorted(int(t) for t in (done or ())),
        "frames": {},
    }
    if extra:
        journal.update(extra)
    return journal


def mark_done(journal: dict, path, t: int, info: Optional[dict] = None) -> None:
    """Record timepoint *t* as written and persist the journal immediately.

    Persisted per frame, before anything else: the failure this guards against gives
    no warning and no chance to flush later.
    """
    t = int(t)
    if t not in journal["done"]:
        journal["done"].append(t)
        journal["done"].sort()
    if info:
        journal.setdefault("frames", {})[str(t)] = info
    write_journal(path, journal)


def check_fingerprint(journal: Optional[dict], fingerprint: str, label: str) -> None:
    """Raise when *journal* vouches for frames made under different settings.

    A no-op when *journal* is ``None`` - an array with no journal predates this
    bookkeeping and has no recorded settings to disagree with.
    """
    if journal is None:
        return
    if journal.get("fingerprint") == fingerprint:
        return
    raise RuntimeError(
        f"Cannot resume {label}: the segmentation settings have changed since the "
        f"existing run (classifier, channels, strategy or post-processing "
        f"parameters).\nResuming would mix frames segmented under two different "
        f"settings in one array.\nEither restore the previous settings, or choose "
        f"'Overwrite' to re-segment from scratch."
    )


def journal_complete(journal: Optional[dict], requested_timepoints=None) -> bool:
    """True when *journal* vouches for every timepoint that would be written.

    With *requested_timepoints* given (an explicit list), every one must be in
    ``done``. Without it, fall back to "every index ``0..t_total-1`` is present"
    from the recorded shape - the check preflight can make before it knows the
    range.

    The fallback checks the actual indices, not just ``len(done) >= t_total``: a
    malformed journal whose ``done`` carries an out-of-range or spurious value
    could hit that count while still missing a real frame. In keeping with this
    module's stance, such a journal must read as *partial* (recompute) and never
    as *complete* (skip).

    A complete array is never resumed into, so a fingerprint disagreement about
    it is moot: nothing new is written, nothing is mixed.
    """
    if not journal:
        return False
    done = {int(t) for t in (journal.get("done") or ())}
    req = list(requested_timepoints or ())
    if req:
        return all(int(t) in done for t in req)
    shape = journal.get("shape") or []
    t_total = int(shape[0]) if shape else 0
    if not t_total or len(done) < t_total:
        return False
    return all(t in done for t in range(t_total))


def journal_state(out_zarr):
    """Summarise how complete the array at *out_zarr* is, for the overwrite dialog.

    Reads only the sidecar JSON - never the raw image and never the zarr - so listing
    dozens of outputs in a modal dialog stays instant.

    Returns ``(state, n_done, t_total)`` where *state* is:

    ``"complete"``
        The journal vouches for every timepoint index ``0..t_total-1``.
    ``"partial"``
        An interrupted run, one whose settings changed part-way through a range,
        or a journal whose ``done`` does not cleanly cover every index.
    ``"unknown"``
        No usable journal. The array may be perfectly fine - it simply predates this
        bookkeeping, so nothing can be asserted about it.
    ``"missing"``
        No array on disk at all.
    """
    out_zarr = Path(out_zarr)
    if not out_zarr.exists():
        return ("missing", 0, 0)
    journal = read_journal(journal_path(out_zarr))
    if journal is None:
        return ("unknown", 0, 0)
    shape = journal.get("shape") or []
    t_total = int(shape[0]) if shape else 0
    n_done = len(journal.get("done") or ())
    if t_total and journal_complete(journal):
        return ("complete", n_done, t_total)
    return ("partial", n_done, t_total)


def describe_state(state: str, n_done: int, t_total: int) -> str:
    """Human-readable suffix for one entry in the overwrite prompt."""
    if state == "complete":
        return f"complete ({n_done}/{t_total})"
    if state == "partial":
        return f"INCOMPLETE ({n_done}/{t_total})"
    if state == "unknown":
        return "unknown (no journal)"
    return "not present"


def discard_journal(out_zarr) -> None:
    """Remove the sidecar for *out_zarr*, if any.

    Called whenever the array itself is recreated: a journal that outlives the data
    it describes is exactly the failure mode this module exists to end.
    """
    try:
        journal_path(out_zarr).unlink(missing_ok=True)
    except Exception:
        pass


def should_recreate(path, shape, dtype, overwrite_existing: bool,
                    timepoint_range) -> bool:
    """Decide whether an output array must be wiped and re-allocated.

    Recreate when it is absent, mis-shaped, or a full-range overwrite was requested.
    **Never** recreate for a sub-range run: doing so would zero every timepoint
    outside the range while leaving the journal claiming they are done.
    """
    import numpy as np

    path = Path(path)
    if not path.exists():
        return True
    try:
        import zarr

        existing = zarr.open(str(path), mode="r")
        if (tuple(existing.shape) != tuple(shape)
                or np.dtype(existing.dtype) != np.dtype(dtype)):
            return True
    except Exception:
        # An array we cannot even open is not one we can resume into.
        return True
    return bool(overwrite_existing) and timepoint_range is None


class ArrayPlan(NamedTuple):
    """What a run has decided to do about one output array.

    ``recreate``
        Wipe and re-allocate the array (and drop its journal).
    ``done``
        Timepoints the journal vouches for under the current fingerprint. Empty
        whenever the journal cannot be trusted or is being overwritten.
    ``complete``
        Nothing left to do for the requested timepoints - the caller can skip this
        array entirely.
    ``legacy``
        The array exists but has no journal: it predates this bookkeeping.
    """

    recreate: bool
    done: set
    complete: bool
    legacy: bool


def plan_output(path, shape, dtype, fingerprint: str, label: str, *,
                overwrite_existing: bool, timepoint_range,
                requested_timepoints) -> ArrayPlan:
    """Decide how to treat one output array before a run touches it.

    Centralises the three rules every engine needs to agree on:

    * A sub-range run never wipes the array (that would zero the frames outside
      the range while the journal still claimed them).
    * Under ``overwrite_existing`` the journal is not consulted at all - the user
      asked for a redo, and a fingerprint disagreement is not an error.
    * Otherwise a fingerprint disagreement *is* an error - but only for a
      **partial** array, the one case a run would actually resume into and so mix
      frames from two settings. A complete array is skipped whatever its
      fingerprint says. See :func:`check_fingerprint` and :func:`journal_complete`.

    An array with no journal is treated as complete for a full-range run. It may
    well be fine - it simply predates journalling - and recomputing every existing
    dataset the first time this code runs would be a costly surprise. It is reported
    as ``unknown`` in the overwrite dialog so the lack of proof stays visible.
    """
    path = Path(path)
    recreate = should_recreate(path, shape, dtype, overwrite_existing, timepoint_range)
    journal = None if recreate else read_journal(journal_path(path))
    legacy = (not recreate) and journal is None and path.exists()

    if journal is not None and overwrite_existing:
        journal = None

    done = set(journal.get("done") or ()) if journal else set()

    if overwrite_existing:
        complete = False
    elif legacy:
        complete = timepoint_range is None
    else:
        complete = journal_complete(journal, requested_timepoints)

    # A complete array is never written into, so it cannot mix two settings -
    # skip it without consulting the fingerprint. Only a partial array that this
    # run would actually resume into needs the settings to still match.
    if journal is not None and not overwrite_existing and not complete:
        check_fingerprint(journal, fingerprint, label)

    return ArrayPlan(recreate=recreate, done=done, complete=complete, legacy=legacy)


def preflight_conflicts(entries, *, overwrite_existing: bool, timepoint_range,
                        requested_timepoints=None) -> list:
    """Find every array whose journal disagrees with the settings about to be used.

    *entries* is an iterable of ``(path, shape, dtype, fingerprint, label)``.
    *requested_timepoints*, when given, lets an already-complete array be skipped
    outright: nothing new is written to it, so its fingerprint cannot matter.

    :func:`plan_output` raises the moment it meets a mismatch, which mid-batch would
    abort a run after some samples had already been rewritten. A fingerprint depends
    only on settings and classifier files - never on image data - so every mismatch is
    knowable up front. Engines call this before touching anything and raise once,
    naming all of them.

    Returns a list of ``(label, reason)``; empty means nothing conflicts.
    """
    if overwrite_existing:
        # The arrays are being redone regardless; a stale journal is not an error.
        return []
    conflicts = []
    for path, shape, dtype, fingerprint, label in entries:
        if should_recreate(path, shape, dtype, overwrite_existing, timepoint_range):
            continue
        journal = read_journal(journal_path(path))
        if journal is None:
            continue
        if journal_complete(journal, requested_timepoints):
            continue  # never resumed into - a fingerprint disagreement is moot
        if journal.get("fingerprint") != fingerprint:
            conflicts.append((label, "settings changed since it was written"))
    return conflicts


def raise_for_conflicts(conflicts, engine: str) -> None:
    """Raise a single actionable error for everything :func:`preflight_conflicts` found."""
    if not conflicts:
        return
    listing = "\n".join(f"  - {label}: {reason}" for label, reason in conflicts)
    raise RuntimeError(
        f"Cannot resume this {engine} run: the segmentation settings have changed "
        f"since these outputs were written.\n{listing}\n"
        f"Resuming would mix frames segmented under two different settings in one "
        f"array.\nEither restore the previous settings, or choose 'Overwrite' to "
        f"re-segment from scratch. Nothing has been modified."
    )


def open_output_zarr(path, dtype, shape, recreate: bool):
    """Open a zarr array for output, creating/recreating as needed.

    When *recreate* is set the array **and** its journal are removed first, so the
    two can never disagree about what is on disk.
    """
    import zarr

    path = Path(path)
    if recreate and path.exists():
        shutil.rmtree(path)
        discard_journal(path)
    if path.exists():
        return zarr.open(str(path), mode="r+")
    discard_journal(path)
    return zarr.open(
        str(path), mode="w", shape=shape, dtype=dtype,
        chunks=(1,) + tuple(shape[1:]),  # one chunk per timepoint
    )
