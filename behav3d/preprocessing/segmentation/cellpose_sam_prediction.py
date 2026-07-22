"""Cellpose-SAM (Cellpose v4) segmentation.

Zero-shot instance segmentation: unlike the Cellpose v3 path in
:mod:`behav3d.preprocessing.segmentation.cellpose_prediction`, no model file,
annotations or fine-tuning are required.

Cellpose v4 cannot share an interpreter with the ``cellpose==3.1.1.2`` that the
v3 path pins, so inference happens in a sidecar environment (see
:mod:`behav3d.preprocessing.segmentation.cpsam_env`) via
:mod:`behav3d.preprocessing.segmentation._cpsam_worker`. This module is the
parent half: it resolves channels/anisotropy from metadata, pre-allocates the
output zarr, drives the worker, and syncs the metadata CSV.

Output matches the pipeline contract exactly - a ``(T, Z, Y, X)`` ``uint16``
label zarr at ``{output_dir}/images/{sample}/{sample}_{cell_type}_segments.zarr``
plus the ``{prefix}_{cell_type}_segments_image_path`` metadata column - so every
tracker consumes it with no changes.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import tempfile
import threading
from collections import deque
from pathlib import Path
from typing import Callable, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm

from behav3d.io.images import load_image
from behav3d.io.formats.zarr import write_zarr_parallel
from behav3d.core.metadata import (
    detect_organoid_types_from_metadata,
    detect_immune_cell_types_from_metadata,
    detect_other_cell_types_from_metadata,
)
from behav3d.preprocessing.segmentation.cpsam_env import (
    build_worker_env,
    cpsam_env_status,
    find_cpsam_python,
)

_WORKER = Path(__file__).with_name("_cpsam_worker.py")

#: Available zero-shot backbones. ``cpsam`` is the original Cellpose-SAM
#: (April 2025); ``cpsam_v2`` (June 2026) is cellpose's current default;
#: ``cpdino`` uses a DINOv3 backbone and needs an extra install, so it is not
#: offered by default.
CPSAM_MODELS = ("cpsam", "cpsam_v2")

#: BEHAV3D defaults for Cellpose-SAM. These intentionally differ from cellpose's
#: own defaults in two places, both of which would otherwise destroy data:
#:
#: * ``min_size``: cellpose defaults to 15 and filters *inside* eval, which is
#:   irreversible and cannot be previewed. We keep everything and apply
#:   BEHAV3D's own previewable size filter afterwards.
#: * ``max_size_fraction``: cellpose defaults to 0.4, silently deleting any
#:   object covering >40% of the frame. That routinely eats large organoids.
DEFAULT_SAM_PARAMS = {
    "diameter": None,            # blank => auto; Cellpose-SAM is size-robust
    "flow_threshold": 0.4,       # ignored by cellpose when do_3D=True
    "cellprob_threshold": 0.0,
    "niter": 0,                  # 0 => proportional to diameter
    "do_3D": True,
    "stitch_threshold": 0.0,     # only used when do_3D=False
    "flow3D_smooth": 0.0,
    "batch_size": 8,             # lower this first if the GPU runs out of memory
    "min_size": 1,
    "max_size_fraction": 1.0,
    "norm_percentile_low": 1.0,
    "norm_percentile_high": 99.0,
    "norm3D": True,
    "sharpen_radius": 0.0,
    "smooth_radius": 0.0,
    "tile_norm_blocksize": 0.0,
    "tile_norm_smooth3D": 0.0,
    "invert": False,
    # BEHAV3D post-process, not a cellpose setting: drop labels that occupy a
    # single slice along Z, Y or X. Cellpose leaves flat fragments on the first
    # and last Z slice (especially in 2D+stitch mode), and they are otherwise
    # indistinguishable from real objects to the trackers.
    "drop_2d_segments": True,
}

#: Interactive size filter, in native voxels. ``None``/0 disables a bound.
DEFAULT_SIZE_FILTER = {"size_min": 0, "size_max": 0}

#: Version tag for the on-disk progress journal, so a future format change can be
#: detected instead of silently misreading an old file as "nothing done yet".
_JOURNAL_VERSION = 1


def power_safe_settings(n_cpus: Optional[int] = None) -> dict:
    """Return the conservative overrides for power-constrained machines.

    A full Cellpose-SAM run is one of the rare workloads that pins the GPU *and*
    every CPU core for minutes at a time. Laptops in particular size their power
    delivery for one rail at a time, and the shortfall on the other is covered by
    the battery; once that battery ages, the transient sag can trip the machine's
    over-current protection and cut power outright - no traceback, no dump, just
    off. These settings trade throughput for a much flatter power profile.
    """
    if n_cpus is None:
        n_cpus = os.cpu_count() or 4
    return {
        "batch_size": 2,                     # smaller GPU bursts
        "n_threads": max(1, int(n_cpus) // 4),  # leave most cores idle
        "cooldown_s": 5.0,                   # recovery window between frames
    }


class CellposeSAMEnvironmentError(RuntimeError):
    """Raised when the Cellpose-SAM sidecar environment is missing or broken."""


def _resolve_category(cell_type, organoid_types, immune_types, other_types):
    """Return ``(prefix, category_name)`` for *cell_type*, mirroring the v3 path."""
    if cell_type in organoid_types:
        return "or", "organoid"
    if cell_type in immune_types:
        return "im", "immune"
    if cell_type in other_types:
        return "ot", "other"
    raise ValueError(
        f"Cell type '{cell_type}' not found in metadata. Available types - "
        f"Organoids: {organoid_types}, Immune: {immune_types}, Other: {other_types}"
    )


def _channels_for_cell_type(config: Optional[dict], cell_type: str) -> list:
    """Return the raw channel indices selected for *cell_type*.

    Reads ``cellpose_sam.channel_selection`` (set on the Cellpose-SAM page's
    channel selection panel) - a plain ``{cell_type: [channel_idx, ...]}``
    map, uniform across samples (like ``cellpose_sam.size_filter``). This is
    intentionally separate from classic Cellpose v3's ``cellpose.channel_labels``
    block; the two backends no longer share channel configuration.
    """
    channel_selection = ((config or {}).get("cellpose_sam") or {}).get("channel_selection") or {}
    channels = channel_selection.get(cell_type) or []
    if not channels:
        raise ValueError(
            f"No channels selected for '{cell_type}'. Configure it on the "
            "Cellpose-SAM page's channel selection panel."
        )
    return sorted(int(c) for c in channels)


def _ensure_sidecar(config: dict | None = None) -> Path:
    """Return the sidecar interpreter, or raise with actionable guidance."""
    status = cpsam_env_status(config)
    if not status["available"]:
        raise CellposeSAMEnvironmentError(
            f"Cellpose-SAM environment is not ready: {status['error']}\n"
            "Use 'Set up Cellpose-SAM environment' in the segmentation panel, or run "
            "behav3d.preprocessing.segmentation.cpsam_env.create_cpsam_env()."
        )
    python = find_cpsam_python(config)
    if python is None:  # pragma: no cover - guarded by status above
        raise CellposeSAMEnvironmentError("Cellpose-SAM interpreter disappeared.")
    return python


def _describe_device(ready_event: dict) -> str:
    """Render a worker ``ready`` event as a human-checkable device line.

    ``cuda:0`` alone is not verifiable by the user: it is an index into whatever
    ``CUDA_VISIBLE_DEVICES`` exposed, not a physical slot. Naming the card, its
    VRAM and where the weights actually landed makes a silent CPU fallback or a
    wrong-GPU selection visible in the log instead of only as a slow run.
    """
    name = ready_event.get("gpu_name")
    if not name:
        return "device=CPU"
    parts = [f"device={ready_event.get('device')} ({name}"]
    if ready_event.get("gpu_total_mb"):
        parts.append(f", {ready_event['gpu_total_mb']} MB")
    parts.append(f", CUDA_VISIBLE_DEVICES={ready_event.get('cuda_visible_devices')})")
    line = "".join(parts)
    weights = ready_event.get("weights_device")
    if weights and not str(weights).startswith("cuda"):
        line += f"  [WARNING] weights are on {weights} - running on CPU!"
    return line


def _stream_stderr(proc: subprocess.Popen, tail: deque, lock: threading.Lock) -> None:
    """Echo the worker's stderr to the real stderr live, in a background thread.

    The worker redirects cellpose's own logger/tqdm chatter to stderr. Reading
    it only after ``proc.wait()`` (the previous approach) means a single long
    ``model.eval()`` call gives zero terminal feedback until it's already
    finished. Streaming it concurrently surfaces that chatter - including
    cellpose's own progress bars - live while the call is still running.
    """
    if proc.stderr is None:
        return
    try:
        for line in proc.stderr:
            sys.stderr.write(line)
            sys.stderr.flush()
            with lock:
                tail.append(line.rstrip("\n"))
    except Exception:
        pass


def _run_worker(
    job: dict,
    python: Path,
    device: str,
    force_cpu: bool,
    on_event: Optional[Callable[[dict], None]] = None,
    should_cancel: Optional[Callable[[], bool]] = None,
    n_threads: Optional[int] = None,
) -> None:
    """Spawn the sidecar worker for one job and stream its JSON Lines events."""
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False, encoding="utf-8") as fh:
        json.dump(job, fh)
        job_path = Path(fh.name)

    stderr_tail: deque = deque(maxlen=20)
    stderr_lock = threading.Lock()
    try:
        proc = subprocess.Popen(
            [str(python), str(_WORKER), str(job_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
            env=build_worker_env(device=device, force_cpu=force_cpu, n_threads=n_threads),
        )

        stderr_thread = threading.Thread(
            target=_stream_stderr, args=(proc, stderr_tail, stderr_lock), daemon=True,
        )
        stderr_thread.start()

        error: Optional[dict] = None
        assert proc.stdout is not None
        for line in proc.stdout:
            line = line.strip()
            if not line:
                continue
            try:
                event = json.loads(line)
            except json.JSONDecodeError:
                # Anything non-JSON on stdout is a bug in the worker's stream
                # hygiene; surface it rather than silently dropping it.
                with stderr_lock:
                    stderr_tail.append(line)
                continue
            if event.get("event") == "error":
                error = event
            elif on_event is not None:
                on_event(event)

            if should_cancel is not None and should_cancel():
                proc.terminate()
                raise RuntimeError("Cellpose-SAM run cancelled.")

        proc.wait()
        stderr_thread.join(timeout=5)

        if error is not None:
            raise RuntimeError(
                f"Cellpose-SAM worker failed: {error.get('message')}\n{error.get('traceback', '')}"
            )
        if proc.returncode != 0:
            with stderr_lock:
                tail_msg = "\n".join(stderr_tail)
            raise RuntimeError(
                f"Cellpose-SAM worker exited with code {proc.returncode}:\n{tail_msg}"
            )
    finally:
        job_path.unlink(missing_ok=True)


def _timepoints_for(t_total: int, timepoint_range: Optional[Tuple[int, int]]) -> list:
    """Clamp *timepoint_range* to the data and return explicit global T indices."""
    if timepoint_range is None:
        return list(range(t_total))
    start, end = timepoint_range
    start = max(0, min(int(start), t_total - 1))
    end = max(start, min(int(end), t_total - 1))
    return list(range(start, end + 1))


def _effective_params(sam_params=None, size_filter=None) -> Tuple[dict, dict, bool]:
    """Merge caller overrides onto the BEHAV3D defaults.

    Split out of :func:`_build_job` so the resume fingerprint is computed from the
    same fully-resolved values the worker actually runs with - comparing the raw
    overrides would call two runs identical whenever one of them relied on a
    default that later changed.
    """
    sam = dict(DEFAULT_SAM_PARAMS)
    sam.update(sam_params or {})
    resolved_filter = dict(DEFAULT_SIZE_FILTER)
    resolved_filter.update(size_filter or {})
    # drop_2d_segments is a BEHAV3D post-process rather than a cellpose eval
    # kwarg, so it is lifted out of the sam block for the worker.
    drop_2d = bool(sam.pop("drop_2d_segments", True))
    return sam, resolved_filter, drop_2d


def _build_job(**kwargs) -> dict:
    """Assemble a worker job spec, filling in BEHAV3D defaults."""
    sam, size_filter, drop_2d = _effective_params(
        kwargs.pop("sam_params", None), kwargs.pop("size_filter", None)
    )
    return {
        "sam": sam,
        "size_filter": size_filter,
        "drop_2d_segments": drop_2d,
        **kwargs,
    }


# ── resume journal ──────────────────────────────────────────────────────────
#
# A power loss or a kill mid-run leaves a partially filled label zarr. The array
# is pre-allocated and written one timepoint at a time, so the frames already in
# it are perfectly good - the only thing missing is a record of *which* ones, and
# a promise that they were produced with the settings now being requested. That
# is all this journal is.


def _journal_path(out_zarr: Path) -> Path:
    """Return the journal path beside *out_zarr* (``..._segments.progress.json``)."""
    return out_zarr.with_suffix(".progress.json")


def _params_fingerprint(model_name, channels, anisotropy, sam, size_filter, drop_2d) -> str:
    """Hash everything that changes what a segmented frame looks like.

    Resuming is only sound when the frames already on disk would have come out the
    same way today. Anything that feeds the worker's ``eval`` call or its
    post-processing therefore belongs in here; run bookkeeping (timepoint range,
    device, thread count, cooldown) deliberately does not, because none of it
    changes the resulting labels.
    """
    payload = json.dumps(
        {
            "model_name": str(model_name),
            "channels": list(channels),
            "anisotropy": round(float(anisotropy), 6),
            "sam": sam,
            "size_filter": size_filter,
            "drop_2d_segments": bool(drop_2d),
        },
        sort_keys=True,
        default=str,
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


def _read_journal(path: Path) -> Optional[dict]:
    """Return the journal at *path*, or ``None`` if absent/unreadable/stale."""
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        # A truncated journal is the expected casualty of a power cut. Treating it
        # as "no journal" costs a recompute; trusting it could corrupt the output.
        return None
    if not isinstance(data, dict) or data.get("version") != _JOURNAL_VERSION:
        return None
    return data


def _write_journal(path: Path, data: dict) -> None:
    """Write *data* to *path* atomically.

    Written on every frame of a run whose whole purpose is surviving abrupt power
    loss, so a half-written journal is a real possibility rather than a theoretical
    one: fsync the temp file, then rename over the target.
    """
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


def run_cellpose_sam_segmentation(
    output_dir: str | Path,
    metadata,
    label_name: Optional[str] = None,
    only_cell_types: Optional[Sequence[str]] = None,
    model_name: str = "cpsam",
    timepoint_range: Optional[Tuple[int, int]] = None,
    device: str = "auto",
    force_cpu: bool = False,
    sam_params: Optional[dict] = None,
    size_filter: Optional[dict] = None,
    overwrite_existing: bool = True,
    skip_existing: bool = False,
    resume: bool = True,
    n_threads: Optional[int] = None,
    cooldown_s: float = 0.0,
    progress_cb: Optional[Callable] = None,
    log_callback: Optional[Callable[[str], None]] = None,
    should_cancel: Optional[Callable[[], bool]] = None,
    config: Optional[dict] = None,
    **_kwargs,
):
    """Run Cellpose-SAM for one or many cell types across every sample.

    Parameters
    ----------
    label_name:
        Single cell type to segment. Ignored when *only_cell_types* is given.
    only_cell_types:
        Run these cell types in one batch. Each one uses the channels selected
        for it on the Cellpose-SAM channel selection panel.
    timepoint_range:
        ``(start_t, end_t)`` inclusive, in global indices. Frames outside the
        range keep whatever they already contained.
    device:
        ``"auto"``, ``"cpu"``, or ``"cuda:<n>"``. *force_cpu* overrides it.
    sam_params / size_filter:
        Overrides for :data:`DEFAULT_SAM_PARAMS` / :data:`DEFAULT_SIZE_FILTER`.
    resume:
        Continue an interrupted run instead of starting over: keep the timepoints
        already recorded in the sidecar journal and segment only the rest. Raises
        if the settings no longer match the ones those frames were produced with.
    n_threads / cooldown_s:
        CPU thread cap and per-frame idle time. See :func:`power_safe_settings` -
        both exist to keep a machine with marginal power delivery alive through a
        long run, at the cost of wall-clock time.

    Returns
    -------
    (metadata, summary) where summary is ``{"processed": [...], "skipped": [...]}``.
    """
    assert isinstance(metadata, pd.DataFrame), "metadata must be a pandas DataFrame"
    log = log_callback or print
    output_dir = Path(output_dir)

    python = _ensure_sidecar(config)

    cell_types = list(only_cell_types) if only_cell_types else ([label_name] if label_name else [])
    if not cell_types:
        raise ValueError("Provide either label_name or only_cell_types.")

    organoid_types = detect_organoid_types_from_metadata(metadata)
    immune_types = detect_immune_cell_types_from_metadata(metadata)
    other_types = detect_other_cell_types_from_metadata(metadata)

    summary = {"processed": [], "skipped": []}
    units = [(ct, idx, sample) for ct in cell_types for idx, sample in metadata.iterrows()]

    # Precompute how many timepoints will actually run, across every unit, so
    # progress can be reported per-timepoint (smoothly, within a single long
    # sample) rather than only once per cell_type x sample unit.
    total_ticks = 0
    for cell_type, _idx, sample in units:
        sample_name = sample["sample_name"]
        sample_dir = output_dir / "images" / str(sample_name)
        raw_zarr = sample_dir / f"{sample_name}.zarr"
        out_zarr = sample_dir / f"{sample_name}_{cell_type}_segments.zarr"
        if not raw_zarr.exists() or (skip_existing and out_zarr.exists()):
            continue
        try:
            t_total = load_image(raw_zarr).shape[0]
        except Exception:
            continue
        wanted = _timepoints_for(t_total, timepoint_range)
        if resume:
            # Best-effort only: this drives the progress bar's denominator. The
            # authoritative (fingerprint-checked) filtering happens per unit below,
            # so a journal that turns out to be unusable just makes the bar jump.
            journal = _read_journal(_journal_path(out_zarr))
            if journal:
                done = set(journal.get("done") or ())
                wanted = [t for t in wanted if t not in done]
        total_ticks += len(wanted)

    completed_ticks = 0
    pbar = tqdm(total=total_ticks, desc="Cellpose-SAM", unit="frame",
                file=sys.stdout, dynamic_ncols=True)
    try:
        for unit_i, (cell_type, idx, sample) in enumerate(units):
            prefix, category = _resolve_category(cell_type, organoid_types, immune_types, other_types)
            path_col = f"{prefix}_{cell_type}_segments_image_path"
            if path_col not in metadata.columns:
                metadata[path_col] = pd.NA
            metadata[path_col] = metadata[path_col].astype("object")

            sample_name = sample["sample_name"]
            tag = f"{cell_type} / {sample_name}"
            pbar.set_description(f"Cellpose-SAM {tag}")
            if progress_cb is not None:
                try:
                    progress_cb(completed_ticks, total_ticks, f"Cellpose-SAM {tag}")
                except Exception:
                    pass

            sample_dir = output_dir / "images" / str(sample_name)
            sample_dir.mkdir(parents=True, exist_ok=True)
            raw_zarr = sample_dir / f"{sample_name}.zarr"
            if not raw_zarr.exists():
                log(f"  [SKIP] Zarr file missing for {sample_name}, run image preprocessing first.")
                summary["skipped"].append(f"{tag} (no zarr)")
                continue

            out_zarr = sample_dir / f"{sample_name}_{cell_type}_segments.zarr"
            if skip_existing and out_zarr.exists():
                log(f"  [SKIP] {tag}: segments already exist.")
                summary["skipped"].append(f"{tag} (exists)")
                metadata.at[idx, path_col] = str(out_zarr)
                continue

            images = load_image(raw_zarr)
            t_total, _c, z, y, x = images.shape
            timepoints = _timepoints_for(t_total, timepoint_range)

            channels = _channels_for_cell_type(config, cell_type)

            anisotropy = float(sample["pixel_distance_z"]) / float(sample["pixel_distance_xy"])

            eff_sam, eff_size_filter, eff_drop_2d = _effective_params(sam_params, size_filter)
            fingerprint = _params_fingerprint(
                model_name, channels, anisotropy, eff_sam, eff_size_filter, eff_drop_2d
            )

            # Recreate the array only when it is absent, mis-shaped, or a full-range
            # overwrite was requested. Never recreate for a sub-range run: doing so
            # would zero every timepoint outside the range. Never recreate when
            # resuming either - that is the whole point.
            recreate = overwrite_existing and timepoint_range is None and not resume
            if not out_zarr.exists():
                recreate = True
            else:
                try:
                    import zarr

                    existing = zarr.open(str(out_zarr), mode="r")
                    if tuple(existing.shape) != (t_total, z, y, x) or np.dtype(existing.dtype) != np.uint16:
                        recreate = True
                except Exception:
                    recreate = True

            journal_file = _journal_path(out_zarr)
            journal = None if recreate else _read_journal(journal_file)
            # A journal only vouches for frames produced under the settings it
            # records; once those differ it says nothing about what is on disk.
            if journal is not None and journal.get("fingerprint") != fingerprint:
                if resume:
                    raise RuntimeError(
                        f"Cannot resume {tag}: the segmentation settings have changed since "
                        f"the interrupted run (channels, model, anisotropy, cellpose params "
                        f"or size filter).\nResuming would mix frames segmented under two "
                        f"different settings in one array.\nEither restore the previous "
                        f"settings, or untick 'Resume interrupted run' to re-segment from "
                        f"scratch."
                    )
                journal = None

            already_done = set(journal.get("done") or ()) if journal else set()
            if resume and already_done:
                remaining = [t for t in timepoints if t not in already_done]
                n_skip = len(timepoints) - len(remaining)
                if not remaining:
                    log(f"  [SKIP] {tag}: all {len(timepoints)} timepoints already segmented.")
                    summary["skipped"].append(f"{tag} (complete)")
                    metadata.at[idx, path_col] = str(out_zarr)
                    continue
                log(f"  [RESUME] {tag}: {n_skip} timepoints already done, {len(remaining)} to go.")
                timepoints = remaining

            write_zarr_parallel(
                out_zarr,
                shape=(t_total, z, y, x),
                dtype=np.uint16,
                chunks=(1, z, y, x),
                overwrite=recreate,
            )

            journal = {
                "version": _JOURNAL_VERSION,
                "fingerprint": fingerprint,
                "shape": [t_total, z, y, x],
                "dtype": "uint16",
                "channels": list(channels),
                "model_name": str(model_name),
                "done": sorted(already_done),
                "frames": (journal.get("frames") if journal else None) or {},
            }
            _write_journal(journal_file, journal)

            log(
                f"Cellpose-SAM ({model_name}) '{cell_type}' [{category}] on {sample_name}: "
                f"channel{'s' if len(channels) > 1 else ''} {channels}, "
                f"{len(timepoints)}/{t_total} timepoints, anisotropy {anisotropy:.2f}"
            )

            job = _build_job(
                mode="run",
                raw_zarr=str(raw_zarr),
                out_zarr=str(out_zarr),
                channels=channels,
                timepoints=timepoints,
                model_name=model_name,
                force_cpu=force_cpu,
                anisotropy=anisotropy,
                sam_params=sam_params,
                size_filter=size_filter,
                n_threads=n_threads,
                cooldown_s=cooldown_s,
            )

            def _on_event(event, _tag=tag, _journal=journal, _journal_file=journal_file):
                nonlocal completed_ticks
                if event.get("event") == "progress":
                    completed_ticks += 1
                    pbar.update(1)
                    # Persisted per frame, before anything else: the failure this
                    # guards against gives no warning and no chance to flush later.
                    t_done = int(event["t"])
                    if t_done not in _journal["done"]:
                        _journal["done"].append(t_done)
                        _journal["done"].sort()
                    _journal["frames"][str(t_done)] = {
                        "n_objects": event.get("n_objects"),
                        "seconds": event.get("seconds"),
                        "peak_vram_mb": event.get("peak_vram_mb"),
                    }
                    _write_journal(_journal_file, _journal)
                    if progress_cb is not None:
                        try:
                            progress_cb(
                                completed_ticks, total_ticks,
                                f"Cellpose-SAM {_tag}: t={event['t']} "
                                f"({event['n_done']}/{event['total']}, {event['n_objects']} objects)",
                            )
                        except Exception:
                            pass
                elif event.get("event") == "log":
                    log(f"  {event.get('msg')}")
                elif event.get("event") == "ready":
                    log(f"  {_describe_device(event)} cellpose={event.get('cellpose_version')}")

            _run_worker(
                job, python, device, force_cpu,
                on_event=_on_event, should_cancel=should_cancel, n_threads=n_threads,
            )

            log(f"  [SAVED] {out_zarr}")
            metadata.at[idx, path_col] = str(out_zarr)
            summary["processed"].append(tag)

        if progress_cb is not None:
            try:
                progress_cb(total_ticks, total_ticks, "Cellpose-SAM done")
            except Exception:
                pass
    finally:
        pbar.close()
    return metadata, summary


def run_cellpose_sam_and_sync_metadata(
    output_dir,
    metadata_loader,
    label_name=None,
    only_cell_types=None,
    **kwargs,
):
    """Run Cellpose-SAM and persist the updated metadata CSV.

    Mirrors :func:`behav3d.preprocessing.segmentation.cellpose_prediction.run_cellpose_and_sync_metadata`
    so both GUI surfaces and the queue can call it the same way.
    """
    updated_metadata, summary = run_cellpose_sam_segmentation(
        output_dir=output_dir,
        metadata=metadata_loader.metadata,
        label_name=label_name,
        only_cell_types=only_cell_types,
        config=getattr(metadata_loader, "behav3d_parameters", None),
        **kwargs,
    )
    metadata_loader.metadata = updated_metadata
    updated_metadata.to_csv(metadata_loader.metadata_csv_path, index=False)
    return updated_metadata, summary


def preview_cellpose_sam(
    output_dir: str | Path,
    metadata,
    sample_name: str,
    cell_type: str,
    timepoint: int = 0,
    model_name: str = "cpsam",
    device: str = "auto",
    force_cpu: bool = False,
    sam_params: Optional[dict] = None,
    log_callback: Optional[Callable[[str], None]] = None,
    config: Optional[dict] = None,
) -> np.ndarray:
    """Segment a single ``(Z, Y, X)`` volume and return the **unfiltered** labels.

    Size filtering is deliberately not applied: the caller filters the returned
    array in-memory so the interactive preview can explore thresholds without
    re-running inference.
    """
    log = log_callback or print
    output_dir = Path(output_dir)
    python = _ensure_sidecar(config)

    row = metadata.loc[metadata["sample_name"] == sample_name]
    if row.empty:
        raise ValueError(f"Sample '{sample_name}' not found in metadata.")
    sample = row.iloc[0]

    raw_zarr = output_dir / "images" / str(sample_name) / f"{sample_name}.zarr"
    if not raw_zarr.exists():
        raise FileNotFoundError(f"Converted zarr missing for {sample_name}: {raw_zarr}")

    channels = _channels_for_cell_type(config, cell_type)

    t_total = load_image(raw_zarr).shape[0]
    timepoint = max(0, min(int(timepoint), t_total - 1))
    anisotropy = float(sample["pixel_distance_z"]) / float(sample["pixel_distance_xy"])

    with tempfile.TemporaryDirectory() as tmp:
        preview_out = Path(tmp) / "preview.npy"
        job = _build_job(
            mode="preview",
            raw_zarr=str(raw_zarr),
            out_zarr=None,
            preview_out=str(preview_out),
            channels=channels,
            timepoints=[timepoint],
            model_name=model_name,
            force_cpu=force_cpu,
            anisotropy=anisotropy,
            sam_params=sam_params,
            size_filter={},
        )
        log(f"Cellpose-SAM preview: {sample_name} / {cell_type} / t={timepoint}, channels {channels}")

        def _on_event(e):
            if e.get("event") == "ready":
                log(f"  {_describe_device(e)} cellpose={e.get('cellpose_version')}")
            elif e.get("event") == "progress":
                vram = f", peak VRAM {e['peak_vram_mb']} MB" if e.get("peak_vram_mb") else ""
                log(f"  t={e['t']}: {e['n_objects']} objects in {e['seconds']}s{vram}")

        _run_worker(job, python, device, force_cpu, on_event=_on_event)
        return np.load(preview_out)
