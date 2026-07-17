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

import json
import subprocess
import tempfile
from pathlib import Path
from typing import Callable, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import yaml

from behav3d.io.images import load_image
from behav3d.io.formats.zarr import write_zarr_parallel
from behav3d.core.metadata import (
    detect_organoid_types_from_metadata,
    detect_immune_cell_types_from_metadata,
    detect_other_cell_types_from_metadata,
)
from behav3d.preprocessing.segmentation.cellpose_prediction import (
    _label_to_channel_from_stored_map,
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


def _load_channel_config(output_dir: Path, channel_labels_config: Optional[dict]) -> dict:
    """Read the Cellpose channel-label map, reusing the existing `cellpose` config block."""
    if channel_labels_config is not None:
        return channel_labels_config
    config_path = Path(output_dir) / "behav3d_parameters.yml"
    if not config_path.exists():
        return {}
    with open(config_path, "r") as fh:
        config = yaml.safe_load(fh) or {}
    # Channel labels are a property of the images, not of the segmentation
    # backend, so Cellpose-SAM deliberately shares the `cellpose` block rather
    # than making users configure the same thing twice.
    return config.get("cellpose", {}) or {}


def _label_to_channel_for_sample(sample_name, channel_labels_config) -> dict:
    """Resolve the label -> channel-index map for one sample."""
    labels_mode = channel_labels_config.get("labels_mode", "same_for_all")
    global_labels = channel_labels_config.get("channel_labels", {})
    per_sample = channel_labels_config.get("per_sample_channel_labels", {})

    if labels_mode == "per_sample":
        if sample_name not in per_sample:
            raise ValueError(
                f"Sample '{sample_name}' not found in per_sample_channel_labels config. "
                f"Available samples: {list(per_sample.keys())}. Configure it in the "
                "'Configure Channel Labels' panel."
            )
        return _label_to_channel_from_stored_map(per_sample[sample_name])
    if global_labels:
        return _label_to_channel_from_stored_map(global_labels)
    raise ValueError(
        "No channel labels configured. Please configure channel labels in the "
        "'Configure Channel Labels' panel before running segmentation."
    )


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


def _run_worker(
    job: dict,
    python: Path,
    device: str,
    force_cpu: bool,
    on_event: Optional[Callable[[dict], None]] = None,
    should_cancel: Optional[Callable[[], bool]] = None,
) -> None:
    """Spawn the sidecar worker for one job and stream its JSON Lines events."""
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False, encoding="utf-8") as fh:
        json.dump(job, fh)
        job_path = Path(fh.name)

    stderr_tail: list[str] = []
    try:
        proc = subprocess.Popen(
            [str(python), str(_WORKER), str(job_path)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
            env=build_worker_env(device=device, force_cpu=force_cpu),
        )

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
        if proc.stderr is not None:
            stderr_tail.extend((proc.stderr.read() or "").splitlines()[-20:])

        if error is not None:
            raise RuntimeError(
                f"Cellpose-SAM worker failed: {error.get('message')}\n{error.get('traceback', '')}"
            )
        if proc.returncode != 0:
            raise RuntimeError(
                "Cellpose-SAM worker exited with code "
                f"{proc.returncode}:\n" + "\n".join(stderr_tail[-20:])
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


def _build_job(**kwargs) -> dict:
    """Assemble a worker job spec, filling in BEHAV3D defaults."""
    sam = dict(DEFAULT_SAM_PARAMS)
    sam.update(kwargs.pop("sam_params", None) or {})
    size_filter = dict(DEFAULT_SIZE_FILTER)
    size_filter.update(kwargs.pop("size_filter", None) or {})
    # drop_2d_segments is a BEHAV3D post-process rather than a cellpose eval
    # kwarg, so it is lifted out of the sam block for the worker.
    drop_2d = bool(sam.pop("drop_2d_segments", True))
    return {
        "sam": sam,
        "size_filter": size_filter,
        "drop_2d_segments": drop_2d,
        **kwargs,
    }


def run_cellpose_sam_segmentation(
    output_dir: str | Path,
    metadata,
    label_name: Optional[str] = None,
    only_cell_types: Optional[Sequence[str]] = None,
    model_name: str = "cpsam",
    timepoint_range: Optional[Tuple[int, int]] = None,
    channel_labels_config: Optional[dict] = None,
    device: str = "auto",
    force_cpu: bool = False,
    sam_params: Optional[dict] = None,
    size_filter: Optional[dict] = None,
    overwrite_existing: bool = True,
    skip_existing: bool = False,
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
        Run these cell types in one batch. Each one uses the channel assigned to
        it in the Cellpose channel-label configuration.
    timepoint_range:
        ``(start_t, end_t)`` inclusive, in global indices. Frames outside the
        range keep whatever they already contained.
    device:
        ``"auto"``, ``"cpu"``, or ``"cuda:<n>"``. *force_cpu* overrides it.
    sam_params / size_filter:
        Overrides for :data:`DEFAULT_SAM_PARAMS` / :data:`DEFAULT_SIZE_FILTER`.

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

    channel_cfg = _load_channel_config(output_dir, channel_labels_config)
    organoid_types = detect_organoid_types_from_metadata(metadata)
    immune_types = detect_immune_cell_types_from_metadata(metadata)
    other_types = detect_other_cell_types_from_metadata(metadata)

    summary = {"processed": [], "skipped": []}
    units = [(ct, idx, sample) for ct in cell_types for idx, sample in metadata.iterrows()]
    total_units = len(units)

    for unit_i, (cell_type, idx, sample) in enumerate(units):
        prefix, category = _resolve_category(cell_type, organoid_types, immune_types, other_types)
        path_col = f"{prefix}_{cell_type}_segments_image_path"
        if path_col not in metadata.columns:
            metadata[path_col] = pd.NA
        metadata[path_col] = metadata[path_col].astype("object")

        sample_name = sample["sample_name"]
        tag = f"{cell_type} / {sample_name}"
        if progress_cb is not None:
            try:
                progress_cb(unit_i, total_units, f"Cellpose-SAM {tag}")
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

        label_to_channel = _label_to_channel_for_sample(sample_name, channel_cfg)
        if cell_type not in label_to_channel:
            raise ValueError(
                f"Channel label '{cell_type}' not configured for sample '{sample_name}'. "
                f"Available labels: {list(label_to_channel.keys())}."
            )
        channels = [label_to_channel[cell_type]]

        anisotropy = float(sample["pixel_distance_z"]) / float(sample["pixel_distance_xy"])

        # Recreate the array only when it is absent, mis-shaped, or a full-range
        # overwrite was requested. Never recreate for a sub-range run: doing so
        # would zero every timepoint outside the range.
        recreate = overwrite_existing and timepoint_range is None
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

        write_zarr_parallel(
            out_zarr,
            shape=(t_total, z, y, x),
            dtype=np.uint16,
            chunks=(1, z, y, x),
            overwrite=recreate,
        )

        log(
            f"Cellpose-SAM ({model_name}) '{cell_type}' [{category}] on {sample_name}: "
            f"channel {channels[0]}, {len(timepoints)}/{t_total} timepoints, anisotropy {anisotropy:.2f}"
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
        )

        def _on_event(event, _tag=tag, _i=unit_i):
            if event.get("event") == "progress":
                if progress_cb is not None:
                    try:
                        progress_cb(
                            _i, total_units,
                            f"Cellpose-SAM {_tag}: t={event['t']} "
                            f"({event['n_done']}/{event['total']}, {event['n_objects']} objects)",
                        )
                    except Exception:
                        pass
            elif event.get("event") == "log":
                log(f"  {event.get('msg')}")
            elif event.get("event") == "ready":
                log(f"  device={event.get('device')} cellpose={event.get('cellpose_version')}")

        _run_worker(job, python, device, force_cpu, on_event=_on_event, should_cancel=should_cancel)

        log(f"  [SAVED] {out_zarr}")
        metadata.at[idx, path_col] = str(out_zarr)
        summary["processed"].append(tag)

    if progress_cb is not None:
        try:
            progress_cb(total_units, total_units, "Cellpose-SAM done")
        except Exception:
            pass
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
    channel_labels_config: Optional[dict] = None,
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

    channel_cfg = _load_channel_config(output_dir, channel_labels_config)
    label_to_channel = _label_to_channel_for_sample(sample_name, channel_cfg)
    if cell_type not in label_to_channel:
        raise ValueError(
            f"Channel label '{cell_type}' not configured for sample '{sample_name}'."
        )

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
            channels=[label_to_channel[cell_type]],
            timepoints=[timepoint],
            model_name=model_name,
            force_cpu=force_cpu,
            anisotropy=anisotropy,
            sam_params=sam_params,
            size_filter={},
        )
        log(f"Cellpose-SAM preview: {sample_name} / {cell_type} / t={timepoint}")
        _run_worker(
            job, python, device, force_cpu,
            on_event=lambda e: (
                log(f"  t={e['t']}: {e['n_objects']} objects in {e['seconds']}s")
                if e.get("event") == "progress" else None
            ),
        )
        return np.load(preview_out)
