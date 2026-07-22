"""Cellpose-SAM inference worker.

This script runs inside the Cellpose-SAM sidecar environment (see
:mod:`behav3d.preprocessing.segmentation.cpsam_env`), where ``cellpose>=4``
shadows the ``cellpose==3.1.1.2`` that the rest of BEHAV3D depends on. It is
never imported by the main process; it is spawned as a subprocess by
:mod:`behav3d.preprocessing.segmentation.cellpose_sam_prediction`.

Protocol
--------
``python _cpsam_worker.py <job.json>``

The job file describes one sample. Machine-readable events are emitted as JSON
Lines on **stdout**; all human/library chatter goes to **stderr**. Events::

    {"event": "ready",    "device": "cuda:0", "model": "cpsam", ...}
    {"event": "progress", "t": 3, "n_done": 1, "total": 10, "n_objects": 42, ...}
    {"event": "log",      "msg": "..."}
    {"event": "error",    "message": "...", "traceback": "..."}
    {"event": "done",     "n_done": 10}

Masks are written directly into a **pre-allocated** output zarr opened ``mode="r+"``
so that a timepoint sub-range only touches its own frames (unlike the v3 path,
which deletes the whole array first and so wipes out-of-range timepoints).
"""
# ---------------------------------------------------------------------------
# Import order matters. torchvision bundles its own zlib.dll / libjpeg.dll and
# injects its directory into the Windows DLL search path when its C extension
# loads. If PIL is imported *after* torchvision, PIL's _imaging resolves against
# torchvision's zlib and dies with "DLL load failed ... not a valid Win32
# application". Importing PIL first pins the correct DLLs. cellpose imports
# torchvision via segment_anything, so this must precede it.
from PIL import Image  # noqa: F401  (import for side effect; do not move)
# ---------------------------------------------------------------------------

import argparse
import json
import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np

# Machine-readable events go on the real stdout; everything else (cellpose's
# logger, tqdm bars, warnings) is redirected to stderr so it cannot corrupt the
# JSON Lines stream.
_EVENTS = sys.stdout
sys.stdout = sys.stderr


def _emit(event: str, **payload) -> None:
    """Write one JSON Lines event to the real stdout."""
    _EVENTS.write(json.dumps({"event": event, **payload}) + "\n")
    _EVENTS.flush()


def _build_normalize(sam: dict) -> dict:
    """Translate BEHAV3D's flat params into cellpose's ``normalize`` dict.

    Key names mirror ``cellpose.transforms.normalize_img``. Only non-default
    filter keys are included so cellpose keeps its own defaults for the rest.
    """
    norm = {
        "percentile": [
            float(sam.get("norm_percentile_low", 1.0)),
            float(sam.get("norm_percentile_high", 99.0)),
        ],
        "norm3D": bool(sam.get("norm3D", True)),
    }
    for key in ("sharpen_radius", "smooth_radius", "tile_norm_blocksize", "tile_norm_smooth3D"):
        value = float(sam.get(key, 0) or 0)
        if value > 0:
            norm[key] = value
    if sam.get("invert"):
        norm["invert"] = True
    return norm


def _resolve_device(force_cpu: bool):
    """Pick the torch device.

    The parent already restricted visibility via ``CUDA_VISIBLE_DEVICES``, so a
    requested GPU is always index 0 here, and force-CPU arrives as an empty
    ``CUDA_VISIBLE_DEVICES`` (no visible devices at all).

    Returns ``(device, using_gpu, info)`` where *info* names the *physical* card.
    ``cuda:0`` is a visibility-local index, so on its own it cannot tell the user
    which GPU actually ran; the name and VRAM can.
    """
    import torch

    if force_cpu or not torch.cuda.is_available():
        return torch.device("cpu"), False, {"gpu_name": None}
    props = torch.cuda.get_device_properties(0)
    return torch.device("cuda:0"), True, {
        "gpu_name": props.name,
        "gpu_total_mb": round(props.total_memory / 2**20),
        # Which physical index the parent exposed as cuda:0 (unset => all visible).
        "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES", "<all>"),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Cellpose-SAM inference worker")
    parser.add_argument("job", help="Path to the JSON job specification")
    args = parser.parse_args(argv)

    job = json.loads(Path(args.job).read_text(encoding="utf-8"))

    mode = job.get("mode", "run")
    sam = job.get("sam", {}) or {}
    size_filter = job.get("size_filter", {}) or {}
    channels = list(job.get("channels") or [])
    timepoints = list(job.get("timepoints") or [])
    force_cpu = bool(job.get("force_cpu", False))

    n_threads = job.get("n_threads")
    cooldown_s = float(job.get("cooldown_s", 0) or 0)

    import zarr
    from cellpose import models

    # The *_NUM_THREADS env vars set by build_worker_env are read at import time by
    # MKL/OpenBLAS, but torch's own ATen pool has to be set through the API, and
    # only after torch is imported. Both halves are needed for the cap to hold.
    if n_threads:
        import torch

        torch.set_num_threads(max(1, int(n_threads)))

    # behav3d is normally importable via the parent env's editable install; fall
    # back to the repo root so the worker also runs from an uninstalled checkout.
    try:
        from behav3d.preprocessing.segmentation import segment_size_filter, segment_2d_filter
    except ImportError:  # pragma: no cover - source checkout fallback
        sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
        from behav3d.preprocessing.segmentation import segment_size_filter, segment_2d_filter

    device, using_gpu, device_info = _resolve_device(force_cpu)

    model = models.CellposeModel(
        gpu=using_gpu,
        pretrained_model=str(job.get("model_name", "cpsam")),
        device=device,
    )

    # Read the placement back off the loaded weights rather than trusting the
    # requested device: cellpose silently falls back to CPU if it cannot use the
    # GPU, and a run that is 30x slower for that reason should say so.
    try:
        device_info["weights_device"] = str(next(model.net.parameters()).device)
    except Exception:
        pass

    _emit(
        "ready",
        device=str(device),
        gpu=using_gpu,
        model=str(job.get("model_name", "cpsam")),
        cellpose_version=__import__("cellpose").version,
        **device_info,
    )

    raw = zarr.open(str(job["raw_zarr"]), mode="r")

    # Cellpose-SAM consumes at most 3 channels and is invariant to their order.
    if len(channels) > 3:
        _emit("log", msg=f"Cellpose-SAM uses at most 3 channels; truncating {channels} to {channels[:3]}.")
        channels = channels[:3]

    do_3D = bool(sam.get("do_3D", True))
    stitch_threshold = float(sam.get("stitch_threshold", 0.0) or 0.0)
    normalize = _build_normalize(sam)

    diameter = sam.get("diameter")
    diameter = float(diameter) if diameter not in (None, "", 0) else None
    niter = int(sam.get("niter", 0) or 0) or None

    # Cellpose's own eval() only treats the input as a real 3D/stitched volume
    # (and forwards z_axis to convert_image) when do_3D or stitch_threshold>0
    # - see models.py: `do_3D=(do_3D or stitch_threshold > 0)`. Outside that,
    # it requires z_axis=None and instead infers a "batch of independent 2D
    # planes" layout with Z on axis 0 and channels elsewhere - the opposite of
    # our channel-first (C, Z, Y, X) raw layout - so this is the one case
    # where the per-timepoint array needs transposing before eval().
    batched_2d = not (do_3D or stitch_threshold > 0)

    eval_kwargs = dict(
        batch_size=int(sam.get("batch_size", 8)),
        normalize=normalize,
        diameter=diameter,
        flow_threshold=float(sam.get("flow_threshold", 0.4)),
        cellprob_threshold=float(sam.get("cellprob_threshold", 0.0)),
        do_3D=do_3D,
        anisotropy=job.get("anisotropy") if do_3D else None,
        flow3D_smooth=float(sam.get("flow3D_smooth", 0) or 0),
        stitch_threshold=0.0 if do_3D else stitch_threshold,
        # We deliberately let cellpose keep everything and filter afterwards, so
        # the interactive size preview can explore thresholds without re-running
        # inference. max_size_fraction defaults to 0.4 upstream, which silently
        # deletes any object covering >40% of the frame (e.g. large organoids).
        min_size=int(sam.get("min_size", 1)),
        max_size_fraction=float(sam.get("max_size_fraction", 1.0)),
        niter=niter,
        channel_axis=1 if batched_2d else 0,
        z_axis=None if batched_2d else 1,
    )

    out = None
    if mode == "run":
        out = zarr.open(str(job["out_zarr"]), mode="r+")

    size_min = size_filter.get("size_min")
    size_max = size_filter.get("size_max")
    size_min = int(size_min) if size_min else None
    size_max = int(size_max) if size_max else None
    # Caller-controlled, not derived from do_3D. 2D+stitch is precisely the mode
    # that produces flat fragments (unstitched masks on the first/last Z slice),
    # so tying this to do_3D would disable it exactly where it is most needed.
    drop_2d = bool(job.get("drop_2d_segments", True))

    total = len(timepoints)
    preview_volume = None

    for n_done, t in enumerate(timepoints, start=1):
        started = time.time()

        img_t = np.asarray(raw[t])           # (C, Z, Y, X)
        img_t = img_t[channels]              # keep only the requested channels
        if batched_2d:
            img_t = np.moveaxis(img_t, 0, 1)  # (C, Z, Y, X) -> (Z, C, Y, X)

        mask, *_ = model.eval(img_t, **eval_kwargs)
        mask = np.asarray(mask)

        # Removes labels that are flat along Z, Y or X. Cheap way to kill the
        # single-slice fragments cellpose leaves at the top/bottom of a stack.
        # Turn off for genuinely thin objects that occupy one slice.
        if drop_2d:
            mask = segment_2d_filter(mask)

        if size_min or size_max:
            mask = segment_size_filter(mask, size_min=size_min, size_max=size_max)

        mask = mask.astype(np.uint16)
        n_objects = int(len(np.unique(mask)) - 1) if mask.any() else 0

        if mode == "run":
            out[t] = mask
        else:
            preview_volume = mask

        peak_mb = None
        if using_gpu:
            import torch

            peak_mb = round(torch.cuda.max_memory_allocated() / 2**20)
            # A batch run is many independent eval() calls in one process. Without
            # this the caching allocator keeps every frame's peak reserved, and a
            # later, larger frame OOMs on memory that is free but fragmented.
            torch.cuda.reset_peak_memory_stats()
            torch.cuda.empty_cache()

        _emit(
            "progress",
            t=int(t),
            n_done=n_done,
            total=total,
            n_objects=n_objects,
            seconds=round(time.time() - started, 2),
            peak_vram_mb=peak_mb,
        )

        # Idle between frames so the CPU/GPU voltage regulators (and, on a laptop,
        # the battery that covers the gap between draw and adapter rating) get a
        # recovery window. Emitted first so the parent records the frame as done
        # before we go quiet - a machine that dies during the cooldown must still
        # be able to resume from here.
        if cooldown_s and n_done < total:
            time.sleep(cooldown_s)

    if mode == "preview" and preview_volume is not None:
        np.save(str(job["preview_out"]), preview_volume)

    _emit("done", n_done=total)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:  # pragma: no cover - reported to the parent
        _emit("error", message=str(exc), traceback=traceback.format_exc())
        sys.exit(1)
