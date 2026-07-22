"""Sidecar environment manager for Cellpose-SAM (Cellpose v4).

BEHAV3D pins ``cellpose==3.1.1.2`` because the bundled fine-tuned models
(``cellpose_organoids__channel0-organoid``, ``cellpose_tcells__channel0-tcell``)
are v3 U-Net weights. Cellpose v4 is a SAM ViT-L and refuses to load them, so the
two versions cannot share an interpreter.

This module manages a *sidecar* virtualenv created with ``--system-site-packages``
from the active BEHAV3D environment. The sidecar's own ``site-packages`` precedes
the parent's on ``sys.path``, so ``cellpose>=4`` shadows ``cellpose==3.1.1.2``
inside the sidecar only, while torch/numpy/scipy/zarr/behav3d all resolve to the
parent environment. Installing with ``--no-deps`` keeps the sidecar at ~13 MB
instead of duplicating the parent's ~4 GB CUDA-enabled torch.

Cellpose-SAM itself runs in :mod:`behav3d.preprocessing.segmentation._cpsam_worker`,
which is executed by this environment's interpreter.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Callable, Optional

#: Directory name of the sidecar, created inside the parent environment so it is
#: implicitly keyed to that environment (a sidecar built against one env's torch
#: must never be reused by another) and removed when the parent env is deleted.
CPSAM_ENV_DIRNAME = "cpsam_env"

#: Minimum cellpose version providing Cellpose-SAM.
CPSAM_REQUIREMENT = "cellpose>=4.0"

#: cellpose>=4 needs ``segment_anything``; every other v4 dependency is already
#: satisfied by a standard BEHAV3D environment. Installed explicitly because we
#: use ``--no-deps`` to stop pip re-resolving (and re-downloading) torch.
CPSAM_EXTRA_REQUIREMENTS = ("segment_anything",)

_status_cache: Optional[dict] = None


def cpsam_env_dir(base_prefix: str | Path | None = None) -> Path:
    """Return the sidecar directory for *base_prefix* (default: the running env)."""
    return Path(base_prefix or sys.prefix) / CPSAM_ENV_DIRNAME


def _venv_python(env_dir: Path) -> Path:
    """Return the interpreter path inside *env_dir* for the current platform."""
    if os.name == "nt":
        return env_dir / "Scripts" / "python.exe"
    return env_dir / "bin" / "python"


def find_cpsam_python(config: dict | None = None) -> Optional[Path]:
    """Locate the sidecar interpreter.

    Resolution order: an explicit ``cellpose_sam.python_path`` in the BEHAV3D
    config, then the default sidecar inside the running environment. Returns
    ``None`` when neither exists.
    """
    if config:
        override = (config.get("cellpose_sam") or {}).get("python_path")
        if override:
            candidate = Path(str(override))
            if candidate.exists():
                return candidate

    candidate = _venv_python(cpsam_env_dir())
    return candidate if candidate.exists() else None


#: Environment variables read by the numeric libraries the worker pulls in, in the
#: order they matter: torch's own ATen pool, MKL (numpy/scipy on conda), OpenBLAS
#: (numpy on pip wheels), and numexpr. All of them default to "one thread per
#: logical CPU", so they must all be capped together or the uncapped one wins.
_THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)


def build_worker_env(
    device: str = "auto",
    force_cpu: bool = False,
    n_threads: int | None = None,
) -> dict:
    """Build the environment mapping for a sidecar subprocess.

    Handles two Windows-specific hazards discovered while validating the sidecar:

    1. **conda DLL directories.** A conda environment's own ``python.exe``
       resolves its DLLs relative to itself, but a venv's ``Scripts/python.exe``
       does not, so binary packages inherited from the parent (Pillow, pyclesperanto)
       fail to load with ``ImportError: DLL load failed``. The parent's DLL
       directories are prepended to ``PATH`` explicitly rather than relying on the
       caller having activated conda.
    2. **Device selection.** ``CUDA_VISIBLE_DEVICES`` is used instead of a torch
       device string so the child physically cannot see other GPUs; ``""`` hides
       every GPU, which is the force-CPU path.

    *n_threads* caps the CPU-side thread pools. Left unset, every library below
    spawns one thread per logical CPU, so the mask post-processing and zstd
    compression between GPU bursts saturate the whole package. On machines whose
    power delivery cannot sustain full CPU *and* full GPU at once that is enough to
    trip a hard shutdown, so being able to run narrower - slower, but to completion
    - matters more than peak throughput.
    """
    env = os.environ.copy()

    if n_threads is not None:
        n_threads = max(1, int(n_threads))
        for var in _THREAD_ENV_VARS:
            env[var] = str(n_threads)

    base = Path(sys.prefix)
    dll_dirs = [
        base,
        base / "Library" / "mingw-w64" / "bin",
        base / "Library" / "usr" / "bin",
        base / "Library" / "bin",
        base / "Scripts",
        base / "bin",
    ]
    existing = [str(d) for d in dll_dirs if d.is_dir()]
    if existing:
        env["PATH"] = os.pathsep.join(existing + [env.get("PATH", "")])

    if force_cpu or device == "cpu":
        # No visible CUDA device => torch falls back to CPU. This is the
        # equivalent of simply not passing --use-gpu to the backend.
        env["CUDA_VISIBLE_DEVICES"] = ""
    elif device and device.startswith("cuda:"):
        # Expose only the requested GPU; the worker then always uses cuda:0.
        env["CUDA_VISIBLE_DEVICES"] = device.split(":", 1)[1]
    # device == "auto" => inherit whatever the parent can see.

    return env


def cpsam_env_status(config: dict | None = None, refresh: bool = False) -> dict:
    """Probe the sidecar and return ``{available, python, cellpose_version, error}``.

    Result is cached for the session because it costs an interpreter start-up;
    pass ``refresh=True`` after creating the environment.
    """
    global _status_cache
    if _status_cache is not None and not refresh:
        return _status_cache

    python = find_cpsam_python(config)
    if python is None:
        _status_cache = {
            "available": False,
            "python": None,
            "cellpose_version": None,
            "error": "Cellpose-SAM environment not set up yet.",
        }
        return _status_cache

    probe = (
        "import json\n"
        "from PIL import Image  # must precede torchvision, see _cpsam_worker\n"
        "import cellpose, torch\n"
        "gpus = [torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count())]\n"
        "print(json.dumps({'cellpose': cellpose.version, 'torch': torch.__version__,\n"
        "                  'cuda': torch.cuda.is_available(), 'gpus': gpus}))\n"
    )
    try:
        proc = subprocess.run(
            [str(python), "-c", probe],
            capture_output=True,
            text=True,
            timeout=180,
            env=build_worker_env(),
        )
    except Exception as exc:  # pragma: no cover - environment specific
        _status_cache = {
            "available": False,
            "python": str(python),
            "cellpose_version": None,
            "error": f"Could not run the Cellpose-SAM interpreter: {exc}",
        }
        return _status_cache

    if proc.returncode != 0:
        _status_cache = {
            "available": False,
            "python": str(python),
            "cellpose_version": None,
            "error": (proc.stderr or "").strip()[-800:] or "Unknown import failure.",
        }
        return _status_cache

    try:
        info = json.loads((proc.stdout or "").strip().splitlines()[-1])
    except Exception:
        _status_cache = {
            "available": False,
            "python": str(python),
            "cellpose_version": None,
            "error": f"Unexpected probe output: {(proc.stdout or '')[:400]}",
        }
        return _status_cache

    version = str(info.get("cellpose", ""))
    if not version.startswith("4"):
        _status_cache = {
            "available": False,
            "python": str(python),
            "cellpose_version": version,
            "error": (
                f"Sidecar resolves cellpose {version}, but Cellpose-SAM needs 4.x. "
                "The sidecar is not shadowing the parent environment correctly."
            ),
        }
        return _status_cache

    _status_cache = {
        "available": True,
        "python": str(python),
        "cellpose_version": version,
        "torch_version": info.get("torch"),
        "cuda": bool(info.get("cuda")),
        "gpus": list(info.get("gpus") or []),
        "error": None,
    }
    return _status_cache


def create_cpsam_env(
    progress_cb: Optional[Callable[[str], None]] = None,
    base_prefix: str | Path | None = None,
) -> Path:
    """Create (or repair) the Cellpose-SAM sidecar. Idempotent.

    Returns the sidecar interpreter path. Raises ``RuntimeError`` on failure.
    """
    log = progress_cb or (lambda msg: None)
    env_dir = cpsam_env_dir(base_prefix)
    python = _venv_python(env_dir)

    if not python.exists():
        log(f"Creating Cellpose-SAM environment at {env_dir} ...")
        # --system-site-packages is what keeps this ~13 MB instead of ~4 GB: the
        # sidecar reuses the parent's torch/numpy/scipy/zarr and only adds cellpose 4.
        proc = subprocess.run(
            [sys.executable, "-m", "venv", "--system-site-packages", str(env_dir)],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0 or not python.exists():
            raise RuntimeError(
                f"Failed to create the Cellpose-SAM environment:\n{proc.stderr or proc.stdout}"
            )
    else:
        log(f"Reusing existing Cellpose-SAM environment at {env_dir}")

    log(f"Installing {CPSAM_REQUIREMENT} (this does not touch your cellpose 3 install) ...")
    # --no-deps is essential: every cellpose 4 dependency except segment_anything is
    # already provided by the parent env, and without it pip would try to re-resolve
    # torch and download several GB.
    proc = subprocess.run(
        [
            str(python), "-m", "pip", "install", "--no-deps",
            "--disable-pip-version-check",
            CPSAM_REQUIREMENT, *CPSAM_EXTRA_REQUIREMENTS,
        ],
        capture_output=True,
        text=True,
        env=build_worker_env(),
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"Failed to install cellpose>=4 into the sidecar:\n{proc.stderr or proc.stdout}"
        )

    status = cpsam_env_status(refresh=True)
    if not status["available"]:
        raise RuntimeError(
            f"Cellpose-SAM environment created but unusable: {status['error']}"
        )

    log(f"Cellpose-SAM ready (cellpose {status['cellpose_version']}, torch {status.get('torch_version')}).")
    return python
