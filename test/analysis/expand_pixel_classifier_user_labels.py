#!/usr/bin/env python3
from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import zarr

# Fill these in, then run:
# /Users/s.deblank-3/miniforge3/envs/behav3d/bin/python test/expand_pixel_classifier_user_labels.py
SRC_PATH = Path("/Users/s.deblank-3/Downloads/afasad/PixelClassifier_UserTcell2Labels.zarr")
OUTPUT_DIR = Path("/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/runs/NatureBriefComm/LowDensity_MultiColor/images/PixelClassification")
OLD_PER_SAMPLE = 3
NEW_PER_SAMPLE = 7
OVERWRITE = True


def _open_zarr_array(path: Path, mode: str = "r"):
    if path.suffix == ".zip":
        store = zarr.storage.ZipStore(path, mode=mode)
        return zarr.open(store, mode=mode)
    return zarr.open(str(path), mode=mode)


def _validated_mapping(old_per_sample: int, new_per_sample: int) -> np.ndarray:
    if old_per_sample < 2:
        raise ValueError("old_per_sample must be at least 2")
    if new_per_sample < old_per_sample:
        raise ValueError(
            f"new_per_sample ({new_per_sample}) must be >= old_per_sample ({old_per_sample})"
        )

    numerator = new_per_sample - 1
    denominator = old_per_sample - 1
    if numerator % denominator != 0:
        raise ValueError(
            "new_per_sample does not preserve the existing evenly spaced label indices exactly. "
            f"For old_per_sample={old_per_sample}, new_per_sample must satisfy "
            f"(new_per_sample - 1) % (old_per_sample - 1) == 0. "
            f"Got new_per_sample={new_per_sample}."
        )

    step = numerator // denominator
    mapping = np.arange(old_per_sample, dtype=int) * step
    if mapping[0] != 0 or mapping[-1] != new_per_sample - 1:
        raise ValueError(
            f"Computed mapping {mapping.tolist()} does not span the expected first/middle/last positions."
        )
    return mapping


def expand_user_labels(
    src_path: Path,
    dst_path: Path,
    old_per_sample: int = 3,
    new_per_sample: int = 5,
    overwrite: bool = False,
) -> None:
    mapping = _validated_mapping(old_per_sample, new_per_sample)

    src = _open_zarr_array(src_path, mode="r")
    if src.ndim < 1:
        raise ValueError(f"Expected at least 1 dimension, got shape {src.shape}")

    total_frames = int(src.shape[0])
    if total_frames % old_per_sample != 0:
        raise ValueError(
            f"Source frame count {total_frames} is not divisible by old_per_sample={old_per_sample}"
        )

    n_samples = total_frames // old_per_sample
    dst_shape = (n_samples * new_per_sample,) + tuple(int(v) for v in src.shape[1:])
    dst_chunks = (1,) + tuple(int(v) for v in src.chunks[1:]) if getattr(src, "chunks", None) else None

    if dst_path.exists():
        if not overwrite:
            raise FileExistsError(
                f"Destination already exists: {dst_path}. Use --overwrite to replace it."
            )
        if dst_path.is_dir():
            shutil.rmtree(dst_path)
        else:
            dst_path.unlink()

    if dst_path.suffix == ".zip":
        work_path = dst_path.with_suffix("")
    else:
        work_path = dst_path

    if work_path.exists():
        if work_path.is_dir():
            shutil.rmtree(work_path)
        else:
            work_path.unlink()

    dst = zarr.open(
        str(work_path),
        mode="w",
        shape=dst_shape,
        chunks=dst_chunks,
        dtype=src.dtype,
    )
    dst[:] = 0

    for sample_idx in range(n_samples):
        src_start = sample_idx * old_per_sample
        dst_start = sample_idx * new_per_sample
        for old_idx, new_idx in enumerate(mapping):
            dst[dst_start + int(new_idx)] = src[src_start + old_idx]

    try:
        dst.attrs.update(dict(src.attrs))
    except Exception:
        pass
    dst.attrs["expanded_from"] = str(src_path)
    dst.attrs["old_per_sample"] = int(old_per_sample)
    dst.attrs["new_per_sample"] = int(new_per_sample)
    dst.attrs["sample_frame_mapping"] = mapping.tolist()

    if dst_path.suffix == ".zip":
        shutil.make_archive(str(work_path), "zip", str(work_path))
        shutil.rmtree(work_path)

    print(f"Source: {src_path}")
    print(f"Destination: {dst_path}")
    print(f"Input shape: {src.shape}")
    print(f"Output shape: {dst_shape}")
    print(f"Samples detected: {n_samples}")
    print(f"Frame mapping per sample: {list(range(old_per_sample))} -> {mapping.tolist()}")


def main() -> None:
    src_path = SRC_PATH.expanduser().resolve()
    if not src_path.exists():
        raise FileNotFoundError(f"Source path does not exist: {src_path}")

    output_dir = OUTPUT_DIR.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    dst_path = output_dir / src_path.name

    expand_user_labels(
        src_path=src_path,
        dst_path=dst_path,
        old_per_sample=OLD_PER_SAMPLE,
        new_per_sample=NEW_PER_SAMPLE,
        overwrite=OVERWRITE,
    )


if __name__ == "__main__":
    main()
