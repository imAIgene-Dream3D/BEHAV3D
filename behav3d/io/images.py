import shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from tqdm import tqdm

from .formats.czi import load_czi, load_czi_metadata, load_elsizes_czi, get_czi_shape, load_czi_timepoint_czyx
from .formats.h5 import load_h5
from .formats.ims import load_ims, load_ims_metadata, load_ims_timepoint_czyx, get_ims_dimension_order, get_ims_shape
from .formats.liff import load_liff, get_liff_shape, load_liff_metadata, load_elsizes_liff, load_liff_timepoint_czyx, get_liff_dimension_order
from .formats.tiff import load_tiff, get_tiff_shape, get_tiff_dimension_order
from .formats.zarr import load_zarr, load_zarr_timepoint, append_to_zarr, save_as_zarr, write_zarr_parallel

_FORMATS_WITH_TP_LOADER = {".czi", ".lif", ".liff", ".ims", ".zarr"}

RAW_TARGET_AXIS_ORDER = "TCZYX"
LABEL_TARGET_AXIS_ORDER = "TZYX"


def get_filepath_stem(path):
    path=Path(path)
    if path.suffix == ".zip" or path.suffix ==".bz2":
        path = Path(path.stem)
    return(path.stem)

def _validate_axis_string(axis_order: str, *, allow: str = "TCZYX") -> str:
    if not isinstance(axis_order, str) or not axis_order:
        raise ValueError("axis_order must be a non-empty string")
    axis_order = axis_order.upper()
    bad = [ax for ax in axis_order if ax not in allow]
    if bad:
        raise ValueError(f"axis_order contains unsupported axes {bad}. Allowed axes: {list(allow)}")
    if len(set(axis_order)) != len(axis_order):
        raise ValueError(f"axis_order has duplicate axes: '{axis_order}'")
    return axis_order

def _expand_missing_axes(arr: np.ndarray, current_order: str, target_order: str) -> tuple[np.ndarray, str]:
    """
    Insert singleton dimensions for axes in target_order missing from current_order.
    Axes are inserted in the correct semantic position per target_order.
    """
    out = arr
    order = str(current_order)
    for ax in target_order:
        if ax not in order:
            insert_at = target_order.index(ax)
            out = np.expand_dims(out, axis=insert_at)
            order = order[:insert_at] + ax + order[insert_at:]
    return out, order

def _reorder_to_target(arr: np.ndarray, current_order: str, target_order: str) -> np.ndarray:
    if current_order == target_order:
        return arr
    perm = [current_order.index(ax) for ax in target_order]
    return np.transpose(arr, axes=perm)

def normalize_raw_to_tczyx(arr: np.ndarray, source_order: str) -> np.ndarray:
    """
    Normalize raw intensity arrays to 5D TCZYX by:
    - reordering existing axes by label
    - inserting missing T/C/Z as singleton axes in the correct position
    """
    source_order = _validate_axis_string(source_order)
    if len(source_order) != arr.ndim:
        raise ValueError(f"source_order '{source_order}' does not match array ndim={arr.ndim}")

    present = "".join(ax for ax in RAW_TARGET_AXIS_ORDER if ax in source_order)
    if present != source_order:
        arr = _reorder_to_target(arr, source_order, present)
    arr, present = _expand_missing_axes(arr, present, RAW_TARGET_AXIS_ORDER)
    return _reorder_to_target(arr, present, RAW_TARGET_AXIS_ORDER)

def normalize_label_to_tzyx(arr: np.ndarray, source_order: str) -> np.ndarray:
    """
    Normalize label arrays to TZYX by:
    - rejecting multi-channel labels (C>1)
    - squeezing C only when axis_order explicitly contains C and C==1
    - reordering existing axes by label
    - inserting missing T/Z as singletons in the correct position
    """
    source_order = _validate_axis_string(source_order)
    if len(source_order) != arr.ndim:
        raise ValueError(f"source_order '{source_order}' does not match array ndim={arr.ndim}")

    order = source_order
    out = arr

    if "C" in order:
        c_idx = order.index("C")
        c_size = int(out.shape[c_idx])
        if c_size == 1:
            out = np.squeeze(out, axis=c_idx)
            order = order.replace("C", "")
        else:
            raise ValueError(
                f"Label images must not have multiple channels. Got C={c_size} for axis_order='{source_order}'. "
                "Provide a single label volume per timepoint (no channel axis), e.g. axis_order='TZYX'."
            )

    present = "".join(ax for ax in LABEL_TARGET_AXIS_ORDER if ax in order)
    if present != order:
        out = _reorder_to_target(out, order, present)
    out, present = _expand_missing_axes(out, present, LABEL_TARGET_AXIS_ORDER)
    return _reorder_to_target(out, present, LABEL_TARGET_AXIS_ORDER)


def _cast_uint32_for_zarr(arr: np.ndarray) -> np.ndarray:
    """
    Zarr v3 does not support uint32. For labels and similar arrays, downcast to
    uint16 when values fit, otherwise int32 (no data loss).
    """
    if arr.dtype != np.uint32:
        return arr
    max_val = int(np.max(np.asarray(arr)))
    if max_val <= np.iinfo(np.uint16).max:
        print(f"  Casting uint32 → uint16 (max label value {max_val} fits in uint16)")
        return arr.astype(np.uint16)
    print(f"  Casting uint32 → int32 (max label value {max_val} exceeds uint16)")
    return arr.astype(np.int32)


def reorder_axes(img, source_order, target_order):
    """
    Reorder image axes from source_order to target_order via np.transpose.

    Both orders must be strings of the same 5 axis labels (e.g. "TZCYX" → "TCZYX").
    Returns the image unchanged if the orders already match.
    """
    if source_order == target_order:
        return img
    perm = [source_order.index(ax) for ax in target_order]
    return np.transpose(img, axes=perm)

def load_image(path, group=None, mode="r", *, warn_if_not_zarr=True, **kwargs):
    """
    Load an image from any supported format.
    Returns the array in its native/constructed axis order (no reordering).

    Parameters
    ----------
    warn_if_not_zarr : bool
        If True (default), raise an error when the path is not already a
        .zarr (intended for pipeline steps that expect converted data).
        Set False when loading a source file that is about to be converted
        (e.g. external segmentation/tracking import).
    """
    path = Path(path)
    if warn_if_not_zarr:
        _ensure_zarr(path, label="Image")
    if path.suffix==".czi":
        img = load_czi(path, **kwargs)
    elif path.suffix==".h5":
        img = load_h5(path, substruct=group)
    elif path.suffix==".ims":
        img = load_ims(path)
    elif path.suffix==".zarr" or str(path).endswith(".zarr.zip"):
        img = load_zarr(path, group=group, mode=mode)
    elif path.suffix==".tif" or path.suffix==".tiff":
        img = load_tiff(path)
    elif path.suffix==".lif" or path.suffix==".liff":
        img = load_liff(path)
    else:
        raise ValueError(f"Unsupported file format: {path.suffix}")
    return(img)

def load_image_timepoint(path, t):
    """Load a single timepoint as a numpy array, optimized for minimal memory in parallel workers."""
    path = Path(path)
    if path.suffix == ".zarr" or str(path).endswith(".zarr.zip"):
        return load_zarr_timepoint(path, t)
    return np.asarray(load_image(path)[t])

def load_image_metadata(path):
    path = Path(path)
    if path.suffix==".czi":
        metadata = load_czi_metadata(path)
    elif path.suffix==".h5":
        pass
    elif path.suffix==".ims":
        metadata = load_ims_metadata(path)
    elif path.suffix==".zarr" or str(path).endswith(".zarr.zip"):
        pass
    elif path.suffix==".tif" or path.suffix==".tiff":
        pass
    elif path.suffix==".lif" or path.suffix==".liff":
        metadata = load_liff_metadata(path)
    else:
        raise ValueError(f"Unsupported file format: {path.suffix}")
    return(metadata)

def _has_tp_loader(path):
    """Return True if the format supports reading a single timepoint from disk."""
    path = Path(path)
    if path.suffix in _FORMATS_WITH_TP_LOADER:
        return True
    if str(path).endswith(".zarr.zip"):
        return True
    return False

def _load_czyx_for_timepoint(path, t, source_order):
    """
    Load a single timepoint from a format that supports per-timepoint I/O.
    Returns a numpy array in the natural subset of CZYX axes present in the
    source (e.g. ZYX if no C axis, CZYX if both exist).

    Axes are reordered to standard hierarchy (C < Z < Y < X) but NO
    singleton dimensions are added for missing axes.

    Supported: CZI, LIF/LIFF, IMS, Zarr.
    TIFF and H5 are NOT supported here (they require a full load handled
    separately in convert_file_to_zarr).
    """
    path = Path(path)

    if path.suffix == ".czi":
        data = load_czi_timepoint_czyx(path, t, source_order)
        remaining = source_order.replace("T", "")

    elif path.suffix in (".lif", ".liff"):
        data = load_liff_timepoint_czyx(path, t)
        remaining = "CZYX"

    elif path.suffix == ".ims":
        data = load_ims_timepoint_czyx(path, t)
        remaining = "CZYX"

    elif path.suffix == ".zarr" or str(path).endswith(".zarr.zip"):
        import dask.array as da
        dask_arr = load_zarr(path)
        remaining = source_order

        # Natural subset: reorder to TCZYX hierarchy but only for axes that exist
        target_full = "".join(ax for ax in RAW_TARGET_AXIS_ORDER if ax in remaining)
        if remaining != target_full:
            axes = [remaining.index(ax) for ax in target_full]
            dask_arr = da.transpose(dask_arr, axes)
            remaining = target_full

        # Extract timepoint T if it exists
        if 'T' in remaining:
            t_idx = remaining.index('T')
            tp_slice = np.asarray(dask_arr.take(t, axis=t_idx))
            return tp_slice  # Natural subset minus T
        else:
            return np.asarray(dask_arr)

    else:
        raise ValueError(
            f"Format '{path.suffix}' does not support per-timepoint loading. "
            "Use the full-load path in convert_file_to_zarr instead."
        )

    # --- Centralized reorder to natural TCZYX subset ---
    # We do NOT add singletons here. We just ensure the axes are in the 
    # correct relative order: T < C < Z < Y < X.
    target_order = "".join(ax for ax in RAW_TARGET_AXIS_ORDER if ax in source_order)
    
    # Remove 'T' from target because this function loads ONE timepoint
    target_czyx = target_order.replace("T", "")
    
    if remaining != target_czyx:
        perm = [remaining.index(ax) for ax in target_czyx]
        data = np.transpose(data, perm)

    return data  # (Natural subset of CZYX)


def convert_file_to_zarr(
    path,
    outpath=None,
    axis_order=None,
    chunks=None,
    overwrite=False,
    t_start=None,
    t_end=None,
    n_workers=1,
    ):
    """
    Convert a single image file to .zarr in TCZYX axis order.

    For formats with per-timepoint loaders (CZI, LIF, IMS, Zarr) each worker
    reads one timepoint from disk and writes it — only one CZYX volume per
    worker is held in RAM.

    For TIFF / H5 the full image is loaded once, then timepoint slices are
    distributed to parallel workers for writing.

    Parameters
    ----------
    path : Path
        Source image path.
    outpath : Path
        Destination .zarr path.
    axis_order : str
        The axis order of the source image (e.g. "TZCYX").
        Must be explicitly provided — no default is assumed.
    chunks : tuple, optional
        Chunk shape for the zarr output.
    overwrite : bool
        If True, overwrite existing output.
    t_start, t_end : int, optional
        Inclusive timepoint range to convert. None means full range.
    n_workers : int
        Number of parallel workers for writing.
    """
    path = Path(path)
    outpath = Path(outpath)
    n_workers = max(1, int(n_workers or 1))

    if (
        ((path.suffix == ".zarr" or str(path).endswith(".zarr.zip")) and path.exists()) or
        (outpath.exists() and not overwrite)
    ):
        print("Skipping conversion to zarr, as file already exists")
        return

    if overwrite and outpath.exists():
        if outpath.is_dir():
            shutil.rmtree(outpath)
        else:
            outpath.unlink()

    print(f"Converting {path.name} → {outpath}")

    if axis_order is None:
        raise ValueError(
            f"axis_order must be specified for '{path.name}'.\n"
            "Cannot assume a default axis order for any format. "
            "Ensure dimension_order is set in the metadata (e.g. via the Preprocessing widget)."
        )
    source_order = _validate_axis_string(axis_order)

    source_shape = get_image_shape(path)
    print(f"  Source shape: {source_shape}, axis_order: {source_order}")

    if len(source_shape) != len(source_order):
        raise ValueError(
            f"axis_order '{source_order}' has {len(source_order)} axes but the "
            f"source image has {len(source_shape)} dimensions {source_shape}. "
            f"Provide an axis_order that matches the actual number of dimensions "
            f"(e.g., 'ZYX' for 3D, 'CZYX' for 4D, 'TCZYX' for 5D)."
        )

    # Raw contract: always write TCZYX (insert missing T/C/Z as singletons)
    if "Y" not in source_order or "X" not in source_order:
        raise ValueError(
            f"axis_order must include at least 'Y' and 'X' for image conversion. Got '{source_order}'."
        )
    target_order = RAW_TARGET_AXIS_ORDER
    dim_size = dict(zip(source_order, source_shape))
    target_shape = tuple(int(dim_size.get(ax, 1)) for ax in target_order)

    n_t_total = dim_size.get('T', 1)
    t_s = max(0, t_start) if t_start is not None else 0
    t_e = min(n_t_total - 1, t_end) if t_end is not None else n_t_total - 1
    if t_s > t_e:
        print(f"  Warning: invalid time range [{t_s}, {t_e}], using full range.")
        t_s, t_e = 0, n_t_total - 1
    n_t_out = t_e - t_s + 1

    # Final output shape in Zarr
    out_shape = list(target_shape)
    out_shape[target_order.index('T')] = n_t_out
    out_shape = tuple(out_shape)

    if chunks is None:
        chunks = [s for s in out_shape]
        chunks[target_order.index('T')] = 1
        chunks = tuple(chunks)

    outpath.parent.mkdir(parents=True, exist_ok=True)

    import os
    n_cpu = os.cpu_count() or 1
    cpu_cap = max(1, n_cpu // 2)
    n_workers_req = int(n_workers) if n_workers is not None else 1
    n_workers_eff = max(1, n_workers_req)

    # Effective workers capped by output length (if T exists) and CPU pool
    n_workers_eff = min(n_workers_eff, n_t_out, cpu_cap)

    print(f"  Parallel writing workers: requested={n_workers_req} cpu_cap={cpu_cap} using={n_workers_eff}")

    if _has_tp_loader(path):
        _convert_with_tp_loader(
            path, outpath, source_order, target_order, out_shape,
            chunks, t_s, t_e, n_t_out, n_workers_eff
        )
    else:
        _convert_full_load(
            path, outpath, source_order, target_order, out_shape,
            chunks, t_s, t_e, n_t_out, n_workers_eff
        )

    print(f"  Done: {outpath.name}")


def _convert_with_tp_loader(path, outpath, source_order, target_order, out_shape,
                            chunks, t_s, t_e, n_t_out, n_workers):
    """Formats that support per-timepoint loading from disk (CZI/LIF/IMS/Zarr)."""
    # Load first timepoint and normalize to CZYX (raw contract: C and/or Z may be missing => singleton)
    first = _load_czyx_for_timepoint(path, t_s, source_order)
    present_czyx = "".join(ax for ax in "CZYX" if ax in source_order)
    first, present_czyx = _expand_missing_axes(np.asarray(first), present_czyx, "CZYX")
    first = _reorder_to_target(first, present_czyx, "CZYX")
    dtype = first.dtype

    # Zarr v3 does not support uint32.  Probe the first timepoint and
    # choose a safe output dtype (same logic as _convert_full_load).
    if dtype == np.uint32:
        max_val = int(np.max(np.asarray(first)))
        if max_val <= np.iinfo(np.uint16).max:
            print(f"  Casting {dtype} → uint16 (max value {max_val} fits in uint16)")
            dtype = np.dtype(np.uint16)
        else:
            print(f"  Casting {dtype} → int32 (max value {max_val} fits in int32)")
            dtype = np.dtype(np.int32)
        first = first.astype(dtype)

    write_zarr_parallel(
        outpath=outpath, shape=out_shape, dtype=dtype,
        chunks=chunks, overwrite=True,
    )

    print(f"  Writing {n_t_out} timepoints as {target_order} {out_shape}, dtype={dtype} ...")

    all_items = list(enumerate(range(t_s, t_e + 1)))
    _first_cache = {0: np.asarray(first)}
    del first

    def _worker(args):
        t_out, t_in = args
        if t_out in _first_cache:
            data = _first_cache.pop(t_out)
        else:
            data = np.asarray(_load_czyx_for_timepoint(path, t_in, source_order))
            # Normalize each timepoint to CZYX (insert missing C/Z as singletons in correct position)
            present = "".join(ax for ax in "CZYX" if ax in source_order)
            data, present = _expand_missing_axes(np.asarray(data), present, "CZYX")
            data = _reorder_to_target(data, present, "CZYX")
        if data.dtype != dtype:
            data = data.astype(dtype)
        
        write_zarr_parallel(outpath=outpath, index=t_out, data=data)

    if n_workers > 1:
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            for _ in tqdm(executor.map(_worker, all_items),
                          total=n_t_out, desc="Converting timepoints"):
                pass
    else:
        for item in tqdm(all_items, desc="Converting timepoints"):
            _worker(item)


def _convert_full_load(path, outpath, source_order, target_order, out_shape,
                       chunks, t_s, t_e, n_t_out, n_workers):
    """TIFF / H5: load full image once, reorder, slice range, write tp-by-tp."""
    print(f"  Loading full image (format does not support per-timepoint reading)...")
    if path.suffix in (".tif", ".tiff"):
        img = load_tiff(path)
    elif path.suffix == ".h5":
        img = load_h5(path)
    else:
        raise ValueError(
            f"_convert_full_load does not support '{path.suffix}'. "
            "Only TIFF and H5 use the full-load path."
        )

    if source_order != target_order:
        # Check if we need to expand anything before transposing
        # (Only if source is missing some axis that IS in target_order)
        expanded_order = source_order
        for ax in target_order:
            if ax not in expanded_order:
                img = np.expand_dims(img, axis=0)
                expanded_order = ax + expanded_order
        
        img = reorder_axes(img, expanded_order, target_order)

    if 'T' in target_order:
        t_idx = target_order.index('T')
        # Simple slice if T is first, otherwise use np.take
        if t_idx == 0:
            img = img[t_s : t_e + 1]
        else:
            img = np.take(img, range(t_s, t_e + 1), axis=t_idx)
    dtype = img.dtype

    # Zarr v3 does not support uint32.  Downcast to uint16 when values
    # fit (typical for label/segmentation images), otherwise promote to
    # int32 so we never lose data.
    if dtype == np.uint32:
        max_val = int(np.max(np.asarray(img)))
        if max_val <= np.iinfo(np.uint16).max:
            print(f"  Casting {dtype} → uint16 (max value {max_val} fits in uint16)")
            img = img.astype(np.uint16)
        else:
            print(f"  Casting {dtype} → int32 (max value {max_val} fits in int32)")
            img = img.astype(np.int32)
        dtype = img.dtype

    write_zarr_parallel(
        outpath=outpath, shape=out_shape, dtype=dtype,
        chunks=chunks, overwrite=True,
    )

    print(f"  Writing {n_t_out} timepoints as {target_order} {out_shape}, dtype={dtype} ...")

    def _worker(t):
        if 'T' in target_order:
            t_idx = target_order.index('T')
            # Extract one T slice
            if t_idx == 0:
                data = np.asarray(img[t])
            else:
                data = np.asarray(np.take(img, t, axis=t_idx))
            write_zarr_parallel(outpath=outpath, index=t, data=data)
        else:
            write_zarr_parallel(outpath=outpath, data=np.asarray(img))

    if n_workers > 1:
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            for _ in tqdm(executor.map(_worker, range(n_t_out)),
                          total=n_t_out, desc="Writing timepoints"):
                pass
    else:
        for t in tqdm(range(n_t_out), desc="Writing timepoints"):
            _worker(t)


def _ensure_zarr(path, label="Image"):
    """Raise a clear error if a path does not point to a .zarr file."""
    path = Path(path)
    if path.suffix != ".zarr" and not str(path).endswith(".zarr.zip"):
        raise ValueError(
            f"{label} '{path.name}' is not in .zarr format ({path.suffix}).\n"
            "Please convert all images to Zarr in the Preprocessing step first."
        )
    return path


def convert_label_file_to_zarr(
    *,
    path,
    outpath,
    axis_order: str,
    overwrite: bool = False,
    chunks=None,
):
    """
    Convert a label image (segments/tracks/masks) to Zarr in TZYX order.

    Rules:
    - Never introduce a channel axis.
    - If axis_order explicitly contains C:
      - C==1 is squeezed/dropped.
      - C>1 raises (multi-channel labels unsupported).
    """
    path = Path(path)
    outpath = Path(outpath)

    # If source is already zarr and target points to the same location,
    # never attempt reconversion (protects against accidental self-overwrite).
    if (path.suffix == ".zarr" or str(path).endswith(".zarr.zip")) and path.exists():
        if path.resolve() == outpath.resolve():
            print("Skipping label conversion to zarr, source is already .zarr")
            return

    if outpath.exists() and not overwrite:
        print("Skipping label conversion to zarr, as file already exists")
        return

    if overwrite and outpath.exists():
        if outpath.is_dir():
            shutil.rmtree(outpath)
        else:
            outpath.unlink()

    axis_order = _validate_axis_string(axis_order)
    img = load_image(path, warn_if_not_zarr=False)
    img_np = np.asarray(img)
    if img_np.ndim != len(axis_order):
        raise ValueError(
            f"axis_order '{axis_order}' has {len(axis_order)} axes but the source image has "
            f"{img_np.ndim} dimensions {img_np.shape}. Provide the correct axis_order for this label image."
        )

    img_norm = normalize_label_to_tzyx(img_np, axis_order)
    img_norm = _cast_uint32_for_zarr(img_norm)
    if chunks is None:
        chunks = (1,) + tuple(img_norm.shape[1:])
    save_as_zarr(img_norm, outpath, chunks=chunks)


def convert_raw_file_to_zarr(
    *,
    path,
    outpath,
    axis_order: str,
    chunks=None,
    overwrite: bool = False,
    t_start=None,
    t_end=None,
    n_workers: int = 1,
):
    """
    Convert a raw intensity image to Zarr enforcing the TCZYX contract.
    This is a thin wrapper around convert_file_to_zarr().
    """
    return convert_file_to_zarr(
        path=path,
        outpath=outpath,
        axis_order=axis_order,
        chunks=chunks,
        overwrite=overwrite,
        t_start=t_start,
        t_end=t_end,
        n_workers=n_workers,
    )


def clip_zarr_timepoints(path, t_start=None, t_end=None):
    """
    Clip an existing TCZYX .zarr in place to the inclusive timepoint range
    [t_start, t_end], clamped to the array's own T extent.

    Writes the clipped data to a temporary sibling path first, then replaces
    the original, so a failed/interrupted clip never leaves a partially
    written zarr behind at the final path.
    """
    path = Path(path)
    img = load_zarr(path)
    n_t_total = img.shape[0]

    t_s = max(0, t_start) if t_start is not None else 0
    t_e = min(n_t_total - 1, t_end) if t_end is not None else n_t_total - 1
    if t_s > t_e:
        print(f"  Warning: invalid time range [{t_s}, {t_e}] for {path.name}, using full range.")
        t_s, t_e = 0, n_t_total - 1

    if t_s == 0 and t_e == n_t_total - 1:
        print(f"  Skipping clip for {path.name}: already within range [0, {n_t_total - 1}]")
        return

    clipped = img[t_s:t_e + 1]

    tmp_path = path.with_name(path.stem + "_clip_tmp" + path.suffix)
    if tmp_path.exists():
        if tmp_path.is_dir():
            shutil.rmtree(tmp_path)
        else:
            tmp_path.unlink()

    save_as_zarr(clipped, tmp_path)

    if path.is_dir():
        shutil.rmtree(path)
    elif path.exists():
        path.unlink()
    tmp_path.rename(path)

    print(f"  Clipped {path.name} to timepoints [{t_s}, {t_e}] ({t_e - t_s + 1} of {n_t_total})")


def load_elsizes(path):
    path = Path(path)
    if path.suffix==".czi":
        elsizes = load_elsizes_czi(path)
    elif path.suffix==".h5":
        pass
    elif path.suffix==".ims":
        pass
    elif path.suffix==".zarr" or str(path).endswith(".zarr.zip"):
        pass
    elif path.suffix==".tif" or path.suffix==".tiff":
        pass
    elif path.suffix==".lif" or path.suffix==".liff":
        elsizes = load_elsizes_liff(path)
    else:
        raise ValueError(f"Unsupported file format: {path.suffix}")
    
    return(elsizes)

def get_image_shape(path):
    """
    Get the shape of an image without fully loading it into memory.
    Supports: .zarr, .zarr.zip, .czi, .lif/.liff, .tif/.tiff, .ims, .h5
    """
    path = Path(path)
    if path.suffix == ".zarr" or str(path).endswith(".zarr.zip"):
        img = load_zarr(path)
        shape = img.shape
    elif path.suffix == ".czi":
        shape = get_czi_shape(path)
    elif path.suffix == ".lif" or path.suffix == ".liff":
        shape = get_liff_shape(path)
    elif path.suffix == ".tif" or path.suffix == ".tiff":
        shape = get_tiff_shape(path)
    elif path.suffix == ".ims":
        shape = get_ims_shape(path)
    elif path.suffix == ".h5":
        import h5py
        with h5py.File(str(path), "r") as f:
            keys = list(f.keys())
            if keys:
                shape = f[keys[0]].shape
            else:
                shape = ()
    else:
        raise ValueError(f"get_image_shape: unsupported format {path.suffix}")
    return tuple(shape)

def get_image_dimension_order(path):
    """
    Try to read the axis order directly from the **image file's** metadata.

    Only returns a value when the axis order is explicitly stored in the file
    (e.g. CZI dimension tags, OME-TIFF XML axes). For every other format the
    axis order is unknown and the user must select it manually.

    Returns
    -------
    str or None
        Axis order string (e.g. "TCZYX", "CZYX", "ZYX") read from the
        image file, or None if the format does not contain axis metadata.
    """
    path = Path(path)
    try:
        if path.suffix == ".czi":
            from aicspylibczi import CziFile
            czifile = CziFile(path)
            raw_dims = czifile.dims
            dim_order = "".join(ax for ax in raw_dims if ax in "TCZYX")
            if dim_order and len(dim_order) == len(set(dim_order)):
                return dim_order
            return None

        elif path.suffix == ".tif" or path.suffix == ".tiff":
            return get_tiff_dimension_order(path)  # None if not OME/ImageJ

        elif path.suffix in (".lif", ".liff"):
            return get_liff_dimension_order(path)

        elif path.suffix == ".ims":
            return get_ims_dimension_order(path)

        else:
            # H5, Zarr — no axis order metadata in the file.
            return None

    except Exception:
        return None


def get_czi_shape_and_dimension_order(path, take_dims="TCZYX"):
    """
    Like calling get_image_shape() and get_image_dimension_order() on a
    .czi path, but opens the CziFile only once instead of twice.
    """
    from aicspylibczi import CziFile

    path = Path(path)
    czifile = CziFile(path)
    raw_dims = czifile.dims

    dim_size = {ax: int(sz[1]) for ax, sz in czifile.get_dims_shape()[0].items()}
    shape = tuple(dim_size.get(ax, 1) for ax in raw_dims if ax in take_dims)

    dim_order = "".join(ax for ax in raw_dims if ax in "TCZYX")
    if not (dim_order and len(dim_order) == len(set(dim_order))):
        dim_order = None

    return shape, dim_order
