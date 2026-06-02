from concurrent.futures import ThreadPoolExecutor, as_completed
import shutil
from pathlib import Path

import dask.array as da
import numpy as np
import pandas as pd
from scipy import ndimage
from tqdm import tqdm

from behav3d.core.metadata import (
    detect_immune_cell_types_from_metadata,
    detect_organoid_types_from_metadata,
    detect_other_cell_types_from_metadata,
    expand_multicolor_celltype_names,
    load_behav3d_metadata,
)
from behav3d.io.images import load_zarr, write_zarr_parallel
from behav3d.preprocessing import erode_mask


def _normalize_segmentation_labels_for_view(label_arr, raw_img_shape, layer_name="labels"):
    """
    Normalize a label array to the raw-image viewer shape (T, Z, Y, X).
    """
    raw_shape = tuple(int(v) for v in raw_img_shape)
    label_shape = tuple(int(v) for v in label_arr.shape)

    if len(raw_shape) != 5:
        raise ValueError(f"{layer_name}: expected raw image shape (T,C,Z,Y,X), got {raw_shape}.")

    target_shape = (raw_shape[0],) + tuple(raw_shape[-3:])

    if label_shape == target_shape:
        return label_arr

    if len(label_shape) == 5:
        if label_shape[0] != target_shape[0] or tuple(label_shape[-3:]) != tuple(target_shape[1:]):
            raise ValueError(
                f"{layer_name}: expected label shape {target_shape} or "
                f"{target_shape[:1] + (1,) + target_shape[1:]}, got {label_shape}."
            )
        if int(label_shape[1]) != 1:
            raise ValueError(
                f"{layer_name}: cannot align multi-channel labels with shape {label_shape} "
                f"to viewer shape {target_shape}."
            )
        return label_arr[:, 0, ...]

    if len(label_shape) != 4:
        raise ValueError(
            f"{layer_name}: expected 4D/5D label array compatible with viewer shape {target_shape}, "
            f"got {label_shape}."
        )

    raise ValueError(f"{layer_name}: expected label shape {target_shape}, got {label_shape}.")


def _build_lazy_overlap_mask(segment_paths):
    occupancies = [load_zarr(path) > 0 for path in segment_paths]
    if not occupancies:
        raise ValueError("segment_paths must contain at least one zarr path.")
    overlap_degree = da.stack(occupancies, axis=0).sum(axis=0)
    return (overlap_degree >= 2).astype(np.uint8)


def _load_segment_array(path):
    path = Path(path)
    try:
        return load_zarr(path)
    except Exception as exc:
        raise ValueError(f"Could not load segmentation zarr '{path}': {exc}") from exc


def _validate_segment_inputs(segment_paths):
    if not segment_paths:
        raise ValueError("segment_paths must contain at least one zarr path.")

    arrays = []
    shape = None
    for path in segment_paths:
        arr = _load_segment_array(path)
        if arr.ndim != 4:
            raise ValueError(
                f"Segmentation input must be TZYX (4D). Got shape {tuple(arr.shape)} for '{path}'."
            )
        if shape is None:
            shape = tuple(int(v) for v in arr.shape)
        elif tuple(int(v) for v in arr.shape) != shape:
            raise ValueError(
                f"All segmentation inputs must share the same TZYX shape. "
                f"Expected {shape}, got {tuple(arr.shape)} for '{path}'."
            )
        arrays.append(arr)
    return arrays, shape


def _default_output_path(path):
    path = Path(path)
    if str(path).endswith(".zarr.zip"):
        stem = path.name[:-9]
        return path.with_name(f"{stem}_multicolor_processed.zarr")
    if path.suffix == ".zarr":
        return path.with_name(f"{path.stem}_multicolor_processed.zarr")
    return path.with_name(f"{path.name}_multicolor_processed.zarr")


def _infer_segment_name(path):
    path = Path(path)
    stem = path.stem
    if stem.endswith(".zarr"):
        stem = Path(stem).stem
    for suffix in ("_multicolor_processed", "_segments"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
    if "_" in stem:
        return stem.split("_")[-1]
    return stem




def _resolve_output_paths(segment_paths, output_paths, overwrite):
    segment_paths = [Path(p) for p in segment_paths]
    if output_paths is not None:
        if len(output_paths) != len(segment_paths):
            raise ValueError(
                f"output_paths must match segment_paths length ({len(segment_paths)}), got {len(output_paths)}."
            )
        return [Path(p) for p in output_paths]
    if overwrite:
        return segment_paths
    return [_default_output_path(path) for path in segment_paths]


def _remove_existing_output(path):
    path = Path(path)
    if not path.exists():
        return
    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()


def _initialize_output_paths(input_arrays, input_paths, output_paths, overwrite):
    output_paths = [Path(path) for path in output_paths]
    for input_arr, input_path, output_path in zip(input_arrays, input_paths, output_paths):
        output_path = Path(output_path)
        same_path = output_path.resolve() == Path(input_path).resolve()

        if same_path:
            if output_path.suffix != ".zarr":
                raise ValueError(
                    f"In-place overwrite requires a directory-backed .zarr path, got '{output_path}'."
                )
            continue

        if output_path.suffix != ".zarr":
            raise ValueError(
                f"Output path '{output_path}' must be a directory-backed .zarr to reuse BEHAV3D zarr writers."
            )
        if output_path.exists():
            if not overwrite:
                raise FileExistsError(
                    f"Output already exists: {output_path}. Use overwrite=True to replace it."
                )
            _remove_existing_output(output_path)

        output_path.parent.mkdir(parents=True, exist_ok=True)
        raw_chunks = getattr(input_arr, "chunks", None)
        if raw_chunks:
            chunks = tuple(int(axis_chunks[0]) for axis_chunks in raw_chunks)
        else:
            chunks = (1,) + tuple(int(v) for v in input_arr.shape[1:])
        write_zarr_parallel(
            outpath=output_path,
            shape=tuple(int(v) for v in input_arr.shape),
            dtype=input_arr.dtype,
            chunks=chunks,
            overwrite=True,
        )
    return output_paths


def _node_key(node):
    return int(node[0]), int(node[1])


def _build_conflict_groups(volumes):
    parent = {}

    def find(node):
        parent.setdefault(node, node)
        if parent[node] != node:
            parent[node] = find(parent[node])
        return parent[node]

    def union(a, b):
        ra = find(a)
        rb = find(b)
        if ra != rb:
            parent[rb] = ra

    n_images = len(volumes)
    for idx_a in range(n_images):
        for idx_b in range(idx_a + 1, n_images):
            overlap_mask = (volumes[idx_a] > 0) & (volumes[idx_b] > 0)
            if not np.any(overlap_mask):
                continue
            pairs = np.stack(
                [volumes[idx_a][overlap_mask], volumes[idx_b][overlap_mask]],
                axis=1,
            )
            for label_a, label_b in np.unique(pairs, axis=0):
                node_a = (idx_a, int(label_a))
                node_b = (idx_b, int(label_b))
                union(node_a, node_b)

    groups = {}
    for node in parent:
        root = find(node)
        groups.setdefault(root, set()).add(node)
    return [sorted(group, key=_node_key) for group in groups.values()]


def _mask_is_subset(mask_a, mask_b):
    return not np.any(mask_a & ~mask_b)


def _drop_contained_nodes(node_masks):
    survivors = dict(node_masks)
    removed = True
    while removed:
        removed = False
        for node in sorted(survivors, key=_node_key):
            mask = survivors[node]
            for other in sorted(survivors, key=_node_key):
                if node == other:
                    continue
                other_mask = survivors[other]
                if not _mask_is_subset(mask, other_mask):
                    continue
                if np.array_equal(mask, other_mask) and _node_key(node) < _node_key(other):
                    continue
                del survivors[node]
                removed = True
                break
            if removed:
                break
    return survivors


def _erode_candidate(mask, erosion_pixels):
    if erosion_pixels is None or int(erosion_pixels) <= 0 or not np.any(mask):
        return mask.astype(bool, copy=False)

    eroded = erode_mask(mask.astype(bool), use_dimensions=3, nr_pixels=int(erosion_pixels))
    return np.asarray(eroded, dtype=bool)


def _partition_union_by_nearest_seeds(candidate_masks, union_mask):
    seed_codes = np.zeros(union_mask.shape, dtype=np.int32)
    code_to_node = {}
    for code, node in enumerate(sorted(candidate_masks, key=_node_key), start=1):
        seed_codes[candidate_masks[node]] = code
        code_to_node[code] = node

    if not np.any(seed_codes):
        return {}

    indices = ndimage.distance_transform_edt(
        seed_codes == 0,
        return_distances=False,
        return_indices=True,
    )
    assigned_codes = seed_codes[tuple(indices)]

    assignments = {}
    for code, node in code_to_node.items():
        node_mask = union_mask & (assigned_codes == code)
        if np.any(node_mask):
            assignments[node] = node_mask
    return assignments


def _apply_group_resolution(volumes, processed_volumes, group_nodes, erosion_pixels, min_size, max_size):
    node_masks = {
        node: volumes[node[0]] == node[1]
        for node in group_nodes
    }
    if not node_masks:
        return

    full_union = np.zeros(volumes[0].shape, dtype=bool)
    for mask in node_masks.values():
        full_union |= mask

    # Keep exactly one copy of shared voxels per conflict group.
    # This removes redundant overlap but does not merge the labels into a
    # combined object or transfer union voxels between channels.
    overlap_degree = np.zeros(full_union.shape, dtype=np.uint8)
    for mask in node_masks.values():
        overlap_degree += mask.astype(np.uint8)
    shared_voxels = overlap_degree >= 2

    if not np.any(shared_voxels):
        assignments = dict(node_masks)
    else:
        sizes = {node: int(mask.sum()) for node, mask in node_masks.items()}
        winner = max(sizes, key=lambda node: (sizes[node], -node[0], -node[1]))

        assignments = {winner: node_masks[winner]}
        for node, mask in node_masks.items():
            if node == winner:
                continue
            cleaned = mask & ~shared_voxels
            cleaned = _erode_candidate(cleaned, erosion_pixels)
            if np.any(cleaned):
                assignments[node] = cleaned

    for image_idx in sorted({node[0] for node in group_nodes}):
        processed_volumes[image_idx][full_union] = 0

    for node, mask in assignments.items():
        processed_volumes[node[0]][mask] = node[1]


def process_multicolor_segments(
    segment_paths,
    output_paths=None,
    overwrite=False,
    erosion_pixels=1,
    min_size=None,
    max_size=None,
    show_progress=True,
    n_workers=1,
):
    """
    Remove redundant overlap across multiple labeled segmentation volumes.

    The cleanup keeps per-channel segmentations separate and only removes the
    shared voxels that would otherwise be duplicated across multicolor channels.
    It does not merge the labels into a single combined segmentation.
    """
    from behav3d.preprocessing.segmentation import segment_size_filter

    segment_paths = [Path(path) for path in segment_paths]
    input_arrays, _shape = _validate_segment_inputs(segment_paths)
    resolved_output_paths = _resolve_output_paths(segment_paths, output_paths, overwrite)
    resolved_output_paths = _initialize_output_paths(
        input_arrays,
        segment_paths,
        resolved_output_paths,
        overwrite,
    )

    n_timepoints = int(input_arrays[0].shape[0])
    n_workers = max(1, int(n_workers or 1))

    def _process_timepoint(t):
        volumes = [np.asarray(arr[t]) for arr in input_arrays]
        processed_volumes = [vol.copy() for vol in volumes]
        for group_nodes in _build_conflict_groups(volumes):
            _apply_group_resolution(
                volumes=volumes,
                processed_volumes=processed_volumes,
                group_nodes=group_nodes,
                erosion_pixels=erosion_pixels,
                min_size=min_size,
                max_size=max_size,
            )
        for idx, processed in enumerate(processed_volumes):
            if min_size is not None or max_size is not None:
                processed_volumes[idx] = segment_size_filter(
                    processed,
                    size_min=min_size,
                    size_max=max_size,
                )
        for input_arr, output_path, processed in zip(input_arrays, resolved_output_paths, processed_volumes):
            write_zarr_parallel(
                outpath=output_path,
                index=t,
                data=processed.astype(input_arr.dtype, copy=False),
            )
        return t

    if n_workers == 1 or n_timepoints <= 1:
        timepoints = tqdm(
            range(n_timepoints),
            desc="Resolving multicolor overlaps",
            unit="tp",
            dynamic_ncols=True,
            disable=(not show_progress) or n_timepoints <= 1,
        )
        for t in timepoints:
            _process_timepoint(t)
    else:
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            futures = {
                executor.submit(_process_timepoint, t): t
                for t in range(n_timepoints)
            }
            progress = tqdm(
                as_completed(futures),
                total=len(futures),
                desc="Resolving multicolor overlaps",
                unit="tp",
                dynamic_ncols=True,
                disable=not show_progress,
            )
            for future in progress:
                try:
                    future.result()
                except Exception:
                    for pending in futures:
                        pending.cancel()
                    raise

    return resolved_output_paths


def calculate_multicolor_overlap(segment_paths, show_progress=True, n_workers=1):
    """
    Calculate pairwise and global overlap counts for multiple segmentation images.
    """
    segment_paths = [Path(path) for path in segment_paths]
    input_arrays, shape = _validate_segment_inputs(segment_paths)
    n_images = len(input_arrays)
    n_timepoints = int(shape[0])
    segment_names = [_infer_segment_name(path) for path in segment_paths]

    pairwise_overlap_by_timepoint = np.zeros((n_timepoints, n_images, n_images), dtype=np.int64)
    pairwise_non_overlap_by_timepoint = np.zeros((n_timepoints, n_images, n_images), dtype=np.int64)
    overlap_any_by_timepoint = np.zeros(n_timepoints, dtype=np.int64)
    overlap_all_by_timepoint = np.zeros(n_timepoints, dtype=np.int64)
    overlap_mask = np.zeros(shape, dtype=np.uint8)
    overlap_degree = np.zeros(shape, dtype=np.uint8)
    n_workers = max(1, int(n_workers or 1))

    def _calculate_timepoint(t):
        occupancies = [np.asarray(arr[t]) > 0 for arr in input_arrays]
        occupancy_count = np.zeros(shape[1:], dtype=np.uint8)
        for occ in occupancies:
            occupancy_count += occ.astype(np.uint8)

        t_overlap_degree = occupancy_count
        t_overlap_mask = (occupancy_count >= 2).astype(np.uint8)
        t_overlap_any = int(np.count_nonzero(t_overlap_mask))
        t_overlap_all = (
            int(np.count_nonzero(occupancy_count == n_images))
            if n_images >= 2
            else 0
        )
        t_pairwise_overlap = np.zeros((n_images, n_images), dtype=np.int64)
        t_pairwise_non_overlap = np.zeros((n_images, n_images), dtype=np.int64)

        for idx_a in range(n_images):
            for idx_b in range(idx_a + 1, n_images):
                overlap_count = int(np.count_nonzero(occupancies[idx_a] & occupancies[idx_b]))
                non_overlap_count = int(np.count_nonzero(occupancies[idx_a] ^ occupancies[idx_b]))
                t_pairwise_overlap[idx_a, idx_b] = overlap_count
                t_pairwise_overlap[idx_b, idx_a] = overlap_count
                t_pairwise_non_overlap[idx_a, idx_b] = non_overlap_count
                t_pairwise_non_overlap[idx_b, idx_a] = non_overlap_count

        return (
            t,
            t_overlap_degree,
            t_overlap_mask,
            t_overlap_any,
            t_overlap_all,
            t_pairwise_overlap,
            t_pairwise_non_overlap,
        )

    if n_workers == 1 or n_timepoints <= 1:
        timepoints = tqdm(
            range(n_timepoints),
            desc="Calculating multicolor overlap",
            unit="tp",
            dynamic_ncols=True,
            disable=(not show_progress) or n_timepoints <= 1,
        )
        results_iter = (_calculate_timepoint(t) for t in timepoints)
    else:
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            futures = {
                executor.submit(_calculate_timepoint, t): t
                for t in range(n_timepoints)
            }
            progress = tqdm(
                as_completed(futures),
                total=len(futures),
                desc="Calculating multicolor overlap",
                unit="tp",
                dynamic_ncols=True,
                disable=not show_progress,
            )
            results = []
            for future in progress:
                try:
                    results.append(future.result())
                except Exception:
                    for pending in futures:
                        pending.cancel()
                    raise
            results_iter = sorted(results, key=lambda item: item[0])

    for (
        t,
        t_overlap_degree,
        t_overlap_mask,
        t_overlap_any,
        t_overlap_all,
        t_pairwise_overlap,
        t_pairwise_non_overlap,
    ) in results_iter:
        overlap_degree[t] = t_overlap_degree
        overlap_mask[t] = t_overlap_mask
        overlap_any_by_timepoint[t] = t_overlap_any
        overlap_all_by_timepoint[t] = t_overlap_all
        pairwise_overlap_by_timepoint[t] = t_pairwise_overlap
        pairwise_non_overlap_by_timepoint[t] = t_pairwise_non_overlap

    pairwise_overlap_matrix = pairwise_overlap_by_timepoint.sum(axis=0)
    pairwise_non_overlap_matrix = pairwise_non_overlap_by_timepoint.sum(axis=0)
    pairwise_named_stats = []
    for idx_a in range(n_images):
        for idx_b in range(idx_a + 1, n_images):
            pairwise_named_stats.append(
                {
                    "name": f"{segment_names[idx_a]}-{segment_names[idx_b]}",
                    "image_a": segment_names[idx_a],
                    "image_b": segment_names[idx_b],
                    "path_a": segment_paths[idx_a],
                    "path_b": segment_paths[idx_b],
                    "overlapping_count": int(pairwise_overlap_matrix[idx_a, idx_b]),
                    "non_overlapping_count": int(pairwise_non_overlap_matrix[idx_a, idx_b]),
                }
            )

    return {
        "segment_paths": segment_paths,
        "segment_names": segment_names,
        "pairwise_overlap_matrix": pairwise_overlap_matrix,
        "pairwise_overlap_by_timepoint": pairwise_overlap_by_timepoint,
        "pairwise_non_overlap_matrix": pairwise_non_overlap_matrix,
        "pairwise_non_overlap_by_timepoint": pairwise_non_overlap_by_timepoint,
        "pairwise_named_stats": pairwise_named_stats,
        "overlap_any_count": int(overlap_any_by_timepoint.sum()),
        "overlap_all_count": int(overlap_all_by_timepoint.sum()),
        "overlap_any_by_timepoint": overlap_any_by_timepoint,
        "overlap_all_by_timepoint": overlap_all_by_timepoint,
        "overlap_mask": overlap_mask,
        "overlap_degree": overlap_degree,
    }


def _coerce_metadata(metadata):
    if metadata is None:
        return None
    if isinstance(metadata, pd.DataFrame):
        return metadata
    return load_behav3d_metadata(metadata)



def apply_multicolor_segment_correction_for_base(
    metadata,
    output_dir,
    base_cell_type,
    n_channels,
    overwrite=False,
    show_progress=True,
    n_workers=1,
):
    """Convenience wrapper for a numbered multicolor cell type family.

    This is a top-level helper so other modules (e.g., the widgets package)
    can import it directly.
    """
    metadata_df = _coerce_metadata(metadata)
    cell_types = expand_multicolor_celltype_names(base_cell_type, n_channels)
    corrected_paths_by_sample = apply_multicolor_segment_correction(
        metadata=metadata,
        output_dir=output_dir,
        cell_types=cell_types,
        overwrite=overwrite,
        n_workers=n_workers,
    )
    
    # Update metadata DataFrame with corrected segment paths
    if metadata_df is not None and corrected_paths_by_sample:
        metadata_df = _update_metadata_with_corrected_paths(
            metadata_df,
            corrected_paths_by_sample,
            cell_types,
        )
    
    return metadata_df if metadata_df is not None else metadata


def _update_metadata_with_corrected_paths(metadata_df, corrected_paths_by_sample, cell_types):
    """Update metadata DataFrame columns with corrected segment paths.
    
    Parameters
    ----------
    metadata_df : pd.DataFrame
        The metadata DataFrame to update.
    corrected_paths_by_sample : dict
        Dict mapping sample_name -> list of corrected paths (one per cell_type, in order).
    cell_types : list
        List of cell types in the same order as the corrected paths.
    
    Returns
    -------
    pd.DataFrame
        Updated metadata DataFrame with corrected segment paths.
    """
    metadata_df = metadata_df.copy()
    sample_name_col = metadata_df["sample_name"].astype("string").str.strip()
    
    for sample_name, corrected_paths in corrected_paths_by_sample.items():
        sample_rows = sample_name_col == sample_name
        if not sample_rows.any():
            continue
        
        # Update each cell type's segment path column
        for cell_type, corrected_path in zip(cell_types, corrected_paths):
            # Find the appropriate column for this cell type
            # Column names follow pattern: <prefix>_<cell_type>_segments_image_path
            for col in metadata_df.columns:
                if col.endswith("_segments_image_path"):
                    col_cell_type = col[col.find("_") + 1 : -len("_segments_image_path")]
                    if col_cell_type == cell_type:
                        metadata_df.loc[sample_rows, col] = str(corrected_path)
                        break
    
    return metadata_df

def _resolve_output_dir(metadata, output_dir=None):
    if output_dir is not None:
        return Path(output_dir).expanduser()
    if isinstance(metadata, (str, Path)):
        return Path(metadata).expanduser().parent
    raise ValueError("output_dir is required when metadata is provided as a DataFrame.")


def _normalize_cell_types(cell_types):
    if isinstance(cell_types, str):
        values = [part.strip() for part in cell_types.split(",")]
    else:
        values = [str(part).strip() for part in cell_types]
    values = [value for value in values if value]
    if not values:
        raise ValueError("cell_types must contain at least one cell type.")
    return values


def _available_cell_types_from_metadata(metadata):
    metadata = _coerce_metadata(metadata)
    if metadata is None:
        return set()

    cell_types = set(detect_organoid_types_from_metadata(metadata))
    cell_types.update(detect_immune_cell_types_from_metadata(metadata))
    cell_types.update(detect_other_cell_types_from_metadata(metadata))

    for col in metadata.columns:
        if col.endswith("_segments_image_path"):
            if col.startswith(("or_", "im_", "ot_")):
                prefix = col[:2]
                cell_type = col[len(prefix) + 1 : -len("_segments_image_path")]
                if cell_type:
                    cell_types.add(cell_type)
            else:
                cell_types.add(col[: -len("_segments_image_path")])
    return cell_types


def _resolve_segment_column(sample, cell_type):
    for prefix in ("or", "im", "ot"):
        col = f"{prefix}_{cell_type}_segments_image_path"
        if col in sample.index:
            return col
    legacy_col = f"{cell_type}_segments_image_path"
    if legacy_col in sample.index:
        return legacy_col
    return None


def _default_segment_path(output_dir, sample_name, cell_type):
    return Path(output_dir, "images", sample_name, f"{sample_name}_{cell_type}_segments.zarr")


def _resolve_sample_segment_path(sample, output_dir, cell_type):
    sample_name = str(sample["sample_name"]).strip()
    default_path = _default_segment_path(output_dir, sample_name, cell_type)
    if default_path.exists():
        return default_path

    col = _resolve_segment_column(sample, cell_type)
    if col is not None:
        value = sample.get(col)
        if pd.notna(value) and str(value).strip():
            return Path(str(value)).expanduser()

    return default_path


def _build_sample_correction_plan(sample, output_dir, cell_types):
    sample_name = str(sample["sample_name"]).strip()
    segment_items = []
    missing_items = []

    for cell_type in cell_types:
        output_path = _resolve_sample_segment_path(sample, output_dir, cell_type)

        if output_path.exists():
            segment_items.append(
                {
                    "cell_type": cell_type,
                    "output_path": output_path,
                }
            )
        else:
            missing_items.append(
                {
                    "sample_name": sample_name,
                    "cell_type": cell_type,
                    "expected_path": output_path,
                }
            )

    return {
        "sample_name": sample_name,
        "segments": segment_items,
        "missing": missing_items,
    }


def _plan_multicolor_segment_correction(metadata, output_dir, cell_types, overwrite):
    sample_names = (
        metadata["sample_name"]
        .dropna()
        .astype("string")
        .str.strip()
        .replace("", pd.NA)
        .dropna()
        .unique()
        .tolist()
    )
    if not sample_names:
        raise ValueError("No valid sample_name values were found in metadata.")

    sample_plans = []
    missing_items = []
    sample_name_col = metadata["sample_name"].astype("string").str.strip()

    for sample_name in sample_names:
        sample_rows = metadata[sample_name_col == sample_name]
        sample = sample_rows.iloc[0]
        sample_plan = _build_sample_correction_plan(sample, output_dir, cell_types)
        sample_plans.append(sample_plan)
        missing_items.extend(sample_plan["missing"])

    if missing_items:
        details = "\n".join(
            f"sample='{item['sample_name']}', cell_type='{item['cell_type']}', expected_path='{item['expected_path']}'"
            for item in missing_items
        )
        print(
            "⚠️ Missing required segmentations for multicolor correction. "
            "Those samples will be skipped:\n"
            f"{details}"
        )

    blocked_samples = {item["sample_name"] for item in missing_items}
    if blocked_samples:
        sample_plans = [
            plan for plan in sample_plans
            if plan["sample_name"] not in blocked_samples
        ]

    return sample_plans





def _resolve_raw_image_path(metadata, output_dir, sample_name):
    from behav3d.analysis.behavior.state.visualization.backprojection import (
        _resolve_raw_image_path as _shared_resolve_raw_image_path,
    )

    output_dir = Path(output_dir)
    shared_path = _shared_resolve_raw_image_path(output_dir=output_dir, sample_name=sample_name, verbose=False)
    if shared_path is not None:
        return shared_path

    metadata_df = _coerce_metadata(metadata)
    if metadata_df is not None and {"sample_name", "raw_image_path"}.issubset(metadata_df.columns):
        rows = metadata_df[metadata_df["sample_name"].astype("string").str.strip() == str(sample_name)]
        for raw_path in rows["raw_image_path"].tolist():
            if pd.isna(raw_path):
                continue
            candidate = Path(str(raw_path)).expanduser()
            if candidate.exists():
                return candidate

    raise FileNotFoundError(
        f"Could not resolve raw image path for sample '{sample_name}' in '{output_dir}'."
    )


def visualize_multicolor_overlap(metadata, output_dir, sample_name, segment_paths, overlap_mask=None, run=True):
    """
    Launch a Napari viewer with split raw channels, segmentations, and overlap pixels.
    """
    import napari

    segment_paths = [Path(path) for path in segment_paths]
    metadata_df = _coerce_metadata(metadata)
    raw_path = _resolve_raw_image_path(metadata_df, output_dir, sample_name)
    raw_img = load_zarr(raw_path)
    raw_shape = tuple(int(v) for v in raw_img.shape)
    if len(raw_shape) != 5:
        raise ValueError(f"Expected raw image shape (T,C,Z,Y,X), got {raw_shape} from '{raw_path}'.")

    if overlap_mask is None:
        overlap_mask = _build_lazy_overlap_mask(segment_paths)
    elif isinstance(overlap_mask, (str, Path)):
        overlap_mask = load_zarr(overlap_mask)
    elif not hasattr(overlap_mask, "shape"):
        overlap_mask = np.asarray(overlap_mask)

    overlap_view = _normalize_segmentation_labels_for_view(
        overlap_mask.astype(np.uint8, copy=False),
        raw_shape,
        layer_name="overlap_pixels",
    )

    viewer = napari.Viewer()
    channel_colors = [
        "cyan", "yellow", "red", "green", "magenta", "blue",
        "gray", "turbo", "viridis", "plasma", "inferno", "twilight",
    ]

    for ch in range(raw_shape[1]):
        viewer.add_image(
            raw_img[:, ch],
            name=f"raw_ch{ch}",
            colormap=channel_colors[ch % len(channel_colors)],
            blending="additive",
            opacity=0.8,
        )

    for segment_path in segment_paths:
        labels = load_zarr(segment_path)
        labels_view = _normalize_segmentation_labels_for_view(
            labels,
            raw_shape,
            layer_name=segment_path.stem,
        )
        viewer.add_labels(
            labels_view,
            name=segment_path.stem,
            opacity=0.6,
            visible=True,
        )

    viewer.add_labels(
        overlap_view.astype(np.uint8, copy=False),
        name="overlap_pixels",
        opacity=0.8,
        visible=True,
    )

    if bool(run):
        napari.run()
    return viewer


def _load_metadata(metadata_path, output_dir):
    if metadata_path:
        return load_behav3d_metadata(metadata_path)
    candidate = Path(output_dir) / "metadata.csv"
    if candidate.exists():
        return load_behav3d_metadata(candidate)
    raise FileNotFoundError(
        "Visualization requested but no metadata was provided and no metadata.csv was found in output_dir."
    )


def apply_multicolor_segment_correction(
    metadata,
    cell_types,
    output_dir=None,
    overwrite=False,
    erosion_pixels=1,
    min_size=None,
    max_size=None,
    n_workers=1,
):
    """
    Apply multicolor overlap correction to all samples in metadata for the requested cell types.

    Segmentations are resolved from the default BEHAV3D output layout when possible:
    ``images/{sample_name}/{sample_name}_{cell_type}_segments.zarr``.
    Corrected output is written directly to the standard path, overwriting the original.
    Samples are processed sequentially. ``n_workers`` controls timepoint-level
    parallelism within each sample.
    
    If no overlap is detected across all samples, returns early without processing.
    """
    metadata_df = _coerce_metadata(metadata)
    if metadata_df is None:
        raise ValueError("metadata is required.")
    if "sample_name" not in metadata_df.columns:
        raise ValueError("metadata must contain a 'sample_name' column.")

    output_dir = _resolve_output_dir(metadata, output_dir=output_dir)
    cell_types = _normalize_cell_types(cell_types)

    available_cell_types = _available_cell_types_from_metadata(metadata_df)
    missing_cell_types = [cell_type for cell_type in cell_types if cell_type not in available_cell_types]
    if missing_cell_types:
        raise ValueError(
            f"Cell types not found in metadata: {missing_cell_types}. "
            f"Available cell types: {sorted(available_cell_types)}"
        )

    sample_plans = _plan_multicolor_segment_correction(
        metadata=metadata_df,
        output_dir=output_dir,
        cell_types=cell_types,
        overwrite=overwrite,
    )
    n_workers = max(1, int(n_workers or 1))

    if not sample_plans:
        print("No samples to process (all missing required segmentations).")
        return {}

    # Quick check: calculate global overlap across all samples before processing
    all_segment_paths = []
    sample_to_paths_idx = {}
    for sample_plan in sample_plans:
        sample_to_paths_idx[sample_plan["sample_name"]] = len(all_segment_paths)
        output_paths = [item["output_path"] for item in sample_plan["segments"]]
        all_segment_paths.extend(output_paths)

    if all_segment_paths:
        global_stats = calculate_multicolor_overlap(
            all_segment_paths,
            show_progress=False,
            n_workers=n_workers,
        )
        if global_stats["overlap_any_count"] == 0:
            print(
                f"ℹ️ No overlap detected across {len(sample_plans)} sample(s). "
                "Skipping multicolor correction."
            )
            corrected_paths_by_sample = {}
            for sample_plan in sample_plans:
                output_paths = [item["output_path"] for item in sample_plan["segments"]]
                corrected_paths_by_sample[sample_plan["sample_name"]] = output_paths
            return corrected_paths_by_sample

    print(
        f"Applying multicolor segment correction to {len(sample_plans)} sample(s) "
        f"sequentially with {n_workers} worker(s) per sample"
    )

    corrected_paths_by_sample = {}
    for sample_plan in sample_plans:
        output_paths = [item["output_path"] for item in sample_plan["segments"]]
        input_stats = calculate_multicolor_overlap(
            output_paths,
            n_workers=n_workers,
        )
        _print_overlap_summary(
            f"Input overlap for sample '{sample_plan['sample_name']}'",
            input_stats,
        )

        corrected_paths = process_multicolor_segments(
            segment_paths=output_paths,
            output_paths=output_paths,
            overwrite=True,
            erosion_pixels=erosion_pixels,
            min_size=min_size,
            max_size=max_size,
            n_workers=n_workers,
        )
        processed_stats = calculate_multicolor_overlap(
            corrected_paths,
            n_workers=n_workers,
        )
        _print_overlap_summary(
            f"Processed overlap for sample '{sample_plan['sample_name']}'",
            processed_stats,
        )
        corrected_paths_by_sample[sample_plan["sample_name"]] = corrected_paths

    return corrected_paths_by_sample


def _print_overlap_summary(label, stats):
    print(f"{label}:")
    print(f"  overlap_any_count = {stats['overlap_any_count']}")
    print(f"  overlap_all_count = {stats['overlap_all_count']}")
    if stats.get("pairwise_named_stats"):
        for pair_stats in stats["pairwise_named_stats"]:
            print(
                f"    {pair_stats['name']} "
                f"({pair_stats['non_overlapping_count']} non-overlapping, "
                f"{pair_stats['overlapping_count']} overlapping)"
            )


def main(
    segment_paths,
    metadata_path=None,
    output_dir=None,
    sample_name=None,
    output_paths=None,
    overwrite=False,
    erosion_pixels=1,
    min_size=None,
    max_size=None,
    visualize=False,
    run=True,
):
    overwrite=False
    erosion_pixels=1
    min_size=100
    max_size=None
    metadata_path=r"F:\BHVD_BEHAV3D\BEHAV3D_python\runs\NatureBriefComm\LowDensity_MultiColor\metadata.csv"
    output_dir=r"F:\BHVD_BEHAV3D\BEHAV3D_python\runs\NatureBriefComm\LowDensity_MultiColor"
    
    segment_paths = [
        r"F:\BHVD_BEHAV3D\BEHAV3D_python\runs\NatureBriefComm\LowDensity_MultiColor\images\BHVD_SB1_Exp012_Img001\BHVD_SB1_Exp012_Img001_tcell1_segments.zarr",
        r"F:\BHVD_BEHAV3D\BEHAV3D_python\runs\NatureBriefComm\LowDensity_MultiColor\images\BHVD_SB1_Exp012_Img001\BHVD_SB1_Exp012_Img001_tcell2_segments.zarr",
        r"F:\BHVD_BEHAV3D\BEHAV3D_python\runs\NatureBriefComm\LowDensity_MultiColor\images\BHVD_SB1_Exp012_Img001\BHVD_SB1_Exp012_Img001_tcell3_segments.zarr"
    ]
    
    sample_name = "BHVD_SB1_Exp012_Img001"
    
    input_paths = [Path(path) for path in segment_paths]
    input_stats = calculate_multicolor_overlap(input_paths)
    _print_overlap_summary("Input overlap", input_stats)

    metadata = _load_metadata(metadata_path, output_dir)
    visualize_multicolor_overlap(
        metadata=metadata,
        output_dir=output_dir,
        sample_name=sample_name,
        segment_paths=processed_paths,
        overlap_mask=processed_stats["overlap_mask"],
    )
        
    processed_paths = process_multicolor_segments(
        segment_paths=input_paths,
        overwrite=overwrite,
        erosion_pixels=erosion_pixels,
        min_size=min_size,
        max_size=max_size,
    )
    print("Processed outputs:")
    for path in processed_paths:
        print(f"  {path}")

    processed_stats = calculate_multicolor_overlap(processed_paths)
    _print_overlap_summary("Processed overlap", processed_stats)

    if visualize:
        if not sample_name or not output_dir:
            raise ValueError("visualize=True requires both sample_name and output_dir.")
        metadata = _load_metadata(metadata_path, output_dir)
        visualize_multicolor_overlap(
            metadata=metadata,
            output_dir=output_dir,
            sample_name=sample_name,
            segment_paths=processed_paths,
            overlap_mask=processed_stats["overlap_mask"],
        )


    metadata = "/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/runs/NatureBriefComm/LowDensity_MultiColor/metadata.csv"
    metadata = pd.read_csv(metadata)
    cell_types = ["tcell1", "tcell2", "tcell3"]
    output_dir = "/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/runs/NatureBriefComm/LowDensity_MultiColor"
    min_size=150
    apply_multicolor_segment_correction(
        metadata=metadata,
        cell_types=cell_types,
        output_dir=output_dir,
        min_size=min_size,
        overwrite=True,
        n_workers=16
    )

    return processed_paths
