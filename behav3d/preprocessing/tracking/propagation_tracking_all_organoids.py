from pathlib import Path
import shutil

import numpy as np
import pandas as pd
from skimage.segmentation import relabel_sequential, watershed
from tqdm import tqdm

from behav3d.core.metadata import detect_organoid_types_from_metadata
from behav3d.io.images import append_to_zarr, load_image
from behav3d.preprocessing.tracking import convert_tracked_image_to_csv
from behav3d.preprocessing.tracking.propagation_tracking import (
    _ensure_metadata_output_columns,
    _resolve_tracking_output_columns,
    _run_single_timepoint_propagation,
)


def _remove_output(path):
    path = Path(path)
    if not path.exists():
        return
    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()


def _resolve_organoid_segments(sample, organoid_types):
    paths = {}
    for organoid_type in organoid_types:
        col_name = f"or_{organoid_type}_segments_image_path"
        segments_path = sample.get(col_name)
        if pd.isna(segments_path) or segments_path is None:
            print(f"Warning: No segmentation data found for {sample['sample_name']}, {organoid_type}. Skipping this organoid type.")
            continue

        segments_path = Path(segments_path)
        if not segments_path.exists():
            print(f"Warning: Segmentation file not found: {segments_path}. Skipping {organoid_type}.")
            continue

        paths[organoid_type] = segments_path
    return paths


def _validate_segment_shapes(sample_name, segment_arrays):
    shapes = {organoid_type: tuple(seg.shape) for organoid_type, seg in segment_arrays.items()}
    unique_shapes = set(shapes.values())
    if len(unique_shapes) > 1:
        raise ValueError(
            f"Organoid segmentations for sample {sample_name} do not share the same shape: {shapes}"
        )


def _pick_seed_coordinate(region_coords, markers):
    if len(region_coords) == 0:
        return None

    centroid = region_coords.mean(axis=0)
    distances = np.sum((region_coords - centroid) ** 2, axis=1)
    for idx in np.argsort(distances):
        coord = tuple(region_coords[idx])
        if markers[coord] == 0:
            return coord
    return tuple(region_coords[np.argmin(distances)])


def _initialize_first_timepoint(t0_segments):
    first_seg = next(iter(t0_segments.values()))
    combined_mask = np.zeros(first_seg.shape, dtype=bool)
    overlap_count = np.zeros(first_seg.shape, dtype=np.uint16)
    combined_labels = np.zeros(first_seg.shape, dtype=np.int32)
    label_coords = {}
    track_to_type = {}
    next_label = 1

    for organoid_type, t0_seg in t0_segments.items():
        relabeled, _, _ = relabel_sequential(np.asarray(t0_seg, dtype=np.int32))
        labels = np.unique(relabeled)
        labels = labels[labels > 0]
        if len(labels) == 0:
            print(f"Warning: No {organoid_type} labels found at timepoint 0. Later appearances will not receive new track IDs.")
            continue

        relabeled[relabeled > 0] += next_label - 1
        current_mask = relabeled > 0
        overlap_count[current_mask] += 1
        combined_mask |= current_mask
        combined_labels[current_mask] = relabeled[current_mask]

        labels = np.unique(relabeled)
        labels = labels[labels > 0]
        for label in labels:
            label = int(label)
            track_to_type[label] = organoid_type
            label_coords[label] = np.argwhere(relabeled == label)
        next_label = int(labels.max()) + 1

    if not track_to_type:
        return np.zeros(first_seg.shape, dtype=np.int32), {}

    markers = combined_labels.copy()
    markers[overlap_count > 1] = 0

    for track_id, region_coords in label_coords.items():
        if np.any(markers == track_id):
            continue
        seed_coord = _pick_seed_coordinate(region_coords, markers)
        if seed_coord is not None:
            markers[seed_coord] = track_id

    tracked_t0 = watershed(combined_mask, markers=markers, mask=combined_mask)
    return np.asarray(tracked_t0, dtype=np.int32), track_to_type


def _merge_organoid_masks(segment_arrays, t):
    merged_mask = None
    for seg in segment_arrays.values():
        current_mask = np.asarray(seg[t]) != 0
        if merged_mask is None:
            merged_mask = current_mask
        else:
            merged_mask |= current_mask
    if merged_mask is None:
        raise ValueError("segment_arrays must contain at least one organoid type")
    return merged_mask.astype(np.uint8)


def _build_output_paths(output_dir, sample_name, organoid_types):
    image_dir = Path(output_dir, "images", sample_name)
    combined_csv_dir = Path(output_dir, "trackdata", sample_name, "all_organoids")
    split = {}
    for organoid_type in organoid_types:
        split[organoid_type] = {
            "img": Path(image_dir, f"{sample_name}_{organoid_type}_tracked.zarr"),
            "csv": Path(output_dir, "trackdata", sample_name, organoid_type, f"{sample_name}_{organoid_type}_tracks.csv"),
        }
    return {
        "image_dir": image_dir,
        "combined_img": Path(image_dir, f"{sample_name}_all_organoids_tracked.zarr"),
        "combined_csv_dir": combined_csv_dir,
        "combined_csv": Path(combined_csv_dir, f"{sample_name}_all_organoids_tracks.csv"),
        "split": split,
    }


def _all_outputs_exist(paths):
    required = [paths["combined_img"], paths["combined_csv"]]
    for organoid_paths in paths["split"].values():
        required.extend([organoid_paths["img"], organoid_paths["csv"]])
    return all(path.exists() for path in required)


def run_propagation_tracking_all_organoids(
    metadata,
    output_dir,
    overwrite=False,
    dilation_nr_pixels=2,
    segment_size_min=100,
    organoid_types=None,
    **kwargs,
):
    organoid_types = list(organoid_types or detect_organoid_types_from_metadata(metadata))
    if not organoid_types:
        print("Warning: No organoid types detected in metadata. Skipping all-organoids propagation tracking.")
        return metadata

    for idx, sample in metadata.iterrows():
        sample_name = sample["sample_name"]
        print(f"Tracking sample with all organoids: {sample_name}")

        segments_paths = _resolve_organoid_segments(sample, organoid_types)
        if not segments_paths:
            print(f"Warning: No organoid segmentations available for {sample_name}. Skipping sample.")
            continue

        available_types = list(segments_paths.keys())
        paths = _build_output_paths(output_dir, sample_name, available_types)
        paths["image_dir"].mkdir(parents=True, exist_ok=True)
        paths["combined_csv_dir"].mkdir(parents=True, exist_ok=True)
        for organoid_type in available_types:
            paths["split"][organoid_type]["csv"].parent.mkdir(parents=True, exist_ok=True)

        if _all_outputs_exist(paths) and not overwrite:
            print("All-organoids tracking already exists... Provide overwrite=True to overwrite... Loading existing tracking data")
        else:
            segment_arrays = {organoid_type: load_image(path) for organoid_type, path in segments_paths.items()}
            _validate_segment_shapes(sample_name, segment_arrays)

            _remove_output(paths["combined_img"])
            for organoid_type in available_types:
                _remove_output(paths["split"][organoid_type]["img"])

            t0_segments = {
                organoid_type: np.asarray(seg[0])
                for organoid_type, seg in segment_arrays.items()
            }
            tracked_t0, track_to_type = _initialize_first_timepoint(t0_segments)
            track_ids_by_type = {
                organoid_type: np.asarray(
                    sorted(track_id for track_id, track_type in track_to_type.items() if track_type == organoid_type),
                    dtype=np.int32,
                )
                for organoid_type in available_types
            }

            seg_prev_tp = tracked_t0.copy()
            n_timepoints = next(iter(segment_arrays.values())).shape[0]
            for t in tqdm(range(n_timepoints), total=n_timepoints):
                if t == 0:
                    tracked_seg = tracked_t0
                else:
                    merged_mask = _merge_organoid_masks(segment_arrays, t)
                    tracked_seg = _run_single_timepoint_propagation(
                        merged_mask,
                        seg_prev_tp.copy(),
                        dilation_nr_pixels=dilation_nr_pixels,
                        segment_size_min=segment_size_min,
                    )

                seg_prev_tp = np.asarray(tracked_seg, dtype=np.int32)
                append_to_zarr(np.expand_dims(seg_prev_tp, axis=0), paths["combined_img"])

                for organoid_type, track_ids in track_ids_by_type.items():
                    split_seg = np.where(np.isin(seg_prev_tp, track_ids), seg_prev_tp, 0)
                    append_to_zarr(
                        np.expand_dims(split_seg.astype(seg_prev_tp.dtype, copy=False), axis=0),
                        paths["split"][organoid_type]["img"],
                    )

            combined_df = convert_tracked_image_to_csv(
                img_path=paths["combined_img"],
                element_size_x=sample["pixel_distance_xy"],
                element_size_y=sample["pixel_distance_xy"],
                element_size_z=sample["pixel_distance_z"],
            )
            if "TrackID" in combined_df.columns:
                combined_df["organoid_type"] = combined_df["TrackID"].map(track_to_type)
            else:
                combined_df["organoid_type"] = pd.Series(dtype=object)
            combined_df.to_csv(paths["combined_csv"], index=False)

            for organoid_type in available_types:
                split_df = combined_df[combined_df["organoid_type"] == organoid_type].copy()
                if "organoid_type" not in split_df.columns:
                    split_df["organoid_type"] = organoid_type
                else:
                    split_df["organoid_type"] = organoid_type
                split_df.to_csv(paths["split"][organoid_type]["csv"], index=False)

        for organoid_type in available_types:
            img_col, csv_col = _resolve_tracking_output_columns(organoid_type, f"or_{organoid_type}_segments_image_path")
            _ensure_metadata_output_columns(metadata, img_col, csv_col)
            metadata.at[idx, img_col] = str(paths["split"][organoid_type]["img"])
            metadata.at[idx, csv_col] = str(paths["split"][organoid_type]["csv"])

    return metadata
