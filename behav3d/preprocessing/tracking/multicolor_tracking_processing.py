from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
from tqdm import tqdm

from behav3d.io.images import load_zarr, write_zarr_parallel
from behav3d.core.metadata import expand_multicolor_celltype_names


def _remove_output(path):
    path = Path(path)
    if not path.exists():
        return
    if path.is_dir():
        shutil.rmtree(path)
    else:
        path.unlink()


def _ensure_metadata_output_columns(metadata, *cols):
    for col in cols:
        if col not in metadata.columns:
            metadata[col] = pd.Series(index=metadata.index, dtype=object)
        elif metadata[col].dtype != object:
            metadata[col] = metadata[col].astype(object)


def _resolve_source_tracking_columns(sample, cell_type):
    img_col = None
    csv_col = None
    for prefix in ("or", "im", "ot"):
        candidate_img = f"{prefix}_{cell_type}_tracks_image_path"
        candidate_csv = f"{prefix}_{cell_type}_tracks_csv_path"
        img_value = sample.get(candidate_img)
        csv_value = sample.get(candidate_csv)
        if pd.notna(img_value) and str(img_value).strip():
            img_col = candidate_img
        if pd.notna(csv_value) and str(csv_value).strip():
            csv_col = candidate_csv
        if img_col is not None and csv_col is not None:
            return img_col, csv_col, prefix

    legacy_img = f"{cell_type}_tracks_image_path"
    legacy_csv = f"{cell_type}_tracks_csv_path"
    if img_col is None and legacy_img in sample.index:
        img_col = legacy_img
    if csv_col is None and legacy_csv in sample.index:
        csv_col = legacy_csv
    return img_col, csv_col, None


def _resolve_source_tracking_paths(sample, output_dir, cell_type):
    sample_name = str(sample["sample_name"]).strip()
    img_col, csv_col, prefix = _resolve_source_tracking_columns(sample, cell_type)

    img_path = None
    csv_path = None
    if img_col is not None:
        img_value = sample.get(img_col)
        if pd.notna(img_value) and str(img_value).strip():
            img_path = Path(str(img_value)).expanduser()
    if csv_col is not None:
        csv_value = sample.get(csv_col)
        if pd.notna(csv_value) and str(csv_value).strip():
            csv_path = Path(str(csv_value)).expanduser()

    if img_path is None:
        img_path = Path(output_dir, "images", sample_name, f"{sample_name}_{cell_type}_tracked.zarr")
    if csv_path is None:
        csv_path = Path(output_dir, "trackdata", sample_name, cell_type, f"{sample_name}_{cell_type}_tracks.csv")
    return img_path, csv_path, prefix


def _build_combined_tracking_paths(output_dir, sample_name, combined_cell_type):
    image_dir = Path(output_dir, "images", sample_name)
    csv_dir = Path(output_dir, "trackdata", sample_name, combined_cell_type)
    return {
        "image_dir": image_dir,
        "csv_dir": csv_dir,
        "img": image_dir / f"{sample_name}_{combined_cell_type}_tracked.zarr",
        "csv": csv_dir / f"{sample_name}_{combined_cell_type}_tracks.csv",
    }


def _plan_multicolor_tracking_sample(sample, output_dir, source_cell_types, combined_cell_type):
    sample_name = str(sample["sample_name"]).strip()
    sources = []
    missing = []
    output_prefix = None
    for cell_type in source_cell_types:
        img_path, csv_path, source_prefix = _resolve_source_tracking_paths(sample, output_dir, cell_type)
        if output_prefix is None and source_prefix in {"or", "im", "ot"}:
            output_prefix = source_prefix
        if not img_path.exists():
            missing.append(
                {
                    "sample_name": sample_name,
                    "cell_type": cell_type,
                    "kind": "tracked image",
                    "expected_path": img_path,
                }
            )
        if not csv_path.exists():
            missing.append(
                {
                    "sample_name": sample_name,
                    "cell_type": cell_type,
                    "kind": "tracks csv",
                    "expected_path": csv_path,
                }
            )
        sources.append(
            {
                "cell_type": cell_type,
                "img": img_path,
                "csv": csv_path,
            }
        )

    outputs = _build_combined_tracking_paths(output_dir, sample_name, combined_cell_type)
    return {
        "sample_name": sample_name,
        "sources": sources,
        "outputs": outputs,
        "missing": missing,
        "output_prefix": output_prefix,
    }


def _plan_multicolor_tracking(metadata, output_dir, source_cell_types, combined_cell_type, overwrite):
    sample_name_col = metadata["sample_name"].astype("string").str.strip()
    sample_names = (
        sample_name_col
        .dropna()
        .replace("", pd.NA)
        .dropna()
        .unique()
        .tolist()
    )
    if not sample_names:
        raise ValueError("No valid sample_name values were found in metadata.")

    sample_plans = []
    missing_items = []
    output_conflicts = []
    for sample_name in sample_names:
        sample = metadata[sample_name_col == sample_name].iloc[0]
        sample_plan = _plan_multicolor_tracking_sample(sample, output_dir, source_cell_types, combined_cell_type)
        sample_plans.append(sample_plan)
        missing_items.extend(sample_plan["missing"])

        if not overwrite:
            for path in (sample_plan["outputs"]["img"], sample_plan["outputs"]["csv"]):
                if path.exists():
                    output_conflicts.append(
                        {
                            "sample_name": sample_name,
                            "output_path": path,
                        }
                    )

    if missing_items:
        details = "\n".join(
            f"sample='{item['sample_name']}', cell_type='{item['cell_type']}', "
            f"kind='{item['kind']}', expected_path='{item['expected_path']}'"
            for item in missing_items
        )
        raise FileNotFoundError(
            "Missing required tracked inputs for multicolor combination:\n"
            f"{details}"
        )

    if output_conflicts:
        details = "\n".join(
            f"sample='{item['sample_name']}', output_path='{item['output_path']}'"
            for item in output_conflicts
        )
        raise FileExistsError(
            "Combined multicolor tracking outputs already exist. "
            "Use overwrite=True to replace them:\n"
            f"{details}"
        )

    return sample_plans


def _build_label_offsets(source_tables):
    offsets = {}
    current_offset = 0
    max_shifted_label = 0
    for item in source_tables:
        df = item["df"]
        max_track = int(df["TrackID"].max()) if "TrackID" in df.columns and not df.empty else 0
        max_segment = int(df["SegmentID"].max()) if "SegmentID" in df.columns and not df.empty else 0
        max_label = max(max_track, max_segment, 0)
        offsets[item["cell_type"]] = current_offset
        max_shifted_label = max(max_shifted_label, current_offset + max_label)
        current_offset += max_label
    return offsets, max_shifted_label


def _validate_tracked_shapes(sample_name, source_arrays):
    shapes = {item["cell_type"]: tuple(int(v) for v in item["arr"].shape) for item in source_arrays}
    unique_shapes = set(shapes.values())
    if len(unique_shapes) != 1:
        raise ValueError(
            f"Tracked images for sample '{sample_name}' do not share the same shape: {shapes}"
        )
    return next(iter(unique_shapes))


def _combine_tracked_timepoint(sample_name, t, source_arrays, offsets):
    first_arr = np.asarray(source_arrays[0]["arr"][t])
    combined = np.zeros(first_arr.shape, dtype=np.int64)
    occupancy = np.zeros(first_arr.shape, dtype=np.uint8)
    active_types = []

    for item in source_arrays:
        cell_type = item["cell_type"]
        vol = np.asarray(item["arr"][t])
        mask = vol > 0
        if not np.any(mask):
            continue
        active_types.append((cell_type, mask))
        occupancy += mask.astype(np.uint8)
        shifted = vol.astype(np.int64, copy=True)
        shifted[mask] += int(offsets[cell_type])
        combined[mask] = shifted[mask]

    if np.any(occupancy > 1):
        for idx_a in range(len(active_types)):
            name_a, mask_a = active_types[idx_a]
            for idx_b in range(idx_a + 1, len(active_types)):
                name_b, mask_b = active_types[idx_b]
                if np.any(mask_a & mask_b):
                    raise ValueError(
                        f"Voxel overlap detected while combining tracked images for sample '{sample_name}' "
                        f"at timepoint {t} between '{name_a}' and '{name_b}'."
                    )
        raise ValueError(
            f"Voxel overlap detected while combining tracked images for sample '{sample_name}' at timepoint {t}."
        )

    return combined


def _count_tracks(df):
    if "TrackID" in df.columns:
        return int(df["TrackID"].nunique())
    if "SegmentID" in df.columns:
        return int(df["SegmentID"].nunique())
    return 0


def combine_multicolor_tracked_outputs(
    metadata,
    output_dir,
    source_cell_types,
    combined_cell_type,
    overwrite=False,
    n_workers=1,
):
    metadata = metadata.copy()
    output_dir = Path(output_dir).expanduser()
    source_cell_types = [str(cell_type).strip() for cell_type in source_cell_types if str(cell_type).strip()]
    combined_cell_type = str(combined_cell_type).strip()

    if not source_cell_types:
        raise ValueError("source_cell_types must contain at least one cell type.")
    if not combined_cell_type:
        raise ValueError("combined_cell_type must be a non-empty string.")
    if "sample_name" not in metadata.columns:
        raise ValueError("metadata must contain a 'sample_name' column.")
    if combined_cell_type in source_cell_types:
        raise ValueError("combined_cell_type must differ from all source_cell_types.")

    n_workers = max(1, int(n_workers or 1))
    sample_plans = _plan_multicolor_tracking(
        metadata=metadata,
        output_dir=output_dir,
        source_cell_types=source_cell_types,
        combined_cell_type=combined_cell_type,
        overwrite=overwrite,
    )

    for sample_plan in sample_plans:
        sample_name = sample_plan["sample_name"]
        print(f"Combining multicolor tracking for sample: {sample_name}")

        source_arrays = []
        source_tables = []
        for source in sample_plan["sources"]:
            source_arrays.append(
                {
                    "cell_type": source["cell_type"],
                    "arr": load_zarr(source["img"]),
                }
            )
            source_tables.append(
                {
                    "cell_type": source["cell_type"],
                    "df": pd.read_csv(source["csv"]),
                }
            )

        for item in source_tables:
            n_tracks = _count_tracks(item["df"])
            print(f"  Found {n_tracks} tracks of {item['cell_type']}")

        shape = _validate_tracked_shapes(sample_name, source_arrays)
        offsets, max_shifted_label = _build_label_offsets(source_tables)
        outputs = sample_plan["outputs"]

        outputs["image_dir"].mkdir(parents=True, exist_ok=True)
        outputs["csv_dir"].mkdir(parents=True, exist_ok=True)
        if overwrite:
            _remove_output(outputs["img"])
            _remove_output(outputs["csv"])

        output_dtype = np.int64 if max_shifted_label > 0 else np.asarray(source_arrays[0]["arr"][0]).dtype
        write_zarr_parallel(
            outpath=outputs["img"],
            shape=shape,
            dtype=output_dtype,
            overwrite=True,
        )

        n_timepoints = int(shape[0])

        def _write_timepoint(t):
            combined = _combine_tracked_timepoint(
                sample_name=sample_name,
                t=t,
                source_arrays=source_arrays,
                offsets=offsets,
            )
            write_zarr_parallel(
                outpath=outputs["img"],
                index=t,
                data=combined,
            )
            return t

        if n_workers == 1 or n_timepoints <= 1:
            for t in tqdm(range(n_timepoints), total=n_timepoints, desc=f"  {sample_name}", unit="tp"):
                _write_timepoint(t)
        else:
            with ThreadPoolExecutor(max_workers=n_workers) as executor:
                for _ in tqdm(
                    executor.map(_write_timepoint, range(n_timepoints)),
                    total=n_timepoints,
                    desc=f"  {sample_name}",
                    unit="tp",
                ):
                    pass

        combined_tables = []
        for item in source_tables:
            df = item["df"].copy()
            offset = int(offsets[item["cell_type"]])
            if "TrackID" in df.columns:
                df["TrackID"] = df["TrackID"].astype(np.int64) + offset
            if "SegmentID" in df.columns:
                df["SegmentID"] = df["SegmentID"].astype(np.int64) + offset
            df["source_cell_type"] = item["cell_type"]
            combined_tables.append(df)

        if combined_tables:
            combined_df = pd.concat(combined_tables, ignore_index=True)
        else:
            combined_df = pd.DataFrame(columns=["source_cell_type"])
        total_tracks = _count_tracks(combined_df)
        print(f"  Combined into a total of {total_tracks} tracks")
        combined_df.to_csv(outputs["csv"], index=False)

        img_col = f"{combined_cell_type}_tracks_image_path"
        csv_col = f"{combined_cell_type}_tracks_csv_path"
        prefixed_img_col = None
        prefixed_csv_col = None
        output_prefix = sample_plan.get("output_prefix")
        if output_prefix in {"or", "im", "ot"}:
            prefixed_img_col = f"{output_prefix}_{combined_cell_type}_tracks_image_path"
            prefixed_csv_col = f"{output_prefix}_{combined_cell_type}_tracks_csv_path"

        cols_to_ensure = [img_col, csv_col]
        if prefixed_img_col is not None and prefixed_csv_col is not None:
            cols_to_ensure.extend([prefixed_img_col, prefixed_csv_col])
        _ensure_metadata_output_columns(metadata, *cols_to_ensure)

        row_idx = metadata.index[metadata["sample_name"].astype("string").str.strip() == sample_name][0]
        metadata.at[row_idx, img_col] = str(outputs["img"])
        metadata.at[row_idx, csv_col] = str(outputs["csv"])
        if prefixed_img_col is not None and prefixed_csv_col is not None:
            metadata.at[row_idx, prefixed_img_col] = str(outputs["img"])
            metadata.at[row_idx, prefixed_csv_col] = str(outputs["csv"])

    return metadata


def combine_multicolor_tracked_outputs_for_base(
    metadata,
    output_dir,
    base_cell_type,
    n_channels,
    overwrite=False,
    n_workers=1,
):
    """Convenience wrapper that combines a numbered multicolor set for one base cell type."""
    source_cell_types = expand_multicolor_celltype_names(base_cell_type, n_channels)
    combined_cell_type = f"{str(base_cell_type).strip()}_merged"
    return combine_multicolor_tracked_outputs(
        metadata=metadata,
        output_dir=output_dir,
        source_cell_types=source_cell_types,
        combined_cell_type=combined_cell_type,
        overwrite=overwrite,
        n_workers=n_workers,
    )
