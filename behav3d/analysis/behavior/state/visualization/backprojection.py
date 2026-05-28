import re
import shutil
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import zarr

from behav3d.io.formats.zarr import append_to_zarr
from behav3d.io.images import load_image
from behav3d.analysis.behavior.utils import _mixed_label_sort_key
from behav3d.analysis.behavior.state.utils import (
    _coerce_hex_color,
    _get_classification_state_colors,
    _normalize_label_color_map,
)


def _resolve_tracked_image_path(output_dir, sample_name, cell_type, verbose=False):
    sample_dir = Path(output_dir, "images", str(sample_name))
    if not sample_dir.exists():
        if verbose:
            print(f"Tracked image resolve: sample folder missing '{sample_dir}'")
        return None

    candidates = [
        sample_dir / f"{sample_name}_{cell_type}_tracked.zarr",
        sample_dir / f"{sample_name}_{cell_type}_tracked.zarr.zip",
    ]
    for candidate in candidates:
        if candidate.exists():
            if verbose:
                print(f"Tracked image resolve: using preferred path '{candidate}'")
            return candidate

    fallback = []
    fallback.extend(sorted(sample_dir.glob("*tracked.zarr")))
    fallback.extend(sorted(sample_dir.glob("*tracked.zarr.zip")))
    token = f"_{cell_type}_"
    for path in fallback:
        if token in path.name:
            if verbose:
                print(f"Tracked image resolve: using fallback path '{path}'")
            return path
    if verbose:
        print(
            "Tracked image resolve: no candidate matched expected patterns "
            f"for sample='{sample_name}', cell_type='{cell_type}'."
        )
    return None


def _behavioral_state_backprojection_dir(output_dir, cell_type):
    return Path(
        output_dir,
        "analysis",
        str(cell_type),
        "behavioral_state_trajectories",
        "backprojection",
    )


def _behavioral_state_backprojection_path(output_dir, sample_name, cell_type):
    return Path(
        _behavioral_state_backprojection_dir(output_dir=output_dir, cell_type=cell_type),
        f"{sample_name}_{cell_type}_behavioral_states.zarr",
    )


def _validate_required_obs_columns(obs, required_cols):
    missing_cols = [c for c in required_cols if c not in obs.columns]
    if len(missing_cols) > 0:
        raise ValueError(
            "Cannot export behavioral state backprojection; missing required adata.obs columns: "
            f"{missing_cols}"
        )


def _build_code_map(obs, state_col):
    labels = (
        pd.Series(obs[state_col])
        .dropna()
        .astype("string")
        .str.strip()
    )
    labels = labels[labels != ""]
    unique_labels = sorted(labels.unique().tolist(), key=_mixed_label_sort_key)
    return {str(label): int(idx + 1) for idx, label in enumerate(unique_labels)}


def _build_state_code_color_map(code_map, state_colors=None):
    if not isinstance(code_map, dict) or len(code_map) == 0:
        return {}
    label_colors = _normalize_label_color_map(code_map.keys(), colors=state_colors)
    return {
        str(int(code)): str(label_colors[str(label)])
        for label, code in code_map.items()
        if str(label) in label_colors
    }


def _write_state_color_attrs_to_zarr(path, code_map, state_colors=None):
    label_colors = _normalize_label_color_map(code_map.keys(), colors=state_colors)
    code_colors = _build_state_code_color_map(code_map, state_colors=state_colors)
    if len(code_colors) == 0:
        return {}
    arr = zarr.open(str(path), mode="a")
    arr.attrs["state_label_colors"] = dict(label_colors)
    arr.attrs["state_code_colors"] = dict(code_colors)
    return code_colors


def _apply_state_code_colors_to_layer(layer, code_colors):
    if layer is None or not isinstance(code_colors, dict) or len(code_colors) == 0:
        return False
    try:
        layer.color = {int(code): str(color) for code, color in code_colors.items()}
        layer.metadata["state_code_colors"] = dict(code_colors)
        return True
    except Exception:
        return False


def _safe_file_mtime(path):
    if path is None:
        return None
    try:
        p = Path(path)
        if p.exists():
            return float(p.stat().st_mtime)
    except Exception:
        return None
    return None


def _infer_sample_name_from_obs(sample_obs, fallback="unknown"):
    if sample_obs is None or not isinstance(sample_obs, pd.DataFrame):
        return str(fallback)
    if "sample_name" not in sample_obs.columns:
        return str(fallback)
    vals = (
        pd.Series(sample_obs["sample_name"])
        .dropna()
        .astype("string")
        .str.strip()
    )
    vals = vals[vals != ""]
    unique_vals = vals.unique().tolist()
    if len(unique_vals) == 0:
        return str(fallback)
    if len(unique_vals) == 1:
        return str(unique_vals[0])
    return str(fallback)


def _validate_sample_time_coverage_against_tracked(
    sample_obs,
    tracked_n_timepoints,
    time_col,
    sample_name,
    policy_hint=None,
    enforce=True,
    max_missing_leading_frames=0,
):
    if int(tracked_n_timepoints) <= 0:
        raise ValueError(f"Tracked image has invalid time dimension: {tracked_n_timepoints}")

    times = pd.to_numeric(sample_obs[time_col], errors="coerce")
    times = times.dropna().astype(np.int64)
    observed_ids = sorted(set(int(v) for v in times.tolist()))
    observed_set = set(observed_ids)

    expected_ids = list(range(int(tracked_n_timepoints)))
    missing_ids = [int(t) for t in expected_ids if int(t) not in observed_set]

    missing_leading_frames = []
    for t in expected_ids:
        if int(t) in observed_set:
            break
        missing_leading_frames.append(int(t))

    diag = {
        "sample_name": str(sample_name),
        "time_col": str(time_col),
        "tracked_n_timepoints": int(tracked_n_timepoints),
        "tracked_frame_min": 0,
        "tracked_frame_max": int(tracked_n_timepoints) - 1,
        "observed_time_min": None if len(observed_ids) == 0 else int(min(observed_ids)),
        "observed_time_max": None if len(observed_ids) == 0 else int(max(observed_ids)),
        "observed_time_count": int(len(observed_ids)),
        "missing_time_count_within_tracked": int(len(missing_ids)),
        "missing_leading_frames": list(missing_leading_frames),
        "missing_time_ids_preview": list(missing_ids[:20]),
        "policy_hint": None if policy_hint is None else str(policy_hint),
    }

    allowed_missing = max(0, int(max_missing_leading_frames))
    coverage_ok = int(len(missing_leading_frames)) <= int(allowed_missing)
    diag["coverage_ok"] = bool(coverage_ok)
    diag["max_missing_leading_frames_allowed"] = int(allowed_missing)

    if bool(enforce) and (not coverage_ok):
        msg = (
            "Classifier rows for early frames are missing; regenerate/apply with correct window policy or clear stale "
            "artifacts. "
            f"sample='{diag['sample_name']}', tracked_frames=0..{diag['tracked_frame_max']}, "
            f"observed_time_min={diag['observed_time_min']}, observed_time_max={diag['observed_time_max']}, "
            f"missing_leading_frames={diag['missing_leading_frames'][:20]}, "
            f"policy_hint={diag['policy_hint']}"
        )
        raise ValueError(msg)

    return diag


def backproject_single_sample_behavioral_states(
    tracked_img_path,
    sample_obs,
    state_col,
    track_col,
    time_col,
    output_path,
    code_map,
    raw_image_path=None,
    background_value=0,
    enforce_time_coverage=True,
    coverage_check_max_missing_leading_frames=0,
    sample_name=None,
    policy_hint=None,
    source_h5ad_path=None,
    state_colors=None,
    require_all_rows_present=False,
    verbose=True,
):
    """
    Replace tracked segment labels by behavioral-state integer codes for one sample.

    Parameters
    ----------
    tracked_img_path : str | pathlib.Path
        Path to sample tracked image (*.zarr or *.zarr.zip).
    sample_obs : pandas.DataFrame
        One-sample subset of adata.obs with at least state/track/time columns.
    state_col, track_col, time_col : str
        Column names in sample_obs.
    output_path : str | pathlib.Path
        Destination path for behavioral-state zarr image.
    code_map : dict[str, int]
        Mapping from state label to integer code (0 reserved for background).
    raw_image_path : str | pathlib.Path | None, default None
        Optional raw image path. If supplied, only compatibility is validated.
        Output is always stored in tracked-image shape.
    background_value : int, default 0
        Output value for unmatched voxels.
    verbose : bool, default True
        Print progress messages.
    require_all_rows_present : bool, default False
        If True, every (time_col, track_col) entry in sample_obs must be present
        in the tracked segmentation at that frame; otherwise raise ValueError.
    """
    tracked_img_path = Path(tracked_img_path)
    output_path = Path(output_path)

    tracked_img = load_image(tracked_img_path)
    tracked_shape = tuple(int(v) for v in tracked_img.shape)
    if len(tracked_shape) < 2:
        raise ValueError(
            f"Tracked image has unsupported shape {tracked_shape} at '{tracked_img_path}'."
        )

    raw_shape = None
    if raw_image_path is not None:
        raw_image_path = Path(raw_image_path)
        raw_img = load_image(raw_image_path)
        raw_shape = tuple(int(v) for v in raw_img.shape)
        if len(raw_shape) not in {4, 5}:
            raise ValueError(
                f"Raw image has unsupported shape {raw_shape} at '{raw_image_path}'."
            )
        if int(raw_shape[0]) != int(tracked_shape[0]):
            raise ValueError(
                "Tracked/raw timepoint mismatch for backprojection: "
                f"tracked_shape={tracked_shape}, raw_shape={raw_shape}."
            )
        tracked_spatial = tuple(tracked_shape[-3:])
        raw_spatial = tuple(raw_shape[-3:])
        if tracked_spatial != raw_spatial:
            raise ValueError(
                "Tracked/raw spatial mismatch for backprojection: "
                f"tracked_shape={tracked_shape}, raw_shape={raw_shape}."
            )

    output_shape = tracked_shape

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        if output_path.is_dir():
            shutil.rmtree(output_path)
        else:
            output_path.unlink()

    all_codes = [int(background_value)] + [int(v) for v in code_map.values()]
    min_code = int(min(all_codes)) if len(all_codes) > 0 else int(background_value)
    max_code = int(max(all_codes)) if len(all_codes) > 0 else int(background_value)
    if min_code >= 0 and max_code <= int(np.iinfo(np.uint16).max):
        output_dtype = np.uint16
    elif min_code >= 0 and max_code <= int(np.iinfo(np.uint32).max):
        output_dtype = np.uint32
    else:
        output_dtype = np.int32

    def _align_mapped_frame_to_output(mapped_t, output_frame_shape, dtype):
        mapped_t = np.asarray(mapped_t)
        output_frame_shape = tuple(output_frame_shape)
        if mapped_t.shape == output_frame_shape:
            return mapped_t.astype(dtype, copy=False)

        # (Z, Y, X) -> (C, Z, Y, X)
        if mapped_t.ndim == 3 and len(output_frame_shape) == 4:
            c = int(output_frame_shape[0])
            tiled = np.repeat(np.expand_dims(mapped_t, axis=0), repeats=c, axis=0)
            if tiled.shape == output_frame_shape:
                return tiled.astype(dtype, copy=False)

        # (C, Z, Y, X) -> (Z, Y, X), only allow C=1 squeeze.
        if mapped_t.ndim == 4 and len(output_frame_shape) == 3 and int(mapped_t.shape[0]) == 1:
            squeezed = mapped_t[0]
            if squeezed.shape == output_frame_shape:
                return squeezed.astype(dtype, copy=False)

        # (C1, Z, Y, X) -> (C2, Z, Y, X), allow C1=1 tiling.
        if mapped_t.ndim == 4 and len(output_frame_shape) == 4 and mapped_t.shape[1:] == output_frame_shape[1:]:
            if int(mapped_t.shape[0]) == int(output_frame_shape[0]):
                return mapped_t.astype(dtype, copy=False)
            if int(mapped_t.shape[0]) == 1:
                tiled = np.repeat(mapped_t, repeats=int(output_frame_shape[0]), axis=0)
                if tiled.shape == output_frame_shape:
                    return tiled.astype(dtype, copy=False)

        raise ValueError(
            f"Could not align mapped frame shape {mapped_t.shape} to output frame shape {output_frame_shape}."
        )

    # Matching contract: one-sample export keyed by (position_t, TrackID).
    # position_t is used directly as frame/time id.
    work = sample_obs[[track_col, time_col, state_col]].copy()
    work["_track_id"] = pd.to_numeric(work[track_col], errors="coerce")
    work["_time_id"] = pd.to_numeric(work[time_col], errors="coerce")
    work["_state_label"] = work[state_col].astype("string").str.strip()
    work["_state_code"] = work["_state_label"].map(code_map)

    work = work.dropna(subset=["_track_id", "_time_id", "_state_code"]).copy()
    work = work[work["_state_label"] != ""].copy()
    work["_track_id"] = work["_track_id"].astype(np.int64)
    work["_time_id"] = work["_time_id"].astype(np.int64)
    work["_state_code"] = work["_state_code"].astype(np.int32)

    work = work.sort_values(["_time_id", "_track_id"]).drop_duplicates(
        subset=["_time_id", "_track_id"],
        keep="last",
    )

    resolved_sample_name = (
        str(sample_name)
        if sample_name is not None and len(str(sample_name).strip()) > 0
        else _infer_sample_name_from_obs(sample_obs=sample_obs, fallback="unknown")
    )

    coverage_diag = _validate_sample_time_coverage_against_tracked(
        sample_obs=work,
        tracked_n_timepoints=tracked_shape[0],
        time_col="_time_id",
        sample_name=resolved_sample_name,
        policy_hint=policy_hint,
        enforce=bool(enforce_time_coverage),
        max_missing_leading_frames=int(coverage_check_max_missing_leading_frames),
    )

    if len(work) > 0:
        observed_time_ids = work["_time_id"].to_numpy(dtype=np.int64, copy=False)
        bad_time_mask = (observed_time_ids < 0) | (observed_time_ids >= int(tracked_shape[0]))
        if bool(np.any(bad_time_mask)):
            bad_times = sorted(set(observed_time_ids[bad_time_mask].tolist()))
            preview = [int(v) for v in bad_times[:20]]
            raise ValueError(
                "Backprojection has rows with timepoints outside tracked image range. "
                f"sample='{resolved_sample_name}', tracked_range=[0,{int(tracked_shape[0]) - 1}], "
                f"bad_time_ids_preview={preview}"
            )

    per_time_maps = {}
    if len(work) > 0:
        for time_id, g in work.groupby("_time_id", observed=False):
            per_time_maps[int(time_id)] = dict(
                zip(
                    g["_track_id"].to_numpy(dtype=np.int64),
                    g["_state_code"].to_numpy(dtype=np.int32),
                )
            )

    n_timepoints = tracked_shape[0]
    rows_requested = int(len(work))
    rows_mapped = 0
    rows_unmatched_ignored = 0

    if verbose:
        print(
            "Backprojecting behavioral states (streaming): "
            f"sample='{resolved_sample_name}', source='{tracked_img_path}', "
            f"dest='{output_path}', timepoints={int(n_timepoints)} | "
            f"match_key=(sample_name, {time_col}, {track_col})"
        )

    progress_stride = max(1, min(50, int(n_timepoints)))
    for t in range(n_timepoints):
        label_t = np.asarray(tracked_img[t]).astype(np.int64, copy=False)
        mapping = per_time_maps.get(int(t), {})

        unique_labels, inverse = np.unique(label_t.ravel(), return_inverse=True)
        if len(mapping) > 0:
            present_track_ids = set(int(v) for v in unique_labels.tolist())
            matched_here = sum(
                1 for track_id in mapping.keys() if int(track_id) in present_track_ids
            )
            unmatched_here = int(len(mapping)) - int(matched_here)
            rows_mapped += int(matched_here)
            rows_unmatched_ignored += int(unmatched_here)
        if bool(require_all_rows_present) and len(mapping) > 0:
            missing_track_ids = sorted(
                int(track_id)
                for track_id in mapping.keys()
                if int(track_id) not in present_track_ids
            )
            if len(missing_track_ids) > 0:
                preview = [int(v) for v in missing_track_ids[:20]]
                raise ValueError(
                    "Backprojection row(s) not found in tracked segmentation frame. "
                    f"sample='{resolved_sample_name}', time_id={int(t)}, "
                    f"missing_track_ids_preview={preview}, n_missing={int(len(missing_track_ids))}, "
                    f"tracked_path='{tracked_img_path}'"
                )

        mapped_unique = np.full(unique_labels.shape, int(background_value), dtype=output_dtype)
        if len(mapping) > 0:
            mapped_unique = np.asarray(
                [
                    int(mapping.get(int(lbl), int(background_value)))
                    if int(lbl) >= 0
                    else int(background_value)
                    for lbl in unique_labels
                ],
                dtype=output_dtype,
            )

        mapped_t = mapped_unique[inverse].reshape(label_t.shape)
        mapped_t = _align_mapped_frame_to_output(
            mapped_t,
            output_shape[1:],
            dtype=output_dtype,
        )
        append_to_zarr(img=mapped_t[np.newaxis, ...], outpath=output_path)

        if verbose and ((t == 0) or ((t + 1) % progress_stride == 0) or (t == (n_timepoints - 1))):
            print(f"  backprojection frame {int(t) + 1}/{int(n_timepoints)}")

    out_arr = zarr.open(str(output_path), mode="a")
    code_to_label = {str(code): str(label) for label, code in code_map.items()}
    out_arr.attrs["state_col"] = str(state_col)
    out_arr.attrs["background_value"] = int(background_value)
    out_arr.attrs["label_map"] = dict(code_to_label)
    out_arr.attrs["label_to_code"] = {str(k): int(v) for k, v in code_map.items()}
    state_code_colors = _build_state_code_color_map(code_map, state_colors=state_colors)
    if len(state_code_colors) > 0:
        out_arr.attrs["state_label_colors"] = _normalize_label_color_map(
            code_map.keys(),
            colors=state_colors,
        )
        out_arr.attrs["state_code_colors"] = dict(state_code_colors)
    out_arr.attrs["track_image_path"] = str(tracked_img_path)
    out_arr.attrs["source_tracked_path"] = str(tracked_img_path)
    out_arr.attrs["source_tracked_mtime"] = _safe_file_mtime(tracked_img_path)
    out_arr.attrs["raw_image_path"] = None if raw_image_path is None else str(raw_image_path)
    out_arr.attrs["source_h5ad_path"] = (
        None if source_h5ad_path is None else str(Path(source_h5ad_path))
    )
    out_arr.attrs["source_h5ad_mtime"] = _safe_file_mtime(source_h5ad_path)
    out_arr.attrs["time_coverage_observed_min"] = coverage_diag.get("observed_time_min", None)
    out_arr.attrs["time_coverage_observed_max"] = coverage_diag.get("observed_time_max", None)
    out_arr.attrs["time_coverage_missing_leading_frames"] = list(
        coverage_diag.get("missing_leading_frames", [])
    )
    out_arr.attrs["backprojection_shape_mode"] = "tracked_shape"

    if verbose:
        summary = (
            f"Saved behavioral-state image: {output_path} | "
            f"rows_requested={int(rows_requested)} | rows_mapped={int(rows_mapped)}"
        )
        if not bool(require_all_rows_present):
            summary = (
                f"{summary} | rows_unmatched_ignored={int(rows_unmatched_ignored)}"
            )
        print(summary)

    return {
        "output_path": str(output_path),
        "tracked_img_path": str(tracked_img_path),
        "n_timepoints": int(n_timepoints),
        "n_obs_rows_used": int(len(work)),
        "rows_requested": int(rows_requested),
        "rows_mapped": int(rows_mapped),
        "rows_unmatched_ignored": int(rows_unmatched_ignored),
    }


def export_behavioral_state_backprojection_zarrs(
    adata,
    output_dir,
    cell_type,
    state_col="full_behavioral_cluster",
    sample_col="sample_name",
    track_col="TrackID",
    time_col="position_t",
    background_value=0,
    enforce_time_coverage=True,
    state_colors=None,
    raise_on_error=True,
    require_all_rows_present=False,
    verbose=True,
):
    """
    Export one behavioral-state label image per sample under:
    output_dir/analysis/<cell_type>/behavioral_state_trajectories/backprojection/<sample>_<cell_type>_behavioral_states.zarr
    """
    if adata is None or not hasattr(adata, "obs"):
        raise ValueError("adata with .obs is required for behavioral-state backprojection.")

    _validate_required_obs_columns(
        adata.obs,
        required_cols=[sample_col, track_col, time_col, state_col],
    )

    output_dir = Path(output_dir)
    code_map = _build_code_map(adata.obs, state_col=state_col)
    if len(code_map) == 0:
        raise ValueError(
            f"Cannot export behavioral states because '{state_col}' has no non-empty labels."
        )
    if state_colors is None:
        state_colors = _get_classification_state_colors(adata, state_col)
    state_color_map = _normalize_label_color_map(code_map.keys(), colors=state_colors)

    obs = adata.obs.copy()
    sample_values = (
        pd.Series(obs[sample_col])
        .dropna()
        .astype("string")
        .str.strip()
    )
    sample_values = sample_values[sample_values != ""]
    sample_names = sorted(sample_values.unique().tolist())

    manifest = {
        "state_col": str(state_col),
        "background_value": int(background_value),
        "enforce_time_coverage": bool(enforce_time_coverage),
        "label_map": {str(v): str(k) for k, v in code_map.items()},
        "label_to_code": {str(k): int(v) for k, v in code_map.items()},
        "state_label_colors": dict(state_color_map),
        "state_code_colors": _build_state_code_color_map(code_map, state_colors=state_color_map),
        "output_paths": {},
        "skipped_samples": [],
        "errors": [],
    }

    for sample_name in sample_names:
        sample_obs = obs[obs[sample_col].astype("string") == str(sample_name)].copy()
        policy_hint = None
        try:
            pre_meta = adata.uns.get("preprocessing", {}) if isinstance(getattr(adata, "uns", {}), dict) else {}
            windowing_meta = pre_meta.get("windowing", {}) if isinstance(pre_meta, dict) else {}
            if isinstance(windowing_meta, dict):
                policy_hint = windowing_meta.get("incomplete_window_policy", None)
        except Exception:
            policy_hint = None

        tracked_img_path = _resolve_tracked_image_path(
            output_dir=output_dir,
            sample_name=sample_name,
            cell_type=cell_type,
            verbose=verbose,
        )
        if tracked_img_path is None:
            reason = (
                "Could not resolve tracked image path from expected names or fallback glob for sample "
                f"'{sample_name}' and cell_type '{cell_type}'."
            )
            manifest["skipped_samples"].append(
                {"sample_name": str(sample_name), "reason": reason}
            )
            warnings.warn(reason, RuntimeWarning)
            continue
        raw_img_path = _resolve_raw_image_path(
            output_dir=output_dir,
            sample_name=sample_name,
            verbose=verbose,
        )
        if raw_img_path is None or not Path(raw_img_path).exists():
            reason = (
                "Could not resolve raw image path from sample folder structure or metadata.csv for sample "
                f"'{sample_name}'."
            )
            manifest["skipped_samples"].append(
                {"sample_name": str(sample_name), "reason": reason}
            )
            warnings.warn(reason, RuntimeWarning)
            continue

        output_path = _behavioral_state_backprojection_path(
            output_dir=output_dir,
            sample_name=str(sample_name),
            cell_type=cell_type,
        )
        try:
            sample_result = backproject_single_sample_behavioral_states(
                tracked_img_path=tracked_img_path,
                sample_obs=sample_obs,
                state_col=state_col,
                track_col=track_col,
                time_col=time_col,
                output_path=output_path,
                code_map=code_map,
                raw_image_path=raw_img_path,
                background_value=background_value,
                enforce_time_coverage=bool(enforce_time_coverage),
                coverage_check_max_missing_leading_frames=0,
                sample_name=str(sample_name),
                policy_hint=policy_hint,
                source_h5ad_path=getattr(adata, "filename", None),
                state_colors=state_color_map,
                require_all_rows_present=bool(require_all_rows_present),
                verbose=verbose,
            )
            manifest["output_paths"][str(sample_name)] = str(sample_result["output_path"])
        except Exception as exc:
            error_entry = {
                "sample_name": str(sample_name),
                "error": str(exc),
            }
            manifest["errors"].append(error_entry)
            warnings.warn(
                "Behavioral-state backprojection failed for sample "
                f"'{sample_name}': {exc}",
                RuntimeWarning,
            )

    if verbose:
        print(
            "Behavioral-state backprojection export finished: "
            f"written={len(manifest['output_paths'])}, "
            f"skipped={len(manifest['skipped_samples'])}, "
            f"errors={len(manifest['errors'])}"
        )

    if bool(raise_on_error) and len(manifest["errors"]) > 0:
        first_err = manifest["errors"][0]
        raise RuntimeError(
            "Behavioral-state backprojection export failed for one or more samples. "
            f"first_error_sample='{first_err.get('sample_name', 'unknown')}', "
            f"first_error='{first_err.get('error', 'unknown')}', "
            f"n_errors={int(len(manifest['errors']))}"
        )

    return manifest


def _resolve_raw_image_path(output_dir, sample_name, verbose=False):
    output_dir = Path(output_dir)
    sample_name = str(sample_name).strip()
    sample_dir = Path(output_dir, "images", sample_name)

    if sample_dir.exists():
        preferred = [
            sample_dir / f"{sample_name}.zarr",
            sample_dir / f"{sample_name}.zarr.zip",
        ]
        for candidate in preferred:
            if candidate.exists():
                if verbose:
                    print(f"Raw image resolve: using preferred path '{candidate}'")
                return candidate

        globbed = sorted(list(sample_dir.glob("*.zarr")) + list(sample_dir.glob("*.zarr.zip")))
        sname = sample_name.lower()
        excluded_tokens = [
            "_tracked",
            "_behavioral_states",
            "_segments",
            "_mask",
            "backprojected",
        ]

        fallback = []
        for path in globbed:
            lname = path.name.lower()
            if sname not in lname:
                continue
            if any(token in lname for token in excluded_tokens):
                continue
            fallback.append(path)

        if len(fallback) > 0:
            chosen = sorted(fallback)[0]
            if len(fallback) > 1:
                warnings.warn(
                    "Multiple raw-image candidates found in sample folder; choosing first in sorted order. "
                    f"chosen='{chosen}', candidates={[str(p) for p in sorted(fallback)]}",
                    RuntimeWarning,
                )
            if verbose:
                print(f"Raw image resolve: using fallback path '{chosen}'")
            return chosen

    metadata_csv_path = Path(output_dir, "metadata.csv")
    if metadata_csv_path.exists():
        try:
            metadata = pd.read_csv(metadata_csv_path, low_memory=False)
            if {"sample_name", "raw_image_path"}.issubset(metadata.columns):
                rows = metadata[
                    metadata["sample_name"].astype("string").str.strip() == sample_name
                ]
                for raw_path in rows["raw_image_path"].tolist():
                    if pd.isna(raw_path):
                        continue
                    p = Path(str(raw_path)).expanduser()
                    if p.exists():
                        if verbose:
                            print(f"Raw image resolve: using metadata.csv path '{p}'")
                        return p
        except Exception as exc:
            warnings.warn(
                f"Could not parse metadata.csv for raw image fallback: {exc}",
                RuntimeWarning,
            )
    if verbose:
        print(
            "Raw image resolve: no candidate found for sample "
            f"'{sample_name}' in '{sample_dir}' or metadata.csv."
        )

    return None


def _resolve_behavioral_state_image_path(output_dir, sample_name, cell_type, verbose=False):
    out_dir = _behavioral_state_backprojection_dir(output_dir=output_dir, cell_type=cell_type)
    candidates = [
        out_dir / f"{sample_name}_{cell_type}_behavioral_states.zarr",
        out_dir / f"{sample_name}_{cell_type}_behavioral_states.zarr.zip",
    ]
    for candidate in candidates:
        if candidate.exists():
            if verbose:
                print(f"Behavioral-state image resolve: using path '{candidate}'")
            return candidate

    if verbose:
        print(
            "Behavioral-state image resolve: no candidate matched expected names "
            f"for sample='{sample_name}', cell_type='{cell_type}'."
        )
    return None


def _resolve_behavioral_states_h5ad_path(output_dir, cell_type):
    return Path(
        output_dir,
        "analysis",
        str(cell_type),
        "behavioral_states",
        f"BEHAV3D_{cell_type}_behavioral_states.h5ad",
    )


def _ensure_behavioral_state_backprojection_for_sample(
    sample_name,
    output_dir,
    cell_type,
    state_col=None,
    track_col="TrackID",
    time_col="position_t",
    sample_col="sample_name",
    background_value=0,
    enforce_time_coverage=True,
    refresh_if_stale=True,
    verbose=True,
):
    output_dir = Path(output_dir)
    sample_name = str(sample_name).strip()
    cell_type = str(cell_type).strip()

    tracked_img_path = _resolve_tracked_image_path(
        output_dir=output_dir,
        sample_name=sample_name,
        cell_type=cell_type,
        verbose=verbose,
    )
    if tracked_img_path is None or not Path(tracked_img_path).exists():
        raise FileNotFoundError(
            "Cannot lazily create behavioral-state backprojection because tracked image is missing for "
            f"sample '{sample_name}', cell_type '{cell_type}'."
        )
    h5ad_path = _resolve_behavioral_states_h5ad_path(output_dir=output_dir, cell_type=cell_type)

    existing_state_path = _resolve_behavioral_state_image_path(
        output_dir=output_dir,
        sample_name=sample_name,
        cell_type=cell_type,
        verbose=verbose,
    )
    if existing_state_path is not None and Path(existing_state_path).exists():
        return Path(existing_state_path)

    raw_img_path = _resolve_raw_image_path(
        output_dir=output_dir,
        sample_name=sample_name,
        verbose=verbose,
    )
    if raw_img_path is None or not Path(raw_img_path).exists():
        raise FileNotFoundError(
            "Cannot lazily create behavioral-state backprojection because raw image could not be resolved for "
            f"sample '{sample_name}'. Expected output_dir/images/<sample>/<sample>.zarr(.zip) "
            "or metadata.csv fallback."
        )

    if not h5ad_path.exists():
        raise FileNotFoundError(
            "Cannot lazily create behavioral-state backprojection because full output h5ad is missing: "
            f"'{h5ad_path}'."
        )

    import scanpy as sc

    adata_full = sc.read_h5ad(h5ad_path)
    if not hasattr(adata_full, "obs"):
        raise ValueError(f"Loaded h5ad has no .obs: '{h5ad_path}'.")

    uns_classification = (
        adata_full.uns.get("classification", {})
        if isinstance(getattr(adata_full, "uns", {}), dict)
        else {}
    )
    uns_preprocessing = (
        adata_full.uns.get("preprocessing", {})
        if isinstance(getattr(adata_full, "uns", {}), dict)
        else {}
    )
    windowing_meta = (
        uns_preprocessing.get("windowing", {})
        if isinstance(uns_preprocessing, dict)
        else {}
    )
    policy_hint = (
        windowing_meta.get("incomplete_window_policy", None)
        if isinstance(windowing_meta, dict)
        else None
    )
    applied_full_output_col = (
        uns_classification.get("applied_full_output_col", None)
        if isinstance(uns_classification, dict)
        else None
    )
    resolved_state_col = (
        str(state_col)
        if state_col is not None
        else (
            str(applied_full_output_col)
            if applied_full_output_col is not None and len(str(applied_full_output_col)) > 0
            else "full_behavioral_cluster"
        )
    )

    _validate_required_obs_columns(
        adata_full.obs,
        required_cols=[sample_col, track_col, time_col, resolved_state_col],
    )

    sample_obs = adata_full.obs[
        adata_full.obs[sample_col].astype("string") == str(sample_name)
    ].copy()
    if len(sample_obs) == 0:
        raise ValueError(
            "Cannot lazily create behavioral-state backprojection because sample "
            f"'{sample_name}' is absent in '{h5ad_path}'."
        )

    code_map = _build_code_map(adata_full.obs, state_col=resolved_state_col)
    if len(code_map) == 0:
        raise ValueError(
            f"Cannot lazily create behavioral-state backprojection because '{resolved_state_col}' has no labels."
        )
    state_color_map = _normalize_label_color_map(
        code_map.keys(),
        colors=_get_classification_state_colors(adata_full, resolved_state_col),
    )

    output_path = _behavioral_state_backprojection_path(
        output_dir=output_dir,
        sample_name=sample_name,
        cell_type=cell_type,
    )
    backproject_single_sample_behavioral_states(
        tracked_img_path=tracked_img_path,
        sample_obs=sample_obs,
        state_col=resolved_state_col,
        track_col=track_col,
        time_col=time_col,
        output_path=output_path,
        code_map=code_map,
        raw_image_path=raw_img_path,
        background_value=background_value,
        enforce_time_coverage=bool(enforce_time_coverage),
        coverage_check_max_missing_leading_frames=0,
        sample_name=sample_name,
        policy_hint=policy_hint,
        source_h5ad_path=h5ad_path,
        state_colors=state_color_map,
        verbose=verbose,
    )
    return output_path


def _read_root_zarr_attrs(path):
    path = Path(path)
    if path.suffix == ".zip":
        store = zarr.storage.ZipStore(path, mode="r")
        try:
            arr = zarr.open(store=store, mode="r")
            return dict(getattr(arr, "attrs", {}))
        finally:
            store.close()
    arr = zarr.open(str(path), mode="r")
    return dict(getattr(arr, "attrs", {}))


def _extract_state_label_map(state_img_path):
    attrs = _read_root_zarr_attrs(state_img_path)
    label_map = attrs.get("label_map", None)
    if isinstance(label_map, dict) and len(label_map) > 0:
        return {str(k): str(v) for k, v in label_map.items()}

    label_to_code = attrs.get("label_to_code", None)
    if isinstance(label_to_code, dict) and len(label_to_code) > 0:
        out = {}
        for label, code in label_to_code.items():
            out[str(code)] = str(label)
        return out
    return {}


def _extract_state_code_color_map(state_img_path):
    attrs = _read_root_zarr_attrs(state_img_path)
    code_colors = attrs.get("state_code_colors", None)
    if isinstance(code_colors, dict) and len(code_colors) > 0:
        return {str(k): str(v) for k, v in code_colors.items()}
    label_colors = attrs.get("state_label_colors", None)
    label_to_code = attrs.get("label_to_code", None)
    if isinstance(label_colors, dict) and isinstance(label_to_code, dict):
        out = {}
        for label, code in label_to_code.items():
            if str(label) in label_colors:
                out[str(code)] = str(label_colors[str(label)])
        return out
    return {}


def _state_code_sort_key(code):
    code_s = str(code)
    if re.fullmatch(r"-?\d+", code_s):
        return (0, int(code_s))
    return (1, code_s)


def _build_state_mapping_text(label_map, code_colors=None):
    if not isinstance(label_map, dict) or len(label_map) == 0:
        return "No mapping metadata found in zarr attrs."

    normalized = {str(k): str(v) for k, v in label_map.items()}
    if "0" not in normalized:
        normalized["0"] = "background"
    code_colors = {} if not isinstance(code_colors, dict) else {str(k): str(v) for k, v in code_colors.items()}

    lines = ["State Class Mapping", "idx - color -> full_name", ""]
    for code in sorted(normalized.keys(), key=_state_code_sort_key):
        color = code_colors.get(str(code), "#000000" if str(code) == "0" else "#808080")
        lines.append(f"{code} - {color} -> {normalized[code]}")
    return "\n".join(lines)


def _add_mapping_dock_widget(viewer, mapping_text=None, title="State Class Mapping", label_map=None, code_colors=None):
    try:
        from qtpy.QtWidgets import QLabel, QPlainTextEdit, QHBoxLayout, QWidget, QVBoxLayout

        widget = QWidget()
        layout = QVBoxLayout()
        if isinstance(label_map, dict) and len(label_map) > 0:
            normalized = {str(k): str(v) for k, v in label_map.items()}
            if "0" not in normalized:
                normalized["0"] = "background"
            code_colors = {} if not isinstance(code_colors, dict) else {str(k): str(v) for k, v in code_colors.items()}
            header = QLabel("<b>pixel_value - color -> name</b>")
            layout.addWidget(header)
            for code in sorted(normalized.keys(), key=_state_code_sort_key):
                color = _coerce_hex_color(
                    code_colors.get(str(code), "#000000" if str(code) == "0" else "#808080")
                )
                row = QWidget()
                row_layout = QHBoxLayout()
                row_layout.setContentsMargins(0, 0, 0, 0)
                code_label = QLabel(str(code))
                code_label.setMinimumWidth(48)
                swatch = QLabel("")
                swatch.setFixedSize(22, 16)
                swatch.setStyleSheet(
                    f"background-color: {color}; border: 1px solid #333;"
                )
                arrow_label = QLabel(f"{color} -> {normalized[code]}")
                row_layout.addWidget(code_label)
                row_layout.addWidget(swatch)
                row_layout.addWidget(arrow_label)
                row.setLayout(row_layout)
                layout.addWidget(row)
        else:
            text_box = QPlainTextEdit()
            text_box.setReadOnly(True)
            text_box.setPlainText(str(mapping_text))
            layout.addWidget(text_box)
        widget.setLayout(layout)
        viewer.window.add_dock_widget(widget, area="right", name=str(title))
        return True
    except Exception as exc:
        warnings.warn(
            f"Could not add mapping dock widget to napari window: {exc}",
            RuntimeWarning,
        )
        return False


def _is_dask_array(arr):
    try:
        import dask.array as da

        return isinstance(arr, da.Array)
    except Exception:
        return False


def _broadcast_to_shape(arr, target_shape):
    target_shape = tuple(int(v) for v in target_shape)
    if _is_dask_array(arr):
        import dask.array as da

        return da.broadcast_to(arr, target_shape)
    return np.broadcast_to(arr, target_shape)


def _align_labels_to_raw_shape_for_view(labels_img, raw_img, layer_name, verbose=False):
    raw_shape = tuple(int(v) for v in raw_img.shape)
    labels_shape = tuple(int(v) for v in labels_img.shape)

    if labels_shape == raw_shape:
        if verbose:
            print(f"{layer_name} view align: already matched to raw shape {raw_shape}")
        return labels_img

    if len(raw_shape) not in {4, 5}:
        raise ValueError(
            f"Unsupported raw image shape {raw_shape} for viewer alignment."
        )
    if len(labels_shape) not in {4, 5}:
        raise ValueError(
            f"Unsupported {layer_name} shape {labels_shape} for viewer alignment."
        )

    if int(labels_shape[0]) != int(raw_shape[0]):
        raise ValueError(
            f"{layer_name}/raw timepoint mismatch: labels_shape={labels_shape}, raw_shape={raw_shape}."
        )
    if tuple(labels_shape[-3:]) != tuple(raw_shape[-3:]):
        raise ValueError(
            f"{layer_name}/raw spatial mismatch: labels_shape={labels_shape}, raw_shape={raw_shape}."
        )

    aligned = labels_img

    # Convert labels TZYX -> TCZYX if raw has channel axis.
    if len(raw_shape) == 5 and len(labels_shape) == 4:
        aligned = aligned[:, None, ...]
        labels_shape = tuple(int(v) for v in aligned.shape)

    # Convert labels TCZYX -> TZYX only for singleton C when raw has no channel axis.
    if len(raw_shape) == 4 and len(labels_shape) == 5:
        if int(labels_shape[1]) != 1:
            raise ValueError(
                f"Cannot align {layer_name} with non-singleton channel axis to raw: "
                f"labels_shape={labels_shape}, raw_shape={raw_shape}."
            )
        aligned = aligned[:, 0, ...]
        labels_shape = tuple(int(v) for v in aligned.shape)

    # Match channel count for TCZYX shapes by broadcasting singleton C.
    if len(raw_shape) == 5 and len(labels_shape) == 5:
        raw_c = int(raw_shape[1])
        labels_c = int(labels_shape[1])
        if labels_c != raw_c:
            if labels_c != 1:
                raise ValueError(
                    f"Cannot align {layer_name} channels to raw: labels_shape={labels_shape}, raw_shape={raw_shape}."
                )
            aligned = _broadcast_to_shape(aligned, raw_shape)
            labels_shape = tuple(int(v) for v in aligned.shape)

    if tuple(labels_shape) != tuple(raw_shape):
        raise ValueError(
            f"Could not align {layer_name} to raw shape: labels_shape={labels_shape}, raw_shape={raw_shape}."
        )

    if verbose:
        print(
            f"{layer_name} view align: aligned from {tuple(int(v) for v in labels_img.shape)} "
            f"to raw shape {raw_shape}"
        )
    return aligned


def show_behavioral_state_backprojection(
    sample_name,
    output_dir,
    cell_type,
    raw_image_path=None,
    tracked_img_path=None,
    state_img_path=None,
    auto_create_if_missing=True,
    refresh_if_stale=True,
    state_col=None,
    state_colors=None,
    run=True,
    verbose=True,
):
    """
    Open a single-sample napari view with raw image, TrackID labels, and behavioral-state labels.
    """
    if sample_name is None or len(str(sample_name).strip()) == 0:
        raise ValueError("sample_name is required")
    if output_dir is None:
        raise ValueError("output_dir is required")
    if cell_type is None or len(str(cell_type).strip()) == 0:
        raise ValueError("cell_type is required")

    sample_name = str(sample_name).strip()
    output_dir = Path(output_dir)
    cell_type = str(cell_type).strip()

    raw_path = Path(raw_image_path) if raw_image_path is not None else _resolve_raw_image_path(
        output_dir=output_dir,
        sample_name=sample_name,
        verbose=verbose,
    )
    if raw_path is None or not Path(raw_path).exists():
        raise FileNotFoundError(
            "Could not find raw image for sample "
            f"'{sample_name}'. Expected '{Path(output_dir, 'images', sample_name, f'{sample_name}.zarr')}' "
            "or '.zarr.zip'."
        )

    tracked_path = Path(tracked_img_path) if tracked_img_path is not None else _resolve_tracked_image_path(
        output_dir=output_dir,
        sample_name=sample_name,
        cell_type=cell_type,
        verbose=verbose,
    )
    if tracked_path is None or not Path(tracked_path).exists():
        raise FileNotFoundError(
            "Could not find tracked image for sample "
            f"'{sample_name}' and cell_type '{cell_type}'."
        )

    if state_img_path is not None:
        state_path = Path(state_img_path)
        if not state_path.exists():
            if bool(auto_create_if_missing):
                state_path = _ensure_behavioral_state_backprojection_for_sample(
                    sample_name=sample_name,
                    output_dir=output_dir,
                    cell_type=cell_type,
                    state_col=state_col,
                    track_col="TrackID",
                    time_col="position_t",
                    sample_col="sample_name",
                    background_value=0,
                    enforce_time_coverage=True,
                    refresh_if_stale=bool(refresh_if_stale),
                    verbose=verbose,
                )
            else:
                raise FileNotFoundError(
                    "Could not find behavioral-state image for sample "
                    f"'{sample_name}' and cell_type '{cell_type}'. Expected "
                    f"'{_behavioral_state_backprojection_path(output_dir, sample_name, cell_type)}'."
                )
    else:
        if bool(auto_create_if_missing):
            state_path = _ensure_behavioral_state_backprojection_for_sample(
                sample_name=sample_name,
                output_dir=output_dir,
                cell_type=cell_type,
                state_col=state_col,
                track_col="TrackID",
                time_col="position_t",
                sample_col="sample_name",
                background_value=0,
                enforce_time_coverage=True,
                refresh_if_stale=bool(refresh_if_stale),
                verbose=verbose,
            )
        else:
            state_path = _resolve_behavioral_state_image_path(
                output_dir=output_dir,
                sample_name=sample_name,
                cell_type=cell_type,
                verbose=verbose,
            )
            if state_path is None or not Path(state_path).exists():
                raise FileNotFoundError(
                    "Could not find behavioral-state image for sample "
                    f"'{sample_name}' and cell_type '{cell_type}'. Expected "
                    f"'{_behavioral_state_backprojection_path(output_dir, sample_name, cell_type)}'."
                )

    raw_img = load_image(raw_path)
    tracked_img = load_image(tracked_path)
    state_img = load_image(state_path)
    tracked_img_view = _align_labels_to_raw_shape_for_view(
        tracked_img,
        raw_img,
        layer_name="TrackID",
        verbose=verbose,
    )
    state_img_view = _align_labels_to_raw_shape_for_view(
        state_img,
        raw_img,
        layer_name="behavioral_state_class",
        verbose=verbose,
    )

    import napari

    viewer = napari.Viewer()
    viewer.add_image(raw_img, name="raw_data")
    viewer.add_labels(tracked_img_view, name="TrackID", visible=False)
    state_layer = viewer.add_labels(state_img_view, name="behavioral_state_class", visible=True)

    label_map = _extract_state_label_map(state_path)
    code_colors = _extract_state_code_color_map(state_path)
    if len(code_colors) == 0 and isinstance(state_colors, dict) and len(state_colors) > 0:
        code_map = {str(label): int(code) for code, label in label_map.items()}
        code_colors = _build_state_code_color_map(code_map, state_colors=state_colors)
    _apply_state_code_colors_to_layer(state_layer, code_colors)
    mapping_text = _build_state_mapping_text(label_map, code_colors=code_colors)
    added_dock = _add_mapping_dock_widget(
        viewer=viewer,
        mapping_text=mapping_text,
        label_map=label_map,
        code_colors=code_colors,
        title="State Class Mapping",
    )
    if (not added_dock) and verbose:
        print(mapping_text)

    if verbose:
        print(
            "Opened backprojection viewer for sample "
            f"'{sample_name}' with raw='{raw_path.name}' shape={tuple(int(v) for v in raw_img.shape)}, "
            f"tracked='{tracked_path.name}' shape={tuple(int(v) for v in tracked_img_view.shape)}, "
            f"states='{state_path.name}' shape={tuple(int(v) for v in state_img_view.shape)}."
        )

    if bool(run):
        napari.run()
    return viewer
