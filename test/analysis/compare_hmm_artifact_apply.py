from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import napari
import numpy as np
import pandas as pd
import scanpy as sc

# from test.state_classification_hmm_widget import (
#     _show_hmm_state_backprojection,
#     _show_intrinsic_hmm_backprojection,
# )


# %%
# Interactive configuration
ORIGINAL_H5AD = r"/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/runs/NatureBriefComm/LowDensity_MultiColor/analysis/tcell/behavioral_states_assigned/processing/BEHAV3D_tcell_behavioral_states_modeldata.h5ad"
APPLIED_H5AD = r"/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/runs/NatureBriefComm/LowDensity_MultiColor/analysis/tcell/behavioral_states/BEHAV3D_tcell_behavioral_states.h5ad"
OUTPUT_DIR = r"/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/runs/NatureBriefComm/LowDensity_MultiColor"
CELL_TYPE = "tcell"
SAMPLE_NAME = "BHVD_SB1_Exp012_Img001"
TRACK_FEATURES_CSV = None
COMPARE_ONLY = False
OPEN_NAPARI = False
VIEW_MODE = "both"  # intrinsic | full | both
COMPARE_COLUMNS = None
NUMERIC_TOLERANCE = 1e-6
MISMATCH_PREVIEW_ROWS = 10
BACKPROJECTION_WORKERS = 4
VERBOSE = True


KEY_COLS = ["sample_name", "TrackID", "position_t"]
DEFAULT_COMPARE_COLS = [
    "hmm_intrinsic_behavioral_state_raw",
    "hmm_intrinsic_behavioral_state",
    "binary_group",
    "behavioral_clusterid",
    "behavioral_state",
    "hmm_intrinsic_behavioral_state_raw_confidence",
    "hmm_intrinsic_behavioral_state_confidence",
]


def _resolve_required_file(path_like, label):
    path = Path(str(path_like)).expanduser()
    if str(path).strip() == "":
        raise ValueError(f"{label} is required.")
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: '{path}'")
    return path


def _resolve_track_features_csv(output_dir, cell_type, override=None):
    if override not in {None, ""}:
        return _resolve_required_file(override, "TRACK_FEATURES_CSV")

    base = Path(output_dir, "analysis", str(cell_type), "track_features")
    filtered_csv = base / f"BEHAV3D_{cell_type}_combined_track_features_filtered.csv"
    combined_csv = base / f"BEHAV3D_{cell_type}_combined_track_features.csv"
    if filtered_csv.exists():
        return filtered_csv
    if combined_csv.exists():
        return combined_csv
    raise FileNotFoundError(
        "Could not resolve track-features CSV automatically. Expected one of: "
        f"'{filtered_csv}' or '{combined_csv}'."
    )


def _load_obs_frame(h5ad_path, sample_name=None):
    adata = sc.read_h5ad(h5ad_path)
    obs = adata.obs.copy()
    missing = [col for col in KEY_COLS if col not in obs.columns]
    if missing:
        raise ValueError(f"Missing required key columns in '{h5ad_path}': {missing}")
    if sample_name not in {None, ""}:
        obs = obs[obs["sample_name"].astype(str) == str(sample_name)].copy()
    if len(obs) == 0:
        raise ValueError(f"No rows found in '{h5ad_path}' for sample_name='{sample_name}'.")
    duplicate_mask = obs.duplicated(subset=KEY_COLS, keep=False)
    if bool(duplicate_mask.any()):
        preview = obs.loc[duplicate_mask, KEY_COLS].head(MISMATCH_PREVIEW_ROWS)
        raise ValueError(
            f"Duplicate comparison keys found in '{h5ad_path}'. Preview:\n{preview.to_string(index=False)}"
        )
    return adata, obs


def _extract_binary_cols(adata_original, adata_applied):
    cols = []
    for adata in [adata_original, adata_applied]:
        clust = getattr(adata, "uns", {}).get("clustering", {})
        if isinstance(clust, dict):
            for col in clust.get("binary_cols_to_merge", []):
                col_s = str(col)
                if col_s not in cols:
                    cols.append(col_s)
    return cols


def _normalize_string_series(series):
    out = pd.Series(series, index=series.index, dtype="string")
    out = out.str.strip()
    return out.fillna("<NA>")


def _series_equal_with_tolerance(left, right, tol):
    left_num = pd.to_numeric(left, errors="coerce")
    right_num = pd.to_numeric(right, errors="coerce")
    both_na = left_num.isna() & right_num.isna()
    close = np.isclose(left_num.to_numpy(dtype=float), right_num.to_numpy(dtype=float), atol=tol, rtol=0.0, equal_nan=True)
    return pd.Series(close | both_na.to_numpy(), index=left.index)


def _print_header(title):
    print("")
    print("=" * len(title))
    print(title)
    print("=" * len(title))


def compare_hmm_h5ad_files(
    *,
    original_h5ad,
    applied_h5ad,
    sample_name=None,
    compare_columns=None,
    numeric_tolerance=1e-6,
):
    adata_original, obs_original = _load_obs_frame(original_h5ad, sample_name=sample_name)
    adata_applied, obs_applied = _load_obs_frame(applied_h5ad, sample_name=sample_name)

    binary_cols = _extract_binary_cols(adata_original, adata_applied)
    if compare_columns is None:
        requested_cols = list(DEFAULT_COMPARE_COLS) + list(binary_cols)
    else:
        requested_cols = [str(col) for col in compare_columns]
    compare_cols = [
        col for col in requested_cols
        if (col in obs_original.columns) and (col in obs_applied.columns)
    ]
    skipped_cols = [col for col in requested_cols if col not in compare_cols]

    original_subset = obs_original[KEY_COLS + compare_cols].copy()
    applied_subset = obs_applied[KEY_COLS + compare_cols].copy()
    original_subset = original_subset.rename(
        columns={col: f"{col}_original" for col in compare_cols}
    )
    applied_subset = applied_subset.rename(
        columns={col: f"{col}_applied" for col in compare_cols}
    )

    merged = original_subset.merge(
        applied_subset,
        on=KEY_COLS,
        how="inner",
        validate="one_to_one",
    )
    if len(merged) == 0:
        raise ValueError("No overlapping rows found between the two h5ad files after key alignment.")

    original_keys = obs_original[KEY_COLS].drop_duplicates()
    applied_keys = obs_applied[KEY_COLS].drop_duplicates()
    only_original = original_keys.merge(applied_keys, on=KEY_COLS, how="left", indicator=True)
    only_original = only_original[only_original["_merge"] == "left_only"].drop(columns=["_merge"])
    only_applied = applied_keys.merge(original_keys, on=KEY_COLS, how="left", indicator=True)
    only_applied = only_applied[only_applied["_merge"] == "left_only"].drop(columns=["_merge"])

    _print_header("HMM Artifact Apply Comparison")
    print(f"original_rows={len(obs_original)}")
    print(f"applied_rows={len(obs_applied)}")
    print(f"aligned_rows={len(merged)}")
    print(f"rows_only_in_original={len(only_original)}")
    print(f"rows_only_in_applied={len(only_applied)}")
    if len(skipped_cols) > 0:
        print(f"skipped_columns_missing_in_one_side={skipped_cols}")

    mismatch_counts = {}
    mismatch_previews = {}

    for col in compare_cols:
        left = merged[f"{col}_original"]
        right = merged[f"{col}_applied"]
        if "confidence" in str(col):
            equal_mask = _series_equal_with_tolerance(left, right, numeric_tolerance)
        else:
            equal_mask = _normalize_string_series(left).eq(_normalize_string_series(right))
        mismatch_mask = ~equal_mask
        mismatch_counts[col] = int(mismatch_mask.sum())
        if bool(mismatch_mask.any()):
            mismatch_previews[col] = merged.loc[
                mismatch_mask,
                KEY_COLS + [f"{col}_original", f"{col}_applied"],
            ].head(MISMATCH_PREVIEW_ROWS)

    print("")
    print("Mismatch counts:")
    for col in compare_cols:
        print(f"  {col}: {mismatch_counts[col]}")

    for col in compare_cols:
        if mismatch_counts[col] == 0:
            continue
        print("")
        print(f"Preview mismatches for '{col}':")
        print(mismatch_previews[col].to_string(index=False))

    summary_groups = {
        "raw_hmm_states": ["hmm_intrinsic_behavioral_state_raw"],
        "intrinsic_labels": ["hmm_intrinsic_behavioral_state"],
        "binary_groups": ["binary_group"] + [col for col in binary_cols if col in compare_cols],
        "full_labels": ["behavioral_clusterid", "behavioral_state"],
        "confidence": [
            "hmm_intrinsic_behavioral_state_raw_confidence",
            "hmm_intrinsic_behavioral_state_confidence",
        ],
    }
    print("")
    print("Summary:")
    for label, cols in summary_groups.items():
        available = [col for col in cols if col in mismatch_counts]
        if len(available) == 0:
            print(f"  {label}: not compared")
            continue
        status = "match" if all(mismatch_counts[col] == 0 for col in available) else "differ"
        print(f"  {label}: {status}")

    if len(only_original) > 0:
        print("")
        print("Rows only in original:")
        print(only_original.head(MISMATCH_PREVIEW_ROWS).to_string(index=False))
    if len(only_applied) > 0:
        print("")
        print("Rows only in applied:")
        print(only_applied.head(MISMATCH_PREVIEW_ROWS).to_string(index=False))

    return {
        "adata_original": adata_original,
        "adata_applied": adata_applied,
        "merged": merged,
        "only_original": only_original,
        "only_applied": only_applied,
        "mismatch_counts": mismatch_counts,
        "compare_cols": compare_cols,
        "binary_cols": binary_cols,
    }


def _copy_layers_from_viewer(target_viewer, source_viewer, *, layer_name_map, visible_default):
    for source_name, target_name in layer_name_map.items():
        if source_name not in [layer.name for layer in source_viewer.layers]:
            continue
        layer = source_viewer.layers[source_name]
        is_visible = bool(visible_default.get(target_name, False))
        if layer.__class__.__name__ == "Image":
            target_viewer.add_image(layer.data, name=target_name, visible=is_visible)
        else:
            target_viewer.add_labels(layer.data, name=target_name, visible=is_visible)


def open_hmm_comparison_viewer(
    *,
    adata_original,
    adata_applied,
    output_dir,
    cell_type,
    sample_name,
    track_features_csv_path,
    view_mode="both",
    n_workers=4,
    verbose=True,
):
    if sample_name in {None, ""}:
        raise ValueError("SAMPLE_NAME is required when OPEN_NAPARI is True.")

    sample_name = str(sample_name).strip()
    view_mode = str(view_mode).strip().lower()
    if view_mode not in {"intrinsic", "full", "both"}:
        raise ValueError("VIEW_MODE must be one of: intrinsic, full, both.")

    viewer_original_intrinsic = None
    viewer_applied_intrinsic = None
    viewer_original_full = None
    viewer_applied_full = None

    if view_mode in {"intrinsic", "both"}:
        viewer_original_intrinsic = _show_intrinsic_hmm_backprojection(
            adata=adata_original,
            sample_name=sample_name,
            output_dir=output_dir,
            cell_type=cell_type,
            track_features_csv_path=track_features_csv_path,
            metadata=None,
            n_workers=n_workers,
            run=False,
            verbose=verbose,
        )
        viewer_applied_intrinsic = _show_intrinsic_hmm_backprojection(
            adata=adata_applied,
            sample_name=sample_name,
            output_dir=output_dir,
            cell_type=cell_type,
            track_features_csv_path=track_features_csv_path,
            metadata=None,
            n_workers=n_workers,
            run=False,
            verbose=verbose,
        )

    if view_mode in {"full", "both"}:
        viewer_original_full = _show_hmm_state_backprojection(
            adata=adata_original,
            sample_name=sample_name,
            output_dir=output_dir,
            cell_type=cell_type,
            state_col="full_behavioral_cluster",
            track_features_csv_path=track_features_csv_path,
            metadata=None,
            n_workers=n_workers,
            run=False,
            verbose=verbose,
        )
        viewer_applied_full = _show_hmm_state_backprojection(
            adata=adata_applied,
            sample_name=sample_name,
            output_dir=output_dir,
            cell_type=cell_type,
            state_col="full_behavioral_cluster",
            track_features_csv_path=track_features_csv_path,
            metadata=None,
            n_workers=n_workers,
            run=False,
            verbose=verbose,
        )

    base_viewer = (
        viewer_original_intrinsic
        or viewer_applied_intrinsic
        or viewer_original_full
        or viewer_applied_full
    )
    if base_viewer is None:
        raise ValueError("No viewer payloads could be created.")

    print("")
    print("Opening Napari comparison viewer...")
    viewer = napari.Viewer()
    base_names = [layer.name for layer in base_viewer.layers]
    if "raw_data" in base_names:
        viewer.add_image(base_viewer.layers["raw_data"].data, name="raw_data", visible=True)
    if "TrackID" in base_names:
        viewer.add_labels(base_viewer.layers["TrackID"].data, name="TrackID", visible=False)

    visible_default = {
        "Original Intrinsic": view_mode == "intrinsic",
        "Applied Intrinsic": view_mode == "intrinsic",
        "Original Full": view_mode in {"full", "both"},
        "Applied Full": view_mode in {"full", "both"},
    }
    if view_mode == "both":
        visible_default["Original Intrinsic"] = False
        visible_default["Applied Intrinsic"] = False

    for source_viewer, target_name in [
        (viewer_original_intrinsic, "Original Intrinsic"),
        (viewer_applied_intrinsic, "Applied Intrinsic"),
        (viewer_original_full, "Original Full"),
        (viewer_applied_full, "Applied Full"),
    ]:
        if source_viewer is None:
            continue
        _copy_layers_from_viewer(
            viewer,
            source_viewer,
            layer_name_map={"ClusterID": target_name},
            visible_default=visible_default,
        )

    for temp_viewer in [viewer_original_intrinsic, viewer_applied_intrinsic, viewer_original_full, viewer_applied_full]:
        if temp_viewer is not None:
            temp_viewer.close()

    print(f"sample_name={sample_name}")
    print(f"view_mode={view_mode}")
    print(f"layers={[layer.name for layer in viewer.layers]}")
    napari.run()
    return viewer


def main():
    original_h5ad = _resolve_required_file(ORIGINAL_H5AD, "ORIGINAL_H5AD")
    applied_h5ad = _resolve_required_file(APPLIED_H5AD, "APPLIED_H5AD")
    output_dir = _resolve_required_file(OUTPUT_DIR, "OUTPUT_DIR")

    result = compare_hmm_h5ad_files(
        original_h5ad=original_h5ad,
        applied_h5ad=applied_h5ad,
        sample_name=SAMPLE_NAME if str(SAMPLE_NAME).strip() != "" else None,
        compare_columns=COMPARE_COLUMNS,
        numeric_tolerance=float(NUMERIC_TOLERANCE),
    )

    if bool(COMPARE_ONLY) or (not bool(OPEN_NAPARI)):
        return result

    track_features_csv_path = _resolve_track_features_csv(
        output_dir=output_dir,
        cell_type=CELL_TYPE,
        override=TRACK_FEATURES_CSV,
    )
    open_hmm_comparison_viewer(
        adata_original=result["adata_original"],
        adata_applied=result["adata_applied"],
        output_dir=output_dir,
        cell_type=CELL_TYPE,
        sample_name=SAMPLE_NAME,
        track_features_csv_path=track_features_csv_path,
        view_mode=VIEW_MODE,
        n_workers=int(BACKPROJECTION_WORKERS),
        verbose=bool(VERBOSE),
    )
    return result


if __name__ == "__main__":
    main()
