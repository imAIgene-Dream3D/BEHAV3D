import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas.api.types import is_numeric_dtype
import json
import copy
import warnings
import time
import shutil
from matplotlib.backends.backend_pdf import PdfPages

from pathlib import Path
import pickle

import anndata as ad
ad.settings.allow_write_nullable_strings = True

import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.metrics import accuracy_score, balanced_accuracy_score, classification_report, confusion_matrix, f1_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from behav3d.core.anndata import df_to_adata
from behav3d.features.rolling_window_features import create_descriptive_track_dataset
from behav3d.analysis.behavior.general import relabel_cluster_ids
from behav3d.analysis.behavior.general.leiden import run_pca, run_leiden_clustering
from behav3d.analysis.behavior.state.visualization.plots.state_composition import (
    save_state_composition_report,
)
from behav3d.analysis.behavior.state.visualization.plots.state_transitions import (
    save_state_transition_report,
)
from behav3d.analysis.behavior.state.visualization.backprojection import (
    export_behavioral_state_backprojection_zarrs,
    show_behavioral_state_backprojection,
)
from behav3d.analysis.behavior.state.visualization.videos.track_max_projection import (
    save_selected_fulltrack_cluster_videos,
)
from behav3d.analysis.behavior.state.utils import (
    A4_LANDSCAPE,
    _apply_log1p_to_feature_matrix,
    _apply_log_scaling_to_continuous_matrix,
    _coerce_log_scaling_params,
    _infer_binary_group_constraints,
    _invert_log_scaling_in_continuous_matrix,
    _mixed_label_sort_key,
    _normalize_log_scale_feature_selectors,
    _rebuild_full_behavioral_cluster_from_intrinsic,
    _require_columns,
    _resolve_log_scale_feature_cols,
    _resolve_positions_csv_path,
    _resolve_state_paths,
    _resolve_state_stage_paths,
    _save_pdf_page_a4,
    _sanitize_filename_token,
    _to_numpy_2d,
    _vdone,
    _vinfo,
    _vsave,
    _vstart,
    cap_values_to_quantile,
)


def _short_reasons(reasons):
    if reasons is None:
        return "none"
    if not isinstance(reasons, (list, tuple)):
        return str(reasons)
    if len(reasons) == 0:
        return "none"
    return f"count={len(reasons)}, first={reasons[0]}"


_PER_SIGNAL_DESCRIPTOR_SUFFIXES = {
    "mean": ("_mean_value",),
    "median": ("_median_value",),
    "range": ("_value_range",),
    "std": ("_standard_deviation",),
    "min": ("_minimum_value",),
    "max": ("_maximum_value",),
    "iqr": ("_interquartile_range",),
    "mad": ("_median_absolute_deviation",),
    "trend": ("_linear_trend_slope_per_time_unit",),
    "lag1_autocorr": ("_lag1_autocorrelation",),
    "mean_abs_first_diff": ("_mean_absolute_first_difference",),
    "fraction_near_bounds": (
        "_fraction_near_lower_bound",
        "_fraction_near_upper_bound",
    ),
    "binary_runs": (
        "_longest_true_length",
        "_mean_true_length",
        "_longest_false_length",
        "_mean_false_length",
    ),
    "transition_rate": ("_transition_rate",),
    "count_dispersion": (
        "_dispersion_index_variance_over_mean",
        "_fraction_of_zeros",
    ),
}

_GLOBAL_DESCRIPTOR_COLUMNS = {
    "summed_displacement": "summed_displacement",
    "net_displacement": "net_displacement",
    "straightness": "straightness",
    "directional_persistence": "directional_persistence",
    "median_turning_angle": "median_turning_angle",
    "fraction_reversed_movement": "fraction_reversed_movement",
    "mean_square_displacement": "mean_square_displacement",
}

_QUANTILE_DESCRIPTOR_PCTS = (10, 25, 75, 90)


def _descriptor_column_mismatch_reasons(feature_cols, *, features, descriptive_features):
    """Return reasons when cached feature columns conflict with descriptor selection."""
    feature_cols = [] if feature_cols is None else list(feature_cols)
    features = [] if features is None else list(features)
    descriptive_features = [] if descriptive_features is None else list(descriptive_features)
    feature_set = {str(c) for c in feature_cols}
    base_features = [str(c) for c in features]
    requested = {str(c) for c in descriptive_features}
    unexpected = []

    for descriptor, suffixes in _PER_SIGNAL_DESCRIPTOR_SUFFIXES.items():
        if descriptor in requested:
            continue
        for feature_name in base_features:
            for suffix in suffixes:
                col = f"{feature_name}{suffix}"
                if col in feature_set:
                    unexpected.append(col)

    if "quantiles" not in requested:
        for feature_name in base_features:
            for pct in _QUANTILE_DESCRIPTOR_PCTS:
                col = f"{feature_name}_quantile_{pct}percent"
                if col in feature_set:
                    unexpected.append(col)

    for descriptor, col in _GLOBAL_DESCRIPTOR_COLUMNS.items():
        if descriptor not in requested and col in feature_set:
            unexpected.append(col)

    if len(unexpected) == 0:
        return []
    unexpected = sorted(set(unexpected))
    return [
        "feature columns conflict with requested descriptive_features "
        f"(unexpected columns: {unexpected[:10]})"
    ]




def build_identity_cluster_mapping(
    adata,
    cluster_col="intrinsic_behavioral_cluster",
):
    """
    Build an identity mapping dict from unique values in a cluster column.

    Example:
        {"dead": "dead", "scanner": "scanner", "static": "static"}
    """
    if not hasattr(adata, "obs"):
        raise ValueError("adata must have an .obs attribute.")
    if cluster_col not in adata.obs.columns:
        raise ValueError(f"Missing '{cluster_col}' in adata.obs.")

    labels = (
        pd.Series(adata.obs[cluster_col])
        .astype("string")
        .dropna()
        .unique()
        .tolist()
    )
    labels = sorted([str(x) for x in labels], key=_mixed_label_sort_key)
    return {label: label for label in labels}




def prepare_state_classification_dataset(
    df_positions,
    features,
    binary_features_to_group,
    window_size=5,
    min_spacing=None,
    descriptive_features=("mean", "median", "std", "net_displacement", "straightness", "mean_square_displacement"),
    lower_quantile_cap=None,
    upper_quantile_cap=0.99,
    log_scale_features=None,
    outfolder=None,
    scale_features=False,
    incomplete_window_policy="partial",
    reuse_prepared_dataset=True,
    save_prepared_dataset=False,
    prepared_dataset_path=None,
    verbose=True,
):
    """Create the windowed descriptive dataset used for clustering/classification.

    Returns:
        sc.AnnData: prepared dataset (scaled if scale_features=True).
    """
    groupby = ["sample_name", "TrackID"]
    non_feature_cols = [
        "sample_name",
        "TrackID",
        "sub_TrackID",
        "position_t",
        "window_start_position_t",
        "window_end_position_t",
        "window_length_frames",
    ]

    if prepared_dataset_path is None and outfolder is not None:
        prepared_dataset_path = Path(outfolder) / "BEHAV3D_behavioral_states.h5ad"
    elif prepared_dataset_path is not None:
        prepared_dataset_path = Path(prepared_dataset_path)

    log_scale_selectors = _normalize_log_scale_feature_selectors(log_scale_features)
    quantile_cap_limits = {}
    log_scaling_params = None
    loaded_from_cache = False
    adata_prepared = None
    cached_pre_meta = {}
    if reuse_prepared_dataset and prepared_dataset_path is not None and prepared_dataset_path.exists():
        try:
            cached_prepared = sc.read_h5ad(prepared_dataset_path)
            pre_match_cached, pre_reasons_cached = _matches_requested_preprocessing_in_adata(
                cached_prepared,
                features=features,
                binary_features_to_group=binary_features_to_group,
                window_size=window_size,
                min_spacing=min_spacing,
                descriptive_features=descriptive_features,
                incomplete_window_policy=incomplete_window_policy,
                lower_quantile_cap=lower_quantile_cap,
                upper_quantile_cap=upper_quantile_cap,
                log_scale_features=log_scale_selectors,
                scale_features=scale_features,
            )
            if pre_match_cached:
                adata_prepared = cached_prepared
                loaded_from_cache = True
                _vinfo(verbose, "state-clustering", f"Loaded prepared window cache from: {prepared_dataset_path}")
            else:
                _vinfo(
                    verbose,
                    "state-clustering",
                    "Prepared cache mismatch; recomputing windows | "
                    f"{_short_reasons(pre_reasons_cached)}",
                )
        except Exception as exc:
            _vinfo(verbose, "state-clustering", f"Prepared cache load failed; recomputing windows ({exc})")

    if not loaded_from_cache:
        safe_features = []
        dropped_features = []
        for col in features:
            if col not in df_positions.columns:
                continue
            series = df_positions[col]
            if pd.api.types.is_numeric_dtype(series):
                safe_features.append(col)
                continue
            coerced = pd.to_numeric(series, errors="coerce")
            if coerced.notna().any():
                df_positions[col] = coerced
                safe_features.append(col)
            else:
                dropped_features.append(col)
        if dropped_features:
            _vinfo(
                verbose,
                "state-clustering",
                f"Dropping non-numeric features from windowing: {dropped_features[:10]}"
                + (" ..." if len(dropped_features) > 10 else ""),
            )

        df_windows_descriptive = create_descriptive_track_dataset(
            df_tracks=df_positions,
            columns_to_summarize=safe_features,
            window_size=window_size,
            step_size=1,
            time_col="position_t",
            id_cols=groupby,
            features_to_compute=list(descriptive_features),
            only_nonbinary=True,
            incomplete_window_policy=incomplete_window_policy,
        )

        descriptive_feature_cols = [
            col for col in df_windows_descriptive.columns
            if (col not in non_feature_cols) and (not col.endswith("_signal_type"))
        ]

        binary_cols_to_merge = [col for col in binary_features_to_group if col in df_positions.columns]
        merge_cols = ["sample_name", "TrackID", "position_t"]
        df_binary = df_positions[merge_cols + binary_cols_to_merge].copy()

        df_prepared = df_windows_descriptive.merge(df_binary, on=merge_cols, how="left", suffixes=("", "_orig"))
        if str(incomplete_window_policy).lower() == "drop":
            df_prepared = df_prepared.dropna(subset=descriptive_feature_cols).copy()

        empty_descriptive_cols = [c for c in descriptive_feature_cols if c in df_prepared.columns and df_prepared[c].isna().all()]
        if len(empty_descriptive_cols) > 0:
            df_prepared = df_prepared.drop(columns=empty_descriptive_cols, errors="ignore")
            descriptive_feature_cols = [c for c in descriptive_feature_cols if c not in empty_descriptive_cols]

        binary_prefixes = tuple(f"{c}_" for c in binary_cols_to_merge)
        kept_features = [
            c for c in descriptive_feature_cols
            if (c not in binary_cols_to_merge) and (not c.startswith(binary_prefixes))
        ]

        if len(kept_features) > 0 and (lower_quantile_cap is not None or upper_quantile_cap is not None):
            df_prepared, quantile_cap_limits = cap_values_to_quantile(
                df_prepared,
                kept_features,
                lower_quantile=lower_quantile_cap,
                upper_quantile=upper_quantile_cap,
                return_limits=True,
            )

        log_scaling_params = _resolve_log_scale_feature_cols(
            kept_features,
            log_scale_selectors,
        )
        if len(log_scaling_params["unresolved_features"]) > 0:
            raise ValueError(
                "log_scale_features could not be resolved to prepared feature columns: "
                f"{log_scaling_params['unresolved_features'][:20]}"
            )
        if len(log_scaling_params["resolved_feature_cols"]) > 0:
            X_log = df_prepared.loc[:, kept_features].to_numpy(dtype=float, copy=True)
            X_log = _apply_log1p_to_feature_matrix(
                X_log,
                feature_cols=kept_features,
                resolved_feature_cols=log_scaling_params["resolved_feature_cols"],
                inplace=True,
            )
            df_prepared.loc[:, kept_features] = X_log
            del X_log

        obs_cols = [c for c in non_feature_cols if c in df_prepared.columns] + [
            c for c in binary_cols_to_merge if c in df_prepared.columns
        ]
        adata_prepared = df_to_adata(df_prepared, kept_features, obs_cols=obs_cols)
        del df_windows_descriptive, df_binary, df_prepared

    else:
        binary_cols_to_merge = [col for col in binary_features_to_group if col in adata_prepared.obs.columns]
        kept_features = [str(c) for c in list(adata_prepared.var_names)]
        cached_pre_meta = (
            adata_prepared.uns.get("preprocessing", {})
            if (hasattr(adata_prepared, "uns") and isinstance(adata_prepared.uns, dict))
            else {}
        )
        if isinstance(cached_pre_meta, dict):
            cached_kept_features = cached_pre_meta.get("kept_features", None)
            if isinstance(cached_kept_features, list) and len(cached_kept_features) > 0:
                valid_kept_features = [str(c) for c in cached_kept_features if str(c) in adata_prepared.var_names]
                if len(valid_kept_features) > 0:
                    kept_features = valid_kept_features
                    adata_prepared = adata_prepared[:, kept_features].copy()

            qcap_meta = cached_pre_meta.get("quantile_capping", {})
            if isinstance(qcap_meta, dict):
                feature_limits = qcap_meta.get("feature_limits", None)
                if isinstance(feature_limits, dict) and len(quantile_cap_limits) == 0:
                    quantile_cap_limits = dict(feature_limits)
                if "lower_quantile" in qcap_meta:
                    lower_quantile_cap = qcap_meta.get("lower_quantile", None)
                if "upper_quantile" in qcap_meta:
                    upper_quantile_cap = qcap_meta.get("upper_quantile", None)
            log_scaling_params = _coerce_log_scaling_params(cached_pre_meta)

    scaler = None
    if scale_features and len(kept_features) > 0:
        duplicate_kept_features = []
        _seen_kept_features = set()
        for col in kept_features:
            scol = str(col)
            if scol in _seen_kept_features and scol not in duplicate_kept_features:
                duplicate_kept_features.append(scol)
            _seen_kept_features.add(scol)
        if len(duplicate_kept_features) > 0:
            raise ValueError(
                "Duplicate kept_features are not supported for deterministic scaling/index assignment. "
                f"Duplicates (first 20): {duplicate_kept_features[:20]}"
            )
        cached_scaler_meta = (
            cached_pre_meta.get("scaler", None)
            if isinstance(cached_pre_meta, dict)
            else None
        )
        if (
            loaded_from_cache
            and isinstance(cached_scaler_meta, dict)
            and ("mean" in cached_scaler_meta)
            and ("scale" in cached_scaler_meta)
        ):
            scaler = StandardScaler()
            scaler.mean_ = np.asarray(cached_scaler_meta["mean"], dtype=float)
            scaler.scale_ = np.asarray(cached_scaler_meta["scale"], dtype=float)
        else:
            X_prepared = _to_numpy_2d(adata_prepared[:, kept_features].X).astype(float, copy=False)
            scaler = StandardScaler().fit(X_prepared)
            scaled_prepared = scaler.transform(X_prepared)
            X_all = _to_numpy_2d(adata_prepared.X).astype(float, copy=True)
            feat_indices = adata_prepared.var_names.get_indexer(kept_features)
            if np.any(feat_indices < 0):
                missing = [kept_features[i] for i, idx in enumerate(feat_indices) if idx < 0]
                raise ValueError(f"Missing kept_features in adata_prepared.var_names: {missing[:20]}")
            X_all[:, feat_indices] = scaled_prepared
            adata_prepared.X = X_all
            del X_prepared, scaled_prepared, X_all

    windowing_params = {
        "features": list(features),
        "binary_features_to_group": list(binary_features_to_group),
        "window_size": int(window_size),
        "min_spacing": None if min_spacing is None else int(min_spacing),
        "descriptive_features": list(descriptive_features),
        "incomplete_window_policy": str(incomplete_window_policy),
        "lower_quantile_cap": lower_quantile_cap,
        "upper_quantile_cap": upper_quantile_cap,
    }
    quantile_capping_params = {
        "lower_quantile": lower_quantile_cap,
        "upper_quantile": upper_quantile_cap,
        "feature_limits": dict(quantile_cap_limits or {}),
    }
    if log_scaling_params is None:
        log_scaling_params = {
            "transform": "log1p",
            "requested_features": list(log_scale_selectors),
            "resolved_feature_cols": [],
        }
    if len(log_scaling_params.get("unresolved_features", [])) > 0:
        raise ValueError(
            "log_scale_features could not be resolved to prepared feature columns: "
            f"{log_scaling_params['unresolved_features'][:20]}"
        )

    adata_prepared.uns["preprocessing"] = {
        "kept_features": list(kept_features),
        "non_feature_cols": list(non_feature_cols),
        "binary_cols_to_merge": list(binary_cols_to_merge),
        "quantile_capping": dict(quantile_capping_params),
        "windowing": dict(windowing_params),
        "features": list(windowing_params.get("features", [])),
        "binary_features_to_group": list(windowing_params.get("binary_features_to_group", [])),
        "window_size": int(windowing_params.get("window_size", 0)),
        "descriptive_features": list(windowing_params.get("descriptive_features", [])),
        "incomplete_window_policy": str(windowing_params.get("incomplete_window_policy", None)),
        "lower_quantile_cap": windowing_params.get("lower_quantile_cap", None),
        "upper_quantile_cap": windowing_params.get("upper_quantile_cap", None),
    }
    if len(log_scaling_params.get("requested_features", [])) > 0 or len(log_scaling_params.get("resolved_feature_cols", [])) > 0:
        adata_prepared.uns["preprocessing"]["log_scaling"] = {
            "transform": "log1p",
            "requested_features": list(log_scaling_params.get("requested_features", [])),
            "resolved_feature_cols": list(log_scaling_params.get("resolved_feature_cols", [])),
        }
    if scaler is not None:
        adata_prepared.uns["preprocessing"]["scaler"] = {
            "mean": scaler.mean_.astype(float),
            "scale": scaler.scale_.astype(float),
        }
    adata_prepared.uns.pop("clustering", None)

    if save_prepared_dataset and prepared_dataset_path is not None:
        prepared_dataset_path.parent.mkdir(parents=True, exist_ok=True)
        adata_prepared.write(prepared_dataset_path, compression="gzip")
    return adata_prepared


def _resolve_legacy_prepared_full_adata_path(state_paths, cell_type):
    return state_paths.state_outdir / f"BEHAV3D_{cell_type}_behavioral_states_prepared_full.h5ad"
























def _coerce_scaler_params(preprocessing_params):
    """Extract scaler params from preprocessing params; supports legacy 'zscore' key."""
    if not isinstance(preprocessing_params, dict):
        return None
    scaler_meta = preprocessing_params.get("scaler", None)
    if not isinstance(scaler_meta, dict):
        scaler_meta = preprocessing_params.get("zscore", None)
    if not isinstance(scaler_meta, dict):
        return None
    mean = scaler_meta.get("mean", None)
    scale = scaler_meta.get("scale", None)
    if mean is None or scale is None:
        return None
    return {
        "mean": np.asarray(mean, dtype=float),
        "scale": np.asarray(scale, dtype=float),
    }


def _build_preprocessing_params_from_uns(preprocessing_meta, feature_cols=None):
    """Build reusable preprocessing parameter dict from adata.uns['preprocessing'] metadata."""
    if not isinstance(preprocessing_meta, dict):
        return None

    kept_features = preprocessing_meta.get("kept_features", feature_cols)
    if kept_features is None:
        kept_features = []
    kept_features = [str(c) for c in kept_features]

    qmeta = preprocessing_meta.get("quantile_capping", {})
    if not isinstance(qmeta, dict):
        qmeta = {}
    feature_limits = qmeta.get("feature_limits", {})
    if not isinstance(feature_limits, dict):
        feature_limits = {}

    params = {
        "continuous_feature_cols": list(kept_features),
        "quantile_capping": {
            "lower_quantile": qmeta.get("lower_quantile", None),
            "upper_quantile": qmeta.get("upper_quantile", None),
            "feature_limits": dict(feature_limits),
        },
        "log_scaling": None,
        "scaler": None,
    }
    log_scaling_params = _coerce_log_scaling_params(preprocessing_meta)
    if log_scaling_params is not None:
        params["log_scaling"] = log_scaling_params
    scaler_params = _coerce_scaler_params(preprocessing_meta)
    if scaler_params is not None:
        params["scaler"] = scaler_params
    return params


def _normalize_preprocessing_params_for_compare(preprocessing_params, windowing_params, feature_cols):
    """Normalize preprocessing + windowing params into deterministic compare payload."""
    feature_cols = [str(c) for c in list(feature_cols or [])]
    pre = preprocessing_params if isinstance(preprocessing_params, dict) else {}
    window = windowing_params if isinstance(windowing_params, dict) else {}

    qmeta = pre.get("quantile_capping", {})
    qmeta = qmeta if isinstance(qmeta, dict) else {}
    feature_limits = qmeta.get("feature_limits", {})
    feature_limits = feature_limits if isinstance(feature_limits, dict) else {}
    feature_limits_norm = {}
    for feat in sorted([str(k) for k in feature_limits.keys()]):
        lim = feature_limits.get(feat, {})
        lim = lim if isinstance(lim, dict) else {}
        lower = lim.get("lower", None)
        upper = lim.get("upper", None)
        feature_limits_norm[feat] = {
            "lower": None if lower is None else float(lower),
            "upper": None if upper is None else float(upper),
        }

    scaler_meta = _coerce_scaler_params(pre)
    scaler_norm = None
    if isinstance(scaler_meta, dict):
        scaler_norm = {
            "mean": np.asarray(scaler_meta.get("mean", []), dtype=float).tolist(),
            "scale": np.asarray(scaler_meta.get("scale", []), dtype=float).tolist(),
        }

    log_scaling_meta = _coerce_log_scaling_params(pre)
    log_scaling_norm = None
    if isinstance(log_scaling_meta, dict):
        log_scaling_norm = {
            "transform": log_scaling_meta.get("transform", None),
            "requested_features": [str(x) for x in list(log_scaling_meta.get("requested_features", []))],
            "resolved_feature_cols": [str(x) for x in list(log_scaling_meta.get("resolved_feature_cols", []))],
        }

    windowing_keys = [
        "features",
        "binary_features_to_group",
        "window_size",
        "min_spacing",
        "descriptive_features",
        "incomplete_window_policy",
        "lower_quantile_cap",
        "upper_quantile_cap",
    ]
    windowing_norm = {}
    for key in windowing_keys:
        val = window.get(key, None)
        if key in {"features", "binary_features_to_group", "descriptive_features"}:
            windowing_norm[key] = [] if val is None else [str(x) for x in list(val)]
        elif key in {"window_size", "min_spacing"}:
            windowing_norm[key] = None if val is None else int(val)
        elif key in {"lower_quantile_cap", "upper_quantile_cap"}:
            windowing_norm[key] = None if val is None else float(val)
        elif key == "incomplete_window_policy":
            windowing_norm[key] = None if val is None else str(val)
        else:
            windowing_norm[key] = val

    return {
        "continuous_feature_cols": list(feature_cols),
        "windowing": windowing_norm,
        "quantile_capping": {
            "lower_quantile": None if qmeta.get("lower_quantile", None) is None else float(qmeta.get("lower_quantile")),
            "upper_quantile": None if qmeta.get("upper_quantile", None) is None else float(qmeta.get("upper_quantile")),
            "feature_limits": feature_limits_norm,
        },
        "log_scaling": log_scaling_norm,
        "scaler": scaler_norm,
    }


def _build_classifier_preprocessing_compare_payload(classifier_artifact):
    """Build normalized preprocessing/windowing compare payload from classifier metadata."""
    if not isinstance(classifier_artifact, dict):
        return None
    feature_cols = classifier_artifact.get("continuous_feature_cols", None)
    if feature_cols is None:
        feature_cols = classifier_artifact.get("feature_cols", [])
    preprocessing_params = classifier_artifact.get("preprocessing", None)
    if not isinstance(preprocessing_params, dict):
        preprocessing_params = {
            "continuous_feature_cols": list(feature_cols),
            "quantile_capping": {"lower_quantile": None, "upper_quantile": None, "feature_limits": {}},
            "log_scaling": None,
            "scaler": None,
        }
    payload_feature_cols = preprocessing_params.get("continuous_feature_cols", feature_cols)
    windowing_params = classifier_artifact.get("windowing", {})
    return _normalize_preprocessing_params_for_compare(
        preprocessing_params=preprocessing_params,
        windowing_params=windowing_params,
        feature_cols=payload_feature_cols,
    )


def _build_adata_preprocessing_compare_payload(adata):
    """Build normalized preprocessing/windowing compare payload from adata.uns."""
    if not (hasattr(adata, "uns") and isinstance(adata.uns, dict)):
        return None
    pre_meta = adata.uns.get("preprocessing", {})
    pre_meta = pre_meta if isinstance(pre_meta, dict) else {}
    preprocessing_params = _build_preprocessing_params_from_uns(
        pre_meta,
        feature_cols=list(getattr(adata, "var_names", [])),
    )
    feature_cols = (
        preprocessing_params.get("continuous_feature_cols", list(getattr(adata, "var_names", [])))
        if isinstance(preprocessing_params, dict)
        else list(getattr(adata, "var_names", []))
    )
    return _normalize_preprocessing_params_for_compare(
        preprocessing_params=preprocessing_params,
        windowing_params=pre_meta.get("windowing", {}),
        feature_cols=feature_cols,
    )


def _compare_preprocessing_params(payload_a, payload_b, rtol=1e-8, atol=1e-12):
    """Compare two normalized preprocessing/windowing payloads and return (is_match, reasons)."""
    reasons = []
    if not isinstance(payload_a, dict):
        return False, ["left preprocessing payload is missing/invalid"]
    if not isinstance(payload_b, dict):
        return False, ["right preprocessing payload is missing/invalid"]

    cols_a = list(payload_a.get("continuous_feature_cols", []))
    cols_b = list(payload_b.get("continuous_feature_cols", []))
    if cols_a != cols_b:
        reasons.append("continuous_feature_cols mismatch")

    win_a = payload_a.get("windowing", {})
    win_b = payload_b.get("windowing", {})
    for key in [
        "features",
        "binary_features_to_group",
        "window_size",
        "min_spacing",
        "descriptive_features",
        "incomplete_window_policy",
        "lower_quantile_cap",
        "upper_quantile_cap",
    ]:
        if win_a.get(key, None) != win_b.get(key, None):
            reasons.append(f"windowing.{key} mismatch")

    q_a = payload_a.get("quantile_capping", {})
    q_b = payload_b.get("quantile_capping", {})
    for key in ["lower_quantile", "upper_quantile"]:
        va = q_a.get(key, None)
        vb = q_b.get(key, None)
        if va is None and vb is None:
            continue
        if (va is None) != (vb is None):
            reasons.append(f"quantile_capping.{key} mismatch")
            continue
        if not np.isclose(float(va), float(vb), rtol=rtol, atol=atol):
            reasons.append(f"quantile_capping.{key} mismatch")

    fl_a = q_a.get("feature_limits", {})
    fl_b = q_b.get("feature_limits", {})
    feats_a = sorted([str(k) for k in fl_a.keys()])
    feats_b = sorted([str(k) for k in fl_b.keys()])
    if feats_a != feats_b:
        reasons.append("quantile_capping.feature_limits keys mismatch")
    else:
        for feat in feats_a:
            lim_a = fl_a.get(feat, {})
            lim_b = fl_b.get(feat, {})
            for bound in ["lower", "upper"]:
                va = lim_a.get(bound, None)
                vb = lim_b.get(bound, None)
                if va is None and vb is None:
                    continue
                if (va is None) != (vb is None):
                    reasons.append(f"feature_limits.{feat}.{bound} mismatch")
                    continue
                if not np.isclose(float(va), float(vb), rtol=rtol, atol=atol):
                    reasons.append(f"feature_limits.{feat}.{bound} mismatch")

    log_a = payload_a.get("log_scaling", None)
    log_b = payload_b.get("log_scaling", None)
    if (log_a is None) != (log_b is None):
        reasons.append("log_scaling presence mismatch")
    elif isinstance(log_a, dict) and isinstance(log_b, dict):
        if log_a.get("transform", None) != log_b.get("transform", None):
            reasons.append("log_scaling.transform mismatch")
        for key in ["requested_features", "resolved_feature_cols"]:
            vals_a = [str(x) for x in list(log_a.get(key, []))]
            vals_b = [str(x) for x in list(log_b.get(key, []))]
            if vals_a != vals_b:
                reasons.append(f"log_scaling.{key} mismatch")

    sc_a = payload_a.get("scaler", None)
    sc_b = payload_b.get("scaler", None)
    if (sc_a is None) != (sc_b is None):
        reasons.append("scaler presence mismatch")
    elif isinstance(sc_a, dict) and isinstance(sc_b, dict):
        for key in ["mean", "scale"]:
            arr_a = np.asarray(sc_a.get(key, []), dtype=float)
            arr_b = np.asarray(sc_b.get(key, []), dtype=float)
            if arr_a.shape != arr_b.shape:
                reasons.append(f"scaler.{key} length mismatch")
                continue
            if not np.allclose(arr_a, arr_b, rtol=rtol, atol=atol):
                reasons.append(f"scaler.{key} values mismatch")

    return len(reasons) == 0, reasons


def _matches_requested_preprocessing_in_adata(
    adata,
    *,
    features,
    binary_features_to_group,
    window_size,
    min_spacing,
    descriptive_features,
    incomplete_window_policy,
    lower_quantile_cap,
    upper_quantile_cap,
    log_scale_features,
    scale_features,
):
    """Check whether adata.uns preprocessing matches requested preprocessing settings.

    Note:
        `min_spacing` is intentionally excluded from preprocessing matching because it
        belongs to model sampling stage, not full window preprocessing stage.
    """
    reasons = []
    if not (hasattr(adata, "uns") and isinstance(adata.uns, dict)):
        return False, ["adata.uns is missing/invalid"]
    pre_meta = adata.uns.get("preprocessing", {})
    if not isinstance(pre_meta, dict):
        return False, ["adata.uns['preprocessing'] is missing/invalid"]

    window = pre_meta.get("windowing", {})
    if not isinstance(window, dict):
        reasons.append("preprocessing.windowing missing/invalid")
        window = {}
    expected_window = {
        "features": [str(x) for x in list(features)],
        "binary_features_to_group": [str(x) for x in list(binary_features_to_group)],
        "window_size": int(window_size),
        "descriptive_features": [str(x) for x in list(descriptive_features)],
        "incomplete_window_policy": str(incomplete_window_policy),
        "lower_quantile_cap": None if lower_quantile_cap is None else float(lower_quantile_cap),
        "upper_quantile_cap": None if upper_quantile_cap is None else float(upper_quantile_cap),
    }
    for key, expected_val in expected_window.items():
        got = window.get(key, None)
        if key in {"features", "binary_features_to_group", "descriptive_features"}:
            got = [] if got is None else [str(x) for x in list(got)]
        elif key in {"window_size"}:
            got = None if got is None else int(got)
        elif key in {"lower_quantile_cap", "upper_quantile_cap"}:
            got = None if got is None else float(got)
        elif key == "incomplete_window_policy":
            got = None if got is None else str(got)
        if got != expected_val:
            reasons.append(f"windowing.{key} mismatch")

    qmeta = pre_meta.get("quantile_capping", {})
    if not isinstance(qmeta, dict):
        reasons.append("preprocessing.quantile_capping missing/invalid")
        qmeta = {}
    got_lower = qmeta.get("lower_quantile", None)
    got_upper = qmeta.get("upper_quantile", None)
    got_lower = None if got_lower is None else float(got_lower)
    got_upper = None if got_upper is None else float(got_upper)
    exp_lower = None if lower_quantile_cap is None else float(lower_quantile_cap)
    exp_upper = None if upper_quantile_cap is None else float(upper_quantile_cap)
    if got_lower != exp_lower:
        reasons.append("quantile_capping.lower_quantile mismatch")
    if got_upper != exp_upper:
        reasons.append("quantile_capping.upper_quantile mismatch")

    expected_log = _resolve_log_scale_feature_cols(
        getattr(adata, "var_names", []),
        log_scale_features,
    )
    if len(expected_log["unresolved_features"]) > 0:
        reasons.append("log_scaling.unresolved_features")
    got_log = _coerce_log_scaling_params(pre_meta)
    got_log_requested = [] if got_log is None else list(got_log.get("requested_features", []))
    got_log_resolved = [] if got_log is None else list(got_log.get("resolved_feature_cols", []))
    if got_log_requested != list(expected_log["requested_features"]):
        reasons.append("log_scaling.requested_features mismatch")
    if got_log_resolved != list(expected_log["resolved_feature_cols"]):
        reasons.append("log_scaling.resolved_feature_cols mismatch")
    got_transform = None if got_log is None else got_log.get("transform", None)
    exp_transform = None if len(expected_log["resolved_feature_cols"]) == 0 else "log1p"
    if got_transform != exp_transform:
        reasons.append("log_scaling.transform mismatch")

    scaler_meta = pre_meta.get("scaler", None)
    has_scaler = isinstance(scaler_meta, dict) and ("mean" in scaler_meta) and ("scale" in scaler_meta)
    if bool(scale_features) != bool(has_scaler):
        reasons.append("scaler presence mismatch")

    kept_features_meta = pre_meta.get("kept_features", None)
    if isinstance(kept_features_meta, list):
        kept_features_meta = [str(x) for x in kept_features_meta]
        var_names = [str(x) for x in list(getattr(adata, "var_names", []))]
        if kept_features_meta != var_names:
            reasons.append("preprocessing.kept_features mismatch with var_names")

    reasons.extend(
        _descriptor_column_mismatch_reasons(
            getattr(adata, "var_names", []),
            features=features,
            descriptive_features=descriptive_features,
        )
    )

    return len(reasons) == 0, reasons














def _apply_preprocessing_to_continuous_matrix(
    X,
    preprocessing_params,
    feature_cols,
    inplace=False,
):
    """Apply saved quantile caps and scaling to a continuous feature matrix."""
    out = np.asarray(X, dtype=float)
    if not inplace:
        out = out.copy()
    feature_cols = list(feature_cols)
    qmeta = (preprocessing_params or {}).get("quantile_capping", {})
    limits = qmeta.get("feature_limits", {}) if isinstance(qmeta, dict) else {}
    if isinstance(limits, dict) and len(limits) > 0:
        for i, feat in enumerate(feature_cols):
            lim = limits.get(feat, None)
            if lim is None:
                continue
            lower = lim.get("lower", None)
            upper = lim.get("upper", None)
            out[:, i] = np.clip(out[:, i], a_min=lower, a_max=upper)

    out = _apply_log_scaling_to_continuous_matrix(
        out,
        preprocessing_params=preprocessing_params,
        feature_cols=feature_cols,
        inplace=True,
    )

    scaler_meta = _coerce_scaler_params(preprocessing_params)
    if scaler_meta is not None:
        mean = np.asarray(scaler_meta["mean"], dtype=float)
        scale = np.asarray(scaler_meta["scale"], dtype=float)
        if out.shape[1] != len(mean) or out.shape[1] != len(scale):
            raise ValueError(
                "Preprocessing scaler parameter size mismatch: "
                f"n_features={out.shape[1]}, len(mean)={len(mean)}, len(scale)={len(scale)}"
            )
        out -= mean
        out /= scale
    return out


def load_state_classifier_artifact(path):
    """Load a saved state classification classifier artifact pickle."""
    with open(path, "rb") as f:
        artifact = pickle.load(f)
    if not isinstance(artifact, dict) or "classifier" not in artifact:
        raise ValueError("Invalid classifier artifact: expected dict with key 'classifier'.")
    return artifact


def _extract_pipeline_metadata_from_classifier_artifact(classifier_artifact, artifact_name):
    """Return strict pipeline metadata payload from a classifier artifact."""
    if not isinstance(classifier_artifact, dict) or "classifier" not in classifier_artifact:
        raise ValueError(
            f"{artifact_name} is not a valid classifier artifact dict."
        )
    pipeline_meta = classifier_artifact.get("pipeline_metadata", None)
    if not isinstance(pipeline_meta, dict):
        raise ValueError(
            f"{artifact_name} is missing required 'pipeline_metadata'. "
            "Re-train classifiers with the updated pipeline."
        )

    pre_meta = pipeline_meta.get("preprocessing", None)
    clust_meta = pipeline_meta.get("clustering", None)
    class_meta = pipeline_meta.get("classification", None)
    if not isinstance(pre_meta, dict):
        raise ValueError(
            f"{artifact_name} pipeline_metadata.preprocessing is missing/invalid. "
            "Re-train classifiers with the updated pipeline."
        )
    if not isinstance(clust_meta, dict):
        raise ValueError(
            f"{artifact_name} pipeline_metadata.clustering is missing/invalid. "
            "Re-train classifiers with the updated pipeline."
        )
    if not isinstance(class_meta, dict):
        raise ValueError(
            f"{artifact_name} pipeline_metadata.classification is missing/invalid. "
            "Re-train classifiers with the updated pipeline."
        )

    return {
        "preprocessing": copy.deepcopy(pre_meta),
        "clustering": copy.deepcopy(clust_meta),
        "classification": copy.deepcopy(class_meta),
    }


def _validate_classifier_artifact_metadata_consistency(
    intrinsic_artifact,
    full_artifact,
):
    """Validate strict metadata consistency when both intrinsic/full artifacts are provided."""
    intrinsic_meta = _extract_pipeline_metadata_from_classifier_artifact(
        intrinsic_artifact,
        artifact_name="continuous_label_classifier_artifact",
    )
    full_meta = _extract_pipeline_metadata_from_classifier_artifact(
        full_artifact,
        artifact_name="full_label_classifier_artifact",
    )

    payload_intrinsic = _build_classifier_preprocessing_compare_payload(intrinsic_artifact)
    payload_full = _build_classifier_preprocessing_compare_payload(full_artifact)
    is_match, mismatch_reasons = _compare_preprocessing_params(
        payload_intrinsic,
        payload_full,
        rtol=1e-8,
        atol=1e-12,
    )
    if not is_match:
        raise ValueError(
            "Intrinsic and full classifier artifacts were trained with different preprocessing/windowing settings: "
            f"{mismatch_reasons}. Re-train classifiers from the same model run."
        )

    intrinsic_bin_cols = intrinsic_meta["clustering"].get("binary_cols_to_merge", None)
    full_bin_cols = full_meta["clustering"].get("binary_cols_to_merge", None)
    if not isinstance(intrinsic_bin_cols, list):
        raise ValueError(
            "continuous_label_classifier_artifact pipeline_metadata.clustering.binary_cols_to_merge is missing/invalid. "
            "Re-train classifiers with the updated pipeline."
        )
    if not isinstance(full_bin_cols, list):
        raise ValueError(
            "full_label_classifier_artifact pipeline_metadata.clustering.binary_cols_to_merge is missing/invalid. "
            "Re-train classifiers with the updated pipeline."
        )
    intrinsic_bin_cols = [str(c) for c in intrinsic_bin_cols]
    full_bin_cols = [str(c) for c in full_bin_cols]
    if intrinsic_bin_cols != full_bin_cols:
        raise ValueError(
            "Intrinsic and full classifier artifacts have different binary_cols_to_merge metadata. "
            f"intrinsic={intrinsic_bin_cols}, full={full_bin_cols}. "
            "Re-train classifiers from the same model run."
        )


def _coerce_classifier_artifact_input(artifact_input, artifact_name):
    if artifact_input is None:
        return None
    if isinstance(artifact_input, (str, Path)):
        return load_state_classifier_artifact(artifact_input)
    if isinstance(artifact_input, dict) and "classifier" in artifact_input:
        return artifact_input
    raise ValueError(
        f"{artifact_name} must be a classifier artifact dict or a path to a pickled classifier artifact, "
        f"got {type(artifact_input)}."
    )


def _resolve_classifier_artifact_variant_input(artifact_input, variant_name, artifact_name):
    """Resolve explicit artifact input (path/artifact/variant-dict) into artifact + source metadata."""
    if artifact_input is None:
        return None, None, None
    if isinstance(artifact_input, (str, Path)):
        return (
            _coerce_classifier_artifact_input(artifact_input, artifact_name),
            "explicit_path",
            str(Path(artifact_input)),
        )
    if isinstance(artifact_input, dict) and "classifier" in artifact_input:
        return artifact_input, "explicit_artifact", None
    if isinstance(artifact_input, dict):
        available = sorted([str(k) for k, v in artifact_input.items() if v is not None])
        if variant_name in artifact_input and artifact_input[variant_name] is not None:
            candidate = artifact_input[variant_name]
            source_path = str(Path(candidate)) if isinstance(candidate, (str, Path)) else None
            return (
                _coerce_classifier_artifact_input(candidate, f"{artifact_name} variant '{variant_name}'"),
                "explicit_variant",
                source_path,
            )
        raise ValueError(
            f"Requested {artifact_name} variant '{variant_name}' not available. "
            f"Available variants: {available}"
        )
    raise ValueError(
        f"{artifact_name} must be a classifier artifact dict, a path, or a variant dict "
        f"(keys: 'full_data'/'train_split' with artifact/path values), got {type(artifact_input)}."
    )


def _resolve_classifier_artifacts_with_defaults(
    *,
    label_classifier_artifact,
    full_label_classifier_artifact,
    continuous_model_variant,
    full_model_variant,
    state_paths,
    verbose=True,
):
    """Resolve intrinsic/full artifacts via explicit input first, then canonical default paths."""

    def _load_default_artifact_if_missing(current_artifact, default_path):
        if current_artifact is not None:
            return current_artifact, None, None
        default_path = Path(default_path)
        if default_path.exists():
            loaded = load_state_classifier_artifact(default_path)
            return loaded, "default_path", str(default_path)
        return None, None, None

    label_classifier_selected, label_source_type, label_source_path = _resolve_classifier_artifact_variant_input(
        label_classifier_artifact,
        continuous_model_variant,
        "continuous_label_classifier_artifact",
    )
    full_label_classifier_selected, full_source_type, full_source_path = _resolve_classifier_artifact_variant_input(
        full_label_classifier_artifact,
        full_model_variant,
        "full_label_classifier_artifact",
    )

    if label_classifier_selected is None:
        label_classifier_selected, default_type, default_path = _load_default_artifact_if_missing(
            label_classifier_selected,
            state_paths.intrinsic_classifier_default_path,
        )
        if default_type is not None:
            label_source_type = default_type
            label_source_path = default_path
            _vinfo(verbose, "state-apply", f"using default intrinsic classifier path: {default_path}")
        else:
            _vinfo(
                verbose,
                "state-apply",
                f"no intrinsic classifier provided and default path not found: {state_paths.intrinsic_classifier_default_path}",
            )

    if full_label_classifier_selected is None:
        full_label_classifier_selected, default_type, default_path = _load_default_artifact_if_missing(
            full_label_classifier_selected,
            state_paths.full_classifier_default_path,
        )
        if default_type is not None:
            full_source_type = default_type
            full_source_path = default_path
            _vinfo(verbose, "state-apply", f"using default full classifier path: {default_path}")
        else:
            _vinfo(
                verbose,
                "state-apply",
                f"no full classifier provided and default path not found: {state_paths.full_classifier_default_path}",
            )

    if label_classifier_selected is None and full_label_classifier_selected is None:
        raise FileNotFoundError(
            "No classifier artifact could be resolved for apply stage. "
            f"Checked intrinsic default path='{state_paths.intrinsic_classifier_default_path}' "
            f"and full default path='{state_paths.full_classifier_default_path}'. "
            "Pass label_classifier_artifact and/or full_label_classifier_artifact, "
            "or train/save classifiers first."
        )

    return {
        "label_classifier_selected": label_classifier_selected,
        "full_label_classifier_selected": full_label_classifier_selected,
        "label_source_type": label_source_type,
        "label_source_path": label_source_path,
        "full_source_type": full_source_type,
        "full_source_path": full_source_path,
    }


def apply_classifier(
    df_positions,
    classifier_artifact_or_path,
    output_col="intrinsic_behavioral_cluster",
    confidence_col="intrinsic_behavioral_cluster_confidence",
    outfolder=None,
    verbose=True,
):
    """
    Apply a saved classifier artifact to raw position data.

    This function:
    1) loads saved windowing/preprocessing/classifier metadata
    2) recreates the windowed descriptive dataset with the same parameters
    3) applies saved preprocessing + classifier
    4) returns predicted AnnData
    """
    artifact = (
        load_state_classifier_artifact(classifier_artifact_or_path)
        if isinstance(classifier_artifact_or_path, (str, Path))
        else classifier_artifact_or_path
    )
    if not isinstance(artifact, dict) or "classifier" not in artifact:
        raise ValueError("classifier_artifact_or_path must be a valid classifier artifact or path to one.")

    windowing = artifact.get("windowing", None)
    if not isinstance(windowing, dict):
        raise ValueError(
            "Classifier artifact is missing 'windowing' settings. "
            "Re-train/save with updated run_state_classification to enable auto window reconstruction."
        )

    required_window_keys = [
        "features",
        "binary_features_to_group",
        "window_size",
        "descriptive_features",
        "incomplete_window_policy",
    ]
    missing_window_keys = [k for k in required_window_keys if k not in windowing]
    if len(missing_window_keys) > 0:
        raise ValueError(f"Classifier artifact windowing settings are incomplete: missing {missing_window_keys}")

    _require_columns(
        df_positions,
        ["sample_name", "TrackID", "position_t"] + list(windowing["features"]),
        context="windowed dataset creation",
    )

    adata_prepared = prepare_state_classification_dataset(
        df_positions=df_positions,
        features=list(windowing["features"]),
        binary_features_to_group=list(windowing["binary_features_to_group"]),
        window_size=int(windowing["window_size"]),
        min_spacing=windowing.get("min_spacing", None),
        descriptive_features=list(windowing["descriptive_features"]),
        lower_quantile_cap=windowing.get("lower_quantile_cap", None),
        upper_quantile_cap=windowing.get("upper_quantile_cap", None),
        outfolder=outfolder,
        scale_features=False,
        incomplete_window_policy=str(windowing.get("incomplete_window_policy", None)),
        reuse_prepared_dataset=False,
        save_prepared_dataset=False,
        verbose=verbose,
    )
    pre_meta = (
        adata_prepared.uns.get("preprocessing", {})
        if isinstance(getattr(adata_prepared, "uns", {}), dict)
        else {}
    )
    default_non_feature_cols = [
        "sample_name",
        "TrackID",
        "sub_TrackID",
        "position_t",
        "window_start_position_t",
        "window_end_position_t",
        "window_length_frames",
    ]
    non_feature_cols = pre_meta.get("non_feature_cols", default_non_feature_cols)
    if not isinstance(non_feature_cols, list):
        non_feature_cols = list(default_non_feature_cols)
    binary_cols_to_merge = pre_meta.get("binary_cols_to_merge", None)
    if not isinstance(binary_cols_to_merge, list):
        binary_cols_to_merge = [
            c for c in list(windowing.get("binary_features_to_group", []))
            if c in adata_prepared.obs.columns
        ]

    if "continuous_feature_cols" in artifact:
        # Full-label classifier artifact path.
        cont_cols = list(artifact["continuous_feature_cols"])
    else:
        # Label-transfer classifier artifact path.
        cont_cols = list(artifact.get("feature_cols", []))
    if len(cont_cols) == 0:
        raise ValueError("Classifier artifact has no stored feature columns ('continuous_feature_cols'/'feature_cols').")

    missing_cont = [c for c in cont_cols if c not in adata_prepared.var_names]
    if len(missing_cont) > 0:
        raise ValueError(
            "Prepared dataset is missing required classifier features: "
            f"{missing_cont[:20]}"
        )

    obs_cols = [c for c in non_feature_cols if c in adata_prepared.obs.columns] + [
        c for c in binary_cols_to_merge if c in adata_prepared.obs.columns
    ]
    adata_query = adata_prepared[:, cont_cols].copy()
    adata_query.obs = adata_query.obs[obs_cols]
    del adata_prepared

    if "continuous_feature_cols" in artifact:
        predict_clusterids_with_full_classifier(
            adata=adata_query,
            classifier_artifact=artifact,
            output_col=output_col,
            confidence_col=confidence_col,
            inplace=True,
        )
    else:
        predict_clusterids_with_classifier(
            adata=adata_query,
            classifier=artifact,
            feature_cols=cont_cols,
            output_col=output_col,
            confidence_col=confidence_col,
            inplace=True,
        )

    return adata_query






def predict_clusterids_with_classifier(
    adata,
    classifier,
    feature_cols=None,
    output_col="intrinsic_behavioral_cluster",
    confidence_col=None,
    inplace=True,
    apply_preprocessing=True,
):
    """Predict cluster labels for an AnnData object using a trained classifier or artifact."""
    target = adata if inplace else adata.copy()
    clf = classifier
    preprocessing_params = None
    if isinstance(classifier, dict) and "classifier" in classifier:
        clf = classifier["classifier"]
        preprocessing_params = classifier.get("preprocessing", None)

    cols = list(target.var_names) if feature_cols is None else list(feature_cols)

    missing = [c for c in cols if c not in target.var_names]
    if len(missing) > 0:
        raise ValueError(f"adata is missing classifier feature columns: {missing[:10]}")

    if apply_preprocessing and preprocessing_params is not None:
        X_query = _to_numpy_2d(target[:, cols].X).astype(float, copy=True)
    else:
        X_query = _to_numpy_2d(target[:, cols].X).astype(float, copy=False)
    if apply_preprocessing and preprocessing_params is not None:
        X_query = _apply_preprocessing_to_continuous_matrix(
            X_query,
            preprocessing_params=preprocessing_params,
            feature_cols=cols,
            inplace=True,
        )
    pred = clf.predict(X_query)
    target.obs[output_col] = pd.Categorical(pd.Series(pred, index=target.obs.index).astype(str))

    if confidence_col is not None and hasattr(clf, "predict_proba"):
        proba = clf.predict_proba(X_query)
        target.obs[confidence_col] = pd.Series(proba.max(axis=1), index=target.obs.index, dtype=float)

    return target
















def _build_binary_feature_dataframe(
    obs_df,
    binary_feature_cols,
    expected_expanded_cols=None,
):
    """Build numeric binary/group feature frame with one-hot encoding for categoricals."""
    pieces = []
    for col in list(binary_feature_cols):
        if col not in obs_df.columns:
            raise ValueError(f"Missing binary/group feature column in adata.obs: '{col}'")
        s = obs_df[col]
        if is_numeric_dtype(s) or pd.api.types.is_bool_dtype(s):
            piece = pd.DataFrame(
                {col: pd.to_numeric(s, errors="coerce").fillna(0.0).astype(float)},
                index=obs_df.index,
            )
        else:
            piece = pd.get_dummies(
                s.astype("string").fillna("missing"),
                prefix=str(col),
                dtype=float,
            )
            piece.index = obs_df.index
        pieces.append(piece)

    if len(pieces) == 0:
        out = pd.DataFrame(index=obs_df.index)
    else:
        out = pd.concat(pieces, axis=1)

    if expected_expanded_cols is not None:
        expected = list(expected_expanded_cols)
        out = out.reindex(columns=expected, fill_value=0.0)
    return out






def predict_clusterids_with_full_classifier(
    adata,
    classifier_artifact,
    output_col="full_behavioral_cluster",
    confidence_col=None,
    inplace=True,
    apply_preprocessing=True,
):
    """Apply a trained full ClusterID classifier artifact to an AnnData object."""
    target = adata if inplace else adata.copy()
    clf = classifier_artifact["classifier"]
    cont_cols = classifier_artifact["continuous_feature_cols"]
    bin_cols = classifier_artifact.get("binary_feature_cols", [])
    bin_expanded_cols = classifier_artifact.get("binary_expanded_feature_cols", None)
    preprocessing_params = classifier_artifact.get("preprocessing", None)

    # Ensure continuous features use identical capping/z-normalization as training.
    if apply_preprocessing and preprocessing_params is not None:
        X_cont = _to_numpy_2d(target[:, cont_cols].X).astype(float, copy=True)
    else:
        X_cont = _to_numpy_2d(target[:, cont_cols].X).astype(float, copy=False)
    if apply_preprocessing and preprocessing_params is not None:
        X_cont = _apply_preprocessing_to_continuous_matrix(
            X_cont,
            preprocessing_params=preprocessing_params,
            feature_cols=cont_cols,
            inplace=True,
        )

    if len(bin_cols) == 0:
        X_query = X_cont
    else:
        bin_df = _build_binary_feature_dataframe(
            obs_df=target.obs,
            binary_feature_cols=bin_cols,
            expected_expanded_cols=bin_expanded_cols,
        )
        X_query = np.hstack([X_cont, bin_df.to_numpy(dtype=float)])
    pred = np.asarray(clf.predict(X_query)).astype(str)
    target.obs[output_col] = pd.Categorical(pd.Series(pred, index=target.obs.index).astype(str))

    if confidence_col is not None and hasattr(clf, "predict_proba"):
        proba = clf.predict_proba(X_query)
        target.obs[confidence_col] = pd.Series(proba.max(axis=1), index=target.obs.index, dtype=float)

    return target












def apply_state_classifiers_to_full_dataset(
    output_dir,
    cell_type,
    label_classifier_artifact=None,
    full_label_classifier_artifact=None,
    continuous_output_col="intrinsic_behavioral_cluster",
    full_output_col="full_behavioral_cluster",
    continuous_model_variant="full_data",
    full_model_variant="full_data",
    combine_binary_with_continuous=True,
    export_state_transitions=True,
    state_transitions_rows_per_page=3,
    state_transitions_include_self_pairs=True,
    state_transitions_min_count=100,
    state_transitions_relative_count=0.0,
    state_paths=None,
    verbose=True,
    export_behavioral_state_images=False,
    behavioral_state_image_col=None,
):
    """Apply trained stage-2 classifier artifacts to full dataset from track-features CSV.

    Confidence columns are always derived as `<output_col>_confidence`.
    If artifacts are omitted, canonical default classifier paths are probed.
    """
    apply_started = _vstart(verbose, "state-apply", "apply state classifiers")
    state_paths = state_paths or _resolve_state_paths(output_dir, cell_type)
    state_outdir = state_paths.state_outdir

    continuous_confidence_col = f"{str(continuous_output_col)}_confidence"
    full_confidence_col = f"{str(full_output_col)}_confidence"

    resolved_artifacts = _resolve_classifier_artifacts_with_defaults(
        label_classifier_artifact=label_classifier_artifact,
        full_label_classifier_artifact=full_label_classifier_artifact,
        continuous_model_variant=continuous_model_variant,
        full_model_variant=full_model_variant,
        state_paths=state_paths,
        verbose=verbose,
    )
    label_classifier_selected = resolved_artifacts["label_classifier_selected"]
    full_label_classifier_selected = resolved_artifacts["full_label_classifier_selected"]
    label_source_type = resolved_artifacts["label_source_type"]
    label_source_path = resolved_artifacts["label_source_path"]
    full_source_type = resolved_artifacts["full_source_type"]
    full_source_path = resolved_artifacts["full_source_path"]

    if (label_classifier_selected is not None) and (full_label_classifier_selected is not None):
        _validate_classifier_artifact_metadata_consistency(
            intrinsic_artifact=label_classifier_selected,
            full_artifact=full_label_classifier_selected,
        )

    reference_classifier_artifact = (
        label_classifier_selected if label_classifier_selected is not None else full_label_classifier_selected
    )
    reference_artifact_name = (
        "continuous_label_classifier_artifact"
        if label_classifier_selected is not None
        else "full_label_classifier_artifact"
    )
    reference_pipeline_meta = _extract_pipeline_metadata_from_classifier_artifact(
        reference_classifier_artifact,
        artifact_name=reference_artifact_name,
    )
    class_meta = dict(reference_pipeline_meta["classification"])
    clust_meta = dict(reference_pipeline_meta["clustering"])
    pre_meta = dict(reference_pipeline_meta["preprocessing"])

    _vinfo(
        verbose,
        "state-apply",
        (
            "Classifier variants: "
            f"continuous={continuous_model_variant if label_classifier_selected is not None else None}, "
            f"full={full_model_variant if full_label_classifier_selected is not None else None}"
        ),
    )
    _vinfo(
        verbose,
        "state-apply",
        (
            "Classifier sources: "
            f"intrinsic=({label_source_type}, {label_source_path}), "
            f"full=({full_source_type}, {full_source_path})"
        ),
    )

    full_predicted_in_primary = False
    adata_full = None
    positions_csv_path = None
    df_positions = None

    def _load_positions_if_needed():
        nonlocal df_positions, positions_csv_path
        if df_positions is None:
            if positions_csv_path is None:
                positions_csv_path = _resolve_positions_csv_path(state_paths=state_paths)
            df_positions = pd.read_csv(positions_csv_path, low_memory=False)
            df_positions = df_positions.sort_values(by=["sample_name", "TrackID", "position_t"])

    binary_cols_to_merge = clust_meta.get("binary_cols_to_merge", None)
    if not isinstance(binary_cols_to_merge, list):
        raise ValueError(
            "Classifier artifact pipeline metadata is missing 'clustering.binary_cols_to_merge'. "
            "Re-train classifiers with the updated pipeline."
        )
    binary_cols_to_merge = [str(c) for c in list(binary_cols_to_merge)]
    binary_group_constraints = clust_meta.get("binary_group_constraints", None)
    if (binary_group_constraints is not None) and (not isinstance(binary_group_constraints, dict)):
        raise ValueError(
            "Classifier artifact pipeline metadata has invalid 'clustering.binary_group_constraints'. "
            "Re-train classifiers with the updated pipeline."
        )
    enforce_binary_group_constraints = (
        isinstance(binary_group_constraints, dict)
        and ("forbidden_binary_combinations" in binary_group_constraints)
    )

    prepared_cache_path = state_paths.prepared_full_adata_path
    legacy_prepared_cache_path = _resolve_legacy_prepared_full_adata_path(
        state_paths=state_paths,
        cell_type=cell_type,
    )
    prepared_cache_read_path = None
    if prepared_cache_path is not None and prepared_cache_path.exists():
        prepared_cache_read_path = prepared_cache_path
    elif legacy_prepared_cache_path.exists():
        prepared_cache_read_path = legacy_prepared_cache_path

    if prepared_cache_read_path is not None:
        cache_load_started = _vstart(verbose, "state-apply", "load prepared full cache")
        cached_adata = sc.read_h5ad(prepared_cache_read_path)
        _vdone(verbose, "state-apply", "load prepared full cache", cache_load_started)
        try:
            cache_payload = _build_adata_preprocessing_compare_payload(cached_adata)
            classifier_payload = _build_classifier_preprocessing_compare_payload(reference_classifier_artifact)
            is_match, mismatch_reasons = _compare_preprocessing_params(
                cache_payload,
                classifier_payload,
                rtol=1e-8,
                atol=1e-12,
            )
            if is_match:
                adata_full = cached_adata
                _vinfo(
                    verbose,
                    "state-apply",
                    f"Using prepared full cache: {prepared_cache_read_path}",
                )
                if prepared_cache_read_path != prepared_cache_path:
                    _vinfo(
                        verbose,
                        "state-apply",
                        f"Legacy cache path detected; future writes use: {prepared_cache_path}",
                    )
            else:
                err_msg = (
                    "ERROR: cached prepared adata preprocessing mismatch; recomputing from input track-features CSV. "
                    f"reasons={list(mismatch_reasons)}"
                )
                warnings.warn(err_msg, RuntimeWarning)
                _vinfo(
                    verbose,
                    "state-apply",
                    (
                        "Prepared full cache mismatch; recomputing "
                        f"(n_reasons={len(mismatch_reasons)}, first_reason={_short_reasons(mismatch_reasons)})"
                    ),
                )
        except Exception as exc:
            warnings.warn(
                f"Could not validate cached prepared adata; recomputing from input track-features CSV ({exc}).",
                RuntimeWarning,
            )
            _vinfo(
                verbose,
                "state-apply",
                f"Prepared full cache validation failed; recomputing ({exc})",
            )
    else:
        _vinfo(
            verbose,
            "state-apply",
            (
                "Prepared full cache missing; recomputing "
                f"(expected_new={prepared_cache_path}, checked_legacy={legacy_prepared_cache_path})"
            ),
        )

    if adata_full is None:
        prep_started = _vstart(verbose, "state-apply", "prepare dataset for classifier apply")
        _load_positions_if_needed()
        if label_classifier_selected is not None:
            adata_full = apply_classifier(
                df_positions=df_positions,
                classifier_artifact_or_path=label_classifier_selected,
                output_col=continuous_output_col,
                confidence_col=continuous_confidence_col,
                outfolder=state_outdir,
                verbose=verbose,
            )
        else:
            adata_full = apply_classifier(
                df_positions=df_positions,
                classifier_artifact_or_path=full_label_classifier_selected,
                output_col=full_output_col,
                confidence_col=full_confidence_col,
                outfolder=state_outdir,
                verbose=verbose,
            )
            full_predicted_in_primary = True
        del df_positions
        _vdone(verbose, "state-apply", "prepare dataset for classifier apply", prep_started)
    else:
        if label_classifier_selected is not None:
            cont_cols = list(
                label_classifier_selected.get(
                    "feature_cols",
                    label_classifier_selected.get("continuous_feature_cols", []),
                )
            )
            missing_cont = [c for c in cont_cols if c not in adata_full.var_names]
            if len(missing_cont) > 0:
                raise ValueError(
                    "Cached prepared adata is missing intrinsic classifier features: "
                    f"{missing_cont[:20]}"
                )
            predict_clusterids_with_classifier(
                adata=adata_full,
                classifier=label_classifier_selected,
                feature_cols=cont_cols,
                output_col=continuous_output_col,
                confidence_col=continuous_confidence_col,
                inplace=True,
                apply_preprocessing=False,
            )
        elif full_label_classifier_selected is not None:
            predict_clusterids_with_full_classifier(
                adata=adata_full,
                classifier_artifact=full_label_classifier_selected,
                output_col=full_output_col,
                confidence_col=full_confidence_col,
                inplace=True,
                apply_preprocessing=False,
            )
            full_predicted_in_primary = True

    if combine_binary_with_continuous and (label_classifier_selected is not None):
        if continuous_output_col != "intrinsic_behavioral_cluster":
            adata_full.obs["intrinsic_behavioral_cluster"] = adata_full.obs[continuous_output_col].astype("category")
        _rebuild_full_behavioral_cluster_from_intrinsic(
            adata=adata_full,
            binary_cols_to_merge=binary_cols_to_merge,
            intrinsic_col="intrinsic_behavioral_cluster",
            binary_group_constraints=binary_group_constraints,
            enforce_binary_group_constraints=enforce_binary_group_constraints,
        )

    if full_label_classifier_selected is not None and (not full_predicted_in_primary):
        predict_clusterids_with_full_classifier(
            adata=adata_full,
            classifier_artifact=full_label_classifier_selected,
            output_col=full_output_col,
            confidence_col=full_confidence_col,
            inplace=True,
            apply_preprocessing=False,
        )

    adata_full.uns["preprocessing"] = dict(pre_meta)
    adata_full.uns["classification"] = dict(class_meta)
    adata_full.uns["classification"]["applied_continuous_output_col"] = str(continuous_output_col)
    adata_full.uns["classification"]["applied_continuous_model_variant"] = (
        None if label_classifier_selected is None else str(continuous_model_variant)
    )
    adata_full.uns["classification"]["applied_full_output_col"] = (
        None if full_label_classifier_selected is None else str(full_output_col)
    )
    adata_full.uns["classification"]["applied_full_model_variant"] = (
        None if full_label_classifier_selected is None else str(full_model_variant)
    )

    state_composition_report_pdf = None
    state_composition_report_pdfs = []
    state_composition_report_auc_csv = None
    state_composition_report_plot_csvs = []
    state_composition_report_error = None

    if full_output_col in adata_full.obs.columns:
        report_dir = state_paths.state_composition_outdir
        report_dir.mkdir(parents=True, exist_ok=True)
        report_token = _sanitize_filename_token(
            full_output_col,
            fallback="full_behavioral_cluster",
        )
        report_pdf_path = report_dir / f"state_composition_report_{report_token}.pdf"
        report_auc_csv_path = report_dir / f"state_composition_report_{report_token}.csv"
        try:
            report_out = save_state_composition_report(
                adata=adata_full,
                output_pdf_path=report_pdf_path,
                output_csv_path=report_auc_csv_path,
                time_col="position_t",
                state_col=str(full_output_col),
                sample_col="sample_name",
                include_pooled_summary=True,
                verbose=verbose,
            )
            pdf_paths = report_out.get("pdf_paths", {})
            if isinstance(pdf_paths, dict) and len(pdf_paths) > 0:
                state_composition_report_pdfs = [str(v) for v in pdf_paths.values()]
                state_composition_report_pdf = str(
                    pdf_paths.get(
                        "stacked_by_sample",
                        state_composition_report_pdfs[0],
                    )
                )
            else:
                state_composition_report_pdf = str(report_out.get("pdf_path", report_pdf_path))
                state_composition_report_pdfs = [state_composition_report_pdf]
            state_composition_report_auc_csv = str(report_out.get("csv_path", report_auc_csv_path))
            plot_csv_paths = report_out.get("plot_data_csv_paths", {})
            if isinstance(plot_csv_paths, dict) and len(plot_csv_paths) > 0:
                state_composition_report_plot_csvs = [str(v) for v in plot_csv_paths.values()]
            else:
                plot_csv_path = report_out.get("plot_data_csv_path", None)
                if plot_csv_path is not None:
                    state_composition_report_plot_csvs = [str(plot_csv_path)]
        except Exception as exc:
            state_composition_report_error = str(exc)
            warnings.warn(
                "Could not generate state composition report after classifier apply: "
                f"{exc}",
                RuntimeWarning,
            )
            _vinfo(verbose, "state-apply", f"State composition report skipped: {exc}")
    else:
        state_composition_report_error = (
            "Could not generate state composition report: "
            f"missing '{full_output_col}' in adata_full.obs."
        )
        warnings.warn(state_composition_report_error, RuntimeWarning)

    adata_full.uns["classification"]["state_composition_report_pdf"] = state_composition_report_pdf
    adata_full.uns["classification"]["state_composition_report_pdfs"] = list(state_composition_report_pdfs)
    adata_full.uns["classification"]["state_composition_report_auc_csv"] = state_composition_report_auc_csv
    adata_full.uns["classification"]["state_composition_report_plot_csvs"] = list(
        state_composition_report_plot_csvs
    )
    adata_full.uns["classification"]["state_composition_report_error"] = state_composition_report_error
    if state_composition_report_pdf is not None:
        _vsave(verbose, "state-apply", "state composition report pdf", state_composition_report_pdf)
    if state_composition_report_auc_csv is not None:
        _vsave(verbose, "state-apply", "state composition report csv", state_composition_report_auc_csv)

    state_transition_report_dir = None
    state_transition_matrix_counts_csv = None
    state_transition_matrix_probs_csv = None
    state_transition_matrix_counts_no_self_csv = None
    state_transition_matrix_probs_no_self_csv = None
    state_transition_ngrams_all_csv = None
    state_transition_ngrams_pooled_csv = None
    state_transition_ngrams_per_end_csv = None
    state_transition_ngrams_per_start_csv = None
    state_transition_sankey_pdf = None
    state_transition_sankey_html_dir = None
    state_transition_report_error = None

    if bool(export_state_transitions):
        if full_output_col in adata_full.obs.columns:
            try:
                transitions_outdir = state_paths.state_transitions_outdir
                transitions_outdir.mkdir(parents=True, exist_ok=True)
                transition_out = save_state_transition_report(
                    adata=adata_full,
                    output_dir=transitions_outdir,
                    state_col=str(full_output_col),
                    id_cols=("sample_name", "TrackID"),
                    time_col="position_t",
                    include_self_pairs=bool(state_transitions_include_self_pairs),
                    rows_per_page=max(1, int(state_transitions_rows_per_page)),
                    sankey_min_count=int(state_transitions_min_count),
                    sankey_relative_count=float(state_transitions_relative_count),
                    verbose=verbose,
                )
                state_transition_report_dir = str(transition_out.get("output_dir", transitions_outdir))
                state_transition_matrix_counts_csv = transition_out.get("transition_matrix_counts_csv", None)
                state_transition_matrix_probs_csv = transition_out.get("transition_matrix_probs_csv", None)
                state_transition_matrix_counts_no_self_csv = transition_out.get(
                    "transition_matrix_counts_no_self_csv", None
                )
                state_transition_matrix_probs_no_self_csv = transition_out.get(
                    "transition_matrix_probs_no_self_csv", None
                )
                state_transition_ngrams_all_csv = transition_out.get("transition_ngrams_all_csv", None)
                state_transition_ngrams_pooled_csv = transition_out.get("transition_ngrams_pooled_csv", None)
                state_transition_ngrams_per_end_csv = transition_out.get("transition_ngrams_per_end_csv", None)
                state_transition_ngrams_per_start_csv = transition_out.get(
                    "transition_ngrams_per_start_csv", None
                )
                state_transition_sankey_pdf = transition_out.get("sankey_pdf_path", None)
                state_transition_sankey_html_dir = transition_out.get("sankey_html_dir", None)
            except Exception as exc:
                state_transition_report_error = str(exc)
                warnings.warn(
                    "Could not generate state transition report after classifier apply: "
                    f"{exc}",
                    RuntimeWarning,
                )
                _vinfo(verbose, "state-apply", f"State transition report skipped: {exc}")
        else:
            state_transition_report_error = (
                "Could not generate state transition report: "
                f"missing '{full_output_col}' in adata_full.obs."
            )
            warnings.warn(state_transition_report_error, RuntimeWarning)

    adata_full.uns["classification"]["state_transition_report_dir"] = state_transition_report_dir
    adata_full.uns["classification"]["state_transition_matrix_counts_csv"] = state_transition_matrix_counts_csv
    adata_full.uns["classification"]["state_transition_matrix_probs_csv"] = state_transition_matrix_probs_csv
    adata_full.uns["classification"]["state_transition_matrix_counts_no_self_csv"] = (
        state_transition_matrix_counts_no_self_csv
    )
    adata_full.uns["classification"]["state_transition_matrix_probs_no_self_csv"] = (
        state_transition_matrix_probs_no_self_csv
    )
    adata_full.uns["classification"]["state_transition_ngrams_all_csv"] = state_transition_ngrams_all_csv
    adata_full.uns["classification"]["state_transition_ngrams_pooled_csv"] = (
        state_transition_ngrams_pooled_csv
    )
    adata_full.uns["classification"]["state_transition_ngrams_per_end_csv"] = (
        state_transition_ngrams_per_end_csv
    )
    adata_full.uns["classification"]["state_transition_ngrams_per_start_csv"] = (
        state_transition_ngrams_per_start_csv
    )
    adata_full.uns["classification"]["state_transition_sankey_pdf"] = state_transition_sankey_pdf
    adata_full.uns["classification"]["state_transition_sankey_html_dir"] = state_transition_sankey_html_dir
    adata_full.uns["classification"]["state_transition_report_error"] = state_transition_report_error
    if state_transition_report_dir is not None:
        _vsave(verbose, "state-apply", "state transition report dir", state_transition_report_dir)
    if state_transition_sankey_pdf is not None:
        _vsave(verbose, "state-apply", "state transition sankey pdf", state_transition_sankey_pdf)

    behavioral_state_backprojection = {
        "enabled": bool(export_behavioral_state_images),
        "state_col": None,
        "output_paths": {},
        "skipped_samples": [],
        "errors": [],
    }
    if bool(export_behavioral_state_images):
        backprojection_started = _vstart(verbose, "state-apply", "export behavioral-state backprojection")
        export_state_col = (
            str(behavioral_state_image_col)
            if behavioral_state_image_col is not None
            else str(full_output_col)
        )
        behavioral_state_backprojection = export_behavioral_state_backprojection_zarrs(
            adata=adata_full,
            output_dir=state_paths.output_dir,
            cell_type=cell_type,
            state_col=export_state_col,
            sample_col="sample_name",
            track_col="TrackID",
            time_col="position_t",
            background_value=0,
            verbose=verbose,
        )
        behavioral_state_backprojection["enabled"] = True
        behavioral_state_backprojection["state_col"] = export_state_col
        _vdone(verbose, "state-apply", "export behavioral-state backprojection", backprojection_started)
    adata_full.uns["classification"]["behavioral_state_backprojection"] = behavioral_state_backprojection

    write_started = _vstart(verbose, "state-apply", "write classified full adata")
    adata_full.write(state_paths.full_output_adata_path, compression="gzip")
    _vdone(verbose, "state-apply", "write classified full adata", write_started)
    _vsave(verbose, "state-apply", "classified full adata", state_paths.full_output_adata_path)
    n_rows = int(adata_full.n_obs)
    n_intrinsic = int(
        adata_full.obs["intrinsic_behavioral_cluster"].astype(str).nunique()
    ) if "intrinsic_behavioral_cluster" in adata_full.obs.columns else 0
    n_full = int(
        adata_full.obs["full_behavioral_cluster"].astype(str).nunique()
    ) if "full_behavioral_cluster" in adata_full.obs.columns else 0
    _vinfo(
        verbose,
        "state-apply",
        f"Apply summary: rows={n_rows}, intrinsic_clusters={n_intrinsic}, full_behavioral_clusters={n_full}",
    )
    _vdone(verbose, "state-apply", "apply state classifiers", apply_started)

    return adata_full


 
if __name__ == "__main__":
    """
    Return lightweight defaults for quick test runs.

    This helper is optional and not used by the normal pipeline.
    """
    from behav3d.core.metadata import load_behav3d_metadata
    ssd_dir = r"/Volumes/T7_Sam/"
    ssd_dir = Path(ssd_dir)
    output_dir = Path(ssd_dir, r"BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE")
    metadata_csv_path = Path(ssd_dir, r"BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE/metadata.csv")
    # metadata_csv_path = Path(ssd_dir, r"BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE/metadata_home.csv")
    metadata = load_behav3d_metadata(metadata_csv_path)
    output_dir = Path("/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE")
    cell_type = "tcell"
    window_size = 5
    max_samples = None
    min_spacing = 10
    n_neighbors = 60
    resolution = 0.2
    features=[
        "percentage_dead_mask",
        # "mean_dead_dye",
        # "nr_dead_mask_pixels",
        # "organoid_contact_pixels",
        # "tcell_contact_pixels",
        # "mean_square_displacement",
        "speed",
        # "directional_persistence",
        "extent",
        "elongation",
        "sphericity",
        "solidity",
        # "oblateness",
        # "prolateness"
        # "surface_to_volume_ratio"
    ]
    
    # Define binary features to use for grouping
    binary_features_to_group = [
        "organoid_contact_pixels",
        "tcell_contact_pixels",
    ]
    
    descriptive_features = (
        "mean",
        "median",
        "std",
        "net_displacement",
        "straightness",
        "mean_square_displacement",
    )
    pca_var_selection = 0.95
    clustering_method = "leiden"
    resolutions = 0.2
    lower_quantile_cap = None
    upper_quantile_cap = 0.99
    incomplete_window_policy = "partial"
    random_state = 12345
    reuse_prepared_dataset = True
    save_prepared_dataset = True
    label_transfer_method = "classifier"
    cell_type = "tcell"
    classifier_backend = "random_forest"
    classifier_n_estimators = 300
    classifier_min_samples_leaf = 2
    classifier_n_jobs = -1
    classifier_max_depth = None
    classifier_min_samples_split = 2
    classifier_max_features = "sqrt"
    classifier_class_weight = None
    classifier_confidence_col = "intrinsic_behavioral_cluster_confidence"
    save_label_classifier = True
    
    
    train_full_label_classifier = True
    verbose=True
    
