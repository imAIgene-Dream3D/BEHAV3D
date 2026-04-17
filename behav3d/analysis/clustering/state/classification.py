import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas.api.types import is_numeric_dtype
import json
import copy
import warnings
import re
import time
import shutil
from itertools import combinations
from matplotlib.backends.backend_pdf import PdfPages
from dataclasses import asdict, dataclass

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
from behav3d.analysis.clustering.general import relabel_cluster_ids
from behav3d.analysis.clustering.general.leiden import run_pca, run_leiden_clustering
from behav3d.analysis.clustering.state.visualization.plots.state_composition import (
    save_state_composition_report,
)
from behav3d.analysis.clustering.state.visualization.plots.state_transitions import (
    save_state_transition_report,
)
from behav3d.analysis.clustering.state.visualization.backprojection import (
    export_behavioral_state_backprojection_zarrs,
    show_behavioral_state_backprojection,
)
from behav3d.analysis.clustering.state.visualization.videos.track_max_projection import (
    save_selected_fulltrack_cluster_videos,
)

A4_PORTRAIT = (8.27, 11.69)
A4_LANDSCAPE = (11.69, 8.27)
STATE_CLASSIFIER_PIPELINE_SCHEMA_VERSION = 3


def _vinfo(verbose, prefix, message):
    if not bool(verbose):
        return
    print(f"[{prefix}] INFO {str(message)}")


def _vstart(verbose, prefix, step_name):
    if bool(verbose):
        print(f"[{prefix}] START {step_name}")
    return time.perf_counter()


def _vdone(verbose, prefix, step_name, t_start):
    if bool(verbose):
        dt = time.perf_counter() - float(t_start)
        print(f"[{prefix}] DONE {step_name} | took {dt:.2f}s")


def _vsave(verbose, prefix, label, path):
    if not bool(verbose):
        return
    print(f"[{prefix}] SAVED {str(label)}: {path}")


def _short_reasons(reasons):
    if reasons is None:
        return "none"
    if not isinstance(reasons, (list, tuple)):
        return str(reasons)
    if len(reasons) == 0:
        return "none"
    return f"count={len(reasons)}, first={reasons[0]}"


def _binary_col_to_group_name(col: str) -> str:
    # Breaking change (schema v3): binary groups are exact selected column names.
    return str(col)


def _binary_value_is_active(val) -> bool:
    if pd.isna(val):
        return False
    try:
        return float(val) == 1.0
    except Exception:
        return False


def _active_binary_groups_from_row(row, binary_cols: list[str]) -> list[str]:
    active = []
    for col in binary_cols:
        if _binary_value_is_active(row[col]):
            active.append(_binary_col_to_group_name(col))
    return active


def _binary_groups_to_label(active_groups: list[str]) -> str:
    if len(active_groups) == 0:
        return "no_contact"
    if len(active_groups) == 1:
        return str(active_groups[0])
    return "_and_".join([str(x) for x in active_groups])


def _infer_binary_group_constraints(df: pd.DataFrame, binary_cols: list[str]) -> dict:
    """Infer allowed binary groups and forbidden binary-group combinations from train-time data."""
    cols = [str(c) for c in list(binary_cols or []) if str(c) in df.columns]
    clean_cols = [_binary_col_to_group_name(c) for c in cols]

    if len(cols) == 0:
        return {
            "inference_rule": "exact_zero_cooccurrence",
            "binary_cols": [],
            "binary_group_names": [],
            "allowed_binary_groups": ["no_contact"],
            "forbidden_binary_combinations": [],
            "support_counts": {},
        }

    support_counts = {name: 0 for name in clean_cols}
    cooccur_counts = {}
    allowed_groups = set()

    for _, row in df[cols].iterrows():
        active = _active_binary_groups_from_row(row, cols)
        allowed_groups.add(_binary_groups_to_label(active))
        for grp in active:
            support_counts[grp] = int(support_counts.get(grp, 0)) + 1
        for size in range(2, len(active) + 1):
            for combo in combinations(active, size):
                combo_key = tuple(sorted([str(x) for x in combo]))
                cooccur_counts[combo_key] = int(cooccur_counts.get(combo_key, 0)) + 1

    forbidden_combinations = []
    supported_groups = [str(name) for name, count in support_counts.items() if int(count) > 0]
    for size in range(2, len(supported_groups) + 1):
        for combo in combinations(supported_groups, size):
            combo_key = tuple(sorted([str(x) for x in combo]))
            if int(cooccur_counts.get(combo_key, 0)) == 0:
                forbidden_combinations.append(list(combo_key))

    return {
        "inference_rule": "exact_zero_cooccurrence",
        "binary_cols": list(cols),
        "binary_group_names": list(clean_cols),
        "allowed_binary_groups": sorted([str(x) for x in allowed_groups], key=_mixed_label_sort_key),
        "forbidden_binary_combinations": sorted(
            [[str(x) for x in combo] for combo in forbidden_combinations],
            key=lambda x: (len(x), tuple(x)),
        ),
        "support_counts": {str(k): int(v) for k, v in support_counts.items()},
    }


def _normalize_binary_group_constraints(binary_group_constraints):
    if binary_group_constraints is None:
        return None
    if not isinstance(binary_group_constraints, dict):
        raise ValueError(
            "binary_group_constraints must be a dict when provided."
        )

    allowed_raw = binary_group_constraints.get("allowed_binary_groups", None)
    forbidden_raw = binary_group_constraints.get("forbidden_binary_combinations", None)
    if allowed_raw is None and forbidden_raw is None:
        raise ValueError(
            "binary_group_constraints['forbidden_binary_combinations'] is required."
        )

    allowed_set = None
    if allowed_raw is not None:
        if not isinstance(allowed_raw, list):
            raise ValueError(
                "binary_group_constraints['allowed_binary_groups'] must be a list."
            )
        allowed_set = set([str(x) for x in allowed_raw])

    forbidden_set = set()
    if forbidden_raw is None:
        raise ValueError(
            "binary_group_constraints['forbidden_binary_combinations'] is required."
        )
    if not isinstance(forbidden_raw, list):
        raise ValueError(
            "binary_group_constraints['forbidden_binary_combinations'] must be a list."
        )
    for item in forbidden_raw:
        if not isinstance(item, (list, tuple)) or len(item) < 2:
            raise ValueError(
                "binary_group_constraints['forbidden_binary_combinations'] entries must be length-2+ lists/tuples."
            )
        combo = tuple(sorted([str(x) for x in item]))
        if len(set(combo)) != len(combo):
            raise ValueError(
                "binary_group_constraints['forbidden_binary_combinations'] entries must not contain duplicate group names."
            )
        forbidden_set.add(combo)

    return {
        "allowed_binary_groups": allowed_set,
        "forbidden_binary_combinations": forbidden_set,
    }


def _assign_binary_group_labels(
    df: pd.DataFrame,
    binary_cols: list[str],
    binary_group_constraints=None,
    enforce_binary_group_constraints: bool = False,
) -> pd.Series:
    """Build a single categorical group label from binary indicator columns."""
    if len(binary_cols) == 0:
        if enforce_binary_group_constraints and binary_group_constraints is not None:
            normalized = _normalize_binary_group_constraints(binary_group_constraints)
            if normalized is not None:
                allowed = normalized.get("allowed_binary_groups", None)
                if allowed is not None and "no_contact" not in allowed:
                    raise ValueError(
                        "Binary group constraints violated: 'no_contact' is required for empty binary_cols."
                    )
        return pd.Series(["no_contact"] * len(df), index=df.index, dtype="object")

    normalized_constraints = None
    if binary_group_constraints is not None:
        normalized_constraints = _normalize_binary_group_constraints(binary_group_constraints)

    labels = []
    forbidden_hits = {}
    unknown_hits = {}
    for _, row in df[binary_cols].iterrows():
        active = _active_binary_groups_from_row(row, binary_cols)
        label = _binary_groups_to_label(active)
        labels.append(label)

        if not enforce_binary_group_constraints or normalized_constraints is None:
            continue

        row_idx = str(row.name)
        forbidden_combinations = normalized_constraints.get("forbidden_binary_combinations", set())
        for size in range(2, len(active) + 1):
            for combo in combinations(active, size):
                combo_key = tuple(sorted([str(x) for x in combo]))
                if combo_key not in forbidden_combinations:
                    continue
                hit = forbidden_hits.setdefault(
                    combo_key,
                    {"count": 0, "first_index": row_idx},
                )
                hit["count"] = int(hit["count"]) + 1

        allowed_groups = normalized_constraints.get("allowed_binary_groups", None)
        if allowed_groups is not None and label not in allowed_groups:
            hit = unknown_hits.setdefault(
                label,
                {"count": 0, "first_index": row_idx},
            )
            hit["count"] = int(hit["count"]) + 1

    if enforce_binary_group_constraints and (len(forbidden_hits) > 0 or len(unknown_hits) > 0):
        parts = []
        if len(forbidden_hits) > 0:
            combo_msgs = []
            for combo, info in sorted(forbidden_hits.items(), key=lambda x: (len(x[0]), x[0])):
                combo_msgs.append(
                    f"{'+'.join(list(combo))} (n={int(info['count'])}, first_index={info['first_index']})"
                )
            parts.append("forbidden_combinations=[" + "; ".join(combo_msgs[:8]) + "]")
        if len(unknown_hits) > 0:
            group_msgs = []
            for group, info in sorted(unknown_hits.items(), key=lambda x: x[0]):
                group_msgs.append(
                    f"{group} (n={int(info['count'])}, first_index={info['first_index']})"
                )
            parts.append("unknown_binary_groups=[" + "; ".join(group_msgs[:8]) + "]")
        raise ValueError(
            "Binary group constraints violated while rebuilding full behavioral clusters. "
            + " | ".join(parts)
            + " | This indicates the supplied data contains binary-group combinations not present in training."
        )

    return pd.Series(labels, index=df.index, dtype="object")


def _add_clean_binary_annotation_columns(df: pd.DataFrame, binary_cols: list[str]) -> pd.DataFrame:
    """Back-compat shim for annotation columns; no-op for schema v3 exact-name grouping."""
    out = df.copy()
    for col in binary_cols:
        if col not in out.columns:
            continue
        clean_name = _binary_col_to_group_name(col)
        if clean_name == col or clean_name in out.columns:
            continue
        vals = pd.to_numeric(out[col], errors="coerce").fillna(0.0)
        out[clean_name] = (vals > 0).astype(np.int8)
    return out


def _rebuild_full_behavioral_cluster_from_intrinsic(
    adata,
    binary_cols_to_merge,
    intrinsic_col="intrinsic_behavioral_cluster",
    binary_group_constraints=None,
    enforce_binary_group_constraints=False,
):
    """Rebuild binary_group/full_behavioral_cluster columns from intrinsic cluster labels."""
    if intrinsic_col not in adata.obs.columns:
        raise ValueError(f"Missing '{intrinsic_col}' in adata.obs.")

    binary_cols = [str(c) for c in list(binary_cols_to_merge or [])]
    adata.obs = _add_clean_binary_annotation_columns(adata.obs, binary_cols)
    adata.obs["behavioral_clusterid"] = adata.obs[intrinsic_col].astype(str)
    adata.obs["binary_group"] = _assign_binary_group_labels(
        adata.obs,
        binary_cols,
        binary_group_constraints=binary_group_constraints,
        enforce_binary_group_constraints=bool(enforce_binary_group_constraints),
    ).astype("category")
    adata.obs["full_behavioral_cluster"] = (
        adata.obs["binary_group"].astype(str) + "_" + adata.obs["behavioral_clusterid"].astype(str)
    ).astype("category")
    adata.obs["behavioral_clusterid"] = adata.obs["behavioral_clusterid"].astype("category")
    return adata


def rename_intrinsic_behavioral_clusters(
    adata,
    mapping,
    binary_cols_to_merge=None,
    keep_unmapped=True,
):
    """
    Relabel intrinsic_behavioral_cluster and rebuild full_behavioral_cluster from binary groups.
    """
    if "intrinsic_behavioral_cluster" not in adata.obs.columns:
        raise ValueError("Missing 'intrinsic_behavioral_cluster' in adata.obs.")

    relabel_cluster_ids(
        adata=adata,
        mapping=mapping,
        cluster_key="intrinsic_behavioral_cluster",
        new_key="intrinsic_behavioral_cluster",
        keep_unmapped=keep_unmapped,
        overwrite_original=True,
    )

    if binary_cols_to_merge is None:
        binary_cols_to_merge = []
        if isinstance(adata.uns, dict):
            binary_cols_to_merge = list(
                adata.uns.get("clustering", {}).get("binary_cols_to_merge", [])
            )
    binary_group_constraints = None
    enforce_binary_group_constraints = False
    if isinstance(adata.uns, dict):
        clustering_meta = adata.uns.get("clustering", {})
        if isinstance(clustering_meta, dict):
            binary_group_constraints = clustering_meta.get("binary_group_constraints", None)
            enforce_binary_group_constraints = (
                isinstance(binary_group_constraints, dict)
                and ("forbidden_binary_combinations" in binary_group_constraints)
            )

    return _rebuild_full_behavioral_cluster_from_intrinsic(
        adata=adata,
        binary_cols_to_merge=binary_cols_to_merge,
        intrinsic_col="intrinsic_behavioral_cluster",
        binary_group_constraints=binary_group_constraints,
        enforce_binary_group_constraints=enforce_binary_group_constraints,
    )


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


def subsample_with_temporal_spacing(
    data,
    id_cols=None,
    time_col="position_t",
    min_spacing=5,
    min_time_from_track_start=None,
    max_samples=None,
    random_state=123,
):
    """Subsample rows with temporal spacing per track for DataFrame or AnnData input."""
    if id_cols is None:
        id_cols = ["sample_name", "TrackID"]

    rng = np.random.RandomState(random_state)

    def _select_track_row_indices(track_df):
        n_rows = len(track_df)
        if n_rows == 0:
            return []

        eligible_positions = np.arange(n_rows, dtype=int)
        if min_time_from_track_start is not None:
            time_values = pd.to_numeric(track_df[time_col], errors="coerce").to_numpy(dtype=float, copy=False)
            finite_mask = np.isfinite(time_values)
            if not np.any(finite_mask):
                return []
            track_start_t = float(np.nanmin(time_values[finite_mask]))
            min_eligible_t = track_start_t + float(min_time_from_track_start)
            eligible_positions = np.flatnonzero(finite_mask & (time_values >= min_eligible_t))

        n_eligible = int(len(eligible_positions))
        if n_eligible == 0:
            return []

        start_offset = int(rng.randint(0, min(int(min_spacing), n_eligible)))
        return eligible_positions[start_offset:: int(min_spacing)].tolist()

    if hasattr(data, "obs") and hasattr(data, "n_obs"):
        obs = data.obs
        missing_cols = [c for c in list(id_cols) + [time_col] if c not in obs.columns]
        if len(missing_cols) > 0:
            raise ValueError(
                "AnnData input is missing required columns for temporal subsampling: "
                f"{missing_cols}"
            )

        obs_work = obs[list(id_cols) + [time_col]].copy()
        obs_work["_obs_idx"] = np.arange(data.n_obs, dtype=int)
        selected_obs_idx = []

        for _, track_obs in obs_work.groupby(id_cols, observed=False):
            track_obs = track_obs.sort_values(time_col).reset_index(drop=True)
            local_selected = _select_track_row_indices(track_obs)
            if len(local_selected) == 0:
                continue
            selected_obs_idx.extend(track_obs.iloc[local_selected]["_obs_idx"].tolist())

        if len(selected_obs_idx) == 0:
            subset_empty = data[:0].copy()
            subset_empty.uns = copy.deepcopy(dict(getattr(data, "uns", {})))
            return subset_empty

        selected_obs_idx = np.asarray(selected_obs_idx, dtype=int)
        if max_samples is not None and len(selected_obs_idx) > int(max_samples):
            selected_obs_idx = rng.choice(selected_obs_idx, size=int(max_samples), replace=False)
        subset = data[selected_obs_idx].copy()
        subset.uns = copy.deepcopy(dict(getattr(data, "uns", {})))
        return subset

    subsampled_rows = []
    for _, track_df in data.groupby(id_cols):
        track_df = track_df.sort_values(time_col).reset_index(drop=True)
        selected_indices = _select_track_row_indices(track_df)
        if len(selected_indices) == 0:
            continue
        subsampled_rows.append(track_df.iloc[selected_indices])

    if len(subsampled_rows) == 0:
        return data.iloc[:0].copy()

    df_subsampled = pd.concat(subsampled_rows, ignore_index=True)
    if max_samples is not None and len(df_subsampled) > int(max_samples):
        df_subsampled = df_subsampled.sample(n=int(max_samples), random_state=random_state).reset_index(drop=True)
    return df_subsampled


def cap_values_to_quantile(df, features, lower_quantile=None, upper_quantile=0.99, return_limits=False):
    """Cap numeric feature values to requested quantiles."""
    if lower_quantile is None and upper_quantile is None:
        if return_limits:
            return df, {}
        return df

    df_capped = df.copy()
    cap_limits = {}

    for feature in features:
        if feature not in df_capped.columns:
            continue
        col = df_capped[feature]
        if pd.api.types.is_bool_dtype(col):
            if return_limits:
                cap_limits[feature] = {"lower": None, "upper": None}
            continue
        series = pd.to_numeric(col, errors="coerce").dropna()
        if len(series) == 0:
            if return_limits:
                cap_limits[feature] = {"lower": None, "upper": None}
            continue
        lower_limit = series.quantile(lower_quantile) if lower_quantile is not None else None
        upper_limit = series.quantile(upper_quantile) if upper_quantile is not None else None
        numeric_col = pd.to_numeric(df_capped[feature], errors="coerce")
        df_capped[feature] = numeric_col.clip(lower=lower_limit, upper=upper_limit)
        if return_limits:
            cap_limits[feature] = {
                "lower": None if lower_limit is None else float(lower_limit),
                "upper": None if upper_limit is None else float(upper_limit),
            }

    if return_limits:
        return df_capped, cap_limits
    return df_capped


def prepare_state_classification_dataset(
    df_positions,
    features,
    binary_features_to_group,
    window_size=5,
    min_spacing=None,
    descriptive_features=("mean", "median", "std", "net_displacement", "straightness", "mean_square_displacement"),
    lower_quantile_cap=None,
    upper_quantile_cap=0.99,
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

    quantile_cap_limits = {}
    loaded_from_cache = False
    adata_prepared = None
    cached_pre_meta = {}
    if reuse_prepared_dataset and prepared_dataset_path is not None and prepared_dataset_path.exists():
        try:
            adata_prepared = sc.read_h5ad(prepared_dataset_path)
            loaded_from_cache = True
            _vinfo(verbose, "state-clustering", f"Loaded prepared window cache from: {prepared_dataset_path}")
        except Exception as exc:
            _vinfo(verbose, "state-clustering", f"Prepared cache load failed; recomputing windows ({exc})")

    if not loaded_from_cache:
        df_windows_descriptive = create_descriptive_track_dataset(
            df_tracks=df_positions,
            columns_to_summarize=features,
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


def _save_pdf_page_a4(pdf, fig, orientation="portrait"):
    fig.set_size_inches(*(A4_PORTRAIT if orientation == "portrait" else A4_LANDSCAPE), forward=True)
    pdf.savefig(fig, bbox_inches="tight")


def _resolve_output_dir(output_dir):
    """Resolve/create the global BEHAV3D output root directory."""
    if output_dir is None:
        raise ValueError("output_dir is required.")
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    return output_dir_path


@dataclass(frozen=True)
class StatePaths:
    output_dir: Path
    analysis_outdir: Path
    feature_outdir: Path
    filtered_tracks_csv: Path
    combined_tracks_csv: Path
    state_outdir: Path
    processing_outdir: Path
    state_clustering_outdir: Path
    model_adata_path: Path
    prepared_full_adata_path: Path
    full_output_adata_path: Path
    intrinsic_outdir: Path
    full_outdir: Path
    intrinsic_qc_outdir: Path
    full_qc_outdir: Path
    intrinsic_classifier_default_path: Path
    full_classifier_default_path: Path
    state_composition_outdir: Path
    state_transitions_outdir: Path


def _resolve_state_outdir(output_dir, cell_type):
    """Return/create state-classification root: output_dir/analysis/<cell_type>/behavioral_states."""
    if cell_type is None or len(str(cell_type).strip()) == 0:
        raise ValueError("cell_type is required.")
    root = _resolve_output_dir(output_dir)
    state_outdir = root / "analysis" / str(cell_type) / "behavioral_states"
    state_outdir.mkdir(parents=True, exist_ok=True)
    return state_outdir


def _resolve_state_paths(output_dir, cell_type):
    """Resolve canonical state-pipeline paths as a typed context."""
    if cell_type is None or len(str(cell_type).strip()) == 0:
        raise ValueError("cell_type is required.")
    root = _resolve_output_dir(output_dir)
    analysis_outdir = root / "analysis" / str(cell_type)
    feature_outdir = analysis_outdir / "track_features"

    state_outdir = root / "analysis" / str(cell_type) / "behavioral_states"
    state_outdir.mkdir(parents=True, exist_ok=True)
    processing_outdir = state_outdir / "processing"
    processing_outdir.mkdir(parents=True, exist_ok=True)

    intrinsic_outdir = processing_outdir / "intrinsic_behavioral_classification"
    full_outdir = processing_outdir / "full_behavioral_classification"
    intrinsic_qc_outdir = intrinsic_outdir / "quality_control"
    full_qc_outdir = full_outdir / "quality_control"

    return StatePaths(
        output_dir=root,
        analysis_outdir=analysis_outdir,
        feature_outdir=feature_outdir,
        filtered_tracks_csv=feature_outdir / f"BEHAV3D_{cell_type}_combined_track_features_filtered.csv",
        combined_tracks_csv=feature_outdir / f"BEHAV3D_{cell_type}_combined_track_features.csv",
        state_outdir=state_outdir,
        processing_outdir=processing_outdir,
        state_clustering_outdir=processing_outdir / "state_clustering",
        model_adata_path=processing_outdir / f"BEHAV3D_{cell_type}_behavioral_states_modeldata.h5ad",
        prepared_full_adata_path=processing_outdir / f"BEHAV3D_{cell_type}_behavioral_states_prepared_full.h5ad",
        full_output_adata_path=state_outdir / f"BEHAV3D_{cell_type}_behavioral_states.h5ad",
        intrinsic_outdir=intrinsic_outdir,
        full_outdir=full_outdir,
        intrinsic_qc_outdir=intrinsic_qc_outdir,
        full_qc_outdir=full_qc_outdir,
        intrinsic_classifier_default_path=intrinsic_outdir / f"intrinsic_state_random_forest_{cell_type}.pkl",
        full_classifier_default_path=full_outdir / f"state_random_forest_{cell_type}.pkl",
        state_composition_outdir=state_outdir / "state_composition",
        state_transitions_outdir=state_outdir / "state_transitions",
    )


def _resolve_analysis_paths(output_dir, cell_type):
    """Resolve standard analysis folders and expected track-feature CSV inputs."""
    if cell_type is None or len(str(cell_type).strip()) == 0:
        raise ValueError("cell_type is required.")
    root = _resolve_output_dir(output_dir)
    analysis_outdir = root / "analysis" / str(cell_type)
    feature_outdir = analysis_outdir / "track_features"
    return {
        "output_dir": root,
        "analysis_outdir": analysis_outdir,
        "feature_outdir": feature_outdir,
        "filtered_tracks_csv": feature_outdir / f"BEHAV3D_{cell_type}_combined_track_features_filtered.csv",
        "combined_tracks_csv": feature_outdir / f"BEHAV3D_{cell_type}_combined_track_features.csv",
    }


def _resolve_legacy_prepared_full_adata_path(state_paths, cell_type):
    return state_paths.state_outdir / f"BEHAV3D_{cell_type}_behavioral_states_prepared_full.h5ad"


def _resolve_positions_csv_path(output_dir=None, cell_type=None, state_paths=None):
    """Resolve the positions/features CSV from standard BEHAV3D analysis layout."""
    if state_paths is None:
        analysis_paths = _resolve_analysis_paths(output_dir, cell_type)
        filtered_path = analysis_paths["filtered_tracks_csv"]
        combined_path = analysis_paths["combined_tracks_csv"]
    else:
        filtered_path = state_paths.filtered_tracks_csv
        combined_path = state_paths.combined_tracks_csv
    if filtered_path.exists():
        return filtered_path
    if combined_path.exists():
        return combined_path
    raise FileNotFoundError(
        "Could not find track-features CSV for state classification. "
        f"Expected one of: '{filtered_path}' or '{combined_path}'. "
        f"Run feature extraction for cell_type='{cell_type}' first, or place the CSV in "
        "output_dir/analysis/<cell_type>/track_features."
    )


def _resolve_processing_outdir(base_outdir):
    """Return/create processing output root under a stage root folder."""
    if base_outdir is None:
        raise ValueError("base_outdir is required to resolve processing output.")
    processing_outdir = Path(base_outdir) / "processing"
    processing_outdir.mkdir(parents=True, exist_ok=True)
    return processing_outdir


def _resolve_processing_subdir(base_outdir, *parts):
    """Return/create a subfolder under processing output root."""
    processing_outdir = _resolve_processing_outdir(base_outdir)
    subdir = processing_outdir.joinpath(*parts)
    subdir.mkdir(parents=True, exist_ok=True)
    return subdir


def _resolve_state_clustering_outdir(base_outdir):
    """Return/create clustering diagnostics/proportion folder under processing root."""
    if base_outdir is None:
        raise ValueError("base_outdir is required to resolve state clustering output.")
    state_clustering_outdir = Path(base_outdir) / "state_clustering"
    state_clustering_outdir.mkdir(parents=True, exist_ok=True)
    return state_clustering_outdir


def _resolve_state_stage_paths(output_dir, cell_type):
    """Compatibility wrapper. Prefer `_resolve_state_paths(...)`."""
    return asdict(_resolve_state_paths(output_dir, cell_type))


def _sanitize_filename_token(value, fallback="value"):
    """Sanitize a value for safe use in filenames."""
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip())
    token = token.strip("._-")
    return token if len(token) > 0 else str(fallback)


def _rmtree_ignore_missing(path):
    """Remove a directory tree while tolerating concurrent missing-file races."""
    path = Path(path)
    if not path.exists():
        return

    def _ignore_missing(func, target, excinfo):
        err = excinfo[1]
        if isinstance(err, FileNotFoundError):
            return
        raise err

    shutil.rmtree(path, onexc=_ignore_missing)


def _select_exemplar_windows_by_cluster(
    adata_windows,
    *,
    cluster_key="intrinsic_behavioral_cluster",
    n_per_cluster=5,
    sample_key="sample_name",
    track_key="TrackID",
    time_key="position_t",
    tmin_key="window_start_position_t",
    tmax_key="window_end_position_t",
    seed=0,
):
    required_cols = [sample_key, track_key, cluster_key, time_key, tmin_key, tmax_key]
    missing = [c for c in required_cols if c not in adata_windows.obs.columns]
    if len(missing) > 0:
        raise ValueError(f"adata_windows.obs missing required columns for exemplar selection: {missing}")

    rng = np.random.default_rng(seed)
    window_df = adata_windows.obs[required_cols].copy()
    window_df = window_df.dropna(subset=[cluster_key, tmin_key, tmax_key]).copy()
    if len(window_df) == 0:
        raise ValueError("No windows available after dropping rows missing cluster or window bounds.")

    cluster_total_counts = window_df.groupby(cluster_key, observed=True).size().to_dict()
    chosen_parts = []

    for _, df_cl in window_df.groupby(cluster_key, sort=False, observed=True):
        if len(df_cl) == 0:
            continue

        df_cl = df_cl.copy()
        df_cl["_draw"] = rng.random(len(df_cl))
        df_cl = df_cl.sort_values(
            ["_draw", str(time_key), str(tmin_key), str(tmax_key)],
            kind="stable",
        )
        df_cl = df_cl.drop_duplicates(
            subset=[sample_key, track_key],
            keep="first",
        ).copy()

        k = min(int(n_per_cluster), len(df_cl))
        if k <= 0:
            continue
        idx = rng.choice(len(df_cl), size=k, replace=False)
        chosen_parts.append(df_cl.iloc[idx].copy())

    if len(chosen_parts) == 0:
        raise ValueError("No exemplar windows were selected.")

    chosen = pd.concat(chosen_parts, axis=0, ignore_index=True)
    chosen = chosen.sort_values(
        [str(cluster_key), str(sample_key), str(track_key), str(tmin_key), str(tmax_key), str(time_key)],
        kind="stable",
    ).reset_index(drop=True)
    chosen["example_rank"] = (
        chosen.groupby(cluster_key, observed=True).cumcount() + 1
    ).astype(int)
    return chosen, cluster_total_counts


def _render_state_cluster_exemplar_videos(
    model_adata,
    *,
    state_paths,
    cell_type,
    state_clustering_outdir,
    cluster_key="intrinsic_behavioral_cluster",
    n_per_cluster=5,
    random_state=123,
    video_fps=6,
    video_dpi=180,
    video_margin=20,
    video_pmin=0.0,
    video_pmax=99.99,
    video_track_color="#63ff33",
    video_show_segments=True,
    video_segment_style="outline",
    video_segment_color="#ffffff",
    verbose=True,
):
    result = {
        "enabled": True,
        "cluster_key": str(cluster_key),
        "selection_csv": None,
        "video_paths_by_cluster": {},
        "global_xy_shape": None,
        "cluster_total_counts": {},
        "segment_outline_errors": {},
        "error": None,
        "render_config": {
            "stage": "after_clustering",
            "cluster_key": str(cluster_key),
            "state_key": str(cluster_key),
            "n_per_cluster": int(n_per_cluster),
            "selection_random_state": int(random_state),
            "layout_mode": "per_cluster",
            "video_fps": int(video_fps),
            "video_dpi": int(video_dpi),
            "video_margin": int(video_margin),
            "video_pmin": float(video_pmin),
            "video_pmax": float(video_pmax),
            "video_track_color": str(video_track_color),
            "video_show_segments": bool(video_show_segments),
            "video_segment_style": str(video_segment_style),
            "video_segment_color": str(video_segment_color),
            "show_state_legend": False,
            "positions_csv_path": None,
        },
    }

    if str(cluster_key) not in model_adata.obs.columns:
        result["error"] = f"Missing clustering exemplar column '{cluster_key}' in model_adata.obs."
        return result

    exemplar_started = _vstart(
        verbose,
        "state-clustering",
        f"render exemplar videos | cluster_key={cluster_key}",
    )
    exemplar_root = (
        Path(state_clustering_outdir)
        / "example_windows"
        / _sanitize_filename_token(cluster_key, fallback="cluster")
    )
    _rmtree_ignore_missing(exemplar_root)
    exemplar_root.mkdir(parents=True, exist_ok=True)

    try:
        positions_csv_path = _resolve_positions_csv_path(state_paths=state_paths)
        result["render_config"]["positions_csv_path"] = str(positions_csv_path)
        df_positions = pd.read_csv(positions_csv_path, low_memory=False)
        df_positions = df_positions.sort_values(["sample_name", "TrackID", "position_t"])
        required_pos_cols = [
            "sample_name",
            "TrackID",
            "position_t",
            "pixel_position_x",
            "pixel_position_y",
            "pixel_position_z",
        ]
        missing_pos_cols = [c for c in required_pos_cols if c not in df_positions.columns]
        if len(missing_pos_cols) > 0:
            raise ValueError(
                "State exemplar videos require pixel-position columns in the positions CSV. "
                f"Missing columns: {missing_pos_cols}. csv='{positions_csv_path}'"
            )

        chosen_exemplars, cluster_total_counts = _select_exemplar_windows_by_cluster(
            adata_windows=model_adata,
            n_per_cluster=int(n_per_cluster),
            cluster_key=str(cluster_key),
            seed=int(random_state),
        )
        result["cluster_total_counts"] = {
            str(k): int(v) for k, v in dict(cluster_total_counts).items()
        }
        chosen_exemplars["cluster_total_windows"] = chosen_exemplars[str(cluster_key)].astype(str).map(
            result["cluster_total_counts"]
        )

        selection_csv = exemplar_root / (
            f"example_window_selection_cluster_{_sanitize_filename_token(cluster_key, fallback='cluster')}.csv"
        )
        chosen_exemplars.to_csv(selection_csv, index=False)
        result["selection_csv"] = str(selection_csv)
        _vsave(verbose, "state-clustering", "example window selection CSV", selection_csv)

        backprojection_video_out = save_selected_fulltrack_cluster_videos(
            chosen_df=chosen_exemplars,
            df_positions=df_positions,
            output_folder=state_paths.output_dir,
            out_dir=exemplar_root,
            cluster_key=str(cluster_key),
            sample_key="sample_name",
            track_key="TrackID",
            tmin_key="window_start_position_t",
            tmax_key="window_end_position_t",
            fps=int(video_fps),
            dpi=int(video_dpi),
            margin=int(video_margin),
            pmin=float(video_pmin),
            pmax=float(video_pmax),
            track_color=str(video_track_color),
            show_segment_outlines=bool(video_show_segments),
            segment_style=str(video_segment_style),
            segment_color=str(video_segment_color),
            cell_type=cell_type,
            verbose=bool(verbose),
        )
        result["video_paths_by_cluster"] = dict(
            backprojection_video_out.get("video_paths_by_cluster", {})
        )
        result["segment_outline_errors"] = dict(
            backprojection_video_out.get("segment_outline_errors", {})
        )
    except Exception as exc:
        result["error"] = str(exc)
        _vinfo(verbose, "state-clustering", f"State exemplar videos skipped: {exc}")

    _vdone(verbose, "state-clustering", "render exemplar videos", exemplar_started)
    return result


def plot_clustering_diagnostics_pdf(
    adata,
    cluster_col,
    feature_cols,
    pdf_path,
    title="Clustering diagnostics",
):

    if adata is None or adata.n_obs == 0:
        return
    if cluster_col not in adata.obs.columns:
        return

    ad = _ensure_umap(adata, n_neighbors=30, random_state=123)

    pdf_path = Path(pdf_path)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(pdf_path) as pdf:
        fig, (ax_umap, ax_corr) = plt.subplots(
            1, 2, figsize=A4_LANDSCAPE, gridspec_kw={"width_ratios": [1.25, 0.75]}
        )
        sc.pl.umap(
            ad,
            color=cluster_col,
            ax=ax_umap,
            show=False,
            title="UMAP",
            legend_loc="on data",
            legend_fontsize=7,
        )
        ax_umap.set_aspect("equal", adjustable="box")

        # Correlation matrix on same page using explicit axis.  
        sc.pl.correlation_matrix(
            ad,
            groupby=cluster_col,
            show_correlation_numbers=True,
            ax=ax_corr,
            show=False,
        )
        ax_corr.set_title("Cluster correlations", fontsize=10)

        fig.suptitle(f"{title} | UMAP + Correlation", fontsize=12, fontweight="bold", y=0.995)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        umap_bbox = ax_umap.get_position().frozen()

        # Keep UMAP unchanged; only move/shrink correlation panel.
        bbox = ax_corr.get_position()
        extra_gap = 0.03
        new_x = min(0.98 - bbox.width, bbox.x0 + extra_gap)
        height_shrink = 0.5
        new_h = bbox.height * height_shrink
        ax_corr.set_position([new_x, bbox.y0 + (bbox.height - new_h) / 2, bbox.width, new_h])
        ax_umap.set_position(umap_bbox)

        _save_pdf_page_a4(pdf, fig, orientation="landscape")
        plt.close(fig)

        use_dendrogram = ad.obs[cluster_col].astype(str).nunique() > 1
        dendrogram_key = f"dendrogram_{cluster_col}"
        dendrogram_arg = False
        if use_dendrogram:
            try:
                if dendrogram_key not in ad.uns:
                    sc.tl.dendrogram(ad, groupby=cluster_col, key_added=dendrogram_key)
                if dendrogram_key in ad.uns:
                    dendrogram_arg = dendrogram_key
            except Exception:
                dendrogram_arg = False

        sc.pl.heatmap(
            ad,
            var_names=list(feature_cols),
            groupby=cluster_col,
            standard_scale="var",
            figsize=A4_LANDSCAPE,
            swap_axes=True,
            dendrogram=dendrogram_arg,
            show_gene_labels=True,
            show=False,
        )
  
        fig_hm = plt.gcf()
        fig_hm.suptitle(f"{title} | Heatmap", y=0.995)
        _save_pdf_page_a4(pdf, fig_hm, orientation="landscape")
        plt.close(fig_hm)


def _export_clustering_diagnostics_csvs(
    adata,
    cluster_col,
    feature_cols,
    outdir,
    filename_prefix="behavioral_states_diagnostics",
):
    """Export CSV companions for clustering diagnostics panels."""
    if adata is None or adata.n_obs == 0:
        return {}
    if cluster_col not in adata.obs.columns:
        return {}

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ad = _ensure_umap(adata, n_neighbors=30, random_state=123)

    diagnostics_paths = {}

    # UMAP panel source data.
    umap = np.asarray(ad.obsm.get("X_umap", np.empty((ad.n_obs, 2))), dtype=float)
    if umap.ndim != 2 or umap.shape[0] != ad.n_obs:
        umap = np.zeros((ad.n_obs, 2), dtype=float)
    if umap.shape[1] < 2:
        umap = np.pad(umap, ((0, 0), (0, 2 - umap.shape[1])), mode="constant")
    umap_df = pd.DataFrame(
        {
            "obs_index": ad.obs.index.astype(str),
            "umap_1": umap[:, 0],
            "umap_2": umap[:, 1],
            "cluster_col": ad.obs[cluster_col].astype(str).to_numpy(),
        }
    )
    for c in ["sample_name", "TrackID", "position_t"]:
        if c in ad.obs.columns:
            umap_df[c] = ad.obs[c].to_numpy()
    umap_csv = outdir / f"{filename_prefix}_umap.csv"
    umap_df.to_csv(umap_csv, index=False)
    diagnostics_paths["umap_csv"] = str(umap_csv)

    valid_features = [str(c) for c in list(feature_cols or []) if str(c) in ad.var_names]
    cluster_order = sorted(ad.obs[cluster_col].astype(str).unique().tolist(), key=_mixed_label_sort_key)
    if len(valid_features) == 0:
        cluster_means = pd.DataFrame(index=cluster_order)
    else:
        X = _to_numpy_2d(ad[:, valid_features].X).astype(float, copy=False)
        feature_df = pd.DataFrame(X, columns=valid_features, index=ad.obs.index)
        feature_df["_cluster"] = ad.obs[cluster_col].astype(str).to_numpy()
        cluster_means = feature_df.groupby("_cluster", observed=True)[valid_features].mean()
        cluster_means = cluster_means.reindex(cluster_order)

    # Correlation panel source data (cluster-by-cluster correlation).
    if cluster_means.shape[0] > 0 and cluster_means.shape[1] > 0:
        corr_df = cluster_means.T.corr()
    else:
        corr_df = pd.DataFrame(index=cluster_means.index, columns=cluster_means.index, dtype=float)
    corr_csv = outdir / f"{filename_prefix}_cluster_correlation.csv"
    corr_df.to_csv(corr_csv)
    diagnostics_paths["cluster_correlation_csv"] = str(corr_csv)

    # Heatmap panel source data (cluster x feature, standard_scale='var' equivalent).
    heatmap_df = cluster_means.copy()
    if heatmap_df.shape[1] > 0:
        for col in heatmap_df.columns:
            series = pd.to_numeric(heatmap_df[col], errors="coerce")
            col_min = series.min()
            col_max = series.max()
            if pd.isna(col_min) or pd.isna(col_max) or float(col_max) == float(col_min):
                heatmap_df[col] = 0.0
            else:
                heatmap_df[col] = (series - col_min) / (col_max - col_min)
    heatmap_csv = outdir / f"{filename_prefix}_heatmap_matrix.csv"
    heatmap_df.to_csv(heatmap_csv)
    diagnostics_paths["heatmap_matrix_csv"] = str(heatmap_csv)

    return diagnostics_paths


def _ensure_umap(adata, n_neighbors=30, random_state=123):
    """Ensure AnnData has a UMAP embedding for plotting."""
    if "X_umap" in adata.obsm:
        return adata

    requested_n_neighbors = int(n_neighbors)
    neighbors_meta = adata.uns.get("neighbors", {}) if isinstance(adata.uns, dict) else {}
    params = neighbors_meta.get("params", {}) if isinstance(neighbors_meta, dict) else {}
    existing_n_neighbors = params.get("n_neighbors", None)
    conn_key = neighbors_meta.get("connectivities_key", "connectivities")
    dist_key = neighbors_meta.get("distances_key", "distances")
    has_graph = (conn_key in adata.obsp) and (dist_key in adata.obsp)

    can_reuse_neighbors = False
    try:
        can_reuse_neighbors = has_graph and (int(existing_n_neighbors) == requested_n_neighbors)
    except (TypeError, ValueError):
        can_reuse_neighbors = False

    if not can_reuse_neighbors:
        if "X_pca" in adata.obsm:
            sc.pp.neighbors(
                adata,
                n_neighbors=requested_n_neighbors,
                use_rep="X_pca",
                random_state=int(random_state),
            )
        else:
            sc.pp.neighbors(adata, n_neighbors=requested_n_neighbors, random_state=int(random_state))
    sc.tl.umap(adata, random_state=int(random_state))
    return adata


def check_model_clusterid_consistency(
    model_adata,
    adata_full,
    model_cluster_col="intrinsic_behavioral_cluster",
    full_cluster_col="intrinsic_behavioral_cluster",
    key_cols=("sample_name", "TrackID", "position_t"),
    max_examples=10,
    verbose=True,
):
    """
    Check whether model rows keep the same cluster label after label transfer into full data.

    Matching priority:
    1) key-based join on key_cols
    2) obs_names intersection (fallback only if keys are unavailable)
    """
    if model_cluster_col not in model_adata.obs.columns:
        raise ValueError(f"Missing '{model_cluster_col}' in model_adata.obs.")
    if full_cluster_col not in adata_full.obs.columns:
        raise ValueError(f"Missing '{full_cluster_col}' in adata_full.obs.")

    model_obs = model_adata.obs.copy()
    full_obs = adata_full.obs.copy()
    model_obs["_cluster_model"] = model_obs[model_cluster_col].astype(str)
    full_obs["_cluster_full"] = full_obs[full_cluster_col].astype(str)

    missing_model_keys = [c for c in key_cols if c not in model_obs.columns]
    missing_full_keys = [c for c in key_cols if c not in full_obs.columns]
    keys_available = len(missing_model_keys) == 0 and len(missing_full_keys) == 0
    used_match_mode = "key_cols" if keys_available else "obs_names"

    if keys_available:
        keys = list(key_cols)
        model_cmp = model_obs[keys + ["_cluster_model"]].copy()
        full_cmp = full_obs[keys + ["_cluster_full"]].copy()
        for c in keys:
            model_cmp[c] = model_cmp[c].astype(str)
            full_cmp[c] = full_cmp[c].astype(str)

        model_dups = int(model_cmp.duplicated(subset=keys, keep=False).sum())
        full_dups = int(full_cmp.duplicated(subset=keys, keep=False).sum())
        if model_dups > 0 or full_dups > 0:
            raise ValueError(
                f"Cannot do one-to-one comparison: duplicate key rows found "
                f"(model duplicates={model_dups}, full duplicates={full_dups}) for key_cols={keys}."
            )

        merged = model_cmp.merge(full_cmp, on=keys, how="left", validate="one_to_one")
    else:
        overlap = model_obs.index.intersection(full_obs.index)
        if len(overlap) == 0:
            raise ValueError(
                "Cannot compare labels: no shared obs_names and key columns are missing for key-based matching. "
                f"missing in model_adata: {missing_model_keys}; missing in adata_full: {missing_full_keys}"
            )
        merged = pd.DataFrame(
            {
                "_cluster_model": model_obs.loc[overlap, "_cluster_model"].to_numpy(),
                "_cluster_full": full_obs.loc[overlap, "_cluster_full"].to_numpy(),
            },
            index=overlap,
        ).reset_index(drop=False).rename(columns={"index": "_row_id"})

    merged["_match"] = merged["_cluster_model"] == merged["_cluster_full"]
    missing_in_full = int(merged["_cluster_full"].isna().sum())
    n_total = int(len(merged))
    n_match = int(merged["_match"].sum())
    n_mismatch = int((~merged["_match"] & merged["_cluster_full"].notna()).sum())

    summary = {
        "match_mode": used_match_mode,
        "n_model_rows_checked": n_total,
        "n_matching_clusterid": n_match,
        "n_mismatched_clusterid": n_mismatch,
        "n_missing_rows_in_adata_full": missing_in_full,
        "match_rate": float(n_match / n_total) if n_total > 0 else np.nan,
    }

    _vinfo(
        verbose,
        "state-apply",
        "clusterid consistency check "
        f"(mode={summary['match_mode']}): "
        f"checked={summary['n_model_rows_checked']}, "
        f"matches={summary['n_matching_clusterid']}, "
        f"mismatches={summary['n_mismatched_clusterid']}, "
        f"missing_in_full={summary['n_missing_rows_in_adata_full']}",
    )
    if n_mismatch > 0:
        _vinfo(verbose, "state-apply", f"mismatch preview suppressed (n_mismatch={n_mismatch})")

    return summary, merged


def _to_numpy_2d(X):
    """Convert dense/sparse matrix-like data to a 2D numpy array."""
    if hasattr(X, "toarray"):
        return X.toarray()
    return np.asarray(X)


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
        "scaler": None,
    }
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

    scaler_meta = pre_meta.get("scaler", None)
    has_scaler = isinstance(scaler_meta, dict) and ("mean" in scaler_meta) and ("scale" in scaler_meta)
    if bool(scale_features) != bool(has_scaler):
        reasons.append("scaler presence mismatch")

    return len(reasons) == 0, reasons


def _matches_model_sampling_cache(
    adata,
    *,
    kept_features,
    min_spacing_used,
    start_cutoff_mode=None,
    start_cutoff_window_size=None,
    max_samples,
    random_state,
    id_cols=("sample_name", "TrackID"),
    time_col="position_t",
    nan_policy="drop_any_nan_in_kept_features",
):
    """Check whether cached model adata matches sampling/cleaning stage settings."""
    reasons = []
    if not (hasattr(adata, "uns") and isinstance(adata.uns, dict)):
        return False, ["model cache uns missing/invalid"]
    pre_meta = adata.uns.get("preprocessing", {})
    if not isinstance(pre_meta, dict):
        return False, ["model cache preprocessing metadata missing/invalid"]
    model_cache = pre_meta.get("model_cache", {})
    if not isinstance(model_cache, dict):
        return False, ["model cache metadata missing"]

    features_meta = model_cache.get("features", {})
    if isinstance(features_meta, dict):
        kept_cached = features_meta.get("kept_features", [])
    else:
        kept_cached = features_meta
    if kept_cached is None:
        kept_cached = []
    kept_cached = [str(x) for x in list(kept_cached)]

    if kept_features is None:
        kept_features = []
    kept_expected = [str(x) for x in list(kept_features)]
    if kept_cached != kept_expected:
        reasons.append("model_cache.features.kept_features mismatch")
    var_names = [str(x) for x in list(getattr(adata, "var_names", []))]
    if var_names != kept_expected:
        reasons.append("model cache var_names mismatch")

    sampling_meta = model_cache.get("sampling", {})
    sampling_meta = sampling_meta if isinstance(sampling_meta, dict) else {}
    try:
        got_spacing = int(sampling_meta.get("min_spacing_used", -1))
    except Exception:
        got_spacing = -1
    if got_spacing != int(min_spacing_used):
        reasons.append("model_cache.sampling.min_spacing_used mismatch")
    got_max_samples = sampling_meta.get("max_samples", None)
    got_max_samples = None if got_max_samples is None else int(got_max_samples)
    exp_max_samples = None if max_samples is None else int(max_samples)
    if got_max_samples != exp_max_samples:
        reasons.append("model_cache.sampling.max_samples mismatch")
    try:
        got_random_state = int(sampling_meta.get("random_state", -1))
    except Exception:
        got_random_state = -1
    if got_random_state != int(random_state):
        reasons.append("model_cache.sampling.random_state mismatch")
    if [str(x) for x in list(sampling_meta.get("id_cols", []))] != [str(x) for x in list(id_cols)]:
        reasons.append("model_cache.sampling.id_cols mismatch")
    if str(sampling_meta.get("time_col", "")) != str(time_col):
        reasons.append("model_cache.sampling.time_col mismatch")
    got_start_cutoff_mode = sampling_meta.get("start_cutoff_mode", None)
    got_start_cutoff_mode = None if got_start_cutoff_mode is None else str(got_start_cutoff_mode)
    exp_start_cutoff_mode = None if start_cutoff_mode is None else str(start_cutoff_mode)
    if got_start_cutoff_mode != exp_start_cutoff_mode:
        reasons.append("model_cache.sampling.start_cutoff_mode mismatch")
    got_start_cutoff_window = sampling_meta.get("start_cutoff_window_size", None)
    got_start_cutoff_window = None if got_start_cutoff_window is None else int(got_start_cutoff_window)
    exp_start_cutoff_window = (
        None if start_cutoff_window_size is None else int(start_cutoff_window_size)
    )
    if got_start_cutoff_window != exp_start_cutoff_window:
        reasons.append("model_cache.sampling.start_cutoff_window_size mismatch")

    cleaning_meta = model_cache.get("cleaning", {})
    cleaning_meta = cleaning_meta if isinstance(cleaning_meta, dict) else {}
    if str(cleaning_meta.get("nan_policy", "")) != str(nan_policy):
        reasons.append("model_cache.cleaning.nan_policy mismatch")

    return len(reasons) == 0, reasons


def _matches_cached_pca_stage(
    adata,
    *,
    pca_var_selection,
    svd_solver,
    random_state,
    ncomps_requested,
):
    """Check whether cached PCA output matches current PCA settings."""
    reasons = []
    if "X_pca" not in adata.obsm:
        return False, ["X_pca missing"]
    if not (hasattr(adata, "uns") and isinstance(adata.uns, dict)):
        return False, ["uns missing/invalid"]
    model_cache = (
        adata.uns.get("preprocessing", {}).get("model_cache", {})
        if isinstance(adata.uns.get("preprocessing", {}), dict)
        else {}
    )
    pca_meta = model_cache.get("pca", {}) if isinstance(model_cache, dict) else {}
    if not isinstance(pca_meta, dict):
        return False, ["model_cache.pca missing"]

    got_pca_var = pca_meta.get("pca_var_selection", np.nan)
    if not np.isclose(float(got_pca_var), float(pca_var_selection), rtol=1e-8, atol=1e-12):
        reasons.append("model_cache.pca.pca_var_selection mismatch")
    if str(pca_meta.get("svd_solver", "")) != str(svd_solver):
        reasons.append("model_cache.pca.svd_solver mismatch")
    try:
        got_random_state = int(pca_meta.get("random_state", -1))
    except Exception:
        got_random_state = -1
    if got_random_state != int(random_state):
        reasons.append("model_cache.pca.random_state mismatch")
    try:
        got_ncomps_requested = int(pca_meta.get("ncomps_requested", -1))
    except Exception:
        got_ncomps_requested = -1
    if got_ncomps_requested != int(ncomps_requested):
        reasons.append("model_cache.pca.ncomps_requested mismatch")
    x_pca = adata.obsm.get("X_pca", None)
    ncomps_effective = -1 if x_pca is None else int(x_pca.shape[1])
    if int(pca_meta.get("ncomps_effective", -1)) != ncomps_effective:
        reasons.append("model_cache.pca.ncomps_effective mismatch")

    return len(reasons) == 0, reasons


def _matches_cached_neighbors_stage(
    adata,
    *,
    n_neighbors,
    random_state,
    use_rep="X_pca",
):
    """Check whether cached neighbors graph matches current neighbor settings."""
    reasons = []
    if not (hasattr(adata, "uns") and isinstance(adata.uns, dict)):
        return False, ["uns missing/invalid"]
    neighbors_uns = adata.uns.get("neighbors", {})
    neighbors_uns = neighbors_uns if isinstance(neighbors_uns, dict) else {}
    params = neighbors_uns.get("params", {})
    params = params if isinstance(params, dict) else {}
    conn_key = neighbors_uns.get("connectivities_key", "connectivities")
    dist_key = neighbors_uns.get("distances_key", "distances")
    has_graph = (conn_key in adata.obsp) and (dist_key in adata.obsp)
    if not has_graph:
        reasons.append("neighbors graph missing")
    try:
        if int(params.get("n_neighbors", -1)) != int(n_neighbors):
            reasons.append("neighbors params.n_neighbors mismatch")
    except Exception:
        reasons.append("neighbors params.n_neighbors invalid")

    model_cache = (
        adata.uns.get("preprocessing", {}).get("model_cache", {})
        if isinstance(adata.uns.get("preprocessing", {}), dict)
        else {}
    )
    neighbors_meta = model_cache.get("neighbors", {}) if isinstance(model_cache, dict) else {}
    if not isinstance(neighbors_meta, dict):
        return False, reasons + ["model_cache.neighbors missing"]
    if int(neighbors_meta.get("n_neighbors", -1)) != int(n_neighbors):
        reasons.append("model_cache.neighbors.n_neighbors mismatch")
    if int(neighbors_meta.get("random_state", -1)) != int(random_state):
        reasons.append("model_cache.neighbors.random_state mismatch")
    if str(neighbors_meta.get("use_rep", "")) != str(use_rep):
        reasons.append("model_cache.neighbors.use_rep mismatch")

    return len(reasons) == 0, reasons


def _matches_cached_umap_stage(
    adata,
    *,
    min_dist,
    random_state,
):
    """Check whether cached UMAP embedding matches current UMAP settings."""
    reasons = []
    if "X_umap" not in adata.obsm:
        return False, ["X_umap missing"]
    if not (hasattr(adata, "uns") and isinstance(adata.uns, dict)):
        return False, ["uns missing/invalid"]
    model_cache = (
        adata.uns.get("preprocessing", {}).get("model_cache", {})
        if isinstance(adata.uns.get("preprocessing", {}), dict)
        else {}
    )
    umap_meta = model_cache.get("umap", {}) if isinstance(model_cache, dict) else {}
    if not isinstance(umap_meta, dict):
        return False, ["model_cache.umap missing"]
    if float(umap_meta.get("min_dist", np.nan)) != float(min_dist):
        reasons.append("model_cache.umap.min_dist mismatch")
    if int(umap_meta.get("random_state", -1)) != int(random_state):
        reasons.append("model_cache.umap.random_state mismatch")
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
    """Load a saved state classification v3 classifier artifact pickle."""
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

    try:
        schema_version = int(pipeline_meta.get("schema_version", -1))
    except Exception as exc:
        raise ValueError(
            f"{artifact_name} has invalid pipeline_metadata.schema_version. "
            "Re-train classifiers with the updated pipeline."
        ) from exc
    if schema_version != int(STATE_CLASSIFIER_PIPELINE_SCHEMA_VERSION):
        raise ValueError(
            f"{artifact_name} has unsupported pipeline_metadata schema_version={schema_version}; "
            f"expected {STATE_CLASSIFIER_PIPELINE_SCHEMA_VERSION}. "
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
        "schema_version": int(schema_version),
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


def _require_columns(df, required_cols, context):
    missing = [c for c in required_cols if c not in df.columns]
    if len(missing) > 0:
        raise ValueError(
            f"Missing required columns for {context}: {missing}. "
            f"Available columns include: {list(df.columns[:20])}"
        )


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


def train_random_forest_classifier(
    X_train,
    y_train,
    random_state=123,
    n_estimators=300,
    min_samples_leaf=2,
    n_jobs=-1,
    max_depth=None,
    min_samples_split=2,
    max_features="sqrt",
    class_weight=None,
):
    """Train a RandomForestClassifier and return (classifier, fit_info)."""
    X_train = np.asarray(X_train, dtype=float)
    y_train = np.asarray(y_train).astype(str)
    if X_train.ndim != 2:
        raise ValueError(f"X_train must be 2D, got shape {X_train.shape}.")
    if len(y_train) == 0:
        raise ValueError("No rows available to train RandomForest classifier.")
    if X_train.shape[0] != len(y_train):
        raise ValueError(
            f"X/y row mismatch for RandomForest training: X rows={X_train.shape[0]}, y rows={len(y_train)}."
        )

    params = {
        "n_estimators": int(n_estimators),
        "min_samples_leaf": int(min_samples_leaf),
        "n_jobs": int(n_jobs),
        "random_state": int(random_state),
        "max_depth": None if max_depth is None else int(max_depth),
        "min_samples_split": int(min_samples_split),
        "max_features": max_features,
        "class_weight": class_weight,
    }

    clf = RandomForestClassifier(**params, oob_score=True)
    clf.fit(X_train, y_train)
    train_acc = float((clf.predict(X_train) == y_train).mean())
    fit_info = {
        "backend": "random_forest",
        "n_train_rows": int(X_train.shape[0]),
        "n_features": int(X_train.shape[1]),
        "classes": [str(c) for c in getattr(clf, "classes_", [])],
        "train_accuracy": train_acc,
        "params": {k: repr(v) for k, v in clf.get_params().items()},
        "oob_score": float(getattr(clf, "oob_score_", np.nan)),
    }
    return clf, fit_info


def train_clusterid_classifier_from_model_adata(
    model_adata,
    cluster_col="intrinsic_behavioral_cluster",
    classifier_backend="random_forest",
    random_state=123,
    n_estimators=300,
    min_samples_leaf=2,
    n_jobs=-1,
    max_depth=None,
    min_samples_split=2,
    max_features="sqrt",
    class_weight=None,
    verbose=True,
):
    """Train a supervised classifier on model_adata features -> cluster labels."""
    if cluster_col not in model_adata.obs.columns:
        raise ValueError(f"Missing '{cluster_col}' in model_adata.obs.")

    feature_cols = list(model_adata.var_names)
    X_train = _to_numpy_2d(model_adata[:, feature_cols].X)
    y_train = model_adata.obs[cluster_col].astype(str).to_numpy()

    backend = str(classifier_backend).lower()
    if backend not in {"random_forest", "rf", "randomforest", "randomforestclassifier"}:
        raise ValueError("Only RandomForestClassifier is currently supported (classifier_backend='random_forest').")

    clf, fit_info = train_random_forest_classifier(
        X_train=X_train,
        y_train=y_train,
        random_state=random_state,
        n_estimators=n_estimators,
        min_samples_leaf=min_samples_leaf,
        n_jobs=n_jobs,
        max_depth=max_depth,
        min_samples_split=min_samples_split,
        max_features=max_features,
        class_weight=class_weight,
    )
    fit_info["cluster_col"] = cluster_col
    _vinfo(
        verbose,
        "state-training",
        (
            "Trained intrinsic classifier "
            f"(backend={fit_info['backend']}, rows={fit_info['n_train_rows']}, "
            f"features={fit_info['n_features']}, train_acc={fit_info['train_accuracy']:.4f})"
        ),
    )
    return clf, fit_info


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


def _mixed_label_sort_key(x):
    s = str(x)
    return (0, int(s)) if s.isdigit() else (1, s)


def _resolve_validation_mode(validation_test_size):
    """Resolve validation mode from user-specified test size."""
    if validation_test_size is None or float(validation_test_size) == 0.0:
        return "oob_only"
    test_size = float(validation_test_size)
    if not (0.0 < test_size < 1.0):
        raise ValueError(
            "validation_test_size must be in (0, 1), or set to 0/None for OOB-only mode."
        )
    return "holdout"


def _validate_stratified_split_feasibility(y, test_size, target_name):
    """Raise ValueError if stratified holdout split cannot be constructed."""
    y = np.asarray(y).astype(str)
    class_counts = pd.Series(y).value_counts(dropna=False)
    n_samples = int(len(y))
    n_classes = int(class_counts.shape[0])

    if n_classes < 2:
        raise ValueError(
            f"Cannot do stratified split for '{target_name}': need at least 2 classes, got {n_classes}."
        )
    if int(class_counts.min()) < 2:
        raise ValueError(
            f"Cannot do stratified split for '{target_name}': each class needs >=2 rows; "
            f"min class count={int(class_counts.min())}. Class counts={class_counts.to_dict()}"
        )

    test_size = float(test_size)
    n_test = int(np.ceil(n_samples * test_size))
    n_train = n_samples - n_test
    if n_test < n_classes:
        raise ValueError(
            f"Cannot do stratified split for '{target_name}': test rows ({n_test}) "
            f"must be >= number of classes ({n_classes}). "
            f"n_samples={n_samples}, test_size={test_size}, class_counts={class_counts.to_dict()}"
        )
    if n_train < n_classes:
        raise ValueError(
            f"Cannot do stratified split for '{target_name}': train rows ({n_train}) "
            f"must be >= number of classes ({n_classes}). "
            f"n_samples={n_samples}, test_size={test_size}, class_counts={class_counts.to_dict()}"
        )


def _build_compact_metric_summary(y_true, y_pred):
    """Build compact scalar classification metrics."""
    y_true = np.asarray(y_true).astype(str)
    y_pred = np.asarray(y_pred).astype(str)
    label_set = pd.Index(np.concatenate([y_true, y_pred])).unique()
    return {
        "n_rows": int(len(y_true)),
        "n_classes": int(len(label_set)),
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "macro_f1": float(f1_score(y_true, y_pred, average="macro", zero_division=0)),
    }


def _evaluate_holdout_predictions(
    *,
    y_true,
    y_pred,
    outfolder=None,
    filename_prefix="classifier_holdout_qc",
    row_normalized_decimals=5,
):
    """Evaluate holdout predictions and return compact metrics + saved artifacts."""
    summary = evaluate_classifier_predictions(
        y_true=y_true,
        y_pred=y_pred,
        outfolder=outfolder,
        filename_prefix=filename_prefix,
        row_normalized_decimals=row_normalized_decimals,
    )
    compact = _build_compact_metric_summary(y_true, y_pred)
    return {
        "metrics": compact,
        "artifacts": dict(summary.get("artifacts", {})),
        "summary": summary,
    }


def evaluate_classifier_predictions(
    y_true,
    y_pred,
    outfolder=None,
    filename_prefix="classifier_qc",
    row_normalized_decimals=5,
):
    """Evaluate predicted labels and optionally save confusion/report artifacts."""
    y_true = np.asarray(y_true).astype(str)
    y_pred = np.asarray(y_pred).astype(str)
    label_set = pd.Index(np.concatenate([y_true, y_pred])).unique().tolist()
    labels = sorted([str(x) for x in label_set], key=_mixed_label_sort_key)

    cm_counts = confusion_matrix(y_true, y_pred, labels=labels)
    cm_true_norm = confusion_matrix(y_true, y_pred, labels=labels, normalize="true")
    report_dict = classification_report(y_true, y_pred, labels=labels, output_dict=True, zero_division=0)

    summary = {
        "n_rows": int(len(y_true)),
        "n_classes": int(len(labels)),
        "classes": [str(c) for c in labels],
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "artifacts": {},
    }

    if outfolder is not None:
        outfolder = Path(outfolder)
        outfolder.mkdir(parents=True, exist_ok=True)

        cm_counts_df = pd.DataFrame(cm_counts, index=labels, columns=labels)
        cm_norm_df = pd.DataFrame(cm_true_norm, index=labels, columns=labels)
        report_df = pd.DataFrame(report_dict).T

        cm_counts_csv = outfolder / f"{filename_prefix}_confusion_counts.csv"
        cm_norm_csv = outfolder / f"{filename_prefix}_confusion_true_normalized.csv"
        report_csv = outfolder / f"{filename_prefix}_classification_report.csv"
        cm_pdf = outfolder / f"{filename_prefix}_confusion_matrices.pdf"

        cm_counts_df.to_csv(cm_counts_csv)
        cm_norm_df.to_csv(cm_norm_csv)
        report_df.to_csv(report_csv)

        fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)

        im0 = axes[0].imshow(cm_counts, cmap="Blues", aspect="auto")
        axes[0].set_title("Confusion Matrix (Counts)")
        axes[0].set_xlabel("Predicted label")
        axes[0].set_ylabel("True label")
        axes[0].set_xticks(np.arange(len(labels)))
        axes[0].set_yticks(np.arange(len(labels)))
        axes[0].set_xticklabels(labels, rotation=90)
        axes[0].set_yticklabels(labels)
        fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

        im1 = axes[1].imshow(cm_true_norm, cmap="Blues", aspect="auto", vmin=0.0, vmax=1.0)
        axes[1].set_title("Confusion Matrix (Row-normalized)")
        axes[1].set_xlabel("Predicted label")
        axes[1].set_ylabel("True label")
        axes[1].set_xticks(np.arange(len(labels)))
        axes[1].set_yticks(np.arange(len(labels)))
        axes[1].set_xticklabels(labels, rotation=90)
        axes[1].set_yticklabels(labels)
        fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

        if len(labels) <= 20:
            for i in range(len(labels)):
                for j in range(len(labels)):
                    axes[0].text(j, i, str(int(cm_counts[i, j])), ha="center", va="center", fontsize=7, color="black")
                    axes[1].text(
                        j,
                        i,
                        f"{cm_true_norm[i, j]:.{int(row_normalized_decimals)}f}",
                        ha="center",
                        va="center",
                        fontsize=7,
                        color="black",
                    )

        fig.suptitle("Classifier QC", fontsize=12)
        fig.savefig(cm_pdf, dpi=300, bbox_inches="tight")
        plt.close(fig)

        summary["artifacts"] = {
            "confusion_counts_csv": str(cm_counts_csv),
            "confusion_true_normalized_csv": str(cm_norm_csv),
            "classification_report_csv": str(report_csv),
            "confusion_matrices_pdf": str(cm_pdf),
        }

    return summary


def evaluate_clusterid_classifier_on_model_adata(
    model_adata,
    classifier,
    cluster_col="intrinsic_behavioral_cluster",
    outfolder=None,
    filename_prefix="behavioral_states_classifier_qc",
    row_normalized_decimals=5,
    verbose=True,
):
    """Evaluate classifier quality on model_adata and optionally save QC outputs."""
    if cluster_col not in model_adata.obs.columns:
        raise ValueError(f"Missing '{cluster_col}' in model_adata.obs.")

    feature_cols = list(model_adata.var_names)
    X_eval = _to_numpy_2d(model_adata[:, feature_cols].X)
    y_true = model_adata.obs[cluster_col].astype(str).to_numpy()
    y_pred = np.asarray(classifier.predict(X_eval)).astype(str)

    summary = evaluate_classifier_predictions(
        y_true=y_true,
        y_pred=y_pred,
        outfolder=outfolder,
        filename_prefix=filename_prefix,
        row_normalized_decimals=row_normalized_decimals,
    )

    _vinfo(
        verbose,
        "state-training",
        (
            "Intrinsic model-data QC: "
            f"n={summary['n_rows']}, classes={summary['n_classes']}, "
            f"accuracy={summary['accuracy']:.4f}, balanced_accuracy={summary['balanced_accuracy']:.4f}"
        ),
    )

    return summary


def _resolve_binary_classifier_feature_cols(adata, binary_cols_to_merge):
    """Resolve binary classifier features strictly from supplied binary columns."""
    cols = []
    for col in list(binary_cols_to_merge):
        if col in adata.obs.columns:
            cols.append(col)
    # Preserve order while removing duplicates.
    return list(dict.fromkeys(cols))


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


def _build_classifier_matrix_from_adata(
    adata,
    continuous_feature_cols,
    binary_feature_cols,
    binary_expanded_feature_cols=None,
    return_binary_feature_names=False,
    row_indices=None,
):
    """Build a numeric feature matrix from adata.X (continuous) and obs (binary)."""
    cont_cols = list(continuous_feature_cols)
    bin_cols = list(binary_feature_cols)

    missing_cont = [c for c in cont_cols if c not in adata.var_names]
    if len(missing_cont) > 0:
        raise ValueError(f"Missing continuous feature columns in adata.var_names: {missing_cont[:10]}")

    if row_indices is None:
        adata_cont = adata[:, cont_cols]
        obs_df = adata.obs
    else:
        row_indices = np.asarray(row_indices)
        adata_cont = adata[row_indices, cont_cols]
        if np.issubdtype(row_indices.dtype, np.integer):
            obs_df = adata.obs.iloc[row_indices]
        else:
            obs_df = adata.obs.loc[row_indices]

    X_cont = _to_numpy_2d(adata_cont.X).astype(float, copy=False)
    if len(bin_cols) == 0:
        if return_binary_feature_names:
            return X_cont, []
        return X_cont

    bin_df = _build_binary_feature_dataframe(
        obs_df=obs_df,
        binary_feature_cols=bin_cols,
        expected_expanded_cols=binary_expanded_feature_cols,
    )
    X_bin = bin_df.to_numpy(dtype=float)
    X = np.hstack([X_cont, X_bin])
    if return_binary_feature_names:
        return X, list(bin_df.columns)
    return X


def train_full_clusterid_classifier_from_adata(
    adata,
    target_col="full_behavioral_cluster",
    binary_feature_cols=None,
    preprocessing_artifact=None,
    random_state=123,
    n_estimators=500,
    min_samples_leaf=2,
    n_jobs=-1,
    max_depth=None,
    min_samples_split=2,
    max_features="sqrt",
    class_weight=None,
    verbose=True,
):
    """
    Train final random-forest classifier using continuous + binary features.

    Returns classifier artifact and fit metadata.
    """
    if target_col not in adata.obs.columns:
        raise ValueError(f"Missing target column '{target_col}' in adata.obs.")

    continuous_feature_cols = list(adata.var_names)
    binary_cols = [] if binary_feature_cols is None else list(binary_feature_cols)
    X_train, expanded_binary_cols = _build_classifier_matrix_from_adata(
        adata=adata,
        continuous_feature_cols=continuous_feature_cols,
        binary_feature_cols=binary_cols,
        return_binary_feature_names=True,
    )
    y_train = adata.obs[target_col].astype(str).to_numpy()
    if len(y_train) == 0:
        raise ValueError("No rows available to train full ClusterID classifier.")

    clf, fit_info = train_random_forest_classifier(
        X_train=X_train,
        y_train=y_train,
        random_state=random_state,
        n_estimators=n_estimators,
        min_samples_leaf=min_samples_leaf,
        n_jobs=n_jobs,
        max_depth=max_depth,
        min_samples_split=min_samples_split,
        max_features=max_features,
        class_weight=class_weight,
    )

    artifact = {
        "classifier": clf,
        "continuous_feature_cols": list(continuous_feature_cols),
        "binary_feature_cols": list(binary_cols),
        "binary_expanded_feature_cols": list(expanded_binary_cols),
        "target_col": str(target_col),
        "preprocessing": preprocessing_artifact,
    }
    fit_info["target_col"] = str(target_col)
    fit_info["continuous_feature_cols"] = list(continuous_feature_cols)
    fit_info["binary_feature_cols"] = list(binary_cols)
    fit_info["n_features_total"] = int(X_train.shape[1])
    fit_info["n_features_continuous"] = int(len(continuous_feature_cols))
    fit_info["n_features_binary"] = int(len(binary_cols))
    fit_info["n_features_binary_expanded"] = int(len(expanded_binary_cols))
    _vinfo(
        verbose,
        "state-training",
        (
            "Trained full classifier "
            f"(rows={fit_info['n_train_rows']}, total_features={fit_info['n_features_total']}, "
            f"train_acc={fit_info['train_accuracy']:.4f})"
        ),
    )
    return artifact, fit_info


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


def train_full_classifier_on_labeled_adata(
    adata,
    binary_features_to_group=None,
    preprocessing_artifact=None,
    random_state=123,
    n_estimators=500,
    min_samples_leaf=2,
    n_jobs=-1,
    max_depth=None,
    min_samples_split=2,
    max_features="sqrt",
    class_weight=None,
    output_col="full_behavioral_cluster",
    confidence_col="full_behavioral_cluster_confidence",
    outfolder=None,
    qc_filename_prefix="behavioral_states_full_classifier_qc_random_forest",
    windowing_artifact=None,
    verbose=True,
):
    """
    Train/apply/evaluate the full classifier (continuous + binary) on a labeled AnnData.

    This lets you run full-classifier training separately from the continuous-only stage.
    """
    preprocessing_params = preprocessing_artifact
    if preprocessing_params is None:
        preprocessing_params = _build_preprocessing_params_from_uns(
            adata.uns.get("preprocessing", {}) if isinstance(adata.uns, dict) else None,
            feature_cols=list(adata.var_names),
        )
    windowing_params = windowing_artifact
    if binary_features_to_group is None:
        binary_features_to_group = []
        if isinstance(adata.uns, dict):
            binary_features_to_group = list(
                adata.uns.get("clustering", {}).get("binary_cols_to_merge", [])
            )

    full_label_classifier_artifact, full_classifier_fit_info = train_full_clusterid_classifier_from_adata(
        adata=adata,
        target_col="full_behavioral_cluster",
        binary_feature_cols=binary_features_to_group,
        preprocessing_artifact=preprocessing_params,
        random_state=random_state,
        n_estimators=n_estimators,
        min_samples_leaf=min_samples_leaf,
        n_jobs=n_jobs,
        max_depth=max_depth,
        min_samples_split=min_samples_split,
        max_features=max_features,
        class_weight=class_weight,
        verbose=verbose,
    )
    if windowing_params is not None:
        full_label_classifier_artifact["windowing"] = dict(windowing_params)

    predict_clusterids_with_full_classifier(
        adata=adata,
        classifier_artifact=full_label_classifier_artifact,
        output_col=output_col,
        confidence_col=confidence_col,
        inplace=True,
        apply_preprocessing=False,
    )
    y_true = adata.obs["full_behavioral_cluster"].astype(str).to_numpy()
    y_pred = adata.obs[output_col].astype(str).to_numpy()
    full_classifier_qc_model_data = evaluate_classifier_predictions(
        y_true=y_true,
        y_pred=y_pred,
        outfolder=(Path(outfolder) if outfolder is not None else None),
        filename_prefix=qc_filename_prefix,
    )
    full_cm_pdf = full_classifier_qc_model_data.get("artifacts", {}).get("confusion_matrices_pdf", None)
    if full_cm_pdf is not None:
        _vsave(verbose, "state-training", "full-classifier confusion matrices", full_cm_pdf)
    _vinfo(
        verbose,
        "state-training",
        (
            "Full model-data QC: "
            f"n={full_classifier_qc_model_data['n_rows']}, "
            f"accuracy={full_classifier_qc_model_data['accuracy']:.4f}, "
            f"balanced_accuracy={full_classifier_qc_model_data['balanced_accuracy']:.4f}"
        ),
    )
    if "classification" not in adata.uns:
        adata.uns["classification"] = {}
    adata.uns["classification"]["full_classifier_fit_info"] = full_classifier_fit_info
    adata.uns["classification"]["full_classifier_qc_model_data"] = full_classifier_qc_model_data
    adata.uns["classification"]["full_classifier_prediction_col"] = str(output_col)
    adata.uns["classification"]["full_classifier_confidence_col"] = (
        None if confidence_col is None else str(confidence_col)
    )
    return full_label_classifier_artifact, full_classifier_fit_info, full_classifier_qc_model_data


def plot_binary_group_behavioral_cluster_grid(
    adata,
    binary_group_col="binary_group",
    behavioral_cluster_col="behavioral_clusterid",
    ncols=1,
    figsize_per_plot=(10.0, 1.8),
    csv_path=None,
    pdf_path=None,
    return_csv=True,
    verbose=True,
):
    """
    Plot one horizontal stacked bar per binary group, showing proportions of behavioral clusters.
    """
    obs = adata.obs.copy()
    if binary_group_col not in obs.columns:
        raise ValueError(f"Missing '{binary_group_col}' in adata.obs.")
    if behavioral_cluster_col not in obs.columns:
        raise ValueError(f"Missing '{behavioral_cluster_col}' in adata.obs.")

    plot_df = obs[[binary_group_col, behavioral_cluster_col]].copy()
    plot_df[binary_group_col] = plot_df[binary_group_col].astype("string").fillna("unassigned")
    plot_df[behavioral_cluster_col] = plot_df[behavioral_cluster_col].astype("string").fillna("unassigned")

    # Keep all behavioral clusters present in the dataset (even if concentrated in
    # one binary group), so no cluster silently disappears in the grid.
    cluster_counts = plot_df[behavioral_cluster_col].value_counts(dropna=False)
    all_clusters = cluster_counts.index.tolist()

    ctab = pd.crosstab(
        plot_df[binary_group_col],
        plot_df[behavioral_cluster_col],
        normalize="index",
        dropna=False,
    ).reindex(columns=all_clusters, fill_value=0.0)

    def _sort_key(x):
        s = str(x)
        return (0, int(s)) if s.isdigit() else (1, s)

    group_names = sorted(ctab.index.tolist(), key=_sort_key)
    cluster_names = sorted(ctab.columns.tolist(), key=_sort_key)

    n_groups = len(group_names)
    if n_groups == 0:
        raise ValueError("No groups available for plotting.")

    nrows = int(np.ceil(n_groups / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(figsize_per_plot[0] * ncols, figsize_per_plot[1] * nrows),
        squeeze=False,
    )
    axes = axes.flatten()

    # Build a larger distinct palette by combining qualitative colormaps.
    palette = []
    for cmap_name in ("tab20", "tab20b", "tab20c"):
        cmap = plt.get_cmap(cmap_name)
        palette.extend([cmap(i) for i in range(cmap.N)])
    if len(cluster_names) > len(palette):
        cmap = plt.get_cmap("hsv")
        palette.extend([cmap(i / max(len(cluster_names), 1)) for i in range(len(cluster_names) - len(palette))])
    colors = {cl: palette[i] for i, cl in enumerate(cluster_names)}

    for i, grp in enumerate(group_names):
        ax = axes[i]
        left = 0.0
        for cl in cluster_names:
            val = float(ctab.loc[grp, cl])
            if val <= 0:
                continue
            ax.barh([0], [val], left=left, color=colors[cl], height=0.92, edgecolor="none", linewidth=0.0)
            left += val
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.7, 0.7)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.grid(False)
        ax.text(
            -0.01,
            0.5,
            str(grp),
            transform=ax.transAxes,
            ha="right",
            va="center",
            fontsize=10,
        )
        grp_n = int((plot_df[binary_group_col] == grp).sum())
        ax.set_title(f"n={grp_n}", fontsize=9, loc="right")

    for j in range(n_groups, len(axes)):
        axes[j].axis("off")

    handles = [
        plt.Line2D([0], [0], marker="s", linestyle="none", color=colors[cl], label=str(cl), markersize=7)
        for cl in cluster_names
    ]
    fig.legend(
        handles=handles,
        labels=[str(c) for c in cluster_names],
        title="behavioral_clusterid",
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        borderaxespad=0.0,
    )
    fig.suptitle("Behavioral Cluster Composition per Binary Group", fontsize=12, y=0.995)
    fig.tight_layout(rect=[0, 0, 0.82, 0.97])

    if verbose and "5" in ctab.columns:
        by_group = (ctab["5"] * 100.0).round(3)
        nonzero = int((by_group > 0).sum())
        _vinfo(
            verbose,
            "state-apply",
            f"Cluster 5 composition summary: total_n={int(cluster_counts.get('5', 0))}, nonzero_groups={nonzero}",
        )

    if csv_path is None and pdf_path is not None:
        csv_path = Path(pdf_path).with_suffix(".csv")
    if csv_path is not None:
        ctab.to_csv(csv_path)
    if pdf_path is not None:
        fig.savefig(pdf_path, dpi=300, bbox_inches="tight")

    if return_csv:
        return fig, ctab
    return fig


def plot_behavioral_cluster_binary_group_grid(
    adata,
    binary_group_col="binary_group",
    behavioral_cluster_col="behavioral_clusterid",
    ncols=1,
    figsize_per_plot=(10.0, 1.8),
    csv_path=None,
    pdf_path=None,
    return_csv=True,
):
    """
    Plot one horizontal stacked bar per behavioral cluster, showing proportions of binary groups.
    """
    obs = adata.obs.copy()
    if binary_group_col not in obs.columns:
        raise ValueError(f"Missing '{binary_group_col}' in adata.obs.")
    if behavioral_cluster_col not in obs.columns:
        raise ValueError(f"Missing '{behavioral_cluster_col}' in adata.obs.")

    plot_df = obs[[binary_group_col, behavioral_cluster_col]].copy()
    plot_df[binary_group_col] = plot_df[binary_group_col].astype("string").fillna("unassigned")
    plot_df[behavioral_cluster_col] = plot_df[behavioral_cluster_col].astype("string").fillna("unassigned")

    # Keep all binary groups present in the dataset so no group disappears.
    group_counts = plot_df[binary_group_col].value_counts(dropna=False)
    all_groups = group_counts.index.tolist()

    ctab = pd.crosstab(
        plot_df[behavioral_cluster_col],
        plot_df[binary_group_col],
        normalize="index",
        dropna=False,
    ).reindex(columns=all_groups, fill_value=0.0)

    def _sort_key(x):
        s = str(x)
        return (0, int(s)) if s.isdigit() else (1, s)

    cluster_names = sorted(ctab.index.tolist(), key=_sort_key)
    group_names = sorted(ctab.columns.tolist(), key=_sort_key)

    n_clusters = len(cluster_names)
    if n_clusters == 0:
        raise ValueError("No behavioral clusters available for plotting.")

    nrows = int(np.ceil(n_clusters / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(figsize_per_plot[0] * ncols, figsize_per_plot[1] * nrows),
        squeeze=False,
    )
    axes = axes.flatten()

    palette = []
    for cmap_name in ("tab20", "tab20b", "tab20c"):
        cmap = plt.get_cmap(cmap_name)
        palette.extend([cmap(i) for i in range(cmap.N)])
    if len(group_names) > len(palette):
        cmap = plt.get_cmap("hsv")
        palette.extend([cmap(i / max(len(group_names), 1)) for i in range(len(group_names) - len(palette))])
    colors = {grp: palette[i] for i, grp in enumerate(group_names)}

    for i, cl in enumerate(cluster_names):
        ax = axes[i]
        left = 0.0
        for grp in group_names:
            val = float(ctab.loc[cl, grp])
            if val <= 0:
                continue
            ax.barh([0], [val], left=left, color=colors[grp], height=0.92, edgecolor="none", linewidth=0.0)
            left += val
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.7, 0.7)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.grid(False)
        ax.text(
            -0.01,
            0.5,
            str(cl),
            transform=ax.transAxes,
            ha="right",
            va="center",
            fontsize=10,
        )
        cl_n = int((plot_df[behavioral_cluster_col] == cl).sum())
        ax.set_title(f"n={cl_n}", fontsize=9, loc="right")

    for j in range(n_clusters, len(axes)):
        axes[j].axis("off")

    handles = [
        plt.Line2D([0], [0], marker="s", linestyle="none", color=colors[grp], label=str(grp), markersize=7)
        for grp in group_names
    ]
    fig.legend(
        handles=handles,
        labels=[str(g) for g in group_names],
        title="binary_group",
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        borderaxespad=0.0,
    )
    fig.suptitle("Binary Group Composition per Behavioral Cluster", fontsize=12, y=0.995)
    fig.tight_layout(rect=[0, 0, 0.82, 0.97])

    if csv_path is None and pdf_path is not None:
        csv_path = Path(pdf_path).with_suffix(".csv")
    if csv_path is not None:
        ctab.to_csv(csv_path)
    if pdf_path is not None:
        fig.savefig(pdf_path, dpi=300, bbox_inches="tight")

    if return_csv:
        return fig, ctab
    return fig


def run_state_clustering(
    features,
    binary_features_to_group,
    output_dir,
    cell_type,
    window_size=5,
    min_spacing=None,
    max_samples=None,
    n_neighbors=60,
    min_dist=0.1,
    resolution=0.2,
    descriptive_features=("mean", "median", "std", "net_displacement", "straightness", "mean_square_displacement"),
    pca_var_selection=0.95,
    clustering_method="leiden",
    lower_quantile_cap=None,
    upper_quantile_cap=0.99,
    incomplete_window_policy="partial",
    random_state=123,
    reuse_prepared_dataset=True,
    plot_exemplar_videos=True,
    exemplar_n_per_cluster=5,
    exemplar_video_fps=6,
    exemplar_video_dpi=180,
    exemplar_video_margin=20,
    exemplar_video_pmin=0.0,
    exemplar_video_pmax=99.99,
    exemplar_video_track_color="#63ff33",
    exemplar_video_show_segments=True,
    exemplar_video_segment_style="outline",
    exemplar_video_segment_color="#ffffff",
    # save_prepared_dataset=True,
    df_positions=None,
    state_paths=None,
    verbose=True,
):
    """Stage 1: load/prepare positions, cluster, return model_adata.

    Stage-1 metadata is stored in:
    - model_adata.uns["preprocessing"]
    - model_adata.uns["clustering"]
    - model_adata.uns["preprocessing"]["prepared_adata_full_*"] cache metadata
    """
    run_started = _vstart(verbose, "state-clustering", "run state clustering")
    state_paths = state_paths or _resolve_state_paths(output_dir, cell_type)
    state_outdir = state_paths.state_outdir
    state_clustering_outdir = _resolve_state_clustering_outdir(state_paths.processing_outdir)
    prepared_adata_full_path = state_paths.prepared_full_adata_path
    legacy_prepared_adata_full_path = _resolve_legacy_prepared_full_adata_path(
        state_paths=state_paths,
        cell_type=cell_type,
    )
    prepared_adata_read_path = None
    if prepared_adata_full_path is not None and prepared_adata_full_path.exists():
        prepared_adata_read_path = prepared_adata_full_path
    elif legacy_prepared_adata_full_path.exists():
        prepared_adata_read_path = legacy_prepared_adata_full_path
    prepared_adata_model_path = state_paths.model_adata_path
    df_tracks_path = None
    stage_prepare_started = time.perf_counter()
    adata_prepared = None
    prepared_cache_used = False
    positions_loaded_from_csv = False

    def _load_positions_if_needed():
        nonlocal df_positions, df_tracks_path, positions_loaded_from_csv
        if df_positions is None:
            df_tracks_path = _resolve_positions_csv_path(state_paths=state_paths)
            df_positions = pd.read_csv(df_tracks_path)
            positions_loaded_from_csv = True
        _require_columns(
            df_positions,
            ["sample_name", "TrackID", "position_t"] + list(features),
            context="state clustering input",
        )
        return df_positions

    if bool(reuse_prepared_dataset) and prepared_adata_read_path is not None:
        try:
            cached_prepared = sc.read_h5ad(prepared_adata_read_path)
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
                scale_features=True,
            )
            if pre_match_cached:
                adata_prepared = cached_prepared
                prepared_cache_used = True
                _vinfo(verbose, "state-clustering", f"reusing prepared full adata cache: {prepared_adata_read_path}")
                if prepared_adata_read_path != prepared_adata_full_path:
                    _vinfo(
                        verbose,
                        "state-clustering",
                        f"using legacy prepared cache for read; future writes use {prepared_adata_full_path}",
                    )
            else:
                _vinfo(
                    verbose,
                    "state-clustering",
                    "prepared full cache mismatch; rebuilding from raw positions | "
                    f"{_short_reasons(pre_reasons_cached)}",
                )
        except Exception as exc:
            _vinfo(verbose, "state-clustering", f"could not load prepared full cache; rebuilding ({exc})")

    if adata_prepared is None:
        _load_positions_if_needed()
        _vinfo(
            verbose,
            "state-clustering",
            "preparing dataset (windowed feature extraction, quantile capping, scaling)",
        )
        if positions_loaded_from_csv:
            _vinfo(verbose, "state-clustering", f"input track-features CSV: {df_tracks_path}")
        else:
            _vinfo(verbose, "state-clustering", "input source: provided df_positions DataFrame")

        prepare_kwargs = {
            "df_positions": df_positions,
            "features": features,
            "binary_features_to_group": binary_features_to_group,
            "window_size": window_size,
            "min_spacing": min_spacing,
            "descriptive_features": list(descriptive_features),
            "lower_quantile_cap": lower_quantile_cap,
            "upper_quantile_cap": upper_quantile_cap,
            "outfolder": state_outdir,
            "scale_features": True,
            "incomplete_window_policy": incomplete_window_policy,
            "reuse_prepared_dataset": reuse_prepared_dataset,
            "prepared_dataset_path": prepared_adata_full_path,
            "verbose": verbose,
        }
        adata_prepared = prepare_state_classification_dataset(**prepare_kwargs)

        pre_match, pre_mismatch_reasons = _matches_requested_preprocessing_in_adata(
            adata_prepared,
            features=features,
            binary_features_to_group=binary_features_to_group,
            window_size=window_size,
            min_spacing=min_spacing,
            descriptive_features=descriptive_features,
            incomplete_window_policy=incomplete_window_policy,
            lower_quantile_cap=lower_quantile_cap,
            upper_quantile_cap=upper_quantile_cap,
            scale_features=True,
        )
        if not pre_match:
            _vinfo(
                verbose,
                "state-clustering",
                "prepared full cache mismatch; recomputing full preprocessing | "
                f"{_short_reasons(pre_mismatch_reasons)}",
            )
            prepare_kwargs["reuse_prepared_dataset"] = False
            adata_prepared = prepare_state_classification_dataset(**prepare_kwargs)
        del df_positions
    else:
        _vinfo(
            verbose,
            "state-clustering",
            "preparing dataset (windowed feature extraction, quantile capping, scaling)",
        )
    _vdone(verbose, "state-clustering", "prepare/reuse full dataset", stage_prepare_started)
    
    pre_meta = (
        adata_prepared.uns.get("preprocessing", {})
        if isinstance(getattr(adata_prepared, "uns", {}), dict)
        else {}
    )
    kept_features = pre_meta.get("kept_features", list(adata_prepared.var_names))
    if not isinstance(kept_features, list):
        kept_features = list(adata_prepared.var_names)
    kept_features = [str(c) for c in kept_features if str(c) in adata_prepared.var_names]
    if len(kept_features) == 0:
        kept_features = [str(c) for c in list(adata_prepared.var_names)]

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
        binary_cols_to_merge = [c for c in list(binary_features_to_group) if c in adata_prepared.obs.columns]
    prepared_preprocessing_meta = (
        dict(adata_prepared.uns.get("preprocessing", {}))
        if isinstance(getattr(adata_prepared, "uns", {}), dict)
        else {}
    )

    if adata_prepared.n_obs < 50:
        raise ValueError("Insufficient rows in full descriptive dataset for clustering.")
    _vinfo(
        verbose,
        "state-clustering",
        "prepared dataset summary | "
        f"rows={adata_prepared.n_obs}, continuous_features={len(kept_features)}, "
        f"binary_group_features={len(binary_cols_to_merge)}",
    )

    if min_spacing is None:
        spacing_to_use = int(window_size)
    else:
        spacing_to_use = int(min_spacing)
    sampling_id_cols = ["sample_name", "TrackID"]
    sampling_time_col = "position_t"
    sampling_start_cutoff_mode = "time_from_track_start"
    sampling_start_cutoff_window_size = None if window_size is None else int(window_size)
    sampling_min_time_from_track_start = (
        None
        if sampling_start_cutoff_window_size is None
        else max(int(sampling_start_cutoff_window_size) - 1, 0)
    )
    nan_policy = "drop_any_nan_in_kept_features"

    stage_sampling_started = _vstart(
        verbose,
        "state-clustering",
        f"sampling + NaN cleanup | min_spacing={spacing_to_use} | max_samples={max_samples}",
    )
    model_adata = None
    model_cache_loaded = False
    if (
        bool(reuse_prepared_dataset)
        and prepared_adata_model_path is not None
        and prepared_adata_model_path.exists()
    ):
        try:
            cached_model_adata = sc.read_h5ad(prepared_adata_model_path)
            pre_match_model, pre_reasons_model = _matches_requested_preprocessing_in_adata(
                cached_model_adata,
                features=features,
                binary_features_to_group=binary_features_to_group,
                window_size=window_size,
                min_spacing=min_spacing,
                descriptive_features=descriptive_features,
                incomplete_window_policy=incomplete_window_policy,
                lower_quantile_cap=lower_quantile_cap,
                upper_quantile_cap=upper_quantile_cap,
                scale_features=True,
            )
            sampling_match, sampling_reasons = _matches_model_sampling_cache(
                cached_model_adata,
                kept_features=kept_features,
                min_spacing_used=spacing_to_use,
                start_cutoff_mode=sampling_start_cutoff_mode,
                start_cutoff_window_size=sampling_start_cutoff_window_size,
                max_samples=max_samples,
                random_state=random_state,
                id_cols=sampling_id_cols,
                time_col=sampling_time_col,
                nan_policy=nan_policy,
            )
            if pre_match_model and sampling_match:
                model_adata = cached_model_adata
                model_cache_loaded = True
                _vinfo(
                    verbose,
                    "state-clustering",
                    f"reusing cached model adata: {prepared_adata_model_path}",
                )
            else:
                mismatch_reasons = []
                if not pre_match_model:
                    mismatch_reasons.extend([f"preprocessing: {r}" for r in pre_reasons_model])
                if not sampling_match:
                    mismatch_reasons.extend([f"sampling: {r}" for r in sampling_reasons])
                _vinfo(
                    verbose,
                    "state-clustering",
                    "cached model adata mismatch; rebuilding sampling/clean stage | "
                    f"{_short_reasons(mismatch_reasons)}",
                )
        except Exception as exc:
            _vinfo(verbose, "state-clustering", f"could not load cached model adata; rebuilding ({exc})")

    subsampled_rows_for_log = None
    rows_after_nan_for_pca = None
    if model_adata is None:
        adata_train = subsample_with_temporal_spacing(
            adata_prepared,
            id_cols=sampling_id_cols,
            time_col=sampling_time_col,
            min_spacing=spacing_to_use,
            min_time_from_track_start=sampling_min_time_from_track_start,
            max_samples=max_samples,
            random_state=random_state,
        )
        
        if adata_train.n_obs < 50:
            raise ValueError("Insufficient rows after global subsampling for clustering.")
        subsampled_rows_for_log = int(adata_train.n_obs)

        X_train = _to_numpy_2d(adata_train[:, kept_features].X).astype(float, copy=False)
        valid_mask = ~np.isnan(X_train).any(axis=1)
        n_valid = int(valid_mask.sum())
        if n_valid < 50:
            raise ValueError("Insufficient rows after dropping NaNs for PCA stage.")
        rows_after_nan_for_pca = int(n_valid)
        if n_valid != adata_train.n_obs:
            _vinfo(verbose, "state-clustering", f"dropped NaN rows before PCA: kept {n_valid}/{adata_train.n_obs}")
        model_adata = adata_train[valid_mask, kept_features].copy()
        del X_train, valid_mask, adata_train
    else:
        if model_adata.n_obs < 50:
            raise ValueError("Cached model dataset has insufficient rows for clustering.")
        subsampled_rows_for_log = int(model_adata.n_obs)
        rows_after_nan_for_pca = int(model_adata.n_obs)

    ncomps_requested = min(len(kept_features), model_adata.n_obs - 1)
    if ncomps_requested < 2:
        raise ValueError("Insufficient rows/features to run PCA stage.")

    _vdone(verbose, "state-clustering", "sampling + NaN cleanup", stage_sampling_started)
    _vinfo(verbose, "state-clustering", f"subsampled rows={subsampled_rows_for_log} (spacing={spacing_to_use})")

    pca_input_shape = (int(model_adata.n_obs), int(len(kept_features)))
    stage_pca_started = _vstart(
        verbose,
        "state-clustering",
        (
            "PCA stage | "
            f"matrix_shape={pca_input_shape}, kept_features={len(kept_features)}, "
            f"rows_after_nan={rows_after_nan_for_pca}, min_var_selection={pca_var_selection}"
        ),
    )
    pca_recomputed = False
    neighbors_recomputed = False

    if model_cache_loaded:
        pca_match, pca_reasons = _matches_cached_pca_stage(
            model_adata,
            pca_var_selection=pca_var_selection,
            svd_solver="full",
            random_state=random_state,
            ncomps_requested=ncomps_requested,
        )
        if pca_match:
            _vinfo(verbose, "state-clustering", "reusing cached PCA stage")
        else:
            _vinfo(
                verbose,
                "state-clustering",
                f"PCA recompute required | reasons={_short_reasons(pca_reasons)}",
            )
            model_adata = run_pca(
                model_adata,
                pca_var_selection=pca_var_selection,
                ncomps=ncomps_requested,
                svd_solver="full",
                random_state=random_state,
            )
            pca_recomputed = True
    else:
        model_adata = run_pca(
            model_adata,
            pca_var_selection=pca_var_selection,
            ncomps=ncomps_requested,
            svd_solver="full",
            random_state=random_state,
        )
        pca_recomputed = True
    _vdone(verbose, "state-clustering", "PCA stage", stage_pca_started)

    stage_neighbors_started = _vstart(
        verbose,
        "state-clustering",
        f"neighbors stage | n_neighbors={n_neighbors}",
    )
    if model_cache_loaded:
        neighbors_match, neighbors_reasons = _matches_cached_neighbors_stage(
            model_adata,
            n_neighbors=n_neighbors,
            random_state=random_state,
            use_rep="X_pca",
        )
        if pca_recomputed:
            neighbors_match = False
            neighbors_reasons = list(neighbors_reasons) + ["upstream PCA recomputed"]
        if neighbors_match:
            _vinfo(verbose, "state-clustering", f"reusing cached neighbors graph (n_neighbors={n_neighbors})")
        else:
            _vinfo(
                verbose,
                "state-clustering",
                f"neighbors recompute required | reasons={_short_reasons(neighbors_reasons)}",
            )
            sc.pp.neighbors(model_adata, n_neighbors=n_neighbors, use_rep="X_pca", random_state=random_state)
            neighbors_recomputed = True
    else:
        sc.pp.neighbors(model_adata, n_neighbors=n_neighbors, use_rep="X_pca", random_state=random_state)
        neighbors_recomputed = True
    _vdone(verbose, "state-clustering", "neighbors stage", stage_neighbors_started)

    stage_umap_started = _vstart(
        verbose,
        "state-clustering",
        f"UMAP stage | min_dist={min_dist}, random_state={random_state}",
    )
    if model_cache_loaded:
        umap_match, umap_reasons = _matches_cached_umap_stage(
            model_adata,
            min_dist=min_dist,
            random_state=random_state,
        )
        if pca_recomputed or neighbors_recomputed:
            umap_match = False
            umap_reasons = list(umap_reasons) + ["upstream representation/graph recomputed"]
        if umap_match:
            _vinfo(verbose, "state-clustering", f"reusing cached UMAP embedding (min_dist={min_dist})")
        else:
            _vinfo(
                verbose,
                "state-clustering",
                f"UMAP recompute required | reasons={_short_reasons(umap_reasons)}",
            )
            sc.tl.umap(
                model_adata, 
                min_dist=min_dist, 
                random_state=random_state,
                n_components=2
                )
    else:
        sc.tl.umap(
            model_adata, 
            min_dist=min_dist, 
            random_state=random_state,
            n_components=2
            )
    _vdone(verbose, "state-clustering", "UMAP stage", stage_umap_started)

    stage_clustering_started = _vstart(
        verbose,
        "state-clustering",
        (
            f"{clustering_method} clustering | "
            f"n_neighbors={n_neighbors}, min_dist={min_dist}, resolution={resolution}"
        ),
    )
    if isinstance(resolution, (list, tuple, np.ndarray, pd.Series)):
        stability_resolutions = [float(r) for r in list(resolution)]
    else:
        stability_resolutions = [float(resolution)]
    method = str(clustering_method).lower()
    if method == "leiden":
        if ("neighbors" not in model_adata.uns) or ("connectivities" not in model_adata.obsp):
            raise ValueError(
                "Leiden clustering requires a precomputed neighbors graph when n_neighbors=None. "
                "Expected adata.uns['neighbors'] and adata.obsp['connectivities']."
            )
        model_adata = run_leiden_clustering(
            model_adata,
            n_neighbors=None,
            resolution=resolution,
            use_rep="X_pca",
            metric="euclidean",
            method="umap",
            key_added="intrinsic_behavioral_cluster",
            random_state=random_state,
            stability_resolutions=stability_resolutions,
        )
        model_adata.obs["intrinsic_behavioral_cluster"] = model_adata.obs["intrinsic_behavioral_cluster"].astype("category")
    elif method == "kmeans":
        k = int(resolution) if isinstance(resolution, (int, float)) and float(resolution) >= 2 else 5
        labels = KMeans(n_clusters=k, random_state=random_state, n_init="auto").fit_predict(model_adata.obsm["X_pca"])
        model_adata.obs["intrinsic_behavioral_cluster"] = pd.Categorical(
            pd.Series(labels.astype(str), index=model_adata.obs.index)
        )
    else:
        raise ValueError("clustering_method must be one of: 'leiden', 'kmeans'.")

    binary_group_constraints = _infer_binary_group_constraints(
        model_adata.obs,
        binary_cols_to_merge,
    )
    _vinfo(
        verbose,
        "state-clustering",
        (
            "binary-group constraints inferred | "
            f"allowed_groups={len(binary_group_constraints.get('allowed_binary_groups', []))}, "
            f"forbidden_combinations={len(binary_group_constraints.get('forbidden_binary_combinations', []))}"
        ),
    )
    _rebuild_full_behavioral_cluster_from_intrinsic(
        adata=model_adata,
        binary_cols_to_merge=binary_cols_to_merge,
        intrinsic_col="intrinsic_behavioral_cluster",
        binary_group_constraints=binary_group_constraints,
        enforce_binary_group_constraints=True,
    )
    _vdone(verbose, "state-clustering", "clustering stage", stage_clustering_started)

    stage_prepared_write_started = time.perf_counter()
    prepared_adata_full = adata_prepared[:, kept_features].copy()
    prepared_adata_full.uns["preprocessing"] = dict(prepared_preprocessing_meta)
    prepared_adata_full.uns["preprocessing"].pop("prepared_adata_full_signature", None)
    prepared_adata_full.uns["preprocessing"]["prepared_adata_full_path"] = str(prepared_adata_full_path)
    prepared_adata_full.write(prepared_adata_full_path, compression="gzip")
    del prepared_adata_full, adata_prepared
    _vdone(verbose, "state-clustering", "write prepared full adata cache", stage_prepared_write_started)

    n_intrinsic = int(model_adata.obs["intrinsic_behavioral_cluster"].astype(str).nunique())
    n_full = int(model_adata.obs["full_behavioral_cluster"].astype(str).nunique())
    _vinfo(
        verbose,
        "state-clustering",
        "stage-1 clustering summary | "
        f"intrinsic_clusters={n_intrinsic}, full_behavioral_clusters={n_full}, model_rows={model_adata.n_obs}",
    )

    diagnostics_csvs = {}
    diagnostics_pdf = None
    if state_clustering_outdir is not None:
        diagnostics_pdf = state_clustering_outdir / "behavioral_clustering_diagnostics.pdf"
        plot_clustering_diagnostics_pdf(
            adata=model_adata,
            cluster_col="intrinsic_behavioral_cluster",
            feature_cols=kept_features,
            pdf_path=diagnostics_pdf,
            title=f"all_data | {clustering_method} (resolution={resolution})",
        )
        diagnostics_csvs = _export_clustering_diagnostics_csvs(
            adata=model_adata,
            cluster_col="intrinsic_behavioral_cluster",
            feature_cols=kept_features,
            outdir=state_clustering_outdir,
            filename_prefix="behavioral_clustering_diagnostics",
        )
        if ("binary_group" in model_adata.obs.columns) and ("behavioral_clusterid" in model_adata.obs.columns):
            pdf1 = state_clustering_outdir / "behavioral_clustering_binary_group_cluster_proportions.pdf"
            fig1, _ = plot_binary_group_behavioral_cluster_grid(
                model_adata,
                pdf_path=pdf1,
                return_csv=True,
                verbose=verbose,
            )
            plt.close(fig1)
            pdf2 = state_clustering_outdir / "behavioral_clustering_cluster_binary_group_proportions.pdf"
            fig2, _ = plot_behavioral_cluster_binary_group_grid(
                model_adata,
                pdf_path=pdf2,
                return_csv=True,
            )
            plt.close(fig2)

    ncomps_effective = (
        int(model_adata.obsm["X_pca"].shape[1])
        if ("X_pca" in model_adata.obsm and model_adata.obsm["X_pca"] is not None)
        else 0
    )
    model_cache_meta = {
        "sampling": {
            "min_spacing_used": int(spacing_to_use),
            "max_samples": None if max_samples is None else int(max_samples),
            "random_state": int(random_state),
            "time_col": str(sampling_time_col),
            "id_cols": [str(c) for c in sampling_id_cols],
            "start_cutoff_mode": str(sampling_start_cutoff_mode),
            "start_cutoff_window_size": (
                None
                if sampling_start_cutoff_window_size is None
                else int(sampling_start_cutoff_window_size)
            ),
        },
        "cleaning": {
            "nan_policy": str(nan_policy),
        },
        "features": {
            "kept_features": [str(c) for c in kept_features],
        },
        "pca": {
            "pca_var_selection": float(pca_var_selection),
            "svd_solver": "full",
            "random_state": int(random_state),
            "ncomps_requested": int(ncomps_requested),
            "ncomps_effective": int(ncomps_effective),
        },
        "neighbors": {
            "n_neighbors": int(n_neighbors),
            "use_rep": "X_pca",
            "random_state": int(random_state),
        },
        "umap": {
            "min_dist": float(min_dist),
            "random_state": int(random_state),
        },
    }
    model_preprocessing_meta = dict(prepared_preprocessing_meta)
    model_preprocessing_meta["model_cache"] = dict(model_cache_meta)
    model_preprocessing_meta.pop("prepared_adata_full_signature", None)
    model_preprocessing_meta["prepared_adata_full_path"] = (
        None if prepared_adata_full_path is None else str(prepared_adata_full_path)
    )
    model_adata.uns["preprocessing"] = model_preprocessing_meta

    model_adata.uns["clustering"] = {
        "clustering_method": clustering_method,
        "resolution": resolution,
        "n_neighbors": int(n_neighbors),
        "random_state": int(random_state),
        "non_feature_cols": list(non_feature_cols),
        "binary_cols_to_merge": list(binary_cols_to_merge),
        "binary_group_constraints": copy.deepcopy(dict(binary_group_constraints)),
        "diagnostics_pdf": None if diagnostics_pdf is None else str(diagnostics_pdf),
        "diagnostics_csvs": dict(diagnostics_csvs),
    }

    exemplar_render_result = {
        "enabled": bool(plot_exemplar_videos),
        "cluster_key": "intrinsic_behavioral_cluster",
        "selection_csv": None,
        "video_paths_by_cluster": {},
        "global_xy_shape": None,
        "cluster_total_counts": {},
        "error": None,
        "render_config": {
            "stage": "after_clustering",
            "cluster_key": "intrinsic_behavioral_cluster",
            "state_key": "intrinsic_behavioral_cluster",
            "n_per_cluster": int(exemplar_n_per_cluster),
            "selection_random_state": int(random_state),
            "layout_mode": "per_cluster",
            "video_fps": int(exemplar_video_fps),
            "video_dpi": int(exemplar_video_dpi),
            "video_margin": int(exemplar_video_margin),
            "video_pmin": float(exemplar_video_pmin),
            "video_pmax": float(exemplar_video_pmax),
            "video_track_color": str(exemplar_video_track_color),
            "video_show_segments": bool(exemplar_video_show_segments),
            "video_segment_style": str(exemplar_video_segment_style),
            "video_segment_color": str(exemplar_video_segment_color),
            "show_state_legend": False,
            "positions_csv_path": None,
        },
    }
    if bool(plot_exemplar_videos):
        exemplar_render_result = _render_state_cluster_exemplar_videos(
            model_adata=model_adata,
            state_paths=state_paths,
            cell_type=cell_type,
            state_clustering_outdir=state_clustering_outdir,
            cluster_key="intrinsic_behavioral_cluster",
            n_per_cluster=int(exemplar_n_per_cluster),
            random_state=int(random_state),
            video_fps=int(exemplar_video_fps),
            video_dpi=int(exemplar_video_dpi),
            video_margin=int(exemplar_video_margin),
            video_pmin=float(exemplar_video_pmin),
            video_pmax=float(exemplar_video_pmax),
            video_track_color=str(exemplar_video_track_color),
            video_show_segments=bool(exemplar_video_show_segments),
            video_segment_style=str(exemplar_video_segment_style),
            video_segment_color=str(exemplar_video_segment_color),
            verbose=verbose,
        )

    model_adata.uns.setdefault("visualization", {})
    model_adata.uns["visualization"].update(
        {
            "exemplar_render_stage": "after_clustering",
            "exemplar_cluster_key": str(exemplar_render_result.get("cluster_key", "intrinsic_behavioral_cluster")),
            "exemplar_selection_csv": exemplar_render_result.get("selection_csv", None),
            "exemplar_video_paths_by_cluster": dict(
                exemplar_render_result.get("video_paths_by_cluster", {})
            ),
            "exemplar_backprojection_video_paths_by_cluster": dict(
                exemplar_render_result.get("video_paths_by_cluster", {})
            ),
            "exemplar_video_canvas_shape": exemplar_render_result.get(
                "global_xy_shape", None
            ),
            "exemplar_cluster_total_counts": dict(
                exemplar_render_result.get("cluster_total_counts", {})
            ),
            "exemplar_render_config": dict(
                exemplar_render_result.get("render_config", {})
            ),
            "exemplar_segment_outline_errors": dict(
                exemplar_render_result.get("segment_outline_errors", {})
            ),
            "exemplar_render_error": exemplar_render_result.get("error", None),
        }
    )

    stage_model_write_started = time.perf_counter()
    model_adata.write(prepared_adata_model_path, compression="gzip")
    _vdone(verbose, "state-clustering", "write model adata cache", stage_model_write_started)
    _vsave(verbose, "state-clustering", "model-stage cache", prepared_adata_model_path)
    _vsave(verbose, "state-clustering", "prepared full adata cache", prepared_adata_full_path)
    _vdone(verbose, "state-clustering", "run state clustering", run_started)
    return model_adata


def train_state_classifiers(
    output_dir,
    cell_type,
    model_adata=None,
    label_transfer_method="classifier",
    classifier_backend="random_forest",
    classifier_n_estimators=300,
    classifier_min_samples_leaf=2,
    classifier_n_jobs=-1,
    classifier_max_depth=None,
    classifier_min_samples_split=2,
    classifier_max_features="sqrt",
    classifier_class_weight=None,
    classifier_confidence_col="intrinsic_behavioral_cluster_confidence",
    save_label_classifier=True,
    train_continuous_classifier=True,
    train_full_classifier=True,
    full_classifier_confidence_col="full_behavioral_cluster_confidence",
    save_full_label_classifier=True,
    df_prepared=None,
    validation_test_size=0.05,
    validation_random_state=None,
    validation_stratify=True,
    state_paths=None,
    verbose=True,
):
    """
    Stage 2: train intrinsic/full classifiers using model_adata (auto-loaded from canonical path if needed).

    Mutates model_adata in-place by adding/updating model_adata.uns["classification"].

    Returns:
        dict: {
            "full_classifier_path": str | None,
            "partial_classifier_path": str | None,
        }
    """
    train_started = _vstart(verbose, "state-training", "train state classifiers")
    state_paths = state_paths or _resolve_state_paths(output_dir, cell_type)
    model_adata_path = state_paths.model_adata_path
    if model_adata is None:
        load_started = _vstart(verbose, "state-training", "load model adata")
        if not model_adata_path.exists():
            raise FileNotFoundError(
                "Could not load model_adata for classifier training. "
                f"Expected '{model_adata_path}'. "
                "Run run_state_clustering(...) first or pass model_adata explicitly."
            )
        model_adata = sc.read_h5ad(model_adata_path)
        _vdone(verbose, "state-training", "load model adata", load_started)
        _vinfo(verbose, "state-training", f"Using model adata with rows={model_adata.n_obs}, features={len(model_adata.var_names)}")
    else:
        _vinfo(verbose, "state-training", f"Using provided model adata with rows={model_adata.n_obs}, features={len(model_adata.var_names)}")

    if isinstance(model_adata, dict):
        raise ValueError(
            "train_state_classifiers now expects model_adata (AnnData), not a stage-1 artifact dict."
        )
    if not (hasattr(model_adata, "uns") and hasattr(model_adata, "var_names")):
        raise ValueError("model_adata must be an AnnData-like object with .uns and .var_names.")
    pre_meta = model_adata.uns.get("preprocessing", None)
    clust_meta = model_adata.uns.get("clustering", None)
    if not isinstance(pre_meta, dict):
        raise ValueError(
            "model_adata.uns['preprocessing'] is missing or invalid. "
            "Run run_state_clustering first."
        )
    if not isinstance(clust_meta, dict):
        raise ValueError(
            "model_adata.uns['clustering'] is missing or invalid. "
            "Run run_state_clustering first."
        )

    windowing_params = pre_meta.get("windowing", None)
    if not isinstance(windowing_params, dict):
        raise ValueError("model_adata.uns['preprocessing']['windowing'] must be a dict.")

    transfer_mode = str(label_transfer_method).lower()
    if transfer_mode != "classifier":
        raise ValueError(
            "train_state_classifiers is training-only and supports only "
            "label_transfer_method='classifier'. Use apply_state_classifiers_to_full_dataset for full-data inference."
        )

    classifier_backend_norm = str(classifier_backend).lower()
    if classifier_backend_norm not in {"random_forest", "rf", "randomforest", "randomforestclassifier"}:
        raise ValueError("Only RandomForestClassifier is currently supported (classifier_backend='random_forest').")
    classifier_backend_name = "random_forest"

    cont_cols = list(model_adata.var_names)
    preprocessing_params = _build_preprocessing_params_from_uns(
        pre_meta,
        feature_cols=cont_cols,
    )
    if preprocessing_params is not None and not isinstance(preprocessing_params, dict):
        raise ValueError("model_adata.uns['preprocessing'] could not be converted to preprocessing params.")

    required_clustering_keys = ["non_feature_cols", "binary_cols_to_merge", "resolution", "n_neighbors", "random_state"]
    missing_clustering_keys = [k for k in required_clustering_keys if k not in clust_meta]
    if len(missing_clustering_keys) > 0:
        raise ValueError(
            "model_adata.uns['clustering'] is missing required keys: "
            f"{missing_clustering_keys}"
        )

    non_feature_cols = clust_meta.get("non_feature_cols", None)
    binary_cols_to_merge = clust_meta.get("binary_cols_to_merge", None)
    if not isinstance(non_feature_cols, list):
        raise ValueError("model_adata.uns['clustering']['non_feature_cols'] must be a list.")
    if not isinstance(binary_cols_to_merge, list):
        raise ValueError("model_adata.uns['clustering']['binary_cols_to_merge'] must be a list.")
    non_feature_cols = [str(c) for c in non_feature_cols]
    binary_cols_to_merge = [str(c) for c in binary_cols_to_merge]
    clustering_method = clust_meta.get("clustering_method", None)
    clustering_resolution = clust_meta.get("resolution", None)
    try:
        clustering_n_neighbors = int(clust_meta.get("n_neighbors", 0))
    except Exception as exc:
        raise ValueError(
            "model_adata.uns['clustering']['n_neighbors'] must be int-castable."
        ) from exc
    try:
        stage1_random_state = int(clust_meta.get("random_state", 123))
    except Exception as exc:
        raise ValueError(
            "model_adata.uns['clustering']['random_state'] must be int-castable."
        ) from exc

    validation_mode = _resolve_validation_mode(validation_test_size)
    if validation_mode == "holdout":
        validation_test_size = float(validation_test_size)
    else:
        validation_test_size = 0.0
    validation_random_state = (
        stage1_random_state if validation_random_state is None else int(validation_random_state)
    )
    validation_stratify = bool(validation_stratify)

    required_intrinsic_cols = ["intrinsic_behavioral_cluster"]
    missing_intrinsic_cols = [c for c in required_intrinsic_cols if c not in model_adata.obs.columns]
    if len(missing_intrinsic_cols) > 0:
        raise ValueError(
            "model_adata is missing required intrinsic label columns for training: "
            f"{missing_intrinsic_cols}. Run run_state_clustering before train_state_classifiers."
        )
    if train_full_classifier:
        required_full_cols = ["full_behavioral_cluster"]
        missing_full_cols = [c for c in required_full_cols if c not in model_adata.obs.columns]
        if len(missing_full_cols) > 0:
            raise ValueError(
                "model_adata is missing required full label columns for training: "
                f"{missing_full_cols}. Run run_state_clustering (and optionally "
                "rename_intrinsic_behavioral_clusters) before train_state_classifiers."
            )

    missing_binary_cols = [c for c in binary_cols_to_merge if c not in model_adata.obs.columns]
    if train_full_classifier and len(missing_binary_cols) > 0:
        raise ValueError(
            "model_adata is missing binary columns required for full-classifier training: "
            f"{missing_binary_cols}. Re-run run_state_clustering so binary columns are in model_adata.obs."
        )

    label_classifier = None  # canonical train-split intrinsic model
    label_classifier_artifact_train_split = None
    classifier_fit_info = None  # canonical train-split fit info
    classifier_fit_info_train_split = None
    classifier_qc_model_data = None
    classifier_model_path = None
    classifier_model_paths = {}
    intrinsic_validation = {
        "target_col": "intrinsic_behavioral_cluster",
        "n_train": 0,
        "n_test": 0,
        "train_metrics": None,
        "test_metrics": None,
        "oob_score_train_split": None,
        "oob_score_full_data": None,
        "holdout_artifacts": {},
    }

    outfolder_path = state_paths.processing_outdir
    outfolder_path.mkdir(parents=True, exist_ok=True)
    intrinsic_outfolder = state_paths.intrinsic_outdir
    full_outfolder = state_paths.full_outdir
    intrinsic_qc_outfolder = state_paths.intrinsic_qc_outdir
    full_qc_outfolder = state_paths.full_qc_outdir
    intrinsic_outfolder.mkdir(parents=True, exist_ok=True)
    full_outfolder.mkdir(parents=True, exist_ok=True)
    intrinsic_qc_outfolder.mkdir(parents=True, exist_ok=True)
    full_qc_outfolder.mkdir(parents=True, exist_ok=True)

    if train_continuous_classifier:
        intrinsic_started = _vstart(verbose, "state-training", "train intrinsic classifier")
        X_intrinsic_all = _to_numpy_2d(model_adata[:, cont_cols].X).astype(float, copy=False)
        y_intrinsic_all = model_adata.obs["intrinsic_behavioral_cluster"].astype(str).to_numpy()
        X_intrinsic_train = X_intrinsic_all
        y_intrinsic_train = y_intrinsic_all
        if validation_mode == "holdout":
            if validation_stratify:
                _validate_stratified_split_feasibility(
                    y_intrinsic_all, validation_test_size, target_name="intrinsic_behavioral_cluster"
                )
            idx_all = np.arange(model_adata.n_obs, dtype=int)
            idx_train, idx_test = train_test_split(
                idx_all,
                test_size=validation_test_size,
                random_state=validation_random_state,
                stratify=(y_intrinsic_all if validation_stratify else None),
            )
            X_intrinsic_train = X_intrinsic_all[idx_train]
            y_intrinsic_train = y_intrinsic_all[idx_train]
            X_intrinsic_test = X_intrinsic_all[idx_test]
            y_intrinsic_test = y_intrinsic_all[idx_test]
            intrinsic_validation["n_train"] = int(len(idx_train))
            intrinsic_validation["n_test"] = int(len(idx_test))
        else:
            intrinsic_validation["n_train"] = int(model_adata.n_obs)
            intrinsic_validation["n_test"] = 0

        label_classifier, classifier_fit_info = train_random_forest_classifier(
            X_train=X_intrinsic_train,
            y_train=y_intrinsic_train,
            random_state=stage1_random_state,
            n_estimators=classifier_n_estimators,
            min_samples_leaf=classifier_min_samples_leaf,
            n_jobs=classifier_n_jobs,
            max_depth=classifier_max_depth,
            min_samples_split=classifier_min_samples_split,
            max_features=classifier_max_features,
            class_weight=classifier_class_weight,
        )
        classifier_fit_info["cluster_col"] = "intrinsic_behavioral_cluster"
        classifier_fit_info_train_split = classifier_fit_info

        label_classifier_artifact_train_split = {
            "classifier": label_classifier,
            "backend": classifier_backend_name,
            "feature_cols": cont_cols,
            "cluster_col": "intrinsic_behavioral_cluster",
            "preprocessing": preprocessing_params,
            "windowing": windowing_params,
        }

        y_pred_intrinsic_train = np.asarray(label_classifier.predict(X_intrinsic_train)).astype(str)
        intrinsic_validation["train_metrics"] = _build_compact_metric_summary(
            y_intrinsic_train, y_pred_intrinsic_train
        )
        intrinsic_validation["oob_score_train_split"] = float(
            classifier_fit_info.get("oob_score", np.nan)
        )

        if validation_mode == "holdout":
            y_pred_intrinsic_test = np.asarray(label_classifier.predict(X_intrinsic_test)).astype(str)
            holdout_intrinsic_eval = _evaluate_holdout_predictions(
                y_true=y_intrinsic_test,
                y_pred=y_pred_intrinsic_test,
                outfolder=intrinsic_qc_outfolder,
                filename_prefix=f"state_classification_intrinsic_classifier_holdout_qc_{classifier_backend_name}",
            )
            intrinsic_validation["test_metrics"] = holdout_intrinsic_eval["metrics"]
            intrinsic_validation["holdout_artifacts"] = holdout_intrinsic_eval["artifacts"]

        classifier_qc_model_data = evaluate_clusterid_classifier_on_model_adata(
            model_adata=model_adata,
            classifier=label_classifier,
            cluster_col="intrinsic_behavioral_cluster",
            outfolder=intrinsic_qc_outfolder,
            filename_prefix=f"state_classification_classifier_qc_{classifier_backend_name}",
            verbose=verbose,
        )
        intrinsic_validation["oob_score_full_data"] = (
            intrinsic_validation["oob_score_train_split"] if validation_mode != "holdout" else None
        )

        if validation_mode == "holdout":
            tm = intrinsic_validation["train_metrics"] or {}
            vm = intrinsic_validation["test_metrics"] or {}
            _vinfo(
                verbose,
                "state-training",
                (
                    "Intrinsic validation_set scores: "
                    f"n_val={intrinsic_validation['n_test']}, "
                    f"val_acc={vm.get('accuracy', np.nan):.4f}, "
                    f"val_bal_acc={vm.get('balanced_accuracy', np.nan):.4f}, "
                    f"val_macro_f1={vm.get('macro_f1', np.nan):.4f}"
                ),
            )
            _vinfo(
                verbose,
                "state-training",
                (
                    "Intrinsic training/OOB scores: "
                    f"n_train={intrinsic_validation['n_train']}, "
                    f"train_acc={tm.get('accuracy', np.nan):.4f}, "
                    f"train_bal_acc={tm.get('balanced_accuracy', np.nan):.4f}, "
                    f"train_macro_f1={tm.get('macro_f1', np.nan):.4f}, "
                    f"oob_train_split={intrinsic_validation.get('oob_score_train_split', np.nan):.4f}"
                ),
            )
        else:
            tm = intrinsic_validation["train_metrics"] or {}
            _vinfo(
                verbose,
                "state-training",
                (
                    "Intrinsic training/OOB scores: "
                    f"n_train={intrinsic_validation['n_train']}, "
                    f"train_acc={tm.get('accuracy', np.nan):.4f}, "
                    f"train_bal_acc={tm.get('balanced_accuracy', np.nan):.4f}, "
                    f"train_macro_f1={tm.get('macro_f1', np.nan):.4f}, "
                    f"oob_train_split={intrinsic_validation.get('oob_score_train_split', np.nan):.4f}"
                ),
            )
        _vdone(verbose, "state-training", "train intrinsic classifier", intrinsic_started)
        del X_intrinsic_all, y_intrinsic_all, X_intrinsic_train, y_intrinsic_train
        if validation_mode == "holdout":
            del X_intrinsic_test, y_intrinsic_test

    full_label_classifier_artifact = None  # canonical train-split full-label artifact
    full_label_classifier_artifact_train_split = None
    full_classifier_fit_info = None  # canonical train-split full fit info
    full_classifier_fit_info_train_split = None
    full_classifier_qc_model_data = None
    full_classifier_model_path = None
    full_classifier_model_paths = {}
    full_validation = {
        "target_col": "full_behavioral_cluster",
        "n_train": 0,
        "n_test": 0,
        "train_metrics": None,
        "test_metrics": None,
        "oob_score_train_split": None,
        "oob_score_full_data": None,
        "holdout_artifacts": {},
    }

    if train_full_classifier:
        full_started = _vstart(verbose, "state-training", "train full classifier")
        y_full_all = model_adata.obs["full_behavioral_cluster"].astype(str).to_numpy()

        if validation_mode == "holdout":
            if validation_stratify:
                _validate_stratified_split_feasibility(
                    y_full_all, validation_test_size, target_name="full_behavioral_cluster"
                )
            idx_all = np.arange(model_adata.n_obs, dtype=int)
            idx_train, idx_test = train_test_split(
                idx_all,
                test_size=validation_test_size,
                random_state=validation_random_state,
                stratify=(y_full_all if validation_stratify else None),
            )

            X_full_train, expanded_binary_train_cols = _build_classifier_matrix_from_adata(
                adata=model_adata,
                continuous_feature_cols=cont_cols,
                binary_feature_cols=binary_cols_to_merge,
                return_binary_feature_names=True,
                row_indices=idx_train,
            )
            X_full_test = _build_classifier_matrix_from_adata(
                adata=model_adata,
                continuous_feature_cols=cont_cols,
                binary_feature_cols=binary_cols_to_merge,
                binary_expanded_feature_cols=expanded_binary_train_cols,
                return_binary_feature_names=False,
                row_indices=idx_test,
            )
            y_full_train = y_full_all[idx_train]
            y_full_test = y_full_all[idx_test]
            full_validation["n_train"] = int(len(idx_train))
            full_validation["n_test"] = int(len(idx_test))
        else:
            X_full_train, expanded_binary_train_cols = _build_classifier_matrix_from_adata(
                adata=model_adata,
                continuous_feature_cols=cont_cols,
                binary_feature_cols=binary_cols_to_merge,
                return_binary_feature_names=True,
            )
            y_full_train = y_full_all
            full_validation["n_train"] = int(model_adata.n_obs)
            full_validation["n_test"] = 0

        clf_full_train_split, fit_info_full_train_split = train_random_forest_classifier(
            X_train=X_full_train,
            y_train=y_full_train,
            random_state=stage1_random_state,
            n_estimators=classifier_n_estimators,
            min_samples_leaf=classifier_min_samples_leaf,
            n_jobs=classifier_n_jobs,
            max_depth=classifier_max_depth,
            min_samples_split=classifier_min_samples_split,
            max_features=classifier_max_features,
            class_weight=classifier_class_weight,
        )
        fit_info_full_train_split["target_col"] = "full_behavioral_cluster"
        fit_info_full_train_split["continuous_feature_cols"] = list(cont_cols)
        fit_info_full_train_split["binary_feature_cols"] = list(binary_cols_to_merge)
        fit_info_full_train_split["n_features_total"] = int(X_full_train.shape[1])
        fit_info_full_train_split["n_features_continuous"] = int(len(cont_cols))
        fit_info_full_train_split["n_features_binary"] = int(len(binary_cols_to_merge))
        fit_info_full_train_split["n_features_binary_expanded"] = int(len(expanded_binary_train_cols))
        full_classifier_fit_info_train_split = fit_info_full_train_split
        full_classifier_fit_info = fit_info_full_train_split

        full_label_classifier_artifact_train_split = {
            "classifier": clf_full_train_split,
            "continuous_feature_cols": list(cont_cols),
            "binary_feature_cols": list(binary_cols_to_merge),
            "binary_expanded_feature_cols": list(expanded_binary_train_cols),
            "target_col": "full_behavioral_cluster",
            "preprocessing": preprocessing_params,
            "windowing": dict(windowing_params),
        }
        full_label_classifier_artifact = full_label_classifier_artifact_train_split

        y_pred_full_train = np.asarray(clf_full_train_split.predict(X_full_train)).astype(str)
        full_validation["train_metrics"] = _build_compact_metric_summary(y_full_train, y_pred_full_train)
        full_validation["oob_score_train_split"] = float(
            fit_info_full_train_split.get("oob_score", np.nan)
        )

        if validation_mode == "holdout":
            y_pred_full_test = np.asarray(clf_full_train_split.predict(X_full_test)).astype(str)
            holdout_full_eval = _evaluate_holdout_predictions(
                y_true=y_full_test,
                y_pred=y_pred_full_test,
                outfolder=full_qc_outfolder,
                filename_prefix=f"state_classification_full_classifier_holdout_qc_{classifier_backend_name}",
            )
            full_validation["test_metrics"] = holdout_full_eval["metrics"]
            full_validation["holdout_artifacts"] = holdout_full_eval["artifacts"]

        if validation_mode == "holdout":
            X_full_all = _build_classifier_matrix_from_adata(
                adata=model_adata,
                continuous_feature_cols=cont_cols,
                binary_feature_cols=binary_cols_to_merge,
                binary_expanded_feature_cols=full_label_classifier_artifact_train_split.get(
                    "binary_expanded_feature_cols", []
                ),
                return_binary_feature_names=False,
            )
        else:
            X_full_all = X_full_train
        y_pred_full_all = np.asarray(
            full_label_classifier_artifact_train_split["classifier"].predict(X_full_all)
        ).astype(str)
        full_classifier_qc_model_data = evaluate_classifier_predictions(
            y_true=y_full_all,
            y_pred=y_pred_full_all,
            outfolder=full_qc_outfolder,
            filename_prefix="state_classification_full_classifier_qc_random_forest",
        )
        full_validation["oob_score_full_data"] = (
            full_validation["oob_score_train_split"] if validation_mode != "holdout" else None
        )
        _vinfo(
            verbose,
            "state-training",
            (
                "Full model-data QC: "
                f"n={full_classifier_qc_model_data['n_rows']}, "
                f"accuracy={full_classifier_qc_model_data['accuracy']:.4f}, "
                f"balanced_accuracy={full_classifier_qc_model_data['balanced_accuracy']:.4f}"
            ),
        )
        if validation_mode == "holdout":
            tm = full_validation["train_metrics"] or {}
            vm = full_validation["test_metrics"] or {}
            _vinfo(
                verbose,
                "state-training",
                (
                    "Full validation_set scores: "
                    f"n_val={full_validation.get('n_test', 0)}, "
                    f"val_acc={vm.get('accuracy', np.nan):.4f}, "
                    f"val_bal_acc={vm.get('balanced_accuracy', np.nan):.4f}, "
                    f"val_macro_f1={vm.get('macro_f1', np.nan):.4f}"
                ),
            )
            _vinfo(
                verbose,
                "state-training",
                (
                    "Full training/OOB scores: "
                    f"n_train={full_validation['n_train']}, "
                    f"train_acc={tm.get('accuracy', np.nan):.4f}, "
                    f"train_bal_acc={tm.get('balanced_accuracy', np.nan):.4f}, "
                    f"train_macro_f1={tm.get('macro_f1', np.nan):.4f}, "
                    f"oob_train_split={full_validation.get('oob_score_train_split', np.nan):.4f}"
                ),
            )
        else:
            tm = full_validation["train_metrics"] or {}
            _vinfo(
                verbose,
                "state-training",
                (
                    "Full training/OOB scores: "
                    f"n_train={full_validation['n_train']}, "
                    f"train_acc={tm.get('accuracy', np.nan):.4f}, "
                    f"train_bal_acc={tm.get('balanced_accuracy', np.nan):.4f}, "
                    f"train_macro_f1={tm.get('macro_f1', np.nan):.4f}, "
                    f"oob_train_split={full_validation.get('oob_score_train_split', np.nan):.4f}"
                ),
            )
        _vdone(verbose, "state-training", "train full classifier", full_started)
        del X_full_train, y_full_train, X_full_all, y_pred_full_train, y_pred_full_all
        if validation_mode == "holdout":
            del X_full_test, y_full_test, y_pred_full_test

    intrinsic_metrics_csv = None
    full_metrics_csv = None
    if intrinsic_qc_outfolder is not None and train_continuous_classifier:
        intrinsic_tm = intrinsic_validation.get("train_metrics") or {}
        intrinsic_vm = intrinsic_validation.get("test_metrics") or {}
        intrinsic_metrics_row = {
            "mode": validation_mode,
            "n_train": intrinsic_validation.get("n_train", 0),
            "n_test": intrinsic_validation.get("n_test", 0),
            "oob_score_train_split": intrinsic_validation.get("oob_score_train_split", np.nan),
            "train_accuracy": intrinsic_tm.get("accuracy", np.nan),
            "train_balanced_accuracy": intrinsic_tm.get("balanced_accuracy", np.nan),
            "train_macro_f1": intrinsic_tm.get("macro_f1", np.nan),
            "test_accuracy": intrinsic_vm.get("accuracy", np.nan),
            "test_balanced_accuracy": intrinsic_vm.get("balanced_accuracy", np.nan),
            "test_macro_f1": intrinsic_vm.get("macro_f1", np.nan),
            "model_data_accuracy": (
                None if classifier_qc_model_data is None else classifier_qc_model_data.get("accuracy", np.nan)
            ),
            "model_data_balanced_accuracy": (
                None if classifier_qc_model_data is None else classifier_qc_model_data.get("balanced_accuracy", np.nan)
            ),
        }
        intrinsic_metrics_csv = (
            intrinsic_qc_outfolder / f"state_classification_intrinsic_classifier_metrics_{classifier_backend_name}.csv"
        )
        pd.DataFrame([intrinsic_metrics_row]).to_csv(intrinsic_metrics_csv, index=False)
        _vsave(verbose, "state-training", "intrinsic metrics csv", intrinsic_metrics_csv)

    if full_qc_outfolder is not None and train_full_classifier:
        full_tm = full_validation.get("train_metrics") or {}
        full_vm = full_validation.get("test_metrics") or {}
        full_metrics_row = {
            "mode": validation_mode,
            "n_train": full_validation.get("n_train", 0),
            "n_test": full_validation.get("n_test", 0),
            "oob_score_train_split": full_validation.get("oob_score_train_split", np.nan),
            "train_accuracy": full_tm.get("accuracy", np.nan),
            "train_balanced_accuracy": full_tm.get("balanced_accuracy", np.nan),
            "train_macro_f1": full_tm.get("macro_f1", np.nan),
            "test_accuracy": full_vm.get("accuracy", np.nan),
            "test_balanced_accuracy": full_vm.get("balanced_accuracy", np.nan),
            "test_macro_f1": full_vm.get("macro_f1", np.nan),
            "model_data_accuracy": (
                None if full_classifier_qc_model_data is None else full_classifier_qc_model_data.get("accuracy", np.nan)
            ),
            "model_data_balanced_accuracy": (
                None if full_classifier_qc_model_data is None else full_classifier_qc_model_data.get("balanced_accuracy", np.nan)
            ),
        }
        full_metrics_csv = (
            full_qc_outfolder / f"state_classification_full_classifier_metrics_{classifier_backend_name}.csv"
        )
        pd.DataFrame([full_metrics_row]).to_csv(full_metrics_csv, index=False)
        _vsave(verbose, "state-training", "full metrics csv", full_metrics_csv)

    model_adata.uns["classification"] = {
        "label_transfer_method": transfer_mode,
        "classifier_backend": classifier_backend_name,
        "classifier_fit_info": classifier_fit_info,
        "classifier_fit_info_train_split": classifier_fit_info_train_split,
        "classifier_qc_model_data": classifier_qc_model_data,
        "classifier_confidence_col": classifier_confidence_col,
        "classifier_model_path": None,
        "classifier_model_paths": {},
        "full_classifier_fit_info": None,
        "full_classifier_fit_info_train_split": None,
        "full_classifier_qc_model_data": None,
        "full_classifier_prediction_col": "full_behavioral_cluster",
        "full_classifier_confidence_col": full_classifier_confidence_col,
        "full_classifier_model_path": None,
        "full_classifier_model_paths": {},
        "model_vs_full_clusterid_check": None,
        "intrinsic_artifacts_dir": None if intrinsic_outfolder is None else str(intrinsic_outfolder),
        "full_artifacts_dir": None if full_outfolder is None else str(full_outfolder),
        "intrinsic_qc_dir": None if intrinsic_qc_outfolder is None else str(intrinsic_qc_outfolder),
        "full_qc_dir": None if full_qc_outfolder is None else str(full_qc_outfolder),
        "intrinsic_metrics_csv": None if intrinsic_metrics_csv is None else str(intrinsic_metrics_csv),
        "full_metrics_csv": None if full_metrics_csv is None else str(full_metrics_csv),
        "validation_config": {
            "mode": validation_mode,
            "test_size": None if validation_mode == "oob_only" else float(validation_test_size),
            "random_state": int(validation_random_state),
            "stratify": bool(validation_stratify),
            "split_unit": "row",
        },
        "intrinsic_validation": intrinsic_validation,
        "full_validation": full_validation,
    }

    if train_full_classifier:
        model_adata.uns["classification"]["full_classifier_fit_info"] = full_classifier_fit_info
        model_adata.uns["classification"]["full_classifier_fit_info_train_split"] = (
            full_classifier_fit_info_train_split
        )
        model_adata.uns["classification"]["full_classifier_qc_model_data"] = full_classifier_qc_model_data

    def _build_pipeline_metadata_payload():
        return {
            "schema_version": int(STATE_CLASSIFIER_PIPELINE_SCHEMA_VERSION),
            "preprocessing": copy.deepcopy(dict(model_adata.uns.get("preprocessing", {}))),
            "clustering": copy.deepcopy(dict(model_adata.uns.get("clustering", {}))),
            "classification": copy.deepcopy(dict(model_adata.uns.get("classification", {}))),
        }

    if label_classifier_artifact_train_split is not None:
        label_classifier_artifact_train_split["pipeline_metadata"] = _build_pipeline_metadata_payload()
    if full_label_classifier_artifact_train_split is not None:
        full_label_classifier_artifact_train_split["pipeline_metadata"] = _build_pipeline_metadata_payload()

    if save_label_classifier and label_classifier_artifact_train_split is not None:
        classifier_model_path = state_paths.intrinsic_classifier_default_path
        classifier_model_path.parent.mkdir(parents=True, exist_ok=True)
        model_adata.uns["classification"]["classifier_model_path"] = str(classifier_model_path)
        classifier_model_paths["train_split"] = str(classifier_model_path)
        model_adata.uns["classification"]["classifier_model_paths"] = dict(classifier_model_paths)
        label_classifier_artifact_train_split["pipeline_metadata"] = _build_pipeline_metadata_payload()
        with open(classifier_model_path, "wb") as f:
            pickle.dump(label_classifier_artifact_train_split, f)
        _vsave(verbose, "state-training", "intrinsic classifier artifact", classifier_model_path)

    if train_full_classifier and save_full_label_classifier and full_label_classifier_artifact_train_split is not None:
        full_classifier_model_path = state_paths.full_classifier_default_path
        full_classifier_model_path.parent.mkdir(parents=True, exist_ok=True)
        model_adata.uns["classification"]["full_classifier_model_path"] = str(full_classifier_model_path)
        full_classifier_model_paths["train_split"] = str(full_classifier_model_path)
        model_adata.uns["classification"]["full_classifier_model_paths"] = dict(full_classifier_model_paths)
        full_label_classifier_artifact_train_split["pipeline_metadata"] = _build_pipeline_metadata_payload()
        with open(full_classifier_model_path, "wb") as f:
            pickle.dump(full_label_classifier_artifact_train_split, f)
        _vsave(verbose, "state-training", "full classifier artifact", full_classifier_model_path)

    write_started = _vstart(verbose, "state-training", "write model adata")
    model_adata.write(state_paths.model_adata_path, compression="gzip")
    _vdone(verbose, "state-training", "write model adata", write_started)
    _vsave(verbose, "state-training", "model adata", state_paths.model_adata_path)
    n_intrinsic = int(model_adata.obs["intrinsic_behavioral_cluster"].astype(str).nunique())
    n_full = int(model_adata.obs["full_behavioral_cluster"].astype(str).nunique()) if "full_behavioral_cluster" in model_adata.obs.columns else 0
    _vinfo(
        verbose,
        "state-training",
        f"Training summary: intrinsic_clusters={n_intrinsic}, full_behavioral_clusters={n_full}",
    )
    _vdone(verbose, "state-training", "train state classifiers", train_started)
    return {
        "full_classifier_path": None if full_classifier_model_path is None else str(full_classifier_model_path),
        "partial_classifier_path": None if classifier_model_path is None else str(classifier_model_path),
    }


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


def test_pipeline():
    model_adata = run_state_clustering(
        features=features,
        binary_features_to_group=binary_features_to_group,
        output_dir=output_dir,
        cell_type=cell_type,
        window_size=window_size,
        max_samples=max_samples,
        min_spacing=min_spacing,
        n_neighbors=n_neighbors,
        resolution=resolution,
        descriptive_features=descriptive_features,
        pca_var_selection=pca_var_selection,
        clustering_method=clustering_method,
        lower_quantile_cap=lower_quantile_cap,
        upper_quantile_cap=upper_quantile_cap,
        incomplete_window_policy=incomplete_window_policy,
        random_state=random_state,
        reuse_prepared_dataset=reuse_prepared_dataset,
        verbose=verbose,
    )
    
    model_adata.obs["intrinsic_behavioral_cluster"]
    model_adata.obs["full_behavioral_cluster"]
    
    mapping = {
        1: "dead",
        2: "static",
        3: "plastic_scanner",
        4: "static",
        5: "round_scanner",
    }
    rename_intrinsic_behavioral_clusters(
        adata=model_adata,
        mapping=mapping,
        binary_cols_to_merge=binary_features_to_group,
    )
    
    build_identity_cluster_mapping(
        model_adata,
        cluster_col="full_behavioral_cluster",
    )
    
    full_mapping = {
        'no_contact_dead': 'dead',
        'no_contact_plastic_scanner': 'plastic_scanner',
        'no_contact_round_scanner': 'round_scanner',
        'no_contact_static': 'static',
        'organoid_contact_pixels_and_tcell_contact_pixels_dead': 'dead',
        'organoid_contact_pixels_and_tcell_contact_pixels_plastic_scanner': 'organoid_contact_pixels',
        'organoid_contact_pixels_and_tcell_contact_pixels_round_scanner': 'organoid_contact_pixels',
        'organoid_contact_pixels_and_tcell_contact_pixels_static': 'organoid_contact_pixels',
        'organoid_contact_pixels_dead': 'dead',
        'organoid_contact_pixels_plastic_scanner': 'organoid_contact_pixels',
        'organoid_contact_pixels_round_scanner': 'organoid_contact_pixels',
        'organoid_contact_pixels_static': 'organoid_contact_pixels',
        'tcell_contact_pixels_dead': 'dead',
        'tcell_contact_pixels_plastic_scanner': 'tcell_contact_pixels',
        'tcell_contact_pixels_round_scanner': 'tcell_contact_pixels',
        'tcell_contact_pixels_static': 'tcell_contact_pixels'
        }
    
    model_adata = relabel_cluster_ids(
        adata=model_adata,
        mapping=full_mapping,
        cluster_key="full_behavioral_cluster",
        # Set overwrite_original=True to make reruns idempotent.
        overwrite_original=True,
    )
    
    classifier_paths = train_state_classifiers(
        output_dir=output_dir,
        cell_type=cell_type,
        model_adata=model_adata,
        label_transfer_method="classifier",
        classifier_backend=classifier_backend,
        classifier_n_estimators=classifier_n_estimators,
        classifier_min_samples_leaf=classifier_min_samples_leaf,
        classifier_n_jobs=classifier_n_jobs,
        classifier_max_depth=classifier_max_depth,
        classifier_min_samples_split=classifier_min_samples_split,
        classifier_max_features=classifier_max_features,
        classifier_class_weight=classifier_class_weight,
        classifier_confidence_col=classifier_confidence_col,
        save_label_classifier=save_label_classifier,
        train_continuous_classifier=True,
        train_full_classifier=True,
        )
    
    adata_full = apply_state_classifiers_to_full_dataset(
        output_dir=output_dir,
        cell_type=cell_type,
        # label_classifier_artifact=classifier_paths["partial_classifier_path"],
        # full_label_classifier_artifact=classifier_paths["full_classifier_path"],
        continuous_output_col="intrinsic_behavioral_cluster",
        full_output_col="full_behavioral_cluster",
        combine_binary_with_continuous=True,
        verbose=verbose,
    )

    if "sample_name" not in adata_full.obs.columns:
        raise ValueError(
            "test_pipeline requires 'sample_name' in adata_full.obs to visualize backprojection."
        )
    sample_names = (
        pd.Series(adata_full.obs["sample_name"])
        .dropna()
        .astype("string")
        .str.strip()
    )
    sample_names = sample_names[sample_names != ""]
    unique_samples = sorted(sample_names.unique().tolist())
    if len(unique_samples) == 0:
        raise ValueError(
            "test_pipeline could not resolve a sample_name from adata_full.obs for backprojection visualization."
        )
    selected_sample_name = "ROCHE_JM1_Exp042-8_Img03_10T_HER2II"
    _vinfo(verbose, "state-apply", f"Opening backprojection viewer for sample '{selected_sample_name}'")
    _viewer = show_behavioral_state_backprojection(
        sample_name=selected_sample_name,
        output_dir=output_dir,
        cell_type=cell_type,
        run=True,
        verbose=verbose,
    )
 
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
    
