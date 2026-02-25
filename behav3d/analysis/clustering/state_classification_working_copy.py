import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas.api.types import is_numeric_dtype
import json
from matplotlib.backends.backend_pdf import PdfPages

from pathlib import Path
import pickle

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

A4_PORTRAIT = (8.27, 11.69)
A4_LANDSCAPE = (11.69, 8.27)


def _binary_col_to_group_name(col: str) -> str:
    name = str(col)
    name = name.replace("_pixels", "")
    name = name.replace("_value", "")
    return name


def _assign_binary_group_labels(df: pd.DataFrame, binary_cols: list[str]) -> pd.Series:
    """Build a single categorical group label from binary indicator columns."""
    if len(binary_cols) == 0:
        return pd.Series(["no_contact"] * len(df), index=df.index, dtype="object")

    labels = []
    for _, row in df[binary_cols].iterrows():
        active = []
        for col in binary_cols:
            val = row[col]
            if pd.notna(val) and float(val) == 1.0:
                active.append(_binary_col_to_group_name(col))
        if len(active) == 0:
            labels.append("no_contact")
        elif len(active) == 1:
            labels.append(active[0])
        else:
            labels.append("_and_".join(active))
    return pd.Series(labels, index=df.index, dtype="object")


def _add_clean_binary_annotation_columns(df: pd.DataFrame, binary_cols: list[str]) -> pd.DataFrame:
    """Add user-facing binary annotation columns (e.g. organoid_contact) to obs."""
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
):
    """Rebuild binary_group/full_behavioral_cluster columns from intrinsic cluster labels."""
    if intrinsic_col not in adata.obs.columns:
        raise ValueError(f"Missing '{intrinsic_col}' in adata.obs.")

    binary_cols = [str(c) for c in list(binary_cols_to_merge or [])]
    adata.obs = _add_clean_binary_annotation_columns(adata.obs, binary_cols)
    adata.obs["behavioral_clusterid"] = adata.obs[intrinsic_col].astype(str)
    adata.obs["binary_group"] = _assign_binary_group_labels(
        adata.obs, binary_cols
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

    return _rebuild_full_behavioral_cluster_from_intrinsic(
        adata=adata,
        binary_cols_to_merge=binary_cols_to_merge,
        intrinsic_col="intrinsic_behavioral_cluster",
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
    df,
    id_cols=None,
    time_col="position_t",
    min_spacing=5,
    max_samples=None,
    random_state=123,
):
    """Subsample rows with temporal spacing per track."""
    if id_cols is None:
        id_cols = ["sample_name", "TrackID"]

    np.random.seed(random_state)
    subsampled_rows = []

    for _, track_df in df.groupby(id_cols):
        track_df = track_df.sort_values(time_col).reset_index(drop=True)
        n = len(track_df)
        if n == 0:
            continue
        start_idx = np.random.randint(0, min(min_spacing, n))
        selected_indices = []
        idx = start_idx
        while idx < n:
            selected_indices.append(idx)
            idx += min_spacing
        subsampled_rows.append(track_df.iloc[selected_indices])

    if len(subsampled_rows) == 0:
        return df.iloc[:0].copy()

    df_subsampled = pd.concat(subsampled_rows, ignore_index=True)
    if max_samples is not None and len(df_subsampled) > max_samples:
        df_subsampled = df_subsampled.sample(n=max_samples, random_state=random_state).reset_index(drop=True)
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
    incomplete_window_policy="drop",
    reuse_prepared_dataset=True,
    save_prepared_dataset=True,
    prepared_dataset_cache_path=None,
    verbose=True,
):
    """Create the windowed descriptive dataset used for clustering/classification.

    Returns:
        tuple[pd.DataFrame, dict]: (df_prepared, info)
        where df_prepared is scaled if scale_features=True.
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

    if prepared_dataset_cache_path is None and outfolder is not None:
        prepared_dataset_cache_path = Path(outfolder) / "state_classification_prepared_windows.csv"
    elif prepared_dataset_cache_path is not None:
        prepared_dataset_cache_path = Path(prepared_dataset_cache_path)

    quantile_limits_cache_path = None
    prepared_info_cache_path = None
    if prepared_dataset_cache_path is not None:
        quantile_limits_cache_path = prepared_dataset_cache_path.with_name(
            prepared_dataset_cache_path.stem + "_quantile_cap_limits.json"
        )
        prepared_info_cache_path = prepared_dataset_cache_path.with_name(
            prepared_dataset_cache_path.stem + "_info.json"
        )

    quantile_cap_limits = {}
    loaded_info = None
    loaded_from_cache = False
    if reuse_prepared_dataset and prepared_dataset_cache_path is not None and prepared_dataset_cache_path.exists():
        try:
            df_prepared = pd.read_csv(prepared_dataset_cache_path)
            if quantile_limits_cache_path is not None and quantile_limits_cache_path.exists():
                with open(quantile_limits_cache_path, "r") as f:
                    quantile_cap_limits = json.load(f)
            if prepared_info_cache_path is not None and prepared_info_cache_path.exists():
                with open(prepared_info_cache_path, "r") as f:
                    loaded_info = json.load(f)
            loaded_from_cache = True
            if verbose:
                print(f"Loaded prepared window dataset from cache: {prepared_dataset_cache_path}")
        except Exception as exc:
            if verbose:
                print(f"Could not load prepared dataset cache ({exc}); recomputing windows.")

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

        if save_prepared_dataset and prepared_dataset_cache_path is not None:
            prepared_dataset_cache_path.parent.mkdir(parents=True, exist_ok=True)
            df_prepared.to_csv(prepared_dataset_cache_path, index=False)
            if quantile_limits_cache_path is not None and len(quantile_cap_limits) > 0:
                with open(quantile_limits_cache_path, "w") as f:
                    json.dump(quantile_cap_limits, f, indent=2)
            if prepared_info_cache_path is not None:
                info_to_cache = {
                    "kept_features": list(kept_features),
                    "non_feature_cols": list(non_feature_cols),
                    "binary_cols_to_merge": list(binary_cols_to_merge),
                    "quantile_cap_limits": dict(quantile_cap_limits),
                    "lower_quantile_cap": lower_quantile_cap,
                    "upper_quantile_cap": upper_quantile_cap,
                    "quantile_limits_cache_path": (
                        None if quantile_limits_cache_path is None else str(quantile_limits_cache_path)
                    ),
                }
                with open(prepared_info_cache_path, "w") as f:
                    json.dump(info_to_cache, f, indent=2)
    else:
        binary_cols_to_merge = [col for col in binary_features_to_group if col in df_prepared.columns]
        descriptive_feature_cols = [
            col for col in df_prepared.columns
            if (col not in non_feature_cols) and (not col.endswith("_signal_type"))
        ]
        binary_prefixes = tuple(f"{c}_" for c in binary_cols_to_merge)
        kept_features = [
            c for c in descriptive_feature_cols
            if (c not in binary_cols_to_merge) and (not c.startswith(binary_prefixes))
        ]

        if isinstance(loaded_info, dict):
            cached_non_feature_cols = loaded_info.get("non_feature_cols", None)
            if isinstance(cached_non_feature_cols, list) and len(cached_non_feature_cols) > 0:
                non_feature_cols = [str(c) for c in cached_non_feature_cols]

            cached_binary_cols = loaded_info.get("binary_cols_to_merge", None)
            if isinstance(cached_binary_cols, list) and len(cached_binary_cols) > 0:
                valid_binary_cols = [str(c) for c in cached_binary_cols if str(c) in df_prepared.columns]
                if len(valid_binary_cols) > 0:
                    binary_cols_to_merge = valid_binary_cols

            cached_kept_features = loaded_info.get("kept_features", None)
            if isinstance(cached_kept_features, list) and len(cached_kept_features) > 0:
                valid_kept_features = [str(c) for c in cached_kept_features if str(c) in df_prepared.columns]
                if len(valid_kept_features) > 0:
                    kept_features = valid_kept_features

            cached_qcap = loaded_info.get("quantile_cap_limits", None)
            if isinstance(cached_qcap, dict) and len(quantile_cap_limits) == 0:
                quantile_cap_limits = dict(cached_qcap)

            if "lower_quantile_cap" in loaded_info:
                lower_quantile_cap = loaded_info.get("lower_quantile_cap", None)
            if "upper_quantile_cap" in loaded_info:
                upper_quantile_cap = loaded_info.get("upper_quantile_cap", None)

    scaler = None
    if scale_features and len(kept_features) > 0:
        scaler = StandardScaler().fit(df_prepared[kept_features])
        df_prepared[kept_features] = scaler.transform(df_prepared[kept_features])

    preprocessing_artifact = build_state_preprocessing_artifact(
        kept_features=kept_features,
        scaler=scaler,
        lower_quantile_cap=lower_quantile_cap,
        upper_quantile_cap=upper_quantile_cap,
        quantile_cap_limits=quantile_cap_limits,
    )
    windowing_artifact = {
        "features": list(features),
        "binary_features_to_group": list(binary_features_to_group),
        "window_size": int(window_size),
        "min_spacing": None if min_spacing is None else int(min_spacing),
        "descriptive_features": list(descriptive_features),
        "incomplete_window_policy": str(incomplete_window_policy),
        "lower_quantile_cap": lower_quantile_cap,
        "upper_quantile_cap": upper_quantile_cap,
    }

    info = {
        "kept_features": kept_features,
        "non_feature_cols": non_feature_cols,
        "binary_cols_to_merge": binary_cols_to_merge,
        "scaler": scaler,
        "preprocessing_artifact": preprocessing_artifact,
        "windowing_artifact": windowing_artifact,
        "quantile_cap_limits": quantile_cap_limits,
        "lower_quantile_cap": lower_quantile_cap,
        "upper_quantile_cap": upper_quantile_cap,
        "quantile_limits_cache_path": None if quantile_limits_cache_path is None else str(quantile_limits_cache_path),
        "prepared_info_cache_path": None if prepared_info_cache_path is None else str(prepared_info_cache_path),
    }
    return df_prepared, info


def _save_pdf_page_a4(pdf, fig, orientation="portrait"):
    fig.set_size_inches(*(A4_PORTRAIT if orientation == "portrait" else A4_LANDSCAPE), forward=True)
    pdf.savefig(fig, bbox_inches="tight")


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

    ad = adata.copy()
    ad = _ensure_umap(ad, n_neighbors=30, random_state=123)

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
        if use_dendrogram:
            try:
                sc.tl.dendrogram(ad, groupby=cluster_col)
            except Exception:
                use_dendrogram = False

        sc.pl.heatmap(
            ad,
            var_names=list(feature_cols),
            groupby=cluster_col,
            standard_scale="var",
            figsize=A4_LANDSCAPE,
            swap_axes=True,
            dendrogram=use_dendrogram,
            show_gene_labels=True,
            show=False,
        )
  
        fig_hm = plt.gcf()
        fig_hm.suptitle(f"{title} | Heatmap", y=0.995)
        _save_pdf_page_a4(pdf, fig_hm, orientation="landscape")
        plt.close(fig_hm)


def _ensure_umap(adata, n_neighbors=30, random_state=123):
    """Ensure AnnData has a UMAP embedding for plotting."""
    if "X_umap" in adata.obsm:
        return adata
    if "X_pca" in adata.obsm:
        sc.pp.neighbors(adata, n_neighbors=int(n_neighbors), use_rep="X_pca", random_state=int(random_state))
    else:
        sc.pp.neighbors(adata, n_neighbors=int(n_neighbors), random_state=int(random_state))
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

    if verbose:
        print(
            "ClusterID consistency check "
            f"(mode={summary['match_mode']}): "
            f"checked={summary['n_model_rows_checked']}, "
            f"matches={summary['n_matching_clusterid']}, "
            f"mismatches={summary['n_mismatched_clusterid']}, "
            f"missing_in_full={summary['n_missing_rows_in_adata_full']}"
        )
        if n_mismatch > 0:
            cols = [c for c in merged.columns if c not in {"_match"}]
            print("First mismatches:")
            print(merged.loc[(~merged["_match"]) & merged["_cluster_full"].notna(), cols].head(max_examples))

    return summary, merged


def _to_numpy_2d(X):
    """Convert dense/sparse matrix-like data to a 2D numpy array."""
    if hasattr(X, "toarray"):
        return X.toarray()
    return np.asarray(X)


def build_state_preprocessing_artifact(
    *,
    kept_features,
    scaler,
    lower_quantile_cap=None,
    upper_quantile_cap=None,
    quantile_cap_limits=None,
):
    """Build serializable preprocessing metadata for reproducible inference."""
    artifact = {
        "continuous_feature_cols": list(kept_features),
        "zscore": None,
        "quantile_capping": {
            "lower_quantile": lower_quantile_cap,
            "upper_quantile": upper_quantile_cap,
            "feature_limits": dict(quantile_cap_limits or {}),
        },
    }
    if scaler is not None:
        artifact["zscore"] = {
            "mean": np.asarray(scaler.mean_, dtype=float),
            "scale": np.asarray(scaler.scale_, dtype=float),
        }
    return artifact


def _build_preprocessing_artifact_from_uns(preprocessing_meta, feature_cols=None):
    """Reconstruct preprocessing artifact from model_adata.uns['preprocessing'] metadata."""
    if not isinstance(preprocessing_meta, dict):
        return None

    kept_features = preprocessing_meta.get("kept_features", feature_cols)
    if kept_features is None:
        kept_features = []
    kept_features = [str(c) for c in kept_features]

    artifact = {
        "continuous_feature_cols": list(kept_features),
        "zscore": None,
        "quantile_capping": dict(preprocessing_meta.get("quantile_capping", {})),
    }

    scaler_meta = preprocessing_meta.get("scaler", None)
    if isinstance(scaler_meta, dict) and ("mean" in scaler_meta) and ("scale" in scaler_meta):
        artifact["zscore"] = {
            "mean": np.asarray(scaler_meta["mean"], dtype=float),
            "scale": np.asarray(scaler_meta["scale"], dtype=float),
        }

    return artifact


def _apply_preprocessing_to_continuous_matrix(X, preprocessing_artifact, feature_cols):
    """Apply saved quantile caps and z-normalization to a continuous feature frame."""
    out = np.asarray(X, dtype=float).copy()
    feature_cols = list(feature_cols)
    qmeta = (preprocessing_artifact or {}).get("quantile_capping", {})
    limits = qmeta.get("feature_limits", {}) if isinstance(qmeta, dict) else {}
    if isinstance(limits, dict) and len(limits) > 0:
        for i, feat in enumerate(feature_cols):
            lim = limits.get(feat, None)
            if lim is None:
                continue
            lower = lim.get("lower", None)
            upper = lim.get("upper", None)
            out[:, i] = np.clip(out[:, i], a_min=lower, a_max=upper)

    zmeta = (preprocessing_artifact or {}).get("zscore", None)
    if zmeta is not None:
        mean = np.asarray(zmeta["mean"], dtype=float)
        scale = np.asarray(zmeta["scale"], dtype=float)
        if out.shape[1] != len(mean) or out.shape[1] != len(scale):
            raise ValueError(
                "Preprocessing zscore parameter size mismatch: "
                f"n_features={out.shape[1]}, len(mean)={len(mean)}, len(scale)={len(scale)}"
            )
        out = (out - mean) / scale
    return out


def load_state_classifier_artifact(path):
    """Load a saved state classification v2 classifier artifact pickle."""
    with open(path, "rb") as f:
        artifact = pickle.load(f)
    if not isinstance(artifact, dict) or "classifier" not in artifact:
        raise ValueError("Invalid classifier artifact: expected dict with key 'classifier'.")
    return artifact


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

    df_prepared, prepared_info = prepare_state_classification_dataset(
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
        incomplete_window_policy=str(windowing.get("incomplete_window_policy", "drop")),
        reuse_prepared_dataset=False,
        save_prepared_dataset=False,
        verbose=verbose,
    )

    non_feature_cols = prepared_info["non_feature_cols"]
    binary_cols_to_merge = prepared_info["binary_cols_to_merge"]

    if "continuous_feature_cols" in artifact:
        # Full-label classifier artifact path.
        cont_cols = list(artifact["continuous_feature_cols"])
    else:
        # Label-transfer classifier artifact path.
        cont_cols = list(artifact.get("feature_cols", []))
    if len(cont_cols) == 0:
        raise ValueError("Classifier artifact has no stored feature columns ('continuous_feature_cols'/'feature_cols').")

    missing_cont = [c for c in cont_cols if c not in df_prepared.columns]
    if len(missing_cont) > 0:
        raise ValueError(
            "Prepared dataset is missing required classifier features: "
            f"{missing_cont[:20]}"
        )

    obs_cols = [c for c in non_feature_cols if c in df_prepared.columns] + [
        c for c in binary_cols_to_merge if c in df_prepared.columns
    ]
    adata_query = df_to_adata(df_prepared, cont_cols, obs_cols=obs_cols)

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
    if verbose:
        print(
            "Trained cluster classifier "
            f"(backend={fit_info['backend']}, rows={fit_info['n_train_rows']}, "
            f"features={fit_info['n_features']}, train_acc={fit_info['train_accuracy']:.4f})"
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
    preprocessing_artifact = None
    if isinstance(classifier, dict) and "classifier" in classifier:
        clf = classifier["classifier"]
        preprocessing_artifact = classifier.get("preprocessing", None)

    cols = list(target.var_names) if feature_cols is None else list(feature_cols)

    missing = [c for c in cols if c not in target.var_names]
    if len(missing) > 0:
        raise ValueError(f"adata is missing classifier feature columns: {missing[:10]}")

    X_query = _to_numpy_2d(target[:, cols].X)
    if apply_preprocessing and preprocessing_artifact is not None:
        X_query = _apply_preprocessing_to_continuous_matrix(
            X_query,
            preprocessing_artifact=preprocessing_artifact,
            feature_cols=cols,
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
    filename_prefix="state_classification_classifier_qc",
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

    if verbose:
        print(
            "Classifier QC on model data: "
            f"n={summary['n_rows']}, classes={summary['n_classes']}, "
            f"accuracy={summary['accuracy']:.4f}, balanced_accuracy={summary['balanced_accuracy']:.4f}"
        )

    return summary


def _resolve_binary_classifier_feature_cols(adata, binary_cols_to_merge):
    """Resolve binary classifier features strictly from supplied binary columns."""
    cols = []
    for col in list(binary_cols_to_merge):
        clean_col = _binary_col_to_group_name(col)
        if clean_col in adata.obs.columns:
            cols.append(clean_col)
        elif col in adata.obs.columns:
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
):
    """Build a numeric feature matrix from adata.X (continuous) and obs (binary)."""
    cont_cols = list(continuous_feature_cols)
    bin_cols = list(binary_feature_cols)

    missing_cont = [c for c in cont_cols if c not in adata.var_names]
    if len(missing_cont) > 0:
        raise ValueError(f"Missing continuous feature columns in adata.var_names: {missing_cont[:10]}")

    X_cont = _to_numpy_2d(adata[:, cont_cols].X).astype(float, copy=False)
    if len(bin_cols) == 0:
        if return_binary_feature_names:
            return X_cont, []
        return X_cont

    bin_df = _build_binary_feature_dataframe(
        obs_df=adata.obs,
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
    if verbose:
        print(
            "Trained full ClusterID classifier "
            f"(rows={fit_info['n_train_rows']}, total_features={fit_info['n_features_total']}, "
            f"train_acc={fit_info['train_accuracy']:.4f})"
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
    preprocessing_artifact = classifier_artifact.get("preprocessing", None)

    # Ensure continuous features use identical capping/z-normalization as training.
    X_cont = _to_numpy_2d(target[:, cont_cols].X).astype(float, copy=False)
    if apply_preprocessing and preprocessing_artifact is not None:
        X_cont = _apply_preprocessing_to_continuous_matrix(
            X_cont,
            preprocessing_artifact=preprocessing_artifact,
            feature_cols=cont_cols,
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
    qc_filename_prefix="state_classification_full_classifier_qc_random_forest",
    windowing_artifact=None,
    verbose=True,
):
    """
    Train/apply/evaluate the full classifier (continuous + binary) on a labeled AnnData.

    This lets you run full-classifier training separately from the continuous-only stage.
    """
    if preprocessing_artifact is None:
        preprocessing_artifact = _build_preprocessing_artifact_from_uns(
            adata.uns.get("preprocessing", {}) if isinstance(adata.uns, dict) else None,
            feature_cols=list(adata.var_names),
        )
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
        preprocessing_artifact=preprocessing_artifact,
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
    if windowing_artifact is not None:
        full_label_classifier_artifact["windowing"] = dict(windowing_artifact)

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
    if verbose and full_cm_pdf is not None:
        print(f"Saved full-classifier confusion matrix PDF: {full_cm_pdf}")
    if verbose:
        print(
            "Full classifier QC on adata_full labels: "
            f"n={full_classifier_qc_model_data['n_rows']}, "
            f"accuracy={full_classifier_qc_model_data['accuracy']:.4f}, "
            f"balanced_accuracy={full_classifier_qc_model_data['balanced_accuracy']:.4f}"
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

    # Diagnostic print to quickly verify cluster presence in plotting table.
    if verbose and "5" in ctab.columns:
        by_group = (ctab["5"] * 100.0).round(3)
        print(f"Cluster 5 total n={int(cluster_counts.get('5', 0))}; per-group proportions (%):")
        print(by_group[by_group > 0].sort_values(ascending=False))

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
    df_positions,
    features,
    binary_features_to_group,
    window_size=5,
    min_spacing=None,
    max_samples=None,
    n_neighbors=60,
    min_dist=0.1,
    resolution=0.2,
    descriptive_features=("mean", "median", "std", "net_displacement", "straightness", "mean_square_displacement"),
    pca_var_selection=0.95,
    clustering_method="leiden",
    resolutions=(0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5),
    leiden_seed_tries=10,
    leiden_subsample_tries=10,
    lower_quantile_cap=None,
    upper_quantile_cap=0.99,
    incomplete_window_policy="drop",
    outfolder=None,
    random_state=123,
    reuse_prepared_dataset=True,
    save_prepared_dataset=True,
    verbose=True,
):
    """Stage 1: prepare continuous features, fit clustering model, and return model_adata.

    Stage-1 metadata is stored in:
    - model_adata.uns["preprocessing"]
    - model_adata.uns["clustering"]
    """
    outfolder_path = None if outfolder is None else Path(outfolder)

    if verbose:
        print("# Preparing dataset for state classification and clustering (windowed feature extraction, quantile capping, scaling)...")
    
    df_prepared, prepared_info = prepare_state_classification_dataset(
        df_positions=df_positions,
        features=features,
        binary_features_to_group=binary_features_to_group,
        window_size=window_size,
        min_spacing=min_spacing,
        descriptive_features=list(descriptive_features),
        lower_quantile_cap=lower_quantile_cap,
        upper_quantile_cap=upper_quantile_cap,
        outfolder=outfolder,
        scale_features=True,
        incomplete_window_policy=incomplete_window_policy,
        reuse_prepared_dataset=reuse_prepared_dataset,
        save_prepared_dataset=save_prepared_dataset,
        verbose=verbose,
    )
    kept_features = list(prepared_info["kept_features"])
    non_feature_cols = list(prepared_info["non_feature_cols"])
    binary_cols_to_merge = list(prepared_info["binary_cols_to_merge"])
    preprocessing_artifact = prepared_info["preprocessing_artifact"]
    windowing_artifact = prepared_info["windowing_artifact"]
    scaler = prepared_info.get("scaler", None)

    if len(df_prepared) < 50:
        raise ValueError("Insufficient rows in full descriptive dataset for clustering.")
    if verbose:
        print(
            "Prepared dataset: "
            f"rows={len(df_prepared)}, continuous_features={len(kept_features)}, "
            f"binary_group_features={len(binary_cols_to_merge)}"
        )

    if min_spacing is None:
        spacing_to_use = window_size
    else:
        spacing_to_use = int(min_spacing)
    
    if verbose:
        print(
            f"# Subsampling dataset for clustering with temporal spacing (min_spacing={spacing_to_use}, "
            f"# max_samples={max_samples})..."
        )
    df_train = subsample_with_temporal_spacing(
        df_prepared,
        id_cols=["sample_name", "TrackID"],
        time_col="position_t",
        min_spacing=spacing_to_use,
        max_samples=max_samples,
        random_state=random_state,
    )
    if verbose:
        print(f"Subsampled {len(df_train)} rows for clustering (spacing={spacing_to_use}).")
    if len(df_train) < 50:
        raise ValueError("Insufficient rows after global subsampling for clustering.")

    df_clean = df_train.dropna(subset=kept_features).copy()
    if len(df_clean) < 50:
        raise ValueError("Insufficient rows after dropping NaNs for PCA stage.")
    if verbose and len(df_clean) != len(df_train):
        print(f"Dropped NaNs for PCA input: kept {len(df_clean)} / {len(df_train)} rows.")
    obs_cols = [c for c in non_feature_cols if c in df_clean.columns] + [
        c for c in binary_cols_to_merge if c in df_clean.columns
    ]
    
    model_adata = df_to_adata(df_clean, kept_features, obs_cols=obs_cols)
    
    if verbose:
        print(
            f"# Running PCA for dimensionality reduction before clustering (n_features={len(kept_features)}, "
            f"# n_samples={len(df_clean)}, min_var_selection={pca_var_selection})..."
        )
        
    model_adata = run_pca(
        model_adata,
        pca_var_selection=pca_var_selection,
        ncomps=min(len(kept_features), len(df_clean) - 1),
        svd_solver="full",
        random_state=random_state,
    )
    
    if verbose:
        print(
            f"# Running {clustering_method} clustering with n_neighbors={n_neighbors}, min_dist={min_dist}, "
            f"# resolution={resolution}..."
        )
    sc.pp.neighbors(model_adata, n_neighbors=n_neighbors, use_rep="X_pca", random_state=random_state)
    sc.tl.umap(model_adata, min_dist=min_dist, random_state=random_state)

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
            stability_resolutions=resolutions,
            n_stability_repeats=int(leiden_seed_tries),
            n_subsample_repeats=int(leiden_subsample_tries),
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

    _rebuild_full_behavioral_cluster_from_intrinsic(
        adata=model_adata,
        binary_cols_to_merge=binary_cols_to_merge,
        intrinsic_col="intrinsic_behavioral_cluster",
    )
    if verbose:
        n_intrinsic = int(model_adata.obs["intrinsic_behavioral_cluster"].astype(str).nunique())
        n_full = int(model_adata.obs["full_behavioral_cluster"].astype(str).nunique())
        print(
            "Finished stage-1 clustering: "
            f"intrinsic_clusters={n_intrinsic}, full_behavioral_clusters={n_full}, "
            f"model_rows={model_adata.n_obs}"
        )

    if outfolder_path is not None:
        diagnostics_pdf = Path(outfolder_path) / "state_classification_diagnostics.pdf"
        plot_clustering_diagnostics_pdf(
            adata=model_adata,
            cluster_col="intrinsic_behavioral_cluster",
            feature_cols=kept_features,
            pdf_path=diagnostics_pdf,
            title=f"all_data | {clustering_method} (resolution={resolution})",
        )

    if "preprocessing" not in model_adata.uns:
        model_adata.uns["preprocessing"] = {}
    model_adata.uns["preprocessing"]["kept_features"] = list(kept_features)
    model_adata.uns["preprocessing"]["quantile_capping"] = preprocessing_artifact.get("quantile_capping", {})
    model_adata.uns["preprocessing"]["windowing"] = dict(windowing_artifact)
    model_adata.uns["preprocessing"]["features"] = list(windowing_artifact.get("features", []))
    model_adata.uns["preprocessing"]["binary_features_to_group"] = list(windowing_artifact.get("binary_features_to_group", []))
    model_adata.uns["preprocessing"]["window_size"] = int(windowing_artifact.get("window_size", 0))
    model_adata.uns["preprocessing"]["descriptive_features"] = list(windowing_artifact.get("descriptive_features", []))
    model_adata.uns["preprocessing"]["incomplete_window_policy"] = str(
        windowing_artifact.get("incomplete_window_policy", "drop")
    )
    model_adata.uns["preprocessing"]["lower_quantile_cap"] = windowing_artifact.get("lower_quantile_cap", None)
    model_adata.uns["preprocessing"]["upper_quantile_cap"] = windowing_artifact.get("upper_quantile_cap", None)
    if scaler is not None:
        model_adata.uns["preprocessing"]["scaler"] = {
            "mean": scaler.mean_.astype(float),
            "scale": scaler.scale_.astype(float),
        }

    model_adata.uns["clustering"] = {
        "clustering_method": clustering_method,
        "resolution": resolution,
        "n_neighbors": int(n_neighbors),
        "random_state": int(random_state),
        "non_feature_cols": list(non_feature_cols),
        "binary_cols_to_merge": list(binary_cols_to_merge),
    }
    return model_adata


def train_state_classifiers(
    model_adata,
    outfolder=None,
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
    verbose=True,
):
    """
    Stage 2: train intrinsic/full classifiers on model_adata with optional holdout validation.

    Mutates model_adata in-place by adding/updating model_adata.uns["classification"].
    """
    if isinstance(model_adata, dict):
        raise ValueError(
            "train_state_classifiers now expects model_adata (AnnData), not a stage-1 artifact dict."
        )
    if not (hasattr(model_adata, "uns") and hasattr(model_adata, "var_names")):
        raise ValueError("model_adata must be an AnnData-like object with .uns and .var_names.")
    if verbose:
        print(
            "# Starting classifier training on model_adata: "
            f"# rows={model_adata.n_obs}, features={len(model_adata.var_names)}"
        )

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

    windowing_artifact = pre_meta.get("windowing", None)
    if not isinstance(windowing_artifact, dict):
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
    preprocessing_artifact = _build_preprocessing_artifact_from_uns(
        pre_meta,
        feature_cols=cont_cols,
    )
    if preprocessing_artifact is not None and not isinstance(preprocessing_artifact, dict):
        raise ValueError("model_adata.uns['preprocessing'] could not be converted to preprocessing artifact.")

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

    label_classifier = None  # canonical full-data intrinsic model
    label_classifier_artifact = None  # canonical full-data intrinsic artifact
    label_classifier_artifact_full_data = None
    label_classifier_artifact_train_split = None
    classifier_fit_info = None  # canonical full-data fit info
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

    outfolder_path = None if outfolder is None else Path(outfolder)
    intrinsic_outfolder = None
    full_outfolder = None
    if outfolder_path is not None:
        outfolder_path.mkdir(parents=False, exist_ok=True)
        intrinsic_outfolder = outfolder_path / "intrinsic_behavioral_classification"
        full_outfolder = outfolder_path / "full_behavioral_classification"
        intrinsic_outfolder.mkdir(parents=True, exist_ok=True)
        full_outfolder.mkdir(parents=True, exist_ok=True)

    if train_continuous_classifier:
        X_intrinsic_all = _to_numpy_2d(model_adata[:, cont_cols].X).astype(float, copy=False)
        y_intrinsic_all = model_adata.obs["intrinsic_behavioral_cluster"].astype(str).to_numpy()

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

            clf_intrinsic_train_split, fit_info_intrinsic_train_split = train_random_forest_classifier(
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
            fit_info_intrinsic_train_split["cluster_col"] = "intrinsic_behavioral_cluster"
            classifier_fit_info_train_split = fit_info_intrinsic_train_split

            y_pred_intrinsic_train = np.asarray(
                clf_intrinsic_train_split.predict(X_intrinsic_train)
            ).astype(str)
            y_pred_intrinsic_test = np.asarray(
                clf_intrinsic_train_split.predict(X_intrinsic_test)
            ).astype(str)
            holdout_intrinsic_eval = _evaluate_holdout_predictions(
                y_true=y_intrinsic_test,
                y_pred=y_pred_intrinsic_test,
                outfolder=intrinsic_outfolder,
                filename_prefix=f"state_classification_intrinsic_classifier_holdout_qc_{classifier_backend_name}",
            )

            intrinsic_validation["n_train"] = int(len(idx_train))
            intrinsic_validation["n_test"] = int(len(idx_test))
            intrinsic_validation["train_metrics"] = _build_compact_metric_summary(
                y_intrinsic_train, y_pred_intrinsic_train
            )
            intrinsic_validation["test_metrics"] = holdout_intrinsic_eval["metrics"]
            intrinsic_validation["oob_score_train_split"] = float(
                fit_info_intrinsic_train_split.get("oob_score", np.nan)
            )
            intrinsic_validation["holdout_artifacts"] = holdout_intrinsic_eval["artifacts"]

            label_classifier_artifact_train_split = {
                "classifier": clf_intrinsic_train_split,
                "backend": classifier_backend_name,
                "feature_cols": cont_cols,
                "cluster_col": "intrinsic_behavioral_cluster",
                "preprocessing": preprocessing_artifact,
                "windowing": windowing_artifact,
            }

        label_classifier, classifier_fit_info = train_clusterid_classifier_from_model_adata(
            model_adata=model_adata,
            cluster_col="intrinsic_behavioral_cluster",
            classifier_backend=classifier_backend_name,
            random_state=stage1_random_state,
            n_estimators=classifier_n_estimators,
            min_samples_leaf=classifier_min_samples_leaf,
            n_jobs=classifier_n_jobs,
            max_depth=classifier_max_depth,
            min_samples_split=classifier_min_samples_split,
            max_features=classifier_max_features,
            class_weight=classifier_class_weight,
            verbose=verbose,
        )
        label_classifier_artifact_full_data = {
            "classifier": label_classifier,
            "backend": classifier_backend_name,
            "feature_cols": cont_cols,
            "cluster_col": "intrinsic_behavioral_cluster",
            "preprocessing": preprocessing_artifact,
            "windowing": windowing_artifact,
        }
        label_classifier_artifact = label_classifier_artifact_full_data
        classifier_qc_model_data = evaluate_clusterid_classifier_on_model_adata(
            model_adata=model_adata,
            classifier=label_classifier,
            cluster_col="intrinsic_behavioral_cluster",
            outfolder=intrinsic_outfolder,
            filename_prefix=f"state_classification_classifier_qc_{classifier_backend_name}",
            verbose=verbose,
        )
        intrinsic_validation["oob_score_full_data"] = float(
            classifier_fit_info.get("oob_score", np.nan)
        )
        y_pred_intrinsic_full_data = np.asarray(label_classifier.predict(X_intrinsic_all)).astype(str)
        if validation_mode == "oob_only":
            intrinsic_validation["n_train"] = int(model_adata.n_obs)
            intrinsic_validation["n_test"] = 0
            intrinsic_validation["train_metrics"] = _build_compact_metric_summary(
                y_intrinsic_all,
                y_pred_intrinsic_full_data,
            )

        if verbose and validation_mode == "holdout":
            tm = intrinsic_validation["train_metrics"] or {}
            vm = intrinsic_validation["test_metrics"] or {}
            print(
                "Intrinsic classifier holdout validation: "
                f"n_train={intrinsic_validation['n_train']}, n_test={intrinsic_validation['n_test']}, "
                f"train_acc={tm.get('accuracy', np.nan):.4f}, "
                f"test_acc={vm.get('accuracy', np.nan):.4f}, "
                f"test_bal_acc={vm.get('balanced_accuracy', np.nan):.4f}, "
                f"test_macro_f1={vm.get('macro_f1', np.nan):.4f}, "
                f"oob_train_split={intrinsic_validation.get('oob_score_train_split', np.nan):.4f}, "
                f"oob_full_data={intrinsic_validation.get('oob_score_full_data', np.nan):.4f}"
            )
        elif verbose:
            tm = intrinsic_validation["train_metrics"] or {}
            print(
                "Intrinsic classifier OOB-only validation: "
                f"n_train={intrinsic_validation['n_train']}, "
                f"train_acc={tm.get('accuracy', np.nan):.4f}, "
                f"train_bal_acc={tm.get('balanced_accuracy', np.nan):.4f}, "
                f"train_macro_f1={tm.get('macro_f1', np.nan):.4f}, "
                f"oob_full_data={intrinsic_validation.get('oob_score_full_data', np.nan):.4f}"
            )

    full_label_classifier_artifact = None  # canonical full-data full-label artifact
    full_label_classifier_artifact_full_data = None
    full_label_classifier_artifact_train_split = None
    full_classifier_fit_info = None  # canonical full-data full fit info
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

            adata_train_split = model_adata[idx_train].copy()
            adata_test_split = model_adata[idx_test].copy()
            X_full_train, expanded_binary_train_cols = _build_classifier_matrix_from_adata(
                adata=adata_train_split,
                continuous_feature_cols=cont_cols,
                binary_feature_cols=binary_cols_to_merge,
                return_binary_feature_names=True,
            )
            X_full_test = _build_classifier_matrix_from_adata(
                adata=adata_test_split,
                continuous_feature_cols=cont_cols,
                binary_feature_cols=binary_cols_to_merge,
                binary_expanded_feature_cols=expanded_binary_train_cols,
                return_binary_feature_names=False,
            )
            y_full_train = adata_train_split.obs["full_behavioral_cluster"].astype(str).to_numpy()
            y_full_test = adata_test_split.obs["full_behavioral_cluster"].astype(str).to_numpy()

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

            y_pred_full_train = np.asarray(clf_full_train_split.predict(X_full_train)).astype(str)
            y_pred_full_test = np.asarray(clf_full_train_split.predict(X_full_test)).astype(str)
            holdout_full_eval = _evaluate_holdout_predictions(
                y_true=y_full_test,
                y_pred=y_pred_full_test,
                outfolder=full_outfolder,
                filename_prefix=f"state_classification_full_classifier_holdout_qc_{classifier_backend_name}",
            )

            full_validation["n_train"] = int(len(idx_train))
            full_validation["n_test"] = int(len(idx_test))
            full_validation["train_metrics"] = _build_compact_metric_summary(y_full_train, y_pred_full_train)
            full_validation["test_metrics"] = holdout_full_eval["metrics"]
            full_validation["oob_score_train_split"] = float(
                fit_info_full_train_split.get("oob_score", np.nan)
            )
            full_validation["holdout_artifacts"] = holdout_full_eval["artifacts"]

            full_label_classifier_artifact_train_split = {
                "classifier": clf_full_train_split,
                "continuous_feature_cols": list(cont_cols),
                "binary_feature_cols": list(binary_cols_to_merge),
                "binary_expanded_feature_cols": list(expanded_binary_train_cols),
                "target_col": "full_behavioral_cluster",
                "preprocessing": preprocessing_artifact,
                "windowing": dict(windowing_artifact),
            }

        full_label_classifier_artifact_full_data, full_classifier_fit_info = train_full_clusterid_classifier_from_adata(
            adata=model_adata,
            target_col="full_behavioral_cluster",
            binary_feature_cols=binary_cols_to_merge,
            preprocessing_artifact=preprocessing_artifact,
            random_state=stage1_random_state,
            n_estimators=classifier_n_estimators,
            min_samples_leaf=classifier_min_samples_leaf,
            n_jobs=classifier_n_jobs,
            max_depth=classifier_max_depth,
            min_samples_split=classifier_min_samples_split,
            max_features=classifier_max_features,
            class_weight=classifier_class_weight,
            verbose=verbose,
        )
        full_label_classifier_artifact_full_data["windowing"] = dict(windowing_artifact)
        full_label_classifier_artifact = full_label_classifier_artifact_full_data

        X_full_all = _build_classifier_matrix_from_adata(
            adata=model_adata,
            continuous_feature_cols=cont_cols,
            binary_feature_cols=binary_cols_to_merge,
            binary_expanded_feature_cols=full_label_classifier_artifact_full_data.get(
                "binary_expanded_feature_cols", []
            ),
            return_binary_feature_names=False,
        )
        y_pred_full_all = np.asarray(
            full_label_classifier_artifact_full_data["classifier"].predict(X_full_all)
        ).astype(str)
        full_classifier_qc_model_data = evaluate_classifier_predictions(
            y_true=y_full_all,
            y_pred=y_pred_full_all,
            outfolder=full_outfolder,
            filename_prefix="state_classification_full_classifier_qc_random_forest",
        )
        full_validation["oob_score_full_data"] = float(
            full_classifier_fit_info.get("oob_score", np.nan)
        )
        if validation_mode == "oob_only":
            full_validation["n_train"] = int(model_adata.n_obs)
            full_validation["n_test"] = 0
            full_validation["train_metrics"] = _build_compact_metric_summary(
                y_full_all,
                y_pred_full_all,
            )
        if verbose:
            print(
                "Full classifier QC on model labels: "
                f"n={full_classifier_qc_model_data['n_rows']}, "
                f"accuracy={full_classifier_qc_model_data['accuracy']:.4f}, "
                f"balanced_accuracy={full_classifier_qc_model_data['balanced_accuracy']:.4f}"
            )
        if verbose and validation_mode == "holdout":
            tm = full_validation["train_metrics"] or {}
            vm = full_validation["test_metrics"] or {}
            print(
                "Full classifier holdout validation: "
                f"n_train={full_validation['n_train']}, n_test={full_validation['n_test']}, "
                f"train_acc={tm.get('accuracy', np.nan):.4f}, "
                f"test_acc={vm.get('accuracy', np.nan):.4f}, "
                f"test_bal_acc={vm.get('balanced_accuracy', np.nan):.4f}, "
                f"test_macro_f1={vm.get('macro_f1', np.nan):.4f}, "
                f"oob_train_split={full_validation.get('oob_score_train_split', np.nan):.4f}, "
                f"oob_full_data={full_validation.get('oob_score_full_data', np.nan):.4f}"
            )
        elif verbose:
            tm = full_validation["train_metrics"] or {}
            print(
                "Full classifier OOB-only validation: "
                f"n_train={full_validation['n_train']}, "
                f"train_acc={tm.get('accuracy', np.nan):.4f}, "
                f"train_bal_acc={tm.get('balanced_accuracy', np.nan):.4f}, "
                f"train_macro_f1={tm.get('macro_f1', np.nan):.4f}, "
                f"oob_full_data={full_validation.get('oob_score_full_data', np.nan):.4f}"
            )

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
        "apply_artifact_path": None,
        "intrinsic_artifacts_dir": None if intrinsic_outfolder is None else str(intrinsic_outfolder),
        "full_artifacts_dir": None if full_outfolder is None else str(full_outfolder),
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

    apply_artifact_path = None

    if outfolder_path is not None:
        if save_label_classifier and label_classifier_artifact_full_data is not None:
            classifier_model_path = intrinsic_outfolder / (
                f"state_classification_label_classifier_{classifier_backend_name}_full_data.pkl"
            )
            with open(classifier_model_path, "wb") as f:
                pickle.dump(label_classifier_artifact_full_data, f)
            model_adata.uns["classification"]["classifier_model_path"] = str(classifier_model_path)
            classifier_model_paths["full_data"] = str(classifier_model_path)
            if label_classifier_artifact_train_split is not None:
                classifier_model_path_train_split = intrinsic_outfolder / (
                    f"state_classification_label_classifier_{classifier_backend_name}_train_split.pkl"
                )
                with open(classifier_model_path_train_split, "wb") as f:
                    pickle.dump(label_classifier_artifact_train_split, f)
                classifier_model_paths["train_split"] = str(classifier_model_path_train_split)
            model_adata.uns["classification"]["classifier_model_paths"] = dict(classifier_model_paths)

        if train_full_classifier and save_full_label_classifier and full_label_classifier_artifact_full_data is not None:
            full_classifier_model_path = full_outfolder / "state_classification_full_label_classifier_random_forest_full_data.pkl"
            with open(full_classifier_model_path, "wb") as f:
                pickle.dump(full_label_classifier_artifact_full_data, f)
            model_adata.uns["classification"]["full_classifier_model_path"] = str(full_classifier_model_path)
            full_classifier_model_paths["full_data"] = str(full_classifier_model_path)
            if full_label_classifier_artifact_train_split is not None:
                full_classifier_model_path_train_split = (
                    full_outfolder / "state_classification_full_label_classifier_random_forest_train_split.pkl"
                )
                with open(full_classifier_model_path_train_split, "wb") as f:
                    pickle.dump(full_label_classifier_artifact_train_split, f)
                full_classifier_model_paths["train_split"] = str(full_classifier_model_path_train_split)
            model_adata.uns["classification"]["full_classifier_model_paths"] = dict(full_classifier_model_paths)

        apply_artifact_path = outfolder_path / "state_classification_apply_artifact.pkl"
        apply_artifact = {
            "schema_version": 1,
            "label_transfer_method": transfer_mode,
            "continuous_label_classifier_artifact": {
                "full_data": label_classifier_artifact_full_data,
                "train_split": label_classifier_artifact_train_split,
            },
            "full_label_classifier_artifact": {
                "full_data": full_label_classifier_artifact_full_data,
                "train_split": full_label_classifier_artifact_train_split,
            },
            "default_continuous_model_variant": "full_data",
            "default_full_model_variant": "full_data",
            "preprocessing": preprocessing_artifact,
            "windowing": dict(windowing_artifact),
            "clustering": {
                "clustering_method": clustering_method,
                "resolution": clustering_resolution,
                "n_neighbors": clustering_n_neighbors,
                "random_state": stage1_random_state,
                "non_feature_cols": list(non_feature_cols),
                "binary_cols_to_merge": list(binary_cols_to_merge),
            },
            "output_columns": {
                "clusterid_output_col": "intrinsic_behavioral_cluster",
                "clusterid_confidence_col": classifier_confidence_col,
                "full_classifier_output_col": "full_behavioral_cluster",
                "full_classifier_confidence_col": full_classifier_confidence_col,
            },
        }
        with open(apply_artifact_path, "wb") as f:
            pickle.dump(apply_artifact, f)
        model_adata.uns["classification"]["apply_artifact_path"] = str(apply_artifact_path)

        qcap = (preprocessing_artifact or {}).get("quantile_capping", {})
        feature_limits = qcap.get("feature_limits", {}) if isinstance(qcap, dict) else {}
        if len(feature_limits) > 0:
            qcap_json_path = outfolder_path / "state_classification_quantile_cap_limits.json"
            with open(qcap_json_path, "w") as f:
                json.dump(feature_limits, f, indent=2)
            qcap_df = pd.DataFrame.from_dict(feature_limits, orient="index")
            qcap_df.index.name = "feature"
            qcap_df.to_csv(outfolder_path / "state_classification_quantile_cap_limits.csv")
            model_adata.uns["classification"]["quantile_cap_limits_json"] = str(qcap_json_path)

        model_adata.write(outfolder_path / "adata_state_classification_model.h5ad", compression="gzip")
    if verbose:
        n_intrinsic = int(model_adata.obs["intrinsic_behavioral_cluster"].astype(str).nunique())
        n_full = int(model_adata.obs["full_behavioral_cluster"].astype(str).nunique()) if "full_behavioral_cluster" in model_adata.obs.columns else 0
        print(
            "Completed stage-2 classifier training: "
            f"intrinsic_clusters={n_intrinsic}, full_behavioral_clusters={n_full}, "
            f"saved_artifact={None if apply_artifact_path is None else str(apply_artifact_path)}"
        )
    return
    # return {
    #     "adata_full": None,
    #     "model_adata": model_adata,
    #     "label_classifier": label_classifier,
    #     "label_classifier_artifact": label_classifier_artifact,
    #     "label_classifier_artifact_full_data": label_classifier_artifact_full_data,
    #     "label_classifier_artifact_train_split": label_classifier_artifact_train_split,
    #     "classifier_fit_info": classifier_fit_info,
    #     "classifier_fit_info_train_split": classifier_fit_info_train_split,
    #     "classifier_qc_model_data": classifier_qc_model_data,
    #     "full_label_classifier_artifact": full_label_classifier_artifact,
    #     "full_label_classifier_artifact_full_data": full_label_classifier_artifact_full_data,
    #     "full_label_classifier_artifact_train_split": full_label_classifier_artifact_train_split,
    #     "full_classifier_fit_info": full_classifier_fit_info,
    #     "full_classifier_fit_info_train_split": full_classifier_fit_info_train_split,
    #     "full_classifier_qc_model_data": full_classifier_qc_model_data,
    #     "intrinsic_validation": intrinsic_validation,
    #     "full_validation": full_validation,
    #     "cluster_check_summary": None,
    #     "classifier_model_path": None if classifier_model_path is None else str(classifier_model_path),
    #     "full_classifier_model_path": None if full_classifier_model_path is None else str(full_classifier_model_path),
    #     "apply_artifact_path": None if apply_artifact_path is None else str(apply_artifact_path),
    # }


def apply_state_classifiers_to_full_dataset(
    df_positions,
    model_adata,
    label_classifier_artifact=None,
    full_label_classifier_artifact=None,
    outfolder=None,
    continuous_output_col="intrinsic_behavioral_cluster",
    continuous_confidence_col="intrinsic_behavioral_cluster_confidence",
    full_output_col="full_behavioral_cluster",
    full_confidence_col="full_behavioral_cluster_confidence",
    continuous_model_variant="full_data",
    full_model_variant="full_data",
    combine_binary_with_continuous=True,
    verbose=True,
):
    """Apply trained stage-2 classifier artifacts to the full dataset."""
    if not (hasattr(model_adata, "uns") and hasattr(model_adata, "var_names")):
        raise ValueError("model_adata must be an AnnData-like object with .uns and .var_names.")

    class_meta = model_adata.uns.get("classification", {}) if isinstance(model_adata.uns, dict) else {}
    clust_meta = model_adata.uns.get("clustering", {}) if isinstance(model_adata.uns, dict) else {}

    def _resolve_artifact_variant(artifact_input, variant_name, artifact_name):
        if artifact_input is None:
            return None
        if isinstance(artifact_input, dict) and "classifier" in artifact_input:
            return artifact_input
        if isinstance(artifact_input, dict):
            available = sorted([str(k) for k, v in artifact_input.items() if v is not None])
            if variant_name in artifact_input and artifact_input[variant_name] is not None:
                candidate = artifact_input[variant_name]
                if isinstance(candidate, dict) and "classifier" in candidate:
                    return candidate
                raise ValueError(
                    f"Invalid {artifact_name} variant '{variant_name}': expected classifier artifact dict."
                )
            raise ValueError(
                f"Requested {artifact_name} variant '{variant_name}' not available. "
                f"Available variants: {available}"
            )
        raise ValueError(
            f"{artifact_name} must be a classifier artifact dict or variant dict "
            f"(keys: 'full_data'/'train_split'), got {type(artifact_input)}."
        )

    def _load_artifact_variant_from_meta(paths_key, fallback_key, variant_name, artifact_name):
        paths = class_meta.get(paths_key, None)
        available = []
        if isinstance(paths, dict):
            available.extend([str(k) for k, v in paths.items() if v is not None])
            path = paths.get(variant_name, None)
            if path is not None:
                return load_state_classifier_artifact(path), sorted(set(available))

        fallback_path = class_meta.get(fallback_key, None)
        if fallback_path is not None:
            available = sorted(set(available + ["full_data"]))
            if variant_name == "full_data":
                return load_state_classifier_artifact(fallback_path), available
        return None, sorted(set(available))

    label_classifier_selected = _resolve_artifact_variant(
        label_classifier_artifact, continuous_model_variant, "continuous_label_classifier_artifact"
    )
    full_label_classifier_selected = _resolve_artifact_variant(
        full_label_classifier_artifact, full_model_variant, "full_label_classifier_artifact"
    )

    if label_classifier_selected is None:
        label_classifier_selected, available_continuous_variants = _load_artifact_variant_from_meta(
            paths_key="classifier_model_paths",
            fallback_key="classifier_model_path",
            variant_name=continuous_model_variant,
            artifact_name="continuous_label_classifier_artifact",
        )
    else:
        available_continuous_variants = [str(continuous_model_variant)]

    if full_label_classifier_selected is None:
        full_label_classifier_selected, available_full_variants = _load_artifact_variant_from_meta(
            paths_key="full_classifier_model_paths",
            fallback_key="full_classifier_model_path",
            variant_name=full_model_variant,
            artifact_name="full_label_classifier_artifact",
        )
    else:
        available_full_variants = [str(full_model_variant)]

    if (label_classifier_selected is None) and (continuous_model_variant != "full_data"):
        raise ValueError(
            f"Requested continuous model variant '{continuous_model_variant}' not available. "
            f"Available variants: {available_continuous_variants}"
        )
    if (full_label_classifier_selected is None) and (full_model_variant != "full_data"):
        raise ValueError(
            f"Requested full model variant '{full_model_variant}' not available. "
            f"Available variants: {available_full_variants}"
        )

    if label_classifier_selected is None and full_label_classifier_selected is None:
        raise ValueError(
            "No classifier artifact provided. Supply label/full artifact or train with train_state_classifiers first."
        )
    if verbose:
        print(
            "Applying classifiers to full dataset: "
            f"continuous_variant={continuous_model_variant if label_classifier_selected is not None else None}, "
            f"full_variant={full_model_variant if full_label_classifier_selected is not None else None}"
        )

    if label_classifier_selected is not None:
        adata_full = apply_classifier(
            df_positions=df_positions,
            classifier_artifact_or_path=label_classifier_selected,
            output_col=continuous_output_col,
            confidence_col=continuous_confidence_col,
            outfolder=outfolder,
            verbose=verbose,
        )
    else:
        adata_full = apply_classifier(
            df_positions=df_positions,
            classifier_artifact_or_path=full_label_classifier_selected,
            output_col=full_output_col,
            confidence_col=full_confidence_col,
            outfolder=outfolder,
            verbose=verbose,
        )

    binary_cols_to_merge = list(clust_meta.get("binary_cols_to_merge", []))
    if combine_binary_with_continuous and (label_classifier_selected is not None):
        if continuous_output_col != "intrinsic_behavioral_cluster":
            adata_full.obs["intrinsic_behavioral_cluster"] = adata_full.obs[continuous_output_col].astype("category")
        _rebuild_full_behavioral_cluster_from_intrinsic(
            adata=adata_full,
            binary_cols_to_merge=binary_cols_to_merge,
            intrinsic_col="intrinsic_behavioral_cluster",
        )

    if full_label_classifier_selected is not None:
        predict_clusterids_with_full_classifier(
            adata=adata_full,
            classifier_artifact=full_label_classifier_selected,
            output_col=full_output_col,
            confidence_col=full_confidence_col,
            inplace=True,
            apply_preprocessing=False,
        )

    adata_full.uns["preprocessing"] = dict(model_adata.uns.get("preprocessing", {}))
    adata_full.uns["clustering"] = dict(clust_meta)
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

    if outfolder is not None:
        outfolder = Path(outfolder)
        outfolder.mkdir(parents=False, exist_ok=True)
        adata_full.write(outfolder / "adata_state_classification_full.h5ad", compression="gzip")

        if ("binary_group" in adata_full.obs.columns) and ("behavioral_clusterid" in adata_full.obs.columns):
            pdf1 = outfolder / "state_classification_binary_group_cluster_proportions.pdf"
            fig, _ = plot_binary_group_behavioral_cluster_grid(
                adata_full,
                pdf_path=pdf1,
                return_csv=True,
                verbose=verbose,
            )
            plt.close(fig)
            pdf2 = outfolder / "state_classification_cluster_binary_group_proportions.pdf"
            fig2, _ = plot_behavioral_cluster_binary_group_grid(
                adata_full,
                pdf_path=pdf2,
                return_csv=True,
            )
            plt.close(fig2)

    if verbose:
        n_rows = int(adata_full.n_obs)
        n_intrinsic = int(
            adata_full.obs["intrinsic_behavioral_cluster"].astype(str).nunique()
        ) if "intrinsic_behavioral_cluster" in adata_full.obs.columns else 0
        n_full = int(
            adata_full.obs["full_behavioral_cluster"].astype(str).nunique()
        ) if "full_behavioral_cluster" in adata_full.obs.columns else 0
        print(
            "Completed apply to full dataset: "
            f"rows={n_rows}, intrinsic_clusters={n_intrinsic}, full_behavioral_clusters={n_full}"
        )

    return adata_full


def test_pipeline():
    model_adata = run_state_clustering(
        df_positions,
        features,
        binary_features_to_group, 
        outfolder=outfolder,
        window_size=window_size,
        max_samples=max_samples,
        min_spacing=min_spacing,
        n_neighbors=n_neighbors,
        resolution=resolution,
        descriptive_features=descriptive_features,
        pca_var_selection=pca_var_selection,
        clustering_method=clustering_method,
        resolutions=resolutions,
        leiden_seed_tries=leiden_seed_tries,
        leiden_subsample_tries=leiden_subsample_tries,
        lower_quantile_cap=lower_quantile_cap,
        upper_quantile_cap=upper_quantile_cap,
        incomplete_window_policy=incomplete_window_policy,
        random_state=random_state,
        reuse_prepared_dataset=reuse_prepared_dataset,
        save_prepared_dataset=save_prepared_dataset,
        verbose=verbose,
    )
    
    model_adata.obs["intrinsic_behavioral_cluster"]
    model_adata.obs["full_behavioral_cluster"]
    
    mapping = {
        1: "dead",
        2: "plastic_scanner",
        3: "round_scanner",
        4: "static",
        5: "static",
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
        'organoid_contact_and_tcell_contact_dead': 'dead',
        'organoid_contact_and_tcell_contact_plastic_scanner': 'organoid_contact',
        'organoid_contact_and_tcell_contact_round_scanner': 'organoid_contact',
        'organoid_contact_and_tcell_contact_static': 'organoid_contact',
        'organoid_contact_dead': 'dead',
        'organoid_contact_plastic_scanner': 'organoid_contact',
        'organoid_contact_round_scanner': 'organoid_contact',
        'organoid_contact_static': 'organoid_contact',
        'tcell_contact_dead': 'dead',
        'tcell_contact_plastic_scanner': 'tcell_contact',
        'tcell_contact_round_scanner': 'tcell_contact',
        'tcell_contact_static': 'tcell_contact'
        }
    
    model_adata = relabel_cluster_ids(
        adata=model_adata,
        mapping=full_mapping,
        cluster_key="full_behavioral_cluster",
    )
    
    train_state_classifiers(
        model_adata=model_adata,
        outfolder=outfolder,
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
        df_positions=df_positions,
        model_adata=stage2_artifacts["model_adata"],
        label_classifier_artifact=stage2_artifacts.get("label_classifier_artifact", None),
        full_label_classifier_artifact=(
            stage2_artifacts.get("full_label_classifier_artifact", None) if train_full_label_classifier else None
        ),
        outfolder=outfolder,
        continuous_output_col="intrinsic_behavioral_cluster",
        continuous_confidence_col=classifier_confidence_col,
        full_output_col="full_behavioral_cluster",
        full_confidence_col=full_classifier_confidence_col,
        combine_binary_with_continuous=True,
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
    outfolder = Path(ssd_dir, r"BHVD_BEHAV3D/BEHAV3D_python/rolling_classification")
    metadata = load_behav3d_metadata(metadata_csv_path)
    analysis_outdir = Path(output_dir, "analysis", "tcell")
    feature_outdir = Path(analysis_outdir, "track_features")
    df_tracks_path = Path(feature_outdir, f"BEHAV3D_tcell_combined_track_features_filtered.csv")
    df_positions = pd.read_csv(df_tracks_path)
    df_positions=df_positions.sort_values(by=["sample_name", "TrackID", "position_t"])
    outfolder = Path("/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/rolling_classification/clustering_then_binary_assignment")
    
    window_size = 5
    max_samples = None
    min_spacing = 10
    n_neighbors = 60
    resolution = 0.15
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
    resolutions = 0.15
    leiden_seed_tries = 10
    leiden_subsample_tries = 10
    lower_quantile_cap = None
    upper_quantile_cap = 0.99
    incomplete_window_policy = "drop"
    random_state = 123
    reuse_prepared_dataset = True
    save_prepared_dataset = True
    label_transfer_method = "classifier"
    
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
    
