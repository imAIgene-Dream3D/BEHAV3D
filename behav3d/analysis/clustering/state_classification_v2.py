import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path

import scanpy as sc

from behav3d.core.anndata import df_to_adata
from behav3d.analysis.clustering.state_classification import (
    _prepare_state_classification_dataset,
    cluster_group,
    subsample_with_temporal_spacing,
)


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


def plot_binary_group_behavioral_cluster_grid(
    adata,
    binary_group_col="binary_group",
    behavioral_cluster_col="behavioral_clusterid",
    ncols=1,
    figsize_per_plot=(10.0, 1.8),
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
    fig.legend(handles=handles, labels=[str(c) for c in cluster_names], title="behavioral_clusterid", loc="upper right")
    fig.suptitle("Behavioral Cluster Composition per Binary Group", fontsize=12, y=0.995)
    fig.tight_layout(rect=[0, 0, 0.9, 0.97])

    # Diagnostic print to quickly verify cluster presence in plotting table.
    if "5" in ctab.columns:
        by_group = (ctab["5"] * 100.0).round(3)
        print(f"Cluster 5 total n={int(cluster_counts.get('5', 0))}; per-group proportions (%):")
        print(by_group[by_group > 0].sort_values(ascending=False))

    return fig


def run_state_classification_v2(
    df_positions,
    features,
    binary_features_to_group,
    window_size=5,
    max_samples=None,
    min_spacing=None,
    n_neighbors=15,
    resolution="auto",
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
):
    """
    V2 state classification:
    - cluster all continuous descriptive data jointly (single global model)
    - annotate binary contact info and combined ClusterID afterwards
    """
    prepared = _prepare_state_classification_dataset(
        df_positions=df_positions,
        features=features,
        binary_features_to_group=binary_features_to_group,
        window_size=window_size,
        descriptive_features=list(descriptive_features),
        lower_quantile_cap=lower_quantile_cap,
        upper_quantile_cap=upper_quantile_cap,
        outfolder=outfolder,
        scale_features=False,
        incomplete_window_policy=incomplete_window_policy,
    )
    df_analysis = prepared["df_analysis"]
    kept_features = prepared["kept_features"]
    non_feature_cols = prepared["non_feature_cols"]
    binary_cols_to_merge = prepared["binary_cols_to_merge"]

    scaler = None
    df_scaled = df_analysis.copy()
    if len(kept_features) > 0:
        scaler = StandardScaler().fit(df_scaled[kept_features])
        df_scaled[kept_features] = scaler.transform(df_scaled[kept_features])


    if len(df_scaled) < 50:
        raise ValueError("Insufficient rows in full descriptive dataset for clustering.")

    if min_spacing is None:
        target_samples = max(500, len(df_scaled) // 20)
        spacing_to_use = max(1, len(df_scaled) // target_samples // 10)
    else:
        spacing_to_use = int(min_spacing)

    df_train = subsample_with_temporal_spacing(
        df_scaled,
        id_cols=["sample_name", "TrackID"],
        time_col="position_t",
        min_spacing=spacing_to_use,
        max_samples=max_samples,
        random_state=random_state,
    )
    
    print(f"Subsampled {len(df_train)} rows for clustering (spacing={spacing_to_use}).")
    if len(df_train) < 50:
        raise ValueError("Insufficient rows after global subsampling for clustering.")

    model_adata = cluster_group(
        df_group=df_train,
        feature_cols=kept_features,
        non_feature_cols=non_feature_cols,
        n_neighbors=n_neighbors,
        resolution=resolution,
        pca_var_selection=pca_var_selection,
        outfolder=Path(outfolder) if outfolder else None,
        group_name="all_data",
        clustering_method=clustering_method,
        resolutions=resolutions,
        leiden_seed_tries=leiden_seed_tries,
        leiden_subsample_tries=leiden_subsample_tries,
        random_state=random_state,
    )
    if model_adata is None:
        raise RuntimeError("Global clustering model could not be fit.")

    if "preprocessing" not in model_adata.uns:
        model_adata.uns["preprocessing"] = {}
    model_adata.uns["preprocessing"]["kept_features"] = list(kept_features)
    if scaler is not None:
        model_adata.uns["preprocessing"]["scaler"] = {
            "mean": scaler.mean_.astype(float),
            "scale": scaler.scale_.astype(float),
        }

    # Inference to full dataset
    # df_full = df_analysis.copy()
    # if scaler is not None:
    #     x = df_full[kept_features].to_numpy(dtype=float)
    #     df_full.loc[:, kept_features] = (x - scaler.mean_) / scaler.scale_

    obs_cols = [c for c in non_feature_cols if c in df_scaled.columns] + [c for c in binary_cols_to_merge if c in df_scaled.columns]
    adata_full = df_to_adata(df_scaled, kept_features, obs_cols=obs_cols)
    sc.tl.ingest(adata_full, model_adata, obs="ClusterID", embedding_method="umap")

    adata_full.obs = _add_clean_binary_annotation_columns(adata_full.obs, binary_cols_to_merge)
    adata_full.obs["behavioral_clusterid"] = adata_full.obs["ClusterID"].astype(str)
    adata_full.obs["binary_group"] = _assign_binary_group_labels(adata_full.obs, binary_cols_to_merge).astype("category")
    adata_full.obs["ClusterID"] = (
        adata_full.obs["binary_group"].astype(str) + "_" + adata_full.obs["behavioral_clusterid"].astype(str)
    ).astype("category")
    adata_full.obs["behavioral_clusterid"] = adata_full.obs["behavioral_clusterid"].astype("category")

    adata_full.uns["state_classification_v2"] = {
        "window_size": int(window_size),
        "features": list(features),
        "descriptive_features": list(descriptive_features),
        "binary_features_to_group": list(binary_cols_to_merge),
        "clustering_method": clustering_method,
        "resolution": resolution,
        "n_neighbors": int(n_neighbors),
    }

    if outfolder is not None:
        outfolder = Path(outfolder)
        outfolder.mkdir(parents=True, exist_ok=True)
        adata_full.write(outfolder / "adata_state_classification_v2_full.h5ad", compression="gzip")
        model_adata.write(outfolder / "adata_state_classification_v2_model.h5ad", compression="gzip")

        fig = plot_binary_group_behavioral_cluster_grid(adata_full)
        fig.savefig(outfolder / "state_classification_v2_binary_group_cluster_proportions.pdf", dpi=300, bbox_inches="tight")
        plt.close(fig)
        proportions_csv = pd.crosstab(
            adata_full.obs["binary_group"].astype("string").fillna("unassigned"),
            adata_full.obs["behavioral_clusterid"].astype("string").fillna("unassigned"),
            normalize="index",
            dropna=False,
        )
        proportions_csv.to_csv(outfolder / "state_classification_v2_binary_group_cluster_proportions.csv")

    return adata_full
