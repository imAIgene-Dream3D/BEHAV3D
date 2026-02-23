import time
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import umap
import scanpy as sc

from sklearn.cluster import KMeans, HDBSCAN
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from sklearn.decomposition import PCA
from sklearn.feature_selection import VarianceThreshold

from behav3d.core.metadata import load_behav3d_metadata, check_behav3d_metadata
from behav3d.core.anndata import (
    df_to_adata,
    adata_add_back_to_df,
    merge_pandas_cols_into_obs_anndata,
    relabel_img_from_adata,
)
from behav3d.io.images import load_image
from behav3d.analysis.clustering.general import (
    select_nonbinary_columnnames, 
    relabel_cluster_ids
)

from behav3d.analysis.filtering import subset_timepoints_from_tracks, subset_selection_of_tracks

from behav3d.features.rolling_window_features import create_descriptive_track_dataset, infer_signal_types
from behav3d.features.state_descriptive_features import drop_highly_correlated_features
from behav3d.analysis.clustering.general.leiden import (
    run_pca, 
    run_leiden_clustering, 
    merge_small_clusters
)
from behav3d.analysis.clustering.state.filtering import filter_short_state_runs

from behav3d.analysis.clustering.state.visualization.plots.clustering import (
    plot_exemplar_track_bars
)   
from behav3d.analysis.clustering.general.visualization.plots import (
    plot_per_cluster_proportions, 
    plot_top_ranking_features,
    plot_number_per_clusters
)

from pathlib import Path

seed = 123
random.seed(seed)
np.random.seed(seed)

A4_PORTRAIT = (8.27, 11.69)
A4_LANDSCAPE = (11.69, 8.27)


def _a4_size(orientation="auto", width=None, height=None):
    """Return A4 size in inches, with optional auto orientation from width/height."""
    if orientation == "portrait":
        return A4_PORTRAIT
    if orientation == "landscape":
        return A4_LANDSCAPE
    if width is not None and height is not None and width > height:
        return A4_LANDSCAPE
    return A4_PORTRAIT


def save_pdf_page_a4(pdf, fig, orientation="auto"):
    """Save a matplotlib figure to PDF with A4 page size."""
    width, height = fig.get_size_inches()
    fig.set_size_inches(*_a4_size(orientation=orientation, width=width, height=height), forward=True)
    pdf.savefig(fig, bbox_inches="tight")


def save_figure_a4(fig, path, orientation="auto", dpi=300):
    """Save a standalone figure to an A4-sized PDF/image."""
    width, height = fig.get_size_inches()
    fig.set_size_inches(*_a4_size(orientation=orientation, width=width, height=height), forward=True)
    fig.savefig(path, bbox_inches="tight", dpi=dpi)


def _minmax_scale_columns(df):
    """Approximate scanpy standard_scale='var' on a cluster x feature matrix."""
    if df is None or df.empty:
        return df
    scaled = df.copy().astype(float)
    col_min = scaled.min(axis=0)
    col_max = scaled.max(axis=0)
    denom = (col_max - col_min).replace(0, 1.0)
    return (scaled - col_min) / denom


def identify_binary_features(df, feature_cols):
    """
    Identify which features are binary (only contain 0 and 1 values).
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe
    feature_cols : list
        List of feature column names to check
        
    Returns
    -------
    binary_features : list
        List of binary feature names
    nonbinary_features : list
        List of non-binary feature names
    """
    binary_features = []
    nonbinary_features = []
    
    for col in feature_cols:
        if col in df.columns:
            unique_vals = df[col].dropna().unique()
            # Check if only contains 0 and 1
            if set(unique_vals).issubset({0, 1, 0.0, 1.0}):
                binary_features.append(col)
            else:
                nonbinary_features.append(col)
    
    return binary_features, nonbinary_features


def create_binary_groups(df, binary_features):
    """
    Create groups based on combinations of binary features.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe
    binary_features : list
        List of binary feature column names
        
    Returns
    -------
    groups : dict
        Dictionary mapping group names to boolean masks
    """
    from itertools import combinations
    
    groups = {}
    
    if len(binary_features) == 0:
        return groups
    
    # Create individual groups for each binary feature (only that feature is 1)
    for feature in binary_features:
        # Create mask where ONLY this feature is 1 and all others are 0
        mask = df[feature] == 1
        for other_feature in binary_features:
            if other_feature != feature:
                mask &= (df[other_feature] == 0)
        
        group_name = f"{feature}_only"
        groups[group_name] = mask
    
    # Create combination groups (2 or more features are 1)
    for r in range(2, len(binary_features) + 1):
        for combo in combinations(binary_features, r):
            # All features in combo must be 1
            mask = pd.Series(True, index=df.index)
            for feature in combo:
                mask &= (df[feature] == 1)
            
            # All features NOT in combo must be 0
            for feature in binary_features:
                if feature not in combo:
                    mask &= (df[feature] == 0)
            
            # Create readable group name
            combo_name = "_and_".join(combo)
            groups[combo_name] = mask
    
    # Create a "none" group (all binary features are 0)
    none_mask = pd.Series(True, index=df.index)
    for feature in binary_features:
        none_mask &= (df[feature] == 0)
    groups["no_contact"] = none_mask
    
    return groups


def subsample_with_temporal_spacing(
    df,
    id_cols=None,
    time_col="position_t",
    min_spacing=5,
    max_samples=None,
    random_state=123
):
    """
    Subsample timepoints with minimum temporal spacing to reduce autocorrelation.
    
    This helps prevent consecutive timepoints from the same track clustering together
    purely due to temporal proximity rather than biological state.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe
    id_cols : list
        Columns identifying unique tracks
    time_col : str
        Time column name
    min_spacing : int
        Minimum number of timepoints between selected samples from the same track
    max_samples : int, optional
        Maximum total number of samples after temporal-spacing subsampling.
        If provided and exceeded, samples are randomly downsampled to this size.
    random_state : int
        Random seed
        
    Returns
    -------
    df_subsampled : pd.DataFrame
        Subsampled dataframe with temporal spacing
    """
    if id_cols is None:
        id_cols = ["sample_name", "TrackID"]
    
    np.random.seed(random_state)
    
    subsampled_rows = []
    
    for _, track_df in df.groupby(id_cols):
        track_df = track_df.sort_values(time_col).reset_index(drop=True)
        n = len(track_df)
        
        if n == 0:
            continue
        
        # Start with a random offset
        start_idx = np.random.randint(0, min(min_spacing, n))
        selected_indices = []
        
        idx = start_idx
        while idx < n:
            selected_indices.append(idx)
            idx += min_spacing
            
        subsampled_rows.append(track_df.iloc[selected_indices])
    
    if len(subsampled_rows) == 0:
        return df.iloc[:0].copy()  # Return empty dataframe with same columns

    df_subsampled = pd.concat(subsampled_rows, ignore_index=True)

    # Optional global cap after full per-track temporal-spacing subsampling
    if max_samples is not None and len(df_subsampled) > max_samples:
        df_subsampled = df_subsampled.sample(n=max_samples, random_state=random_state).reset_index(drop=True)

    return df_subsampled




def validate_clustering(X, labels, min_silhouette=0.3, max_davies_bouldin=1.5):
    """
    Validate whether clustering results show clear evidence of distinct subclusters.
    
    Uses silhouette score and Davies-Bouldin index to assess cluster quality.
    Only returns True if there's strong evidence for multiple distinct clusters.
    
    Parameters
    ----------
    X : np.ndarray
        Feature matrix (n_samples, n_features)
    labels : np.ndarray
        Cluster labels
    min_silhouette : float
        Minimum silhouette score to accept clustering (higher = better separation)
        Range: [-1, 1], where >0.3 indicates reasonable structure
    max_davies_bouldin : float
        Maximum Davies-Bouldin index to accept clustering (lower = better separation)
        Range: [0, inf], where <1.5 indicates good separation
        
    Returns
    -------
    is_valid : bool
        True if clustering shows clear evidence of distinct subclusters
    metrics : dict
        Dictionary containing validation metrics
    """
    from sklearn.metrics import silhouette_score, davies_bouldin_score
    
    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)
    
    # If only one cluster, no subclusters exist
    if n_clusters <= 1:
        return False, {"n_clusters": n_clusters, "silhouette": None, "davies_bouldin": None}
    
    # Calculate silhouette score (higher is better, range [-1, 1])
    # >0.7: strong structure, 0.5-0.7: reasonable, 0.25-0.5: weak, <0.25: no structure
    try:
        silhouette = silhouette_score(X, labels)
    except:
        silhouette = -1.0
    
    # Calculate Davies-Bouldin index (lower is better, range [0, inf])
    # <1.0: excellent, 1.0-1.5: good, >1.5: poor separation
    try:
        davies_bouldin = davies_bouldin_score(X, labels)
    except:
        davies_bouldin = float('inf')
    
    metrics = {
        "n_clusters": n_clusters,
        "silhouette": silhouette,
        "davies_bouldin": davies_bouldin
    }
    
    # Require BOTH good silhouette AND good Davies-Bouldin for validation
    is_valid = (silhouette >= min_silhouette) and (davies_bouldin <= max_davies_bouldin)
    
    return is_valid, metrics


def build_cluster_feature_correlation(feature_df, labels):
    """
    Build a cluster-vs-cluster correlation matrix from cluster mean feature profiles.

    Parameters
    ----------
    feature_df : pd.DataFrame
        Feature matrix used for clustering/evaluation
    labels : array-like
        Cluster labels per row in feature_df

    Returns
    -------
    corr : pd.DataFrame or None
        Correlation matrix between clusters based on mean feature vectors.
        Returns None when fewer than 2 clusters are present.
    """
    tmp = feature_df.copy()
    tmp["_cluster_tmp"] = pd.Categorical(labels)
    centroids = tmp.groupby("_cluster_tmp", observed=False).mean(numeric_only=True)
    if centroids.shape[0] < 2:
        return None
    return centroids.T.corr()


def build_cluster_feature_matrix(feature_df, labels):
    """
    Build a cluster-by-feature mean matrix for heatmap plotting.
    """
    tmp = feature_df.copy()
    numeric_tmp = tmp.select_dtypes(include=[np.number]).copy()
    if numeric_tmp.shape[1] == 0:
        numeric_tmp = tmp.apply(pd.to_numeric, errors="coerce")
        numeric_tmp = numeric_tmp.dropna(axis=1, how="all")
    if numeric_tmp.shape[1] == 0:
        return None
    numeric_tmp["_cluster_tmp"] = pd.Categorical(labels.astype(str))
    centroids = numeric_tmp.groupby("_cluster_tmp", observed=False).mean(numeric_only=True)
    if centroids.shape[0] < 1 or centroids.shape[1] < 1:
        return None
    return centroids


def plot_cluster_feature_correlation(
    corr,
    title,
    pdf=None,
    figsize=(6, 5),
):
    """
    Plot a cluster-vs-cluster feature-profile correlation matrix.
    If pdf is provided, writes directly to that PDF.
    """
    if corr is None:
        return

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        corr,
        vmin=-1,
        vmax=1,
        cmap="coolwarm",
        annot=True,
        fmt=".2f",
        square=True,
        cbar_kws={"label": "Correlation"},
        ax=ax,
    )
    ax.set_title(title)
    plt.tight_layout()
    if pdf is not None:
        save_pdf_page_a4(pdf, fig, orientation="portrait")
        plt.close(fig)
    return fig


def plot_cluster_feature_heatmap(pdf, adata, cluster_col, title):
    """
    Plot Scanpy heatmap (with dendrogram) without grouped feature labels.
    """
    if adata is None or adata.n_obs == 0 or adata.n_vars == 0:
        return

    sc.tl.dendrogram(adata, groupby=cluster_col)
    sc.pl.heatmap(
        adata,
        var_names=list(adata.var_names),
        groupby=cluster_col,
        standard_scale="var",
        figsize=A4_LANDSCAPE,
        swap_axes=True,
        dendrogram=True,
        show_gene_labels=True,
        show=False,
    )
    fig = plt.gcf()
    fig.suptitle(title, y=0.995)
    save_pdf_page_a4(pdf, fig, orientation="landscape")
    plt.close(fig)


def plot_scanpy_heatmap_page(pdf, adata, cluster_col, title):
    """
    Save a native Scanpy heatmap as its own A4 PDF page (vector quality).
    """
    if adata is None or adata.n_obs == 0 or adata.n_vars == 0:
        return

    cluster_values = adata.obs[cluster_col].astype(str)
    n_clusters = cluster_values.nunique()
    use_dendrogram = n_clusters > 1

    if use_dendrogram:
        try:
            sc.tl.dendrogram(adata, groupby=cluster_col)
        except Exception as exc:
            print(f"Warning: dendrogram failed for {cluster_col} ({exc}). Plotting heatmap without dendrogram.")
            use_dendrogram = False

    sc.pl.heatmap(
        adata,
        var_names=list(adata.var_names),
        groupby=cluster_col,
        standard_scale="var",
        figsize=A4_LANDSCAPE,
        swap_axes=True,
        dendrogram=use_dendrogram,
        show_gene_labels=True,
        show=False,
    )
    fig = plt.gcf()
    fig.suptitle(title, y=0.995)
    save_pdf_page_a4(pdf, fig, orientation="landscape")
    plt.close(fig)


def plot_result_summary_page(pdf, adata, cluster_col, corr, feature_matrix, title, selected=False):
    """
    One A4 page per tested setting with UMAP + correlation.
    """
    fig, (ax_umap, ax_corr) = plt.subplots(1, 2, figsize=A4_LANDSCAPE, gridspec_kw={"width_ratios": [1.0, 1.0]})

    # Panel 1: UMAP
    sc.pl.umap(
        adata,
        color=cluster_col,
        ax=ax_umap,
        show=False,
        title="UMAP",
        legend_loc="on data",
        legend_fontsize=7,
    )
    ax_umap.set_aspect("equal", adjustable="box")

    # Panel 2: Cluster correlation
    if corr is not None and not corr.empty:
        sns.heatmap(
            corr,
            vmin=-1,
            vmax=1,
            cmap="coolwarm",
            annot=True,
            fmt=".2f",
            square=True,
            cbar_kws={"label": "Correlation"},
            ax=ax_corr,
        )
        ax_corr.set_title("Cluster Correlation")
    else:
        ax_corr.text(0.5, 0.5, "No correlation matrix", ha="center", va="center", transform=ax_corr.transAxes)
        ax_corr.set_axis_off()

    page_title = title + (" [SELECTED]" if selected else "")
    fig.suptitle(page_title, fontsize=12, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    save_pdf_page_a4(pdf, fig, orientation="landscape")
    plt.close(fig)


def _plot_method_umap(pdf, adata, col, title, selected=False):
    fig, ax = plt.subplots(figsize=A4_PORTRAIT)
    sc.pl.umap(
        adata,
        color=col,
        ax=ax,
        show=False,
        title=title,
        legend_loc="on data",
        legend_fontsize=8,
    )
    ax.set_aspect("equal", adjustable="box")
    if selected:
        ax.set_title(ax.get_title(), color="red", fontweight="bold")
    save_pdf_page_a4(pdf, fig, orientation="portrait")
    plt.close(fig)


def run_kmeans_search(
    adata,
    X_pca,
    X_umap,
    feature_df,
    min_silhouette,
    max_davies_bouldin,
    random_state,
):
    print("    Iterative KMeans search (k=1..8) on X_pca, validated on X_umap...")
    results = []
    for k in range(1, 9):
        if k == 1:
            labels = np.zeros(X_umap.shape[0], dtype=int)
            inertia = np.nan
        else:
            kmeans = KMeans(n_clusters=k, random_state=random_state, n_init=10)
            labels = kmeans.fit_predict(X_umap)
            inertia = float(kmeans.inertia_)

        is_valid, metrics = validate_clustering(
            X_umap,
            labels,
            min_silhouette=min_silhouette,
            max_davies_bouldin=max_davies_bouldin
        )
        col_name = f"kmeans_{k}"
        adata.obs[col_name] = pd.Categorical(labels.astype(str))
        corr = build_cluster_feature_correlation(feature_df, labels.astype(str))
        feature_matrix = build_cluster_feature_matrix(feature_df, labels)
        results.append({
            "param_name": "k",
            "param_value": k,
            "col": col_name,
            "labels": labels,
            "metrics": metrics,
            "valid": is_valid,
            "inertia": inertia,
            "corr": corr,
            "feature_matrix": feature_matrix,
        })

        sil_str = "N/A" if metrics["silhouette"] is None else f"{metrics['silhouette']:.3f}"
        db_str = "N/A" if metrics["davies_bouldin"] is None else f"{metrics['davies_bouldin']:.3f}"
        print(f"      k={k}: n_clusters={metrics['n_clusters']}, silhouette={sil_str}, db={db_str}, valid={is_valid}")

    scored = [r for r in results if r["metrics"]["silhouette"] is not None]
    if len(scored) == 0:
        best_result = results[0]
        print("    No valid multi-cluster solution found, defaulting to first tested setting.")
    else:
        best_result = max(scored, key=lambda r: r["metrics"]["silhouette"])
    return results, best_result


def run_leiden_search(
    adata,
    feature_df,
    resolutions,
    n_stability_repeats,
    n_subsample_repeats,
    random_state,
):
    print("    Leiden auto-selection + diagnostics over tested resolutions...")
    auto_col = "leiden_auto_selection_tmp"
    adata = run_leiden_clustering(
        adata,
        n_neighbors=None,
        resolution="auto",
        stability_resolutions=tuple(float(r) for r in resolutions),
        n_stability_repeats=n_stability_repeats,
        n_subsample_repeats=n_subsample_repeats,
        metric="euclidean",
        method="umap",
        use_rep="X_pca",
        key_added=auto_col,
        random_state=random_state,
    )
    best_res_auto = float(adata.uns.get("leiden_stability_best_res", resolutions[0]))
    labels_by_resolution = adata.uns.get("leiden_stability_labels_by_resolution", {})
    print(f"    Auto-selected Leiden resolution from stability search: {best_res_auto}")

    results = []
    for res in resolutions:
        col_name = f"leiden_{str(res).replace('.', '_')}"
        labels = None
        for key_str, key_labels in labels_by_resolution.items():
            if np.isclose(float(key_str), float(res)):
                labels = np.asarray(key_labels).astype(int)
                break
        if labels is None:
            raise ValueError(
                f"Leiden labels for resolution={res} not found in stability search output. "
                "Expected cached labels in adata.uns['leiden_stability_labels_by_resolution']."
            )

        adata.obs[col_name] = pd.Categorical((labels + 1).astype(str))
        metrics = {
            "n_clusters": int(len(np.unique(labels))),
            "silhouette": None,
            "davies_bouldin": None,
        }
        corr = build_cluster_feature_correlation(feature_df, labels.astype(str))
        feature_matrix = build_cluster_feature_matrix(feature_df, labels)
        results.append({
            "param_name": "resolution",
            "param_value": float(res),
            "col": col_name,
            "labels": labels,
            "metrics": metrics,
            "valid": None,
            "inertia": np.nan,
            "corr": corr,
            "feature_matrix": feature_matrix,
        })
        print(f"      res={res}: n_clusters={metrics['n_clusters']}")

    best_result = None
    for res_data in results:
        if np.isclose(float(res_data["param_value"]), best_res_auto):
            best_result = res_data
            break
    if best_result is None:
        best_result = results[0]
        print("    Warning: auto-selected resolution not found in diagnostics list; used first tested resolution.")

    return results, best_result


def plot_kmeans_diagnostics(pdf, adata, results, best_result, min_silhouette, max_davies_bouldin):
    x = [r["param_value"] for r in results]
    sils = [np.nan if r["metrics"]["silhouette"] is None else r["metrics"]["silhouette"] for r in results]
    dbs = [np.nan if r["metrics"]["davies_bouldin"] is None else r["metrics"]["davies_bouldin"] for r in results]

    fig_metrics, axes = plt.subplots(1, 2, figsize=A4_LANDSCAPE)
    axes[0].plot(x, sils, "o-", color="tab:green")
    axes[0].set_title("Silhouette over search")
    axes[0].set_xlabel(best_result["param_name"])
    axes[0].set_ylabel("Silhouette score")
    axes[0].axvline(best_result["param_value"], color="red", linestyle="--", linewidth=1)
    axes[0].axhline(min_silhouette, color="black", linestyle=":", linewidth=1)

    axes[1].plot(x, dbs, "o-", color="tab:red")
    axes[1].set_title("Davies-Bouldin over search")
    axes[1].set_xlabel(best_result["param_name"])
    axes[1].set_ylabel("Davies-Bouldin index")
    axes[1].axvline(best_result["param_value"], color="red", linestyle="--", linewidth=1)
    axes[1].axhline(max_davies_bouldin, color="black", linestyle=":", linewidth=1)
    plt.tight_layout()
    save_pdf_page_a4(pdf, fig_metrics, orientation="landscape")
    plt.close(fig_metrics)

    for res in results:
        selected = np.isclose(float(res["param_value"]), float(best_result["param_value"]))
        title = (
            f"k={res['param_value']} | "
            f"Sil={res['metrics']['silhouette'] if res['metrics']['silhouette'] is not None else 'N/A'} | "
            f"DB={res['metrics']['davies_bouldin'] if res['metrics']['davies_bouldin'] is not None else 'N/A'}"
        )
        plot_result_summary_page(
            pdf=pdf,
            adata=adata,
            cluster_col=res["col"],
            corr=res.get("corr"),
            feature_matrix=res.get("feature_matrix"),
            title=title,
            selected=selected,
        )
        plot_scanpy_heatmap_page(
            pdf=pdf,
            adata=adata,
            cluster_col=res["col"],
            title=f"Scanpy Heatmap (k={res['param_value']})",
        )


def plot_leiden_diagnostics(pdf, adata, results, best_result):
    stability_summary = adata.uns.get("leiden_stability_summary", None)
    if stability_summary is not None and len(stability_summary) > 0:
        tested_resolutions = sorted({float(r["param_value"]) for r in results})
        fig_ari, ax_ari = plt.subplots(figsize=A4_LANDSCAPE)
        ax_ari.plot(
            stability_summary["resolution"],
            stability_summary["mean_ari"],
            "o-",
            color="tab:blue",
            label="Seed mean ARI",
        )
        if "std_ari" in stability_summary.columns:
            ax_ari.fill_between(
                stability_summary["resolution"],
                stability_summary["mean_ari"] - stability_summary["std_ari"],
                stability_summary["mean_ari"] + stability_summary["std_ari"],
                color="tab:blue",
                alpha=0.2,
            )
        if "mean_ari_subsample" in stability_summary.columns:
            ax_ari.plot(
                stability_summary["resolution"],
                stability_summary["mean_ari_subsample"],
                "o-",
                color="tab:orange",
                label="Subsample mean ARI",
            )
        ax_ari.axvline(best_result["param_value"], color="red", linestyle="--", linewidth=1)
        ax_ari.set_title("Leiden Stability ARI over Resolution")
        ax_ari.set_xlabel("Resolution")
        ax_ari.set_ylabel("ARI")
        # Force x-axis to show the exact tested resolutions (e.g., includes 0.05)
        ax_ari.set_xticks(tested_resolutions)
        ax_ari.set_xlim(min(tested_resolutions) - 0.005, max(tested_resolutions) + 0.005)
        ax_ari.legend()
        plt.tight_layout()
        save_pdf_page_a4(pdf, fig_ari, orientation="landscape")
        plt.close(fig_ari)

    for res in results:
        selected = np.isclose(float(res["param_value"]), float(best_result["param_value"]))
        title = (
            f"resolution={res['param_value']} | "
            f"n_clusters={res['metrics']['n_clusters']}"
        )
        plot_result_summary_page(
            pdf=pdf,
            adata=adata,
            cluster_col=res["col"],
            corr=res.get("corr"),
            feature_matrix=res.get("feature_matrix"),
            title=title,
            selected=selected,
        )
        plot_scanpy_heatmap_page(
            pdf=pdf,
            adata=adata,
            cluster_col=res["col"],
            title=f"Scanpy Heatmap (resolution={res['param_value']})",
        )


def cluster_group(
    df_group,
    feature_cols,
    non_feature_cols,
    min_silhouette=0.3,
    max_davies_bouldin=1.5,
    n_neighbors=30,
    resolution="auto",
    pca_var_selection=0.95,
    outfolder=None,
    group_name=None,
    clustering_method="kmeans",
    resolutions=(0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5),
    leiden_seed_tries=10,
    leiden_subsample_tries=10,
    random_state=123
):
    """
    Perform clustering on a single group of data.
    
    Parameters
    ----------
    df_group : pd.DataFrame
        Subset of data for this group
    feature_cols : list
        Feature columns to use for clustering
    non_feature_cols : list
        Non-feature columns (metadata)
    n_neighbors : int
        Number of neighbors for Leiden clustering
    resolution : float
        Resolution parameter for Leiden clustering
    min_cluster_size : int
        Minimum cluster size for merging small clusters
    pca_var_selection : float
        Variance threshold for PCA
    random_state : int
        Random seed
        
    Returns
    -------
    adata : AnnData
        AnnData object with clustering results
    """
    # Drop NaN values
    df_clean = df_group.dropna(subset=feature_cols).copy()
    
    if len(df_clean) < 50:  # Skip if too few samples
        return None
    
    
    # Create AnnData object
    adata = df_to_adata(df_clean, feature_cols, obs_cols=non_feature_cols)
    adata.uns["preprocessing"] = {
        "kept_features": list(feature_cols),
        # "scaler": {
        #     "mean": scaler.mean_.astype(float),
        #     "scale": scaler.scale_.astype(float),
        # }
    }
    
    # Run PCA
    adata = run_pca(
        adata,
        pca_var_selection=pca_var_selection,
        ncomps=min(len(feature_cols), len(df_clean) - 1),
        svd_solver='full', 
        random_state=random_state
    )

    # Build shared UMAP embedding used for validation and visualization
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_pca", random_state=random_state)
    sc.tl.umap(adata, min_dist=0.1, random_state=random_state)
    X_pca = adata.obsm["X_pca"]
    X_umap = adata.obsm["X_umap"]
    feature_df = df_clean[feature_cols].reset_index(drop=True)

    method = clustering_method.lower()
    if method == "kmeans":
        results, best_result = run_kmeans_search(
            adata=adata,
            X_pca=X_pca,
            X_umap=X_umap,
            feature_df=feature_df,
            min_silhouette=min_silhouette,
            max_davies_bouldin=max_davies_bouldin,
            random_state=random_state,
        )
    elif method == "leiden":
        if resolution == "auto":
            results, best_result = run_leiden_search(
                adata=adata,
                feature_df=feature_df,
                resolutions=resolutions,
                n_stability_repeats=leiden_seed_tries,
                n_subsample_repeats=leiden_subsample_tries,
                random_state=random_state,
            )
        else:
            fixed_res = float(resolution)
            col_name = f"leiden_{str(fixed_res)}"
            adata = run_leiden_clustering(
                adata,
                n_neighbors=None,
                resolution=fixed_res,
                metric="euclidean",
                method="umap",
                use_rep="X_pca",
                key_added=col_name,
                random_state=random_state,
            )
            labels = adata.obs[col_name].astype("category").cat.codes.to_numpy()
            metrics = {
                "n_clusters": int(len(np.unique(labels))),
                "silhouette": None,
                "davies_bouldin": None,
            }
            results = [{
                "param_name": "resolution",
                "param_value": fixed_res,
                "col": col_name,
                "labels": labels,
                "metrics": metrics,
                "valid": None,
                "inertia": np.nan,
                "corr": build_cluster_feature_correlation(feature_df, labels.astype(str)),
                "feature_matrix": build_cluster_feature_matrix(feature_df, labels),
            }]
            best_result = results[0]
            print(f"    Fixed Leiden resolution={fixed_res}: n_clusters={metrics['n_clusters']}")
    else:
        raise ValueError(f"Unknown clustering_method '{clustering_method}'. Use 'kmeans' or 'leiden'.")

    adata.obs["ClusterID"] = pd.Categorical(adata.obs[best_result["col"]].astype(str))
    if method == "kmeans":
        best_sil = best_result["metrics"]["silhouette"]
        best_db = best_result["metrics"]["davies_bouldin"]
        print(
            f"    ✓ Final selection: {best_result['param_name']}={best_result['param_value']} "
            f"(Sil={best_sil if best_sil is not None else 'N/A'}, DB={best_db if best_db is not None else 'N/A'})"
        )
    else:
        print(
            f"    ✓ Final selection: {best_result['param_name']}={best_result['param_value']} "
            f"(n_clusters={best_result['metrics']['n_clusters']})"
        )

    # Save per-group diagnostics
    if outfolder is not None and group_name is not None:
        outfolder = Path(outfolder)
        outfolder.mkdir(parents=True, exist_ok=True)
        if method == "leiden":
            stability_summary = adata.uns.get("leiden_stability_summary", None)
            if stability_summary is not None and len(stability_summary) > 0:
                csv_path = outfolder / f"recursive_leiden_stability_summary_{group_name}.csv"
                stability_summary.to_csv(csv_path, index=False)
                print(f"Saved Leiden stability summary CSV to {csv_path}")
        pdf_path = outfolder / f"recursive_{clustering_method}_steps_{group_name}.pdf"
        with PdfPages(pdf_path) as pdf:
            if method == "leiden":
                plot_leiden_diagnostics(
                    pdf=pdf,
                    adata=adata,
                    results=results,
                    best_result=best_result,
                )
            else:
                plot_kmeans_diagnostics(
                    pdf=pdf,
                    adata=adata,
                    results=results,
                    best_result=best_result,
                    min_silhouette=min_silhouette,
                    max_davies_bouldin=max_davies_bouldin,
                )

    # Clean temporary clustering columns
    for res in results:
        if res["col"] in adata.obs.columns and res["col"] != "ClusterID":
            del adata.obs[res["col"]]
    if "leiden_auto_selection_tmp" in adata.obs.columns:
        del adata.obs["leiden_auto_selection_tmp"]

    return adata


def plot_group_umaps(group_results, ncols=3, figsize_per_plot=(5, 5)):
    """
    Plot UMAP results for all groups in a grid.
    
    Parameters
    ----------
    group_results : dict
        Dictionary mapping group names to AnnData objects
    ncols : int
        Number of columns in the grid
    figsize_per_plot : tuple
        Size of each subplot (width, height)
        
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    n_groups = len(group_results)
    nrows = int(np.ceil(n_groups / ncols))
    
    fig, axes = plt.subplots(
        nrows, ncols, 
        figsize=(figsize_per_plot[0] * ncols, figsize_per_plot[1] * nrows)
    )
    
    # Flatten axes for easier iteration
    if n_groups == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    for idx, (group_name, adata) in enumerate(group_results.items()):
        ax = axes[idx]
        
        if adata is None:
            ax.text(0.5, 0.5, f"{group_name}\n(insufficient data)", 
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_xticks([])
            ax.set_yticks([])
            continue
        
        # Get UMAP coordinates
        umap_coords = adata.obsm["X_umap"]
        cluster_labels = adata.obs["ClusterID"].astype(str)
        
        # Plot
        scatter = ax.scatter(
            umap_coords[:, 0], 
            umap_coords[:, 1],
            c=cluster_labels.astype('category').cat.codes,
            cmap='tab20',
            s=2,
            alpha=0.5
        )
        
        ax.set_title(f"{group_name}\n(n={len(adata)})", fontsize=10, fontweight='bold')
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_aspect("equal", adjustable="box")
        
        # Add legend
        unique_clusters = cluster_labels.unique()
        if len(unique_clusters) <= 20:  # Only show legend if not too many clusters
            handles = [plt.Line2D([0], [0], marker='o', color='w', 
                                markerfacecolor=plt.cm.tab20(i % 20), 
                                markersize=8, label=cluster)
                      for i, cluster in enumerate(sorted(unique_clusters))]
            ax.legend(handles=handles, loc='best', fontsize=6, 
                     title='Cluster', title_fontsize=7, framealpha=0.8)
    
    # Hide unused subplots
    for idx in range(n_groups, len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    return fig



def plot_feature_distributions_to_pdf(df, feature_cols, pdf, title_prefix="Feature Distributions"):
    """
    Plot feature distributions and append the figure as one page in a PDF.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe.
    feature_cols : list
        Feature columns to plot.
    pdf : matplotlib.backends.backend_pdf.PdfPages
        Open PDF writer where the page should be saved.
    title_prefix : str
        Title suffix to annotate plot context.
    """
    if len(feature_cols) == 0:
        return

    n_cols = 3
    n_rows = 4
    plots_per_page = n_cols * n_rows

    for start in range(0, len(feature_cols), plots_per_page):
        chunk = feature_cols[start:start + plots_per_page]
        fig, axes = plt.subplots(n_rows, n_cols, figsize=A4_PORTRAIT)
        axes = np.array(axes).flatten()

        for i, feature in enumerate(chunk):
            sns.violinplot(y=df[feature], ax=axes[i])
            axes[i].set_title(f"{feature}\n({title_prefix})", fontsize=9)

        for i in range(len(chunk), len(axes)):
            axes[i].axis("off")

        fig.tight_layout()
        save_pdf_page_a4(pdf, fig, orientation="portrait")
        plt.close(fig)


def cap_values_to_quantile(df, features, lower_quantile=None, upper_quantile=0.99, return_limits=False):
    """
    Cap features to the specified quantiles.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe
    features : list
        List of feature columns to cap
    lower_quantile : float
        Lower quantile threshold (0 to 1). Values below this will be set to the quantile value.
        If None, no lower capping is performed.
    upper_quantile : float
        Upper quantile threshold (0 to 1). Values above this will be set to the quantile value.
        If None, no upper capping is performed.
        
    return_limits : bool
        If True, also return per-feature cap limits used during clipping.

    Returns
    -------
    df_capped : pd.DataFrame
        Dataframe with capped values
    limits : dict, optional
        Returned only when `return_limits=True`. Mapping:
        {feature: {"lower": float|None, "upper": float|None}}
    """
    if lower_quantile is None and upper_quantile is None:
        if return_limits:
            return df, {}
        return df
        
    df_capped = df.copy()
    
    if lower_quantile is not None:
        print(f"Capping features to {lower_quantile*100}th lower quantile...")
    if upper_quantile is not None:
        print(f"Capping features to {upper_quantile*100}th upper quantile...")
    
    skipped_non_numeric = []
    skipped_boolean = []
    cap_limits = {}

    for feature in features:
        if feature in df_capped.columns:
            col = df_capped[feature]

            # Quantile interpolation on boolean data fails in newer numpy/pandas.
            # Binary/boolean features should not be outlier-capped.
            if pd.api.types.is_bool_dtype(col):
                skipped_boolean.append(feature)
                if return_limits:
                    cap_limits[feature] = {"lower": None, "upper": None}
                continue

            # Ensure numeric dtype for quantile computation.
            series = pd.to_numeric(col, errors="coerce").dropna()
            if len(series) > 0:
                lower_limit = None
                upper_limit = None
                
                if lower_quantile is not None:
                    lower_limit = series.quantile(lower_quantile)
                    
                if upper_quantile is not None:
                    upper_limit = series.quantile(upper_quantile)
                
                # Clip values using numeric representation to avoid dtype issues.
                numeric_col = pd.to_numeric(df_capped[feature], errors="coerce")
                df_capped[feature] = numeric_col.clip(lower=lower_limit, upper=upper_limit)
                if return_limits:
                    cap_limits[feature] = {
                        "lower": None if lower_limit is None else float(lower_limit),
                        "upper": None if upper_limit is None else float(upper_limit),
                    }
            else:
                skipped_non_numeric.append(feature)
                if return_limits:
                    cap_limits[feature] = {"lower": None, "upper": None}

    if skipped_boolean:
        print(f"Skipped boolean features for quantile capping: {len(skipped_boolean)}")
    if skipped_non_numeric:
        print(f"Skipped non-numeric features for quantile capping: {len(skipped_non_numeric)}")
            
    if return_limits:
        return df_capped, cap_limits
    return df_capped


def _to_tczyx(arr, name="image"):
    """Coerce image-like array to TCZYX shape."""
    if arr.ndim == 5:
        return arr
    if arr.ndim == 4:
        # Assume TZYX
        return np.expand_dims(arr, axis=1)
    if arr.ndim == 3:
        # Assume ZYX
        return np.expand_dims(np.expand_dims(arr, axis=0), axis=1)
    raise ValueError(f"{name} must be 3D/4D/5D. Got shape={getattr(arr, 'shape', None)}")


def _resolve_track_image_column(df_sample_meta, track_image_col=None):
    """Resolve which metadata column contains tracked-label image path."""
    if track_image_col is not None:
        if track_image_col not in df_sample_meta.columns:
            raise ValueError(f"track_image_col '{track_image_col}' not found in metadata columns.")
        return track_image_col

    candidates = [
        c for c in df_sample_meta.columns
        if c.endswith("_tracks_image_path")
    ]
    candidates = [c for c in candidates if pd.notna(df_sample_meta.iloc[0][c])]
    if len(candidates) == 0:
        raise ValueError(
            "Could not auto-detect track image column. "
            "Pass track_image_col explicitly (e.g. 'tcell_tracks_image_path')."
        )
    if len(candidates) > 1:
        raise ValueError(
            f"Multiple track image path columns found: {candidates}. "
            "Pass track_image_col explicitly."
        )
    return candidates[0]


def open_cluster_backprojection_napari(
    adata,
    metadata,
    outfolder=None,
    sample_name=None,
    cluster_col="ClusterID",
    sample_col="sample_name",
    trackid_col="TrackID",
    time_col="position_t",
    raw_image_col="raw_image_path",
    track_image_col="tcell_tracks_image_path",
    track_channel=0,
    show_trackid_layer=True,
    show_cluster_mapping_widget=True,
    fill_missing_timepoints_per_track=True,
):
    """
    Open napari with raw channels split into separate layers and ClusterID overlay labels.

    Cluster labels are mapped per voxel using (sample_name, TrackID, position_t) from adata.obs.
    """
    try:
        import napari
    except ImportError as exc:
        raise ImportError("napari is required for visualization. Install with `pip install napari`.") from exc

    if isinstance(metadata, (str, Path)):
        metadata = load_behav3d_metadata(metadata)
    if not isinstance(metadata, pd.DataFrame):
        raise TypeError("metadata must be a pandas DataFrame or path to metadata csv.")

    if sample_col not in adata.obs.columns:
        raise ValueError(f"adata.obs is missing required sample column '{sample_col}'.")
    if cluster_col not in adata.obs.columns:
        raise ValueError(f"adata.obs is missing required cluster column '{cluster_col}'.")
    if trackid_col not in adata.obs.columns:
        raise ValueError(f"adata.obs is missing required trackid column '{trackid_col}'.")
    if time_col not in adata.obs.columns:
        raise ValueError(f"adata.obs is missing required time column '{time_col}'.")

    adata_samples = pd.Series(adata.obs[sample_col].astype(str).unique()).sort_values().tolist()
    if sample_name is None:
        if len(adata_samples) == 0:
            raise ValueError("No sample_name values found in adata.obs.")
        sample_name = adata_samples[0]
        print(f"sample_name not provided. Using first sample: {sample_name}")
    sample_name = str(sample_name)

    if sample_name not in adata_samples:
        raise ValueError(f"sample_name '{sample_name}' not found in adata.obs[{sample_col}].")

    if "sample_name" not in metadata.columns:
        raise ValueError("metadata must contain 'sample_name' column.")
    df_sample_meta = metadata[metadata["sample_name"].astype(str) == sample_name]
    if df_sample_meta.empty:
        raise ValueError(f"sample_name '{sample_name}' not found in metadata.")
    df_sample_meta = df_sample_meta.iloc[[0]].copy()

    if raw_image_col not in df_sample_meta.columns:
        raise ValueError(f"raw_image_col '{raw_image_col}' not found in metadata.")

    track_image_col = _resolve_track_image_column(df_sample_meta, track_image_col=track_image_col)

    raw_img_path = Path(df_sample_meta.iloc[0][raw_image_col])
    track_img_path = Path(df_sample_meta.iloc[0][track_image_col])

    raw = _to_tczyx(load_image(raw_img_path), name="raw image")
    track = _to_tczyx(load_image(track_img_path), name="track image")

    if track.shape[1] <= track_channel:
        raise ValueError(
            f"track_channel={track_channel} out of bounds for track image with {track.shape[1]} channels."
        )

    track_tzyx = track[:, track_channel, ...]

    # relabel_img_from_adata currently expects numeric labels.
    # If cluster labels are strings (e.g. "no_contact_1"), encode to ints for display.
    overlay_obs_col = cluster_col
    cluster_id_mapping = None
    if not pd.api.types.is_numeric_dtype(adata.obs[cluster_col]):
        overlay_obs_col = f"__{cluster_col}_napari_codes__"
        adata = adata.copy()
        cluster_values = adata.obs[cluster_col].astype(str)
        categories = sorted(cluster_values.dropna().unique().tolist())
        cluster_id_mapping = {cat: i + 1 for i, cat in enumerate(categories)}
        adata.obs[overlay_obs_col] = cluster_values.map(cluster_id_mapping).astype(np.int64)

    # Align adata time values to image frame index if needed.
    adata = adata.copy()
    sample_mask = adata.obs[sample_col].astype(str) == str(sample_name)
    sample_times = pd.to_numeric(adata.obs.loc[sample_mask, time_col], errors="coerce").dropna().astype(int)
    n_frames = int(track_tzyx.shape[0])
    if len(sample_times) > 0:
        unique_times = np.sort(sample_times.unique())
        expected = np.arange(n_frames)
        if not (len(unique_times) > 0 and np.array_equal(unique_times, expected)):
            if np.array_equal(unique_times, np.arange(1, n_frames + 1)):
                # Common case: timepoints are 1..T instead of 0..T-1
                adata.obs.loc[sample_mask, time_col] = (
                    pd.to_numeric(adata.obs.loc[sample_mask, time_col], errors="coerce") - 1
                )
                print("Adjusted time index from 1..T to 0..T-1 for backprojection.")
            elif len(unique_times) <= n_frames:
                # Robust fallback: map sorted unique time values to sequential frame indices.
                tmap = {t: i for i, t in enumerate(unique_times)}
                mapped = pd.to_numeric(adata.obs.loc[sample_mask, time_col], errors="coerce").map(tmap)
                adata.obs.loc[sample_mask, time_col] = mapped
                print(
                    "Mapped non-sequential time values to frame indices for backprojection "
                    f"(unique_timepoints={len(unique_times)}, n_frames={n_frames})."
                )

    # Quick diagnostics for TrackID overlap.
    img_track_ids = np.unique(track_tzyx.astype(np.int64))
    if hasattr(img_track_ids, "compute"):
        img_track_ids = img_track_ids.compute()
    img_track_ids = np.asarray(img_track_ids)
    img_track_ids = set(img_track_ids[img_track_ids > 0].tolist())
    obs_track_ids = pd.to_numeric(adata.obs.loc[sample_mask, trackid_col], errors="coerce").dropna().astype(int)
    obs_track_ids = set(obs_track_ids.tolist())
    overlap = len(img_track_ids.intersection(obs_track_ids))
    if overlap == 0:
        print(
            "Warning: no TrackID overlap between image labels and adata obs for this sample. "
            "Overlay will be empty."
        )

    # Build overlay via explicit (time, TrackID) mapping so we can optionally fill
    # missing timepoints per track before rasterization.
    df_overlay = adata.obs[[sample_col, trackid_col, time_col, overlay_obs_col]].copy()
    df_overlay = df_overlay[df_overlay[sample_col].astype(str) == str(sample_name)].copy()
    df_overlay[trackid_col] = pd.to_numeric(df_overlay[trackid_col], errors="coerce").astype("Int64")
    df_overlay[time_col] = pd.to_numeric(df_overlay[time_col], errors="coerce").astype("Int64")
    df_overlay[overlay_obs_col] = pd.to_numeric(df_overlay[overlay_obs_col], errors="coerce").astype("Int64")
    df_overlay = df_overlay.dropna(subset=[trackid_col, time_col, overlay_obs_col]).copy()
    df_overlay[trackid_col] = df_overlay[trackid_col].astype(np.int64)
    df_overlay[time_col] = df_overlay[time_col].astype(np.int64)
    df_overlay[overlay_obs_col] = df_overlay[overlay_obs_col].astype(np.int64)

    if fill_missing_timepoints_per_track and len(df_overlay) > 0:
        frames = np.arange(n_frames, dtype=np.int64)
        filled_parts = []
        for tid, g in df_overlay.groupby(trackid_col):
            g = g[[time_col, overlay_obs_col]].drop_duplicates(subset=[time_col], keep="last").set_index(time_col).sort_index()
            g = g.reindex(frames)
            g[overlay_obs_col] = g[overlay_obs_col].ffill().bfill()
            g = g.dropna(subset=[overlay_obs_col])
            if len(g) == 0:
                continue
            g = g.reset_index().rename(columns={"index": time_col})
            g[trackid_col] = int(tid)
            filled_parts.append(g[[trackid_col, time_col, overlay_obs_col]])
        if len(filled_parts) > 0:
            df_overlay = pd.concat(filled_parts, ignore_index=True)

    overlay_by_time = {}
    for t, g in df_overlay.groupby(time_col, observed=False):
        overlay_by_time[int(t)] = (
            g[trackid_col].to_numpy(dtype=np.int64),
            g[overlay_obs_col].to_numpy(dtype=np.int64),
        )

    labels_cluster = np.zeros(track_tzyx.shape, dtype=np.int32)
    for t in range(n_frames):
        labels_t = track_tzyx[t]
        if hasattr(labels_t, "compute"):
            labels_t = labels_t.compute()
        labels_t = np.asarray(labels_t, dtype=np.int64)
        max_label_t = int(labels_t.max()) if labels_t.size else 0
        if max_label_t <= 0:
            continue
        lut = np.zeros(max_label_t + 1, dtype=np.int32)
        tids_vals = overlay_by_time.get(t, None)
        if tids_vals is not None:
            tids, vals = tids_vals
            valid = (tids >= 0) & (tids <= max_label_t)
            if np.any(valid):
                lut[tids[valid]] = vals[valid].astype(np.int32)
        labels_cluster[t] = lut[labels_t]

    scale_tzyx = None
    if all(c in df_sample_meta.columns for c in ["pixel_distance_z", "pixel_distance_xy"]):
        z = float(df_sample_meta.iloc[0]["pixel_distance_z"])
        xy = float(df_sample_meta.iloc[0]["pixel_distance_xy"])
        scale_tzyx = (1.0, z, xy, xy)

    viewer = napari.Viewer(title=f"Cluster backprojection: {sample_name}")

    for c in range(raw.shape[1]):
        viewer.add_image(
            raw[:, c, ...],
            name=f"raw_ch{c}",
            scale=scale_tzyx,
            blending="additive",
        )

    if show_trackid_layer:
        viewer.add_labels(
            track_tzyx.astype(np.int32),
            name="TrackID",
            scale=scale_tzyx,
            opacity=0.25,
        )

    viewer.add_labels(
        labels_cluster.astype(np.int32),
        name=f"{cluster_col}_overlay",
        scale=scale_tzyx,
        opacity=0.5,
    )

    if cluster_id_mapping is not None:
        viewer.layers[f"{cluster_col}_overlay"].metadata["cluster_id_mapping"] = cluster_id_mapping
        print(
            f"Encoded non-numeric '{cluster_col}' for napari overlay. "
            f"Mapping stored in layer metadata under 'cluster_id_mapping'."
        )
        if show_cluster_mapping_widget:
            try:
                from qtpy.QtWidgets import QWidget, QVBoxLayout, QLabel, QTableWidget, QTableWidgetItem

                mapping_widget = QWidget()
                layout = QVBoxLayout(mapping_widget)
                layout.addWidget(QLabel(f"{cluster_col} mapping"))

                table = QTableWidget(len(cluster_id_mapping), 2)
                table.setHorizontalHeaderLabels(["Code", "Cluster Name"])
                table.verticalHeader().setVisible(False)

                for row, (name, code) in enumerate(sorted(cluster_id_mapping.items(), key=lambda x: x[1])):
                    table.setItem(row, 0, QTableWidgetItem(str(code)))
                    table.setItem(row, 1, QTableWidgetItem(str(name)))

                table.resizeColumnsToContents()
                layout.addWidget(table)
                viewer.window.add_dock_widget(mapping_widget, area="right", name=f"{cluster_col} Mapping")
            except Exception as exc:
                print(f"Could not create napari mapping widget ({exc}).")

    return viewer


def _group_to_cluster_prefix(group_name: str) -> str:
    """Convert internal group name to a compact ClusterID prefix."""
    prefix = str(group_name)
    prefix = prefix.replace("_pixels", "")
    prefix = prefix.replace("_only", "")
    prefix = prefix.replace("_and_", "_")
    prefix = prefix.strip("_")
    if prefix == "":
        prefix = "group"
    return prefix


def _resolve_group_leiden_resolution(resolution, group_name):
    """
    Resolve Leiden resolution setting for one group.

    Supported values for `resolution`:
    - "auto"
    - float/int
    - dict[group_name_or_prefix -> float]
    """
    if isinstance(resolution, dict):
        if group_name in resolution:
            return resolution[group_name]
        group_prefix = _group_to_cluster_prefix(group_name)
        if group_prefix in resolution:
            return resolution[group_prefix]
        for k, v in resolution.items():
            k = str(k)
            if group_name.startswith(k) or group_prefix.startswith(k):
                return v
        raise ValueError(
            f"Missing Leiden resolution for group '{group_name}'. "
            f"Provided keys: {list(resolution.keys())}"
        )
    return resolution


def _prepare_state_classification_dataset(
    df_positions,
    features,
    binary_features_to_group,
    window_size=5,
    descriptive_features=["mean", "median", "std", "net_displacement", "straightness", "mean_square_displacement"],
    lower_quantile_cap=None,
    upper_quantile_cap=0.99,
    outfolder=None,
    scale_features=True,
    incomplete_window_policy="drop",
    prepared_dataset_cache_path=None,
    reuse_prepared_dataset=True,
    save_prepared_dataset=True,
):
    """Shared preprocessing for state classification and full-dataset inference."""
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
        if prepared_dataset_cache_path.suffix.lower() != ".csv":
            prepared_dataset_cache_path = prepared_dataset_cache_path.with_suffix(".csv")

    loaded_from_cache = False
    quantile_cap_limits = {}
    if (
        reuse_prepared_dataset
        and prepared_dataset_cache_path is not None
        and prepared_dataset_cache_path.exists()
    ):
        try:
            df_analysis = pd.read_csv(prepared_dataset_cache_path)
            binary_cols_to_merge = [col for col in binary_features_to_group if col in df_analysis.columns]
            descriptive_feature_cols = [
                col for col in df_analysis.columns
                if (col not in non_feature_cols)
                and (not col.endswith("_signal_type"))
            ]
            binary_prefixes = tuple(f"{c}_" for c in binary_cols_to_merge)
            kept_features = [
                c for c in descriptive_feature_cols
                if (c not in binary_cols_to_merge)
                and (not c.startswith(binary_prefixes))
            ]
            df_windows_descriptive = df_analysis.copy()
            loaded_from_cache = True
            print(f"Loaded prepared window dataset from cache: {prepared_dataset_cache_path}")
        except Exception as exc:
            print(f"Could not load prepared dataset cache ({exc}); recomputing windows.")

    if not loaded_from_cache:
        print("Creating descriptive track dataset...")
        df_windows_descriptive = create_descriptive_track_dataset(
            df_tracks=df_positions,
            columns_to_summarize=features,
            window_size=window_size,
            step_size=1,
            time_col="position_t",
            id_cols=groupby,
            features_to_compute=descriptive_features,
            only_nonbinary=True,
            incomplete_window_policy=incomplete_window_policy,
        )

        descriptive_feature_cols = [
            col for col in df_windows_descriptive.columns
            if (col not in non_feature_cols)
            and (not col.endswith("_signal_type"))
        ]

        binary_cols_to_merge = [col for col in binary_features_to_group if col in df_positions.columns]
        merge_cols = ["sample_name", "TrackID", "position_t"]
        df_binary = df_positions[merge_cols + binary_cols_to_merge].copy()

        df_analysis = df_windows_descriptive.merge(
            df_binary,
            on=merge_cols,
            how="left",
            suffixes=("", "_orig"),
        )
        policy = str(incomplete_window_policy).lower()
        if policy == "drop":
            df_analysis = df_analysis.dropna(subset=descriptive_feature_cols).copy()
        elif policy in {"partial", "as_far_as_possible"}:
            df_analysis = df_analysis.copy()
        else:
            raise ValueError(
                "incomplete_window_policy must be one of: "
                "'drop', 'partial' (alias: 'as_far_as_possible')."
            )

        # Drop descriptive feature columns that are entirely empty.
        empty_descriptive_cols = [
            c for c in descriptive_feature_cols
            if c in df_analysis.columns and df_analysis[c].isna().all()
        ]
        if len(empty_descriptive_cols) > 0:
            print(f"Dropping {len(empty_descriptive_cols)} empty descriptive columns.")
            df_analysis = df_analysis.drop(columns=empty_descriptive_cols, errors="ignore")
            df_windows_descriptive = df_windows_descriptive.drop(columns=empty_descriptive_cols, errors="ignore")
            descriptive_feature_cols = [c for c in descriptive_feature_cols if c not in empty_descriptive_cols]

        binary_prefixes = tuple(f"{c}_" for c in binary_cols_to_merge)
        kept_features = [
            c for c in descriptive_feature_cols
            if (c not in binary_cols_to_merge)
            and (not c.startswith(binary_prefixes))
        ]
        if len(kept_features) != len(descriptive_feature_cols):
            removed = sorted(set(descriptive_feature_cols) - set(kept_features))
            print(f"Excluded binary-derived grouping features from clustering feature set: {removed}")

        if len(kept_features) > 0:
            if lower_quantile_cap is not None or upper_quantile_cap is not None:
                if outfolder is not None:
                    outfolder = Path(outfolder)
                    outfolder.mkdir(parents=True, exist_ok=True)
                    capping_pdf_path = outfolder / "feature_distributions_quantile_capping_comparison.pdf"
                    print("Plotting feature distributions before and after quantile capping...")
                    with PdfPages(capping_pdf_path) as pdf:
                        plot_feature_distributions_to_pdf(
                            df=df_analysis,
                            feature_cols=kept_features,
                            pdf=pdf,
                            title_prefix="Before Quantile Capping",
                        )
                        df_analysis, quantile_cap_limits = cap_values_to_quantile(
                            df_analysis,
                            kept_features,
                            lower_quantile=lower_quantile_cap,
                            upper_quantile=upper_quantile_cap,
                            return_limits=True,
                        )
                        plot_feature_distributions_to_pdf(
                            df=df_analysis,
                            feature_cols=kept_features,
                            pdf=pdf,
                            title_prefix="After Quantile Capping",
                        )
                    print(f"Saved quantile capping comparison to {capping_pdf_path}")
                else:
                    df_analysis, quantile_cap_limits = cap_values_to_quantile(
                        df_analysis,
                        kept_features,
                        lower_quantile=lower_quantile_cap,
                        upper_quantile=upper_quantile_cap,
                        return_limits=True,
                    )

        if (
            save_prepared_dataset
            and prepared_dataset_cache_path is not None
        ):
            try:
                prepared_dataset_cache_path.parent.mkdir(parents=True, exist_ok=True)
                df_analysis.to_csv(prepared_dataset_cache_path, index=False)
                print(f"Saved prepared window dataset cache: {prepared_dataset_cache_path}")
            except Exception as exc:
                print(f"Warning: failed to save prepared dataset cache ({exc}).")

    scaler = None
    if scale_features and len(kept_features) > 0:
        print("Scaling features globally...")
        scaler = StandardScaler().fit(df_analysis[kept_features])
        df_analysis[kept_features] = scaler.transform(df_analysis[kept_features])

    return {
        "df_windows_descriptive": df_windows_descriptive,
        "df_analysis": df_analysis,
        "kept_features": kept_features,
        "non_feature_cols": non_feature_cols,
        "binary_cols_to_merge": binary_cols_to_merge,
        "scaler": scaler,
        "quantile_cap_limits": quantile_cap_limits,
        "lower_quantile_cap": lower_quantile_cap,
        "upper_quantile_cap": upper_quantile_cap,
    }


def _infer_full_dataset_from_group_models(
    df_analysis,
    group_reference_adatas,
    binary_cols_to_merge,
    non_feature_cols,
    cluster_col="ClusterID",
):
    """Infer per-timepoint cluster labels on full dataset using trained per-group references."""
    if not isinstance(group_reference_adatas, dict) or len(group_reference_adatas) == 0:
        raise ValueError("group_reference_adatas must be a non-empty dict of {group_name: adata_or_path}.")

    ref_models = {}
    common_feature_order = None
    for group_name, ref in group_reference_adatas.items():
        if ref is None:
            continue
        ref_adata = sc.read_h5ad(ref) if isinstance(ref, (str, Path)) else ref
        if ref_adata is None or ref_adata.n_obs == 0:
            continue
        if cluster_col not in ref_adata.obs.columns:
            raise ValueError(f"Reference model '{group_name}' is missing '{cluster_col}' in obs.")
        feature_order = list(ref_adata.var_names)
        if common_feature_order is None:
            common_feature_order = feature_order
        elif feature_order != common_feature_order:
            raise ValueError(
                "All reference group models must share identical feature columns/order "
                "to build one merged full AnnData."
            )
        ref_models[group_name] = ref_adata

    if len(ref_models) == 0:
        raise ValueError("No usable reference group models provided.")

    missing_features = [c for c in common_feature_order if c not in df_analysis.columns]
    if missing_features:
        raise ValueError(f"Full dataset is missing required model features: {missing_features[:10]}")

    groups = create_binary_groups(df_analysis, binary_cols_to_merge)
    df_out = df_analysis.copy()
    df_out["ClusterID"] = "unassigned"
    df_out["ClusterGroup"] = "unassigned"

    for group_name, mask in groups.items():
        if group_name not in ref_models:
            continue
        ref_adata = ref_models[group_name]
        idx = df_out.index[mask]
        if len(idx) == 0:
            continue

        df_group = df_out.loc[idx].copy()
        feature_order = list(ref_adata.var_names)

        scaler_meta = ref_adata.uns.get("preprocessing", {}).get("scaler", None)
        if scaler_meta is not None:
            means = np.asarray(scaler_meta["mean"], dtype=float)
            scales = np.asarray(scaler_meta["scale"], dtype=float)
            if len(means) != len(feature_order) or len(scales) != len(feature_order):
                raise ValueError(f"Scaler metadata shape mismatch for group '{group_name}'.")
            X = df_group[feature_order].to_numpy(dtype=float)
            df_group.loc[:, feature_order] = (X - means) / scales

        obs_cols = [c for c in non_feature_cols if c in df_group.columns] + [
            c for c in binary_cols_to_merge if c in df_group.columns
        ]
        adata_query = df_to_adata(df_group, feature_order, obs_cols=obs_cols)
        sc.tl.ingest(adata_query, ref_adata, obs=cluster_col, embedding_method="umap")

        prefix = _group_to_cluster_prefix(group_name)
        assigned = adata_query.obs[cluster_col].astype(str).map(lambda x: f"{prefix}_{x}")
        df_out.loc[idx, "ClusterID"] = assigned.to_numpy()
        df_out.loc[idx, "ClusterGroup"] = group_name

    merged_obs_cols = [c for c in df_out.columns if c not in common_feature_order]
    adata_full = df_to_adata(df_out, common_feature_order, obs_cols=merged_obs_cols)
    adata_full.obs["ClusterID"] = adata_full.obs["ClusterID"].astype("category")
    adata_full.obs["ClusterGroup"] = adata_full.obs["ClusterGroup"].astype("category")
    return adata_full


def apply_trained_state_clustering_to_full_dataset(
    df_positions,
    group_reference_adatas,
    features,
    binary_features_to_group,
    window_size=5,
    descriptive_features=["mean", "median", "std", "net_displacement", "straightness", "mean_square_displacement"],
    cluster_col="ClusterID",
    incomplete_window_policy="drop",
    prepared_dataset_cache_path=None,
    reuse_prepared_dataset=True,
    save_prepared_dataset=True,
):
    """
    Apply trained per-group clustering models to the full (non-subsampled) dataset.

    Parameters
    ----------
    df_positions : pd.DataFrame
        Full positions dataframe.
    group_reference_adatas : dict
        Mapping {group_name: AnnData or .h5ad path}. Each reference adata must contain
        trained cluster labels in `cluster_col`.
    features : list
        Base columns to summarize in rolling windows (same as training).
    binary_features_to_group : list
        Binary columns used for grouping.
    window_size : int
        Rolling window size (must match training).
    descriptive_features : list
        Descriptive summary features (must match training).
    cluster_col : str
        Cluster label column in reference adatas.
    incomplete_window_policy : str
        How to handle timepoints without a full descriptive window:
        - "drop": remove those timepoints
        - "partial": keep and fill from neighboring valid rows per track

    Returns
    -------
    adata_full : AnnData
        Full dataset with inferred ClusterID per timepoint.
    """
    prepared = _prepare_state_classification_dataset(
        df_positions=df_positions,
        features=features,
        binary_features_to_group=binary_features_to_group,
        window_size=window_size,
        descriptive_features=descriptive_features,
        lower_quantile_cap=None,
        upper_quantile_cap=None,
        outfolder=None,
        scale_features=True,
        incomplete_window_policy=incomplete_window_policy,
        prepared_dataset_cache_path=prepared_dataset_cache_path,
        reuse_prepared_dataset=reuse_prepared_dataset,
        save_prepared_dataset=save_prepared_dataset,
    )
    adata_full = _infer_full_dataset_from_group_models(
        df_analysis=prepared["df_analysis"],
        group_reference_adatas=group_reference_adatas,
        binary_cols_to_merge=prepared["binary_cols_to_merge"],
        non_feature_cols=prepared["non_feature_cols"],
        cluster_col=cluster_col,
    )
    adata_full.uns["inference"] = {
        "cluster_col": cluster_col,
        "model_groups": sorted([k for k, v in group_reference_adatas.items() if v is not None]),
        "window_size": int(window_size),
        "descriptive_features": list(descriptive_features),
    }
    return adata_full


def run_state_classification(
    df_positions,
    features,
    binary_features_to_group,
    window_size=5,
    chosen_intervals=10,
    max_samples=None,
    min_spacing=None,
    n_neighbors=15,
    resolution="auto",
    descriptive_features=["mean", "median", "std", "net_displacement", "straightness", "mean_square_displacement"],
    pca_var_selection=0.95,
    clustering_method="kmeans",
    resolutions=(0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5),
    leiden_seed_tries=10,
    leiden_subsample_tries=10,
    lower_quantile_cap=None,
    upper_quantile_cap=0.99,
    incomplete_window_policy="drop",
    outfolder=None,
    random_state=123,
    prepared_dataset_cache_path=None,
    reuse_prepared_dataset=True,
    save_prepared_dataset=True,
):
    """
    Main function to run state classification grouped by binary features.
    
    Parameters
    ----------
    df_positions : pd.DataFrame
        Input dataframe with track positions and features
    features : list
        List of all feature column names
    binary_features_to_group : list
        List of binary feature names to use for grouping
    window_size : int
        Window size for rolling features
    chosen_intervals : int
        Interval for subsampling timepoints (deprecated, kept for compatibility)
    max_samples : int, optional
        Maximum total samples per group after temporal-spacing subsampling.
        If None, no global cap is applied.
    min_spacing : int, optional
        Minimum temporal spacing between samples from the same track.
        If None, will be automatically estimated based on group size.
    n_neighbors : int
        Number of neighbors for Leiden clustering
    resolution : str | float | dict
        Leiden resolution mode:
        - "auto": use stability auto-tuning across `resolutions`
        - float/int: fixed resolution for all groups
        - dict: per-group fixed resolution (keys can be full group names or prefixes)
    pca_var_selection : float
        Variance threshold for PCA
    clustering_method : str
        Method to use: "kmeans" or "leiden"
    resolutions : tuple
        Resolutions to test if using Leiden
    leiden_seed_tries : int
        Number of random-seed repeats for Leiden stability scan.
    leiden_subsample_tries : int
        Number of subsample repeats for Leiden stability scan.
    lower_quantile_cap : float, optional
        Lower quantile to cap features at (default None).
    upper_quantile_cap : float, optional
        Upper quantile to cap features at (default 0.99).
    incomplete_window_policy : str
        How to handle timepoints without a full descriptive window:
        - "drop": remove those timepoints
        - "partial": keep and fill from neighboring valid rows per track
    outfolder : Path or str
        Output folder for saving results
    random_state : int
        Random seed
    Returns
    -------
    adata_merged : AnnData
        Merged full AnnData with prefixed ClusterID labels across binary groups.
    """
    
    prepared = _prepare_state_classification_dataset(
        df_positions=df_positions,
        features=features,
        binary_features_to_group=binary_features_to_group,
        window_size=window_size,
        descriptive_features=descriptive_features,
        lower_quantile_cap=lower_quantile_cap,
        upper_quantile_cap=upper_quantile_cap,
        outfolder=outfolder,
        scale_features=True,
        incomplete_window_policy=incomplete_window_policy,
        prepared_dataset_cache_path=prepared_dataset_cache_path,
        reuse_prepared_dataset=reuse_prepared_dataset,
        save_prepared_dataset=save_prepared_dataset,
    )
    df_windows_descriptive = prepared["df_windows_descriptive"]
    df_analysis = prepared["df_analysis"]
    kept_features = prepared["kept_features"]
    non_feature_cols = prepared["non_feature_cols"]
    binary_cols_to_merge = prepared["binary_cols_to_merge"]
    
    print("Creating binary feature groups...")
    groups = create_binary_groups(df_analysis, binary_cols_to_merge)
    
    print(f"Found {len(groups)} groups (before subsampling):")
    for group_name, mask in groups.items():
        print(f"  - {group_name}: {mask.sum()} samples")
    
    # Cluster each group with per-group subsampling
    group_results = {}
    for group_name, mask in groups.items():
        print(f"\nProcessing group: {group_name}")

        df_group_full = df_analysis[mask].copy()
        
        if len(df_group_full) < 50:
            print(f"  Skipping {group_name}: insufficient data ({len(df_group_full)} samples)")
            group_results[group_name] = None
            continue
        
        # Per-group temporal subsampling to reduce autocorrelation
        # Only calculate estimated spacing if not provided by user
        if min_spacing is None:
            target_samples = max(500, len(df_group_full) // 20)  # At least 500 or 5% of data
            spacing_to_use = max(1, len(df_group_full) // target_samples // 10)  # Rough estimate
        else:
            spacing_to_use = min_spacing
        
        print(f"  Subsampling with temporal spacing (min_spacing={spacing_to_use})...")
        df_group = subsample_with_temporal_spacing(
            df_group_full,
            id_cols=["sample_name", "TrackID"],
            time_col="position_t",
            min_spacing=spacing_to_use,
            max_samples=max_samples,
            random_state=random_state
        )
        
        print(f"  After subsampling: {len(df_group)} samples")
        
        if len(df_group) < 50:
            print(f"  Skipping {group_name}: insufficient data after subsampling ({len(df_group)} samples)")
            group_results[group_name] = None
            continue
        
        # Fit preprocessing on the subset used for clustering, then apply
        # that learned transform to define the per-group cluster model.
        scaler_group = None
        df_group_train = df_group.copy()
        if len(kept_features) > 0:
            scaler_group = StandardScaler().fit(df_group_train[kept_features])
            df_group_train[kept_features] = scaler_group.transform(df_group_train[kept_features])
        
        adata = cluster_group(
            df_group=df_group_train,
            feature_cols=kept_features,
            non_feature_cols=non_feature_cols,
            n_neighbors=n_neighbors,
            resolution=_resolve_group_leiden_resolution(resolution, group_name),
            pca_var_selection=pca_var_selection,
            clustering_method=clustering_method,
            resolutions=resolutions,
            leiden_seed_tries=leiden_seed_tries,
            leiden_subsample_tries=leiden_subsample_tries,
            outfolder=Path(outfolder) if outfolder else None,
            group_name=group_name,
            random_state=random_state
        )
        
        group_results[group_name] = adata
        
        if adata is not None:
            if "preprocessing" not in adata.uns:
                adata.uns["preprocessing"] = {}
            adata.uns["preprocessing"]["kept_features"] = list(kept_features)
            if scaler_group is not None:
                adata.uns["preprocessing"]["scaler"] = {
                    "mean": scaler_group.mean_.astype(float),
                    "scale": scaler_group.scale_.astype(float),
                }

            n_clusters = len(adata.obs["ClusterID"].unique())
            print(f"  Found {n_clusters} clusters")
    
    # Plot results
    print("\nPlotting UMAP grid...")
    fig = plot_group_umaps(group_results, ncols=3)
    
    trained_models = {k: v for k, v in group_results.items() if v is not None}
    if len(trained_models) > 0:
        adata_merged = _infer_full_dataset_from_group_models(
            df_analysis=df_analysis,
            group_reference_adatas=trained_models,
            binary_cols_to_merge=binary_cols_to_merge,
            non_feature_cols=non_feature_cols,
            cluster_col="ClusterID",
        )
        adata_merged.uns["inference"] = {
            "cluster_col": "ClusterID",
            "model_groups": sorted(list(trained_models.keys())),
            "window_size": int(window_size),
            "descriptive_features": list(descriptive_features),
        }
        print(f"Inferred full dataset labels: n_obs={adata_merged.n_obs}, n_vars={adata_merged.n_vars}")
    else:
        print("No groups produced cluster results; returning empty merged AnnData.")
        empty_df = df_analysis.iloc[:0].copy()
        if "ClusterID" not in empty_df.columns:
            empty_df["ClusterID"] = pd.Series(dtype="object")
        if "ClusterGroup" not in empty_df.columns:
            empty_df["ClusterGroup"] = pd.Series(dtype="object")
        merged_obs_cols = [c for c in empty_df.columns if c not in kept_features]
        adata_merged = df_to_adata(empty_df, kept_features, obs_cols=merged_obs_cols)

    if outfolder is not None:
        outfolder = Path(outfolder)
        outfolder.mkdir(parents=True, exist_ok=True)
        
        # Save figure
        fig_path = outfolder / "state_classification_umap_grid.pdf"
        save_figure_a4(fig, fig_path, orientation="auto", dpi=300)
        print(f"Saved UMAP grid to {fig_path}")
        
        # Save individual group results
        for group_name, adata in group_results.items():
            if adata is not None:
                group_path = outfolder / f"adata_{group_name}.h5ad"
                adata.write(group_path, compression="gzip")
                print(f"Saved {group_name} results to {group_path}")
        if adata_merged is not None:
            merged_path = outfolder / "adata_merged_group_clusters.h5ad"
            adata_merged.write(merged_path, compression="gzip")
            print(f"Saved merged full results to {merged_path}")
    
    plt.show()
    return adata_merged


def test():
    """
    Test function demonstrating usage of state_classification.
    """
    # Example paths (modify as needed)
    ssd_dir = Path("/Users/s.deblank-3/Downloads/BHVD_SB1_Exp003_62T_lowTcellDensity/behav3d")
    metadata_csv_path = ssd_dir / "metadata.csv"
    
    outfolder = ssd_dir / "state_classification"
    outfolder.mkdir(parents=True, exist_ok=True)
    
    # Load metadata
    metadata = load_behav3d_metadata(metadata_csv_path)
    
    # Load track features
    analysis_outdir = ssd_dir / "analysis" / "tcell"
    feature_outdir = analysis_outdir / "track_features"
    df_tracks_path = feature_outdir / "BEHAV3D_tcell_combined_track_features.csv"
    
    df_positions = pd.read_csv(df_tracks_path)
    df_positions = df_positions.sort_values(by=["sample_name", "TrackID", "position_t"])
    
    # Define features
    features = [
        "percentage_dead_mask",
        "organoid_contact_pixels",
        "tcell_contact_pixels",
        "speed",
    ]
    
    # Define binary features to use for grouping
    binary_features_to_group = [
        "organoid_contact",
        "tcell_contact",
    ]
    
    # Run state classification
    adata_merged = run_state_classification(
        df_positions=df_positions,
        features=features,
        binary_features_to_group=binary_features_to_group,
        window_size=5,
        chosen_intervals=10,
        n_neighbors=15,
        resolution=0.2,
        pca_var_selection=0.95,
        outfolder=outfolder,
        random_state=seed
    )
    
    # Additional analysis for merged output
    if adata_merged is not None and adata_merged.n_obs > 0:
        plot_top_ranking_features(
            adata_merged,
            groupby="ClusterID",
            n_features=20,
        )
        plt.suptitle("Merged - Top Ranking Features", y=1.02)
        plt.tight_layout()
        plt.show()

        plot_number_per_clusters(adata_merged.obs, cluster_col="ClusterID")
        plt.title("Merged - Cluster Sizes")
        plt.show()
    open_cluster_backprojection_napari(
        adata_merged,
        metadata,
        track_image_col="tcell_tracks_image_path",
        sample_name="ROCHE_JM1_Exp042-9_Img04_10T_HER2-I"
    )


if __name__ == "__main__":
    
    ssd_dir = r"/Volumes/T7_Sam/"
    # ssd_dir = r"F:/"
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
    
    binary_features_to_group = [
        "organoid_contact_pixels",
        "tcell_contact_pixels",
    ]

    # Create descriptive features per window
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
    window_size=5

    groupby = ["sample_name", "TrackID"]
    descriptive_features = ["mean", "median", "std", "net_displacement", "straightness", "mean_square_displacement"]
    # max_samples = 20000      # or None
    max_samples=None
    min_spacing = 10          # or None for auto
    n_neighbors = 60
    resolution = 0.2         # currently not used in leiden search path
    min_cluster_size = None  # currently unused
    pca_var_selection = 0.95
    clustering_method = "leiden"   # or "leiden"
    incomplete_window_policy="partial"
    clustering_method
    resolutions = (0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5)
    leiden_seed_tries = 20
    leiden_subsample_tries = 20
    lower_quantile_cap = None      # e.g. 0.01
    upper_quantile_cap = 0.99      # e.g. None to disable upper capping
    outfolder = Path(ssd_dir, "BHVD_BEHAV3D/BEHAV3D_python/rolling_classification")
    random_state = 12345
    reuse_prepared_dataset=True
    save_prepared_dataset=True
    pass
    # test()
    # adata_merged = sc.read_h5ad("/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/rolling_classification/adata_merged_group_clusters.h5ad")
