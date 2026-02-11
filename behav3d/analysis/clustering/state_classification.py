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
from behav3d.core.anndata import df_to_adata, adata_add_back_to_df, merge_pandas_cols_into_obs_anndata
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
    max_samples_per_track=None,
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
    max_samples_per_track : int, optional
        Maximum number of samples to take per track. If None, no limit.
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
            
            if max_samples_per_track is not None and len(selected_indices) >= max_samples_per_track:
                break
        
        subsampled_rows.append(track_df.iloc[selected_indices])
    
    if len(subsampled_rows) == 0:
        return df.iloc[:0].copy()  # Return empty dataframe with same columns
    
    return pd.concat(subsampled_rows, ignore_index=True)




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


def cluster_group(
    df_group,
    feature_cols,
    non_feature_cols,
    n_neighbors=15,
    resolution=0.2,
    min_cluster_size=100,
    pca_var_selection=0.95,
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
    
    # Scale features
    scaler = StandardScaler()
    df_clean[feature_cols] = scaler.fit_transform(df_clean[feature_cols])
    
    # Create AnnData object
    adata = df_to_adata(df_clean, feature_cols, obs_cols=non_feature_cols)
    adata.uns["preprocessing"] = {
        "kept_features": list(feature_cols),
        "scaler": {
            "mean": scaler.mean_.astype(float),
            "scale": scaler.scale_.astype(float),
        }
    }
    
    # Run PCA
    adata = run_pca(
        adata,
        pca_var_selection=pca_var_selection,
        ncomps=min(len(feature_cols), len(df_clean) - 1),
        svd_solver='full', 
        random_state=random_state
    )
    
    
    # Run Leiden clustering
    adata = run_leiden_clustering(
        adata, 
        n_neighbors=n_neighbors,
        resolution=resolution, 
        metric="euclidean",
        method="umap",
        use_rep="X_pca",
        key_added="ClusterID",
        random_state=random_state
    )
    
    # Merge small clusters
    adata = merge_small_clusters(
        adata,
        key="ClusterID",
        min_size=min_cluster_size,
        use_rep="X_pca",
    )
    
    # Validate clustering: only keep multiple clusters if there's clear evidence
    # Use PCA representation for validation
    X_pca = adata.obsm["X_pca"]
    labels = adata.obs["ClusterID"].values
    
    is_valid, metrics = validate_clustering(
        X_pca, 
        labels,
        min_silhouette=0.3,  # Require reasonable cluster separation
        max_davies_bouldin=1.5  # Require good cluster compactness
    )
    
    
    # Format metrics for display
    sil_str = f"{metrics['silhouette']:.3f}" if metrics['silhouette'] is not None else 'N/A'
    db_str = f"{metrics['davies_bouldin']:.3f}" if metrics['davies_bouldin'] is not None else 'N/A'
    
    print(f"    Clustering validation: n_clusters={metrics['n_clusters']}, "
          f"silhouette={sil_str}, davies_bouldin={db_str}")
    
    if not is_valid and metrics['n_clusters'] > 1:
        print("    ⚠ Clustering validation FAILED - collapsing to single cluster")
        print("      (No clear evidence for distinct subclusters)")
        # Collapse all to a single cluster
        adata.obs["ClusterID"] = "0"
    else:
        print(f"    ✓ Clustering validation PASSED - keeping {metrics['n_clusters']} cluster(s)")
    
    # Compute UMAP
    sc.tl.umap(
        adata,
        min_dist=0.1,
        random_state=random_state,
    )
    
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


def run_state_classification(
    df_positions,
    features,
    binary_features_to_group,
    window_size=5,
    chosen_intervals=10,
    min_spacing=None,
    n_neighbors=15,
    resolution=0.2,
    min_cluster_size=100,
    pca_var_selection=0.95,
    outfolder=None,
    random_state=123
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
    min_spacing : int, optional
        Minimum temporal spacing between samples from the same track.
        If None, will be automatically estimated based on group size.
    n_neighbors : int
        Number of neighbors for Leiden clustering
    resolution : float
        Resolution parameter for Leiden clustering
    min_cluster_size : int
        Minimum cluster size for merging
    pca_var_selection : float
        Variance threshold for PCA
    outfolder : Path or str
        Output folder for saving results
    random_state : int
        Random seed
        
    Returns
    -------
    group_results : dict
        Dictionary mapping group names to AnnData objects
    df_windows_descriptive : pd.DataFrame
        Dataframe with windowed descriptive features
    """

    binary_features_to_group = [
        "organoid_contact_pixels",
        "tcell_contact_pixels",
    ]

    # Create descriptive features per window
    features=[
        "percentage_dead_mask",
        # "mean_dead_dye",
        # "nr_dead_mask_pixels",
        "organoid_contact_pixels",
        "tcell_contact_pixels",
        # "mean_square_displacement",
        "speed",
        # "directional_persistence",
        "volume",
        "extent",
        "elongation",
        "sphericity",
        "solidity",
    ]
    window_size=5

    groupby = ["sample_name", "TrackID"]
    descriptive_features = ["mean", "median", "std", "net_displacement", "straightness", "mean_square_displacement"]
    
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
    )
    
    # Define non-feature columns
    non_feature_cols = [
        "sample_name",
        "TrackID",
        "sub_TrackID",
        "position_t",
        "window_start_position_t",
        "window_end_position_t",
        "window_length_frames",
    ]
    
    # Get descriptive feature columns
    descriptive_feature_cols = [
        col for col in df_windows_descriptive.columns
        if (col not in non_feature_cols)
        and (not col.endswith("_signal_type"))
    ]
    
    # Merge binary features from original df_positions into windowed data
    # The position_t in df_windows_descriptive represents the END of the window (current timepoint)
    # so we can merge directly on position_t
    binary_cols_to_merge = [col for col in binary_features_to_group if col in df_positions.columns]
    merge_cols = ["sample_name", "TrackID", "position_t"]
    
    df_binary = df_positions[merge_cols + binary_cols_to_merge].copy()
    df_windows_descriptive = df_windows_descriptive.merge(
        df_binary,
        on=merge_cols,
        how="left",
        suffixes=("", "_orig")
    )
    
    # Drop NaN values
    df_analysis = df_windows_descriptive.dropna(subset=descriptive_feature_cols).copy()
    
    # Reduce redundancy of similar features
    print("Dropping highly correlated features...")
    df_analysis, dropped, report = drop_highly_correlated_features(
        data=df_analysis,
        feature_cols=descriptive_feature_cols,
        threshold=0.95
    )
    print(f"Dropped {len(dropped)} highly correlated features.")
    descriptive_feature_cols = [c for c in descriptive_feature_cols if c not in dropped]
    
    # Remove features with no variance
    selector = VarianceThreshold(threshold=1e-4)
    selector.fit(df_analysis[descriptive_feature_cols])
    
    keep_mask = selector.get_support()
    kept_features = df_analysis[descriptive_feature_cols].columns[keep_mask].tolist()
    dropped_low_var = df_analysis[descriptive_feature_cols].columns[~keep_mask].tolist()
    print(f"Dropped {len(dropped_low_var)} low-variance features.")
    
    # Create groups based on binary features BEFORE subsampling
    # This allows us to subsample each group independently
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
            max_samples_per_track=None,
            random_state=random_state
        )
        
        print(f"  After subsampling: {len(df_group)} samples")
        
        if len(df_group) < 50:
            print(f"  Skipping {group_name}: insufficient data after subsampling ({len(df_group)} samples)")
            group_results[group_name] = None
            continue
        
        adata = cluster_group(
            df_group=df_group,
            feature_cols=kept_features,
            non_feature_cols=non_feature_cols,
            n_neighbors=n_neighbors,
            resolution=resolution,
            min_cluster_size=min_cluster_size,
            pca_var_selection=pca_var_selection,
            random_state=random_state
        )
        
        group_results[group_name] = adata
        
        if adata is not None:
            n_clusters = len(adata.obs["ClusterID"].unique())
            print(f"  Found {n_clusters} clusters")
    
    # Plot results
    print("\nPlotting UMAP grid...")
    fig = plot_group_umaps(group_results, ncols=3)
    
    if outfolder is not None:
        outfolder = Path(outfolder)
        outfolder.mkdir(parents=True, exist_ok=True)
        
        # Save figure
        fig_path = outfolder / "state_classification_umap_grid.pdf"
        fig.savefig(fig_path, bbox_inches='tight', dpi=300)
        print(f"Saved UMAP grid to {fig_path}")
        
        # Save individual group results
        for group_name, adata in group_results.items():
            if adata is not None:
                group_path = outfolder / f"adata_{group_name}.h5ad"
                adata.write(group_path, compression="gzip")
                print(f"Saved {group_name} results to {group_path}")
    
    plt.show()
    
    return group_results, df_windows_descriptive


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
    group_results, df_windows = run_state_classification(
        df_positions=df_positions,
        features=features,
        binary_features_to_group=binary_features_to_group,
        window_size=5,
        chosen_intervals=10,
        n_neighbors=15,
        resolution=0.2,
        min_cluster_size=100,
        pca_var_selection=0.95,
        outfolder=outfolder,
        random_state=seed
    )
    
    # Additional analysis for each group
    for group_name, adata in group_results.items():
        if adata is None:
            continue
            
        print(f"\n=== Analysis for {group_name} ===")
        
        # Plot top ranking features
        plot_top_ranking_features(
            adata,
            groupby="ClusterID",
            n_features=10,
        )
        plt.suptitle(f"{group_name} - Top Ranking Features", y=1.02)
        plt.tight_layout()
        plt.show()
        
        # Plot cluster sizes
        plot_number_per_clusters(adata.obs, cluster_col="ClusterID")
        plt.title(f"{group_name} - Cluster Sizes")
        plt.show()


if __name__ == "__main__":
    pass
    # test()
