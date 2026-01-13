from behav3d.core.metadata import load_behav3d_metadata, check_behav3d_metadata

import pandas as pd
from pathlib import Path

import pandas as pd
import numpy as np

import umap

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import random

from sklearn.cluster import KMeans, HDBSCAN
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from sklearn.decomposition import PCA
from sklearn.feature_selection import VarianceThreshold


from pathlib import Path

import seaborn as sns
# import hdbscan

from behav3d.analysis.classification import *
from behav3d.analysis.classification.clustering.state.features import *
from behav3d.analysis.classification.clustering import *
from behav3d.analysis.classification.clustering.plotting import *
from behav3d.analysis.classification.clustering.state.videos import *
# %matplotlib inline

import time
seed = 123
random.seed(seed)
np.random.seed(seed)

def run_rolling_window_analysis(
    output_dir,
    cell_type="tcell",
    df_tracks_path=None,
    time_col="position_t",
    id_cols=["sample_name", "TrackID"],
    features=[
        "percentage_dead_mask",
        "organoid_contact_pixels",
        "tcell_contact_pixels",
        "speed",
    ],
    # Clustering subsampling
    window_size=10,
    chosen_intervals=10,
    
    # PCA / Leiden / UMAP
    drop_highly_correlated=True,
    drop_low_variance=True,
    low_var_threshold=1e-4,
    pca_var_selection=0.95,
    n_neighbors=50,
    # If leiden_resolution= "auto" or None, use stability-based resolution selection
    leiden_resolution=0.2,
    leiden_stability_resolutions=(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0),
    umap_min_dist=0.1,
    
    # postprocessing
    min_cluster_size=None,
    minimal_state_length=0,
    
    # plotting / saving
    plot_results=True,
    seed=123,
):
    """
    Rolling-window descriptive feature extraction + PCA + Leiden clustering on window subset,
    then ingest cluster labels back to full windowed dataset and apply short-run filtering.

    Returns:
        dict with:
          - adata_full, adata_sub
          - df_analysis (full windowed df with ClusterID added)
          - df_subset (subset df with ClusterID added)
          - kept_features, dropped_low_var
          - paths (written outputs)
    """
    print(f"--------------- Performing {cell_type} rolling-window Leiden analysis ---------------")
    start_time = time.time()

    output_dir = Path(output_dir)
    analysis_outdir = Path(output_dir, "analysis", cell_type)
    analysis_outdir.mkdir(parents=True, exist_ok=True)
    
    feature_outdir = Path(analysis_outdir, "track_features")
    feature_outdir.mkdir(parents=True, exist_ok=True)
    
    outfolder = Path(analysis_outdir, "rolling_classification")
    outfolder.mkdir(parents=True, exist_ok=True)

    if df_tracks_path is None:
        df_tracks_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features_filtered.csv")

    # --- Load & sort positions/tracks ---
    df_positions = pd.read_csv(df_tracks_path, low_memory=False)
    df_positions = df_positions.sort_values(by=[*id_cols, time_col])

    # --- Create rolling-window descriptive features ---
    df_windows_descriptive = create_descriptive_track_dataset(
        df_tracks=df_positions,
        columns_to_summarize=features,
        window_size=window_size,
        time_col=time_col,
        id_cols=id_cols,
    )

    # # Optional: persist & reload (often useful when this is heavy)
    # df_windows_path = Path(outfolder,"df_windows_descriptive.csv")
    # df_windows_descriptive.to_csv(df_windows_path, index=False)

    # --- Prepare feature columns ---
    non_feature_cols = [
        *id_cols,
        "sub_TrackID",
        time_col,
        f"window_start_{time_col}",
        f"window_end_{time_col}",
        "window_length_frames",
    ]
    non_feature_cols = [c for c in non_feature_cols if c in df_windows_descriptive.columns]

    descriptive_feature_cols = [
        c for c in df_windows_descriptive.columns
        if (c not in non_feature_cols) and (not c.endswith("_signal_type"))
    ]

    # --- Convert to AnnData ---
    adata_full = df_to_adata(df_windows_descriptive, descriptive_feature_cols, obs_cols=non_feature_cols)
    
    # Remove NAs (First "window_size" timepoints of each track)
    adata_full = adata_full[~np.isnan(adata_full.X).any(axis=1)].copy()
    
    kept_features = descriptive_feature_cols.copy()
    if drop_highly_correlated:
        # Reduce redundancy of similar features
        adata_full, dropped, report = drop_highly_correlated_features(
            adata_full,
            feature_cols=descriptive_feature_cols,
            threshold=0.95
        )
        print(f"Dropped {len(dropped)} highly correlated features.")
        kept_features = [c for c in descriptive_feature_cols if c not in dropped]
    
    if drop_low_variance:
        # --- Remove low-variance features, then scale ---
        adata_full, kept_features, dropped = drop_low_variance_features(
           adata_full,
            feature_cols=descriptive_feature_cols,
            low_var_threshold=low_var_threshold
        )
        print(f"Dropped {len(dropped)} low-variance features.")

    adata_full.layers["raw"]= adata_full.X.copy()
    
    scaler = StandardScaler().fit(adata_full[:, kept_features].X)
    adata_full.X = scaler.transform(adata_full[:, kept_features].X)

    
    adata_full.uns["preprocessing"] = {
        "kept_features": list(kept_features),
        "dropped_low_variance": drop_low_variance,
        "dropped_highly_correlated": drop_highly_correlated,
        "variance_threshold": float(low_var_threshold) if drop_low_variance else None,
        "scaler": {
            "mean": scaler.mean_.astype(float),
            "scale": scaler.scale_.astype(float),
        },
        "windowing": {
            "window_size": window_size,
        },
    }

    # --- Subset windows for clustering ---
    df_analysis = df_windows_descriptive.copy()
    df_subset = subset_windowed_tracks(
        df_windowed=df_analysis,
        step_size=chosen_intervals,
    )
    adata_sub = df_to_adata(df_subset, kept_features, obs_cols=non_feature_cols)
    adata_sub.uns["preprocessing"] = adata_full.uns["preprocessing"]

    # --- PCA + Leiden on subset ---
    adata_sub = run_pca(
        adata_sub,
        pca_var_selection=pca_var_selection,
        ncomps=len(kept_features),
        svd_solver="full",
        random_state=seed,
    )

    if leiden_resolution is None or leiden_resolution == "auto":
        adata_sub = run_leiden_clustering(
            adata_sub,
            resolution="auto",
            stability_resolutions=leiden_stability_resolutions,
            key_added="ClusterID",
            random_state=seed,
        )
    else:
        adata_sub = run_leiden_clustering(
            adata_sub,
            n_neighbors=n_neighbors,
            resolution=leiden_resolution,
            metric="euclidean",
            method="umap",
            use_rep="X_pca",
            key_added="ClusterID",
            random_state=seed,
        )

    if min_cluster_size is not None:
        adata_sub = merge_small_clusters(
            adata_sub,
            key="ClusterID",
            min_size=min_cluster_size,
            use_rep="X_pca",
        )

    # --- UMAP embedding on subset for viz ---
    sc.tl.umap(adata_sub, min_dist=umap_min_dist, random_state=seed)
    sc.tl.rank_genes_groups(
        adata_sub,
        groupby="ClusterID",      # <-- your cluster column
        method="wilcoxon",     # robust choice
        use_raw=False
        )
    
    
    if plot_results:
        sc.pl.umap(adata_sub, color="ClusterID", size=2, alpha=0.5)
        sc.pl.rank_genes_groups(adata_sub, n_genes=15, swap_axes=True, rotation=90)
        # Group features by prefix for heatmap/matrixplot
        feature_dict = {f: [c for c in kept_features if c.startswith(f + "_")] for f in features}
        known_prefixes = tuple(f + "_" for f in features)
        feature_dict["other"] = [c for c in kept_features if not c.startswith(known_prefixes)]

        sc.tl.dendrogram(adata_sub, groupby="ClusterID")
        sc.pl.heatmap(
            adata_sub,
            var_names=feature_dict,
            var_group_labels=list(feature_dict.keys()),
            var_group_rotation=0,
            groupby="ClusterID",
            standard_scale="var",
            figsize=(25, 20),
            swap_axes=True,
            dendrogram=True,
            show_gene_labels=True,
        )
        sc.pl.matrixplot(
            adata_sub,
            var_names=feature_dict,
            groupby="ClusterID",
            standard_scale="var",
            figsize=(20, 20),
            swap_axes=True,
            dendrogram=True,
        )

    # --- Ingest labels back to full dataset ---
    sc.tl.ingest(
        adata_full,
        adata_sub,
        obs="ClusterID",
        embedding_method="umap",
    )

    # --- Filter short state runs (optional) ---
    if minimal_state_length > 0:
        adata_full.obs["ClusterID_unfiltered"] = adata_full.obs["ClusterID"].copy()
        adata_full = filter_short_state_runs(
            adata_full,
            cluster_key="ClusterID",
            id_cols=list(id_cols),
            time_key=time_col,
            length_removed=minimal_state_length,
            new_key=f"ClusterID",
        )

    # --- Save outputs ---
    adata_sub_outpath = Path(outfolder, "adata_sub.h5ad")
    adata_full_outpath = Path(outfolder, "adata_full.h5ad")
    
    adata_full.write(adata_full_outpath, compression="gzip")
    adata_sub.write(adata_sub_outpath, compression="gzip")

    # --- Downstream analysis plots / transitions ---
    if plot_results:
        compute_cluster_transition_matrix(
            adata=adata_full,
            cluster_key="ClusterID",
            id_cols=list(id_cols),
            plot=True,
        )
        # also run these variants if you want them like in your notebook:
        compute_cluster_transition_matrix(
            adata=adata_full,
            cluster_key="ClusterID",
            id_cols=list(id_cols),
            only_transitions=False,
            plot=True,
        )
        compute_cluster_transition_matrix(
            adata=adata_full,
            cluster_key="ClusterID",
            id_cols=list(id_cols),
            only_transitions=True,
            plot=True,
        )

        plot_exemplar_track_bars(
            adata_full,
            n_tracks=50,
            track_key="TrackID",
            time_key=time_col,
            state_key="ClusterID",
            seed=seed,
            cmap_name="tab20",
            ax=None,
        )

    end_time = time.time()
    elapsed = end_time - start_time
    h = int(elapsed // 3600)
    m = int((elapsed % 3600) // 60)
    s = int(elapsed % 60)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")

    return
