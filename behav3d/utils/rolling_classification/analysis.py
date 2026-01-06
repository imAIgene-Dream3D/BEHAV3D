from behav3d.utils import load_behav3d_metadata, check_behav3d_metadata

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

from behav3d.utils.rolling_classification import *
from behav3d.utils.rolling_classification.features import *
from behav3d.utils.rolling_classification.clustering import *
from behav3d.utils.rolling_classification.plotting import *
from behav3d.utils.rolling_classification.videos import *
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

def run_state_based_analysis(
    output_dir,
    cell_type="tcell",
    
    # Input
    adata_full_path=None,  # if None, will look under output_dir/analysis/<cell_type>/rolling_classification/adata_full.h5ad
    state_col="ClusterID",
    groupby_cols=("sample_name", "TrackID"),
    time_col="position_t",

    # Track filtering / truncation
    min_length=100,
    max_length=100,

    # Feature extraction / normalization
    use_trigrams=True,
    do_block_scaling=True,
    do_l2_normalization=False,

    # Feature selection
    drop_highly_correlated=True,
    corr_threshold=0.95,
    drop_low_variance=False,
    low_var_threshold=1e-4,

    # PCA / Leiden / UMAP
    do_pca=False,
    pca_var_selection=0.95,
    n_neighbors=30,
    leiden_resolution=0.2,
    leiden_metric="euclidean",  # or cosine
    leiden_use_rep="X",   # "X" matches your snippet; could also be "X_pca"
    umap_min_dist=0.1,

    # Plotting
    plot_results=True,
    heatmap_figsize=(25, 20),
    matrixplot_figsize=(20, 40),
    umap_size=1,
    umap_alpha=0.5,

    # Exemplar track plotting
    plot_exemplars=True,
    n_per_cluster=10,
    exemplar_state_keys=("ClusterID", "ClusterID_filt5"),

    # Saving
    save_outputs=True,
    output_subdir_name="state_feature_classification",

    seed=123,
):
    print("--------------- Performing state-based full track analysis ---------------")
    start_time = time.time()

    output_dir = Path(output_dir)
    analysis_outdir = Path(output_dir, "analysis", cell_type)
    analysis_outdir.mkdir(parents=True, exist_ok=True)

    rolling_folder = Path(analysis_outdir, "rolling_classification")
    
    outfolder = Path(analysis_outdir, output_subdir_name)
    outfolder.mkdir(parents=True, exist_ok=True)

    if adata_full_path is None:
        adata_full_path = Path(rolling_folder, "adata_full.h5ad")

    adata_full_path = Path(adata_full_path)
    if not adata_full_path.exists():
        raise FileNotFoundError(f"Could not find adata_full.h5ad at: {adata_full_path}")

    adata_full = sc.read_h5ad(adata_full_path)

    # --------- Filter + truncate tracks ----------
    adata_filt = filter_and_truncate_tracks(
        adata_full,
        groupby_cols=list(groupby_cols),
        time_col=time_col,
        min_length=min_length,
        max_length=max_length,
    )

    # --------- Extract state-describing features ----------
    adata_state_features, blocks = extract_descibing_track_state_features(
        adata_filt,
        state_col=state_col,
        use_trigrams=use_trigrams,
    )

    # --------- (Optional) Block scaling + L2 normalize ----------
    if do_block_scaling:
        adata_state_features = scale_feature_blocks(adata_state_features, blocks=blocks)

    if do_l2_normalization:
        adata_state_features = l2_normalize_features_blocks(adata_state_features, blocks=blocks)

    # --------- Feature selection ----------
    dropped_high_corr = []
    dropped_low_var = []

    if drop_highly_correlated:
        adata_state_features, dropped_high_corr, report = drop_highly_correlated_features(
            data=adata_state_features,
            feature_cols=adata_state_features.var_names,
            threshold=corr_threshold,
        )
        print(f"Dropped {len(dropped_high_corr)} highly correlated features.")

    kept_features = list(adata_state_features.var_names)

    if drop_low_variance:
        adata_state_features, kept_features, dropped_low_var = drop_low_variance_features(
            data=adata_state_features,
            feature_cols=kept_features,
            low_var_threshold=low_var_threshold,
        )
        print(f"Dropped {len(dropped_low_var)} low-variance features.")

    leiden_use_rep="X"
    # --------- PCA ----------
    if do_pca:
        adata_state_features = run_pca(
            adata_state_features,
            ncomps=len(adata_state_features.var_names),
            pca_var=pca_var_selection,  # matches your snippet signature
            random_state=seed,
        )
        leiden_use_rep="X_pca",

    # --------- Leiden ----------
    # Your snippet uses use_rep="X" (not X_pca). Keep it configurable.
    leiden_kwargs = dict(
        n_neighbors=n_neighbors,
        resolution=leiden_resolution,
        method="umap",
        use_rep=leiden_use_rep,
        key_added="ClusterID",
        random_state=seed,
    )

    adata_state_features = run_leiden_clustering(
        adata_state_features, 
         n_neighbors=n_neighbors,
            resolution=leiden_resolution,
            method="umap",
            use_rep=leiden_use_rep,
            key_added="ClusterID",
            random_state=seed,
            metric=leiden_metric
        )

    # --------- UMAP ----------
    sc.tl.umap(
        adata_state_features,
        min_dist=umap_min_dist,
        random_state=seed,
    )

    # --------- Plots ----------
    if plot_results:
        sc.pl.umap(adata_state_features, color="ClusterID", size=umap_size, alpha=umap_alpha)

        sc.tl.dendrogram(adata_state_features, groupby="ClusterID")

        sc.pl.heatmap(
            adata_state_features,
            var_names=adata_state_features.var_names,
            groupby="ClusterID",
            standard_scale="var",
            figsize=heatmap_figsize,
            swap_axes=True,
            dendrogram=True,
            show_gene_labels=True,
        )

        sc.pl.matrixplot(
            adata_state_features,
            var_names=adata_state_features.var_names,
            groupby="ClusterID",
            standard_scale="var",
            figsize=matrixplot_figsize,
            swap_axes=True,
            dendrogram=True,
        )

    # --------- Exemplar tracks by cluster ----------
    if plot_exemplars:
        # This assumes plot_exemplar_tracks_by_cluster signature: (adata_tracks, adata_clusters, n_per_cluster, state_key)
        # where `adata_clusters` contains the clustering in .obs["ClusterID"].
        plot_exemplar_tracks_by_cluster(
            adata_filt,
            adata_state_features,
            n_per_cluster=n_per_cluster,
            state_key="ClusterID",
        )
        for k in exemplar_state_keys:
            if k in adata_filt.obs.columns:
                plot_exemplar_tracks_by_cluster(
                    adata_filt,
                    adata_state_features,
                    n_per_cluster=n_per_cluster,
                    state_key=k,
                )

    # --------- Rank features per cluster ----------
    sc.tl.rank_genes_groups(
        adata_state_features,
        groupby="ClusterID",
        method="wilcoxon",
        use_raw=False,
    )
    if plot_results:
        sc.pl.rank_genes_groups(adata_state_features, n_genes=15)

    # --------- Save outputs ----------
    adata_feat_out = Path(outfolder, "adata_state_features_clustered.h5ad")

    adata_state_features.write(adata_feat_out, compression="gzip")

    elapsed = time.time() - start_time
    h = int(elapsed // 3600)
    m = int((elapsed % 3600) // 60)
    s = int(elapsed % 60)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")

    return