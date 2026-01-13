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
from behav3d.utils.rolling_classification.clustering.state.features import *
from behav3d.utils.rolling_classification.clustering import *
from behav3d.utils.rolling_classification.clustering.plotting import *
from behav3d.utils.rolling_classification.clustering.state.videos import *
# %matplotlib inline

import time
seed = 123
random.seed(seed)
np.random.seed(seed)

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