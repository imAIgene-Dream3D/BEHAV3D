"""Deprecated legacy module kept for historical reference.

This file is not part of the actively maintained clustering pipeline.
"""

from pathlib import Path
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

from behav3d.features.state_descriptive_features import (
    extract_descibing_track_state_features, 
    scale_feature_blocks, 
    l2_normalize_features_blocks,
    drop_highly_correlated_features
)

from behav3d.analysis.clustering.state.visualization.plots.state_transitions import compute_cluster_transition_matrix
from behav3d.analysis.clustering.general import relabel_cluster_ids
from behav3d.analysis.filtering import filter_and_truncate_tracks_anndata
from behav3d.analysis.clustering.general.leiden import run_pca, run_leiden_clustering
from behav3d.analysis.clustering.track.visualization.plots.exemplar_track_per_cluster import plot_exemplar_tracks_by_cluster
from behav3d.analysis.clustering.general.visualization.plots import plot_top_ranking_features

import numpy as np
#%matplotlib inline

if __name__ == "__main__":
    seed = 12345
    random.seed(seed)
    np.random.seed(seed)

    ssd_dir = r"/Volumes/T7_Sam/"
    # ssd_dir = r"F:/"
    ssd_dir = Path(ssd_dir)
    outfolder = Path(ssd_dir, "BHVD_BEHAV3D/BEHAV3D_python/rolling_classification")
    adata_full = sc.read_h5ad(Path(outfolder,"adata_full.h5ad"))

    mapping =  {
            "1": "Dead",
            "2": "T cell interaction",
            "3": "Interacting Scanners",
            "4": "Organoid + T cell aggregation",
            "5": "Static",
            "6": "Non-interacting Movers",
            "7": "Organoid interaction",
        }
    adata_full = relabel_cluster_ids(
        adata_full,
        mapping,
        cluster_key="ClusterID",
    )
    compute_cluster_transition_matrix(
            adata=adata_full,
            cluster_key="ClusterID",
            id_cols = ["sample_name", "TrackID"],
            plot=True
        )
    #HMM STATE MAPPING
    mapping =  {
            "1": "Dead",
            "2": "Static",
            "3": "Scanner",
            "4": "Organoid contact",
            "5": "T cell interaction",
            "6": "Organoid + T cell aggregation",
            "7": "Scanner",
        }
    adata_full = relabel_cluster_ids(
        adata_full,
        mapping,
        cluster_key="hmm_state",
    )

    adata_full = filter_short_state_runs(
            adata_full,
            cluster_key="ClusterID",
            id_cols = ["sample_name", "TrackID"],
            time_key = "position_t",
            length_removed = 5,
            new_key = "ClusterID_filt1"
        )

    min_length = 100
    max_length = 100
    adata_filt = filter_and_truncate_tracks_anndata(
        adata_full,
        groupby_cols=["sample_name", "TrackID"],   # <-- set to your actual column
        time_col="position_t",      # <-- set to your actual column (e.g., "t", "time", "frame"); or None
        min_length=min_length,
        max_length=max_length,
    )

    state_col = 'ClusterID'
    state_col = 'ClusterID_filt1'
    test, blocks = extract_descibing_track_state_features(
        adata_filt,
        state_col=state_col,
        # use_trigrams=False
    )

    """
    Add L2 normalization??
    """
    test = scale_feature_blocks(test, blocks=blocks)

    test = l2_normalize_features_blocks(test, blocks=blocks)

    test, dropped, report = drop_highly_correlated_features(
            data=test,
            feature_cols=test.var_names,
            threshold=0.95
        )
    print(f"Dropped {len(dropped)} highly correlated features.")
    kept_features = [c for c in test.var_names if c not in dropped]

    # --- Remove low-variance features, then scale ---
    # test, kept_features, dropped = drop_low_variance_features(
    #     data=test,
    #     feature_cols=kept_features,
    #     low_var_threshold=1e-4
    # )
    # print(f"Dropped {len(dropped)} low-variance features.")
    
    to_run = test
    to_run = run_pca(
        to_run,
        ncomps=len(to_run.var_names),
        pca_var_selection=0.95, 
        random_state=seed
        )


    # to_run = run_leiden_clustering(
    #     to_run, 
    #     n_neighbors=50,
    #     resolution=0.2, 
    #     metric="cosine",
    #     method="umap",
    #     use_rep="X_pca",
    #     key_added="ClusterID",
    #     random_state=seed
    #     )


    to_run = run_leiden_clustering(
        to_run, 
        n_neighbors=30,
        resolution=0.25, 
        # metric="cosine",
        method="umap",
        use_rep="X",
        key_added="ClusterID",
        random_state=seed
        )

    sc.tl.umap(
        to_run,
        min_dist=0.1,
        random_state=seed,
    )

    sc.pl.umap(to_run, color="ClusterID", size=4, alpha=0.5)


    sc.tl.dendrogram(to_run, groupby="ClusterID")
    sc.pl.heatmap(
        to_run,
        var_names=to_run.var_names,
        # var_group_labels=list(feature_dict.keys()),
        # var_group_rotation=0,
        groupby="ClusterID",
        standard_scale="var",
        figsize=(25, 20),
        swap_axes=True,
        dendrogram=True,
        show_gene_labels=True,
    )

    sc.pl.matrixplot(
        to_run,
        var_names=to_run.var_names,    
        groupby="ClusterID",
        standard_scale="var",
        figsize=(20, 40),
        swap_axes=True,
        dendrogram=True
    )

    plot_exemplar_tracks_by_cluster(
        adata_filt,
        to_run,
        n_per_cluster=10,
        state_key="ClusterID",
    )

    plot_exemplar_tracks_by_cluster(
        adata_filt,
        to_run,
        n_per_cluster=10,
        state_key="ClusterID_filt1",
    )

    sc.tl.rank_genes_groups(
        to_run,
        groupby="ClusterID",      # <-- your cluster column
        method="wilcoxon",     # robust choice
        use_raw=False
    )
    sc.pl.rank_genes_groups(to_run, n_genes=15)

    plot_top_ranking_features(
        to_run,
        groupby="ClusterID",
        n_features=10
    )


    """
    Plot sankey plot from one state to another
    """
    df_paths = paths_between_states(
        adata_full,
        start_state="round_static",
        end_state="organoid_contact",
        state_col="ClusterID",
        collapse_bouts=True,
        mode="next_end",
    )

    plot_paths_by_count(
        df_paths,
        top_n=25,
        min_count=1,
        title="Most common paths from state 1 to state 5",
    )

    colors = [
        "#d62728",
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
    ]

    fig = plot_sankey_diagram_between_states(
        df_paths,
        state_colors=colors,
        min_count=100,
    )
    fig.write_image(
        r"C:/Users/Samde/Downloads/newplot.pdf",
        width=1400,
        height=700,
        scale=2,
    )
    fig.show()

    """
    Plot state composition over time
    """
    df_fig, fig, ax = plot_state_composition_over_time(
        adata_full, 
        time_col="position_t", 
        state_col=cluster_column, 
        relative=False
        )

    df_fig, fig, axes= plot_state_composition_over_time(
        adata_full, 
        time_col="position_t", 
        state_col=cluster_column, 
        relative=True,
        group_by_sample=True
        )
