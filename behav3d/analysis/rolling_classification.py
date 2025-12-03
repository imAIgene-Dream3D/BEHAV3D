from behav3d.utils import load_behav3d_metadata, check_behav3d_metadata

import pandas as pd
from pathlib import Path

import torch
import zarr
import dask.array as da
from tqdm import tqdm

import argparse
from dtaidistance import dtw, dtw_ndim
import pandas as pd
import numpy as np

import umap

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
import random

from sklearn.manifold import TSNE
from sklearn.cluster import KMeans, HDBSCAN
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

from scipy import sparse
import igraph as ig
import leidenalg as la

from pathlib import Path
from behav3d.utils import format_time
from behav3d.utils.filtering import plot_filter_count

import yaml
import time
import seaborn as sns
# import hdbscan

from pandas.api.types import is_numeric_dtype
import math

from behav3d.utils.rolling_classification import *
from behav3d.utils.rolling_classification.plotting import *
%matplotlib inline

seed = 123
random.seed(seed)
np.random.seed(seed)

if __name__ == "__main__":
    # ssd_dir = r"/Volumes/T7_Sam/"
    ssd_dir = r"F:/"
    ssd_dir = Path(ssd_dir)
    
    output_dir = Path(ssd_dir, r"BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE")
    metadata_csv_path = Path(ssd_dir, r"BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE/metadata.csv")

    output_dir = Path(ssd_dir, r"BHVD_BEHAV3D/BEHAV3D_python/runs/CombinedAnalysis_AmberMacrophage")
    metadata_csv_path = Path(ssd_dir, r"BHVD_BEHAV3D/BEHAV3D_python/runs/CombinedAnalysis_AmberMacrophage/metadata.csv")


    downloads_folder = r"C:/Users/Samde/Downloads"
    # downloads_folder = "/Users/s.deblank-3/Downloads"
    # output_dir = r"F:/BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE"
    # metadata_csv_path = r"F:/BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE/metadata.csv"

    metadata = load_behav3d_metadata(metadata_csv_path)


    # from tsfresh import extract_features

    analysis_outdir = Path(output_dir, "analysis", "tcell")
    feature_outdir = Path(analysis_outdir, "track_features")
    df_tracks_path = Path(feature_outdir, f"BEHAV3D_tcell_combined_track_features_filtered.csv")
    # df_tracks_path = Path(feature_outdir, f"BEHAV3D_tcell_combined_track_features.csv")
    df_positions = pd.read_csv(df_tracks_path)
    df_positions=df_positions.sort_values(by=["sample_name", "TrackID", "position_t"])

    # --- Create descriptive features per value ---
    # window_size=100
    # chosen_intervals = 50

    window_size = 25
    # window_size = None

    groupby=["sample_name", "TrackID"]
    features=[
        "percentage_dead_mask",
        "nr_dead_mask_pixels",
        "organoid_contact_pixels",
        "tcell_contact_pixels",
        "mean_square_displacement",
        "speed",
    ]

    df_tracks = df_positions[features+["sample_name", "TrackID", "position_t"]]
    df_windows_descriptive = create_descriptive_track_dataset(
        df_tracks=df_positions,
        columns_to_summarize=features,
        window_size = window_size,
        step_size = 1,
        time_col = "position_t",
        id_cols = ["sample_name", "TrackID"],
    )

    # df_windows_descriptive.to_csv(Path(downloads_folder,r"df_windows_descriptive.csv"), index=False)
    df_windows_descriptive = pd.read_csv(Path(downloads_folder,r"df_windows_descriptive.csv"))
    
    ### Sample the dataset (either to fraction of the tracks)
    non_feature_cols = [
        "sample_name",
        "TrackID",
        "sub_TrackID",
        "position_t",
        "window_start_position_t",
        "window_end_position_t",
        "window_length_frames",
    ]
    
    
    df_analysis = df_windows_descriptive.copy()
    
    # Drop anything that ends with "_signal_type" or other metadata
    descriptive_feature_cols = [
        col for col in df_windows_descriptive.columns
        if (col not in non_feature_cols)
        and (not col.endswith("_signal_type"))
    ]
    df_analysis = df_analysis.dropna(subset=descriptive_feature_cols)
    
    ### Reduce redundancy of similar features
    df_analysis, dropped, report = drop_highly_correlated_features(
        df=df_analysis,
        feature_cols=descriptive_feature_cols,
        threshold=0.95
    )
    print(f"Dropped {len(dropped)} features.")
    descriptive_feature_cols = list(set(descriptive_feature_cols) - set(dropped))
    
    
    # --- 3) Subset windows (custom) ---
    chosen_intervals = 25
    df_analysis_subset = subset_windowed_tracks(
        df_windowed=df_analysis,
        step_size=chosen_intervals,
    )
    
    adata_full = df_to_adata(df_analysis, descriptive_feature_cols, obs_cols=non_feature_cols)
    adata_sub  = df_to_adata(df_analysis_subset, descriptive_feature_cols, obs_cols=non_feature_cols)
    
    sc.pp.scale(adata_full, zero_center=True, max_value=None) #, max_value=10)
    adata_sub.X = adata_full.X[df_analysis_subset.index.values, :]
    
    pca_var=0.95
    sc.pp.pca(adata_sub)
    var_ratio = adata_sub.uns["pca"]["variance_ratio"]
    n_pcs = int(np.searchsorted(np.cumsum(var_ratio), pca_var) + 1)
    adata_sub.obsm["X_pca"] = adata_sub.obsm["X_pca"][:, :n_pcs]
    
    sc.pp.neighbors(
        adata_sub,
        n_neighbors=100,
        metric="cosine",
        method="umap",
        use_rep="X_pca",
    )
    
    sc.tl.leiden(
        adata_sub,
        resolution=0.35,
        random_state=seed,
        key_added="ClusterID",
    )
    
    sc.tl.umap(
        adata_sub,
        min_dist=0.1,
        random_state=seed,
    )
    
    sc.tl.ingest(
        adata_full,
        adata_sub,
        obs="ClusterID",
        embedding_method="umap"  # also transfers UMAP coords into adata_full.obsm["X_umap"]
    )
    
    # write results back
    df_analysis = adata_add_back_to_df(
        df_analysis, adata_full,
        cols_from_obs=["ClusterID"],
    )
    
    df_analysis_subset = adata_add_back_to_df(
        df_analysis_subset, adata_sub,
        cols_from_obs=["ClusterID"],
    )
    
    ###### ANALYSIS VALUES
    plot_number_per_clusters(df_analysis, cluster_col="ClusterID")
    plot_per_cluster_proportions(df_analysis)
    
    sc.pl.umap(adata_sub, color="ClusterID")
    
    plt.figure(figsize=(6,5))
    ax = sns.scatterplot(
    data=df_analysis, x="UMAP1", y="UMAP2",
    hue="cluster_label_leiden", palette="tab20",
    s=2, alpha=0.3, edgecolor=None  
    )
    # bigger legend dots without changing plot dots
    leg = ax.legend(
        bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.,
        markerscale=6,  # <-- enlarge legend markers
        scatterpoints=1 # keep one dot per legend entry
    )
    plt.title("UMAP colored by Leiden clusters (−1 = noise)")
    plt.tight_layout()
    
    
    #######################################
    ####### HMM STATE CLASSIFICATION ######
    #######################################
    
    # Fraction of tracks
    df_analysis_subset, sampled_keys = subset_full_tracks(
        df=df_analysis,
        fraction=0.3,
        random_state=seed,
        id_cols=["sample_name", "TrackID"],
        return_selected_keys=True
    )
    
    df_analysis_subset, hmm_model, selection_df = run_hmm_state_classification(
        df_features=df_analysis_subset,
        feature_cols=descriptive_feature_cols,
        n_states="auto"
    )
    
    # Apply HMM to full dataset
    df_analysis, _, _ = run_hmm_state_classification(
        df_features=df_analysis,
        feature_cols=descriptive_feature_cols,
        model=hmm_model
    )
    
    
    
    
    
    
    
    
    
    
    

    # ---------- 1) PCA & UMAP quick looks ----------
    plt.figure(figsize=(12,5))
    plt.subplot(1,2,1)
    sns.scatterplot(data=df_analysis, x="PC1", y="PC2", s=8, alpha=0.6, edgecolor=None)
    plt.title("PCA (PC1 vs PC2)")
    plt.show()
    
    plt.subplot(1,2,2)
    sns.scatterplot(data=df_analysis, x="UMAP1", y="UMAP2", s=8, alpha=0.6, edgecolor=None)
    plt.title("UMAP (2D)")
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(6,5))
    ax = sns.scatterplot(
    data=df_analysis, x="UMAP1", y="UMAP2",
    hue="cluster_label_leiden", palette="tab20",
    s=2, alpha=0.3, edgecolor=None  
    )
    # bigger legend dots without changing plot dots
    leg = ax.legend(
        bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.,
        markerscale=6,  # <-- enlarge legend markers
        scatterpoints=1 # keep one dot per legend entry
    )
    plt.title("UMAP colored by Leiden clusters (−1 = noise)")
    plt.tight_layout()
    
    plot_feature_cluster_heatmap(
        df_analysis,
        feature_cols=descriptive_feature_cols,
        cluster_col="hmm_state",
        figsize=(8.27, 11.69),
    )

    create_cluster_videos(
        df_analysis,
        df_positions,
        output_folder= output_dir,
        out_dir = downloads_folder,
        # normalize_per_channel: bool = False,
        fps = 6,
        # dpi: int = 200,
        margin = (20, 20, 20),
        pmin = 0.0,
        pmax = 99.99,
        examples_per_cluster = 6,
        # seed: int = 0,
        # figsize_per_row=(12.0, 4.0),
        # traj_pad_frac: float = 0.05,
    )
    
    create_cluster_overview_video(
        df_analysis,
        df_positions,
        output_folder=output_dir,
        out_dir=downloads_folder,
        examples_per_cluster=3,   # 1 example per cluster per row
        fps=6,
        margin = (20, 20, 20),
        pmin = 0.0,
        pmax = 99.99,
        normalize_per_channel=True,
        seed=1234
        # figsize_per_example=(6.0, 3.0),  # smaller per example if many clusters
    )

    create_fulltrack_cluster_videos(
        examples_per_cluster=3,
        clusters=[0],  
        df_windows=df_analysis,
        df_positions=df_positions,
        output_folder=output_dir,  # folder containing images/<sample>/<sample>.zarr
        out_dir=downloads_folder,
        # clusters=None,                # or e.g. [0, 1, 2]
        fps=6,
        margin=20,
        track_color="#63ff33",
        pmin=0.0,
        pmax=99,
        seed=1234,
        normalize_per_channel=False,
        mask_margin=False,
    )
    
    ##################################################
    ##### OLD METHOD WITHOUT SCANPY INGESTION ########
    ###################################################
    df_analysis[descriptive_feature_cols] = StandardScaler().fit_transform(df_analysis[descriptive_feature_cols])
    # Windows
    # chosen_intervals = 25
    # df_analysis_subset = subset_windowed_tracks(
    #     df_windowed=df_analysis,
    #     step_size=chosen_intervals,
    # )
    
   
    ### 2) Scale features and run PCA
    # X_scaled = StandardScaler().fit_transform(X_reduced)
    pca_model = PCA(n_components=0.95, random_state=seed)
    X_pca = pca_model.fit_transform(df_analysis_subset[descriptive_feature_cols])
    df_analysis_subset["PC1"] = X_pca[:, 0]
    df_analysis_subset["PC2"] = X_pca[:, 1]


    ### 3) Embed the data using UMAP
    umap_model = umap.UMAP(
                n_components=2, 
                n_neighbors=100, 
                min_dist=0.1, 
                # init="random", 
                metric = "cosine",
                random_state=seed,
                )
    umap_embedding = umap_model.fit_transform(X_pca)
    df_analysis_subset["UMAP1"] = umap_embedding[:,0]
    df_analysis_subset["UMAP2"] = umap_embedding[:,1]

    ### 4) Run Leiden clustering
    labels_leiden = run_leiden_clustering(
        X=X_pca,                # or umap_embedding, or X_scaled
        n_neighbors=80,
        metric="cosine",
        resolution=0.35,
        random_state=seed,
    )
    df_analysis_subset["cluster_label_leiden"] = labels_leiden
    df_analysis_subset["ClusterID"] = labels_leiden
    df_analysis_subset["ClusterID"].value_counts()
    
    
    ############
    # HDBSCAN
    ############
    hdbscan_clusterer=HDBSCAN(
        min_cluster_size=100,
        min_samples=50,
        metric="euclidean",
        alpha=1.0,
        cluster_selection_method="eom",
        cluster_selection_epsilon=0.0,
        allow_single_cluster=False,
        leaf_size=40,
        algorithm="auto",             
        n_jobs=None,
        copy=False
        )
    cluster_labels = hdbscan_clusterer.fit_predict(umap_embedding)
    df_analysis["cluster_label_hdbscan"] = cluster_labels
    # df_analysis["ClusterID"] = cluster_labels

    n_clusters = 8  # you can tune this
    kmeans_clusterer = KMeans(
        n_clusters=n_clusters, 
        random_state=seed, 
        n_init="auto"
        )
    cluster_labels = kmeans_clusterer.fit_predict(umap_embedding)
    df_analysis["cluster_label_kmeans"] = cluster_labels





    df_analysis["ClusterID"].value_counts()




    plot_all_clusters_window_max_projection(
        df_positions,
        df_analysis,
        max_windows=5,
        cluster_col="cluster_label_hdbscan",
    )
        
    plot_umap_feature_grid(df_analysis, feature_cols=descriptive_feature_cols)

    plot_feature_cluster_heatmap(
        df_analysis,
        feature_cols=descriptive_feature_cols,
        cluster_col="hmm_state",
        figsize=(8.27, 11.69),
    )

    # plot_clustering_feature_heatmap(
    #     df_umap=df_analysis,
    #     info_cols=feature_cols,
    #     sample_cols=non_feature_cols,
    #     outpath=r"/Users/s.deblank-3/Downloads/test.pdf",
    #     rows_per_page = 7,
    #     nr_cols = 2,
    #     figsize = (8.27, 11.69),
    #     plot_results=True,
    #     show_points=False,        # overlay individual samples
    #     point_alpha=0.5,         # transparency for individual samples
    #     point_size=8,            # size for individual points
    #     mean_marker_size=60,     # size for mean markers
    # )
    compare_cluster_distribution(df_analysis, "cluster_label_hdbscan", "cluster_label_kmeans")

    # ---------- 3) UMAP colored by clusters ----------
    plt.figure(figsize=(6,5))
    sns.scatterplot(
        data=df_analysis, x="UMAP1", y="UMAP2",
        hue="cluster_label_kmeans", palette="tab20",
        s=10, alpha=0.8, edgecolor=None
    )
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.title("UMAP colored by KMeans clusters")
    # plt.savefig(Path(r"~/Downloads/umapKmeans.pdf"), bbox_inches="tight")
    plt.show()

    plt.figure(figsize=(6,5))
    sns.scatterplot(
        data=df_analysis, x="UMAP1", y="UMAP2",
        hue="cluster_label_hdbscan", palette="tab20",
        s=10, alpha=0.8, edgecolor=None
    )
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.title("UMAP colored by HDBSCAN clusters (−1 = noise)")
    plt.show()

   


    plt.savefig(r"/Users/s.deblank-3/Downloads/umapLeiden.pdf", bbox_inches="tight")
    plt.show()

    # plt.figure(figsize=(6,5))
    # sns.scatterplot(
    #     data=df_analysis, x="UMAP1", y="UMAP2",
    #     hue="cluster_label_dtw_hdbscan", palette="tab20",
    #     s=10, alpha=0.8, edgecolor=None, legend=False
    # )
    # plt.title("UMAP colored by HDBSCAN clusters (−1 = noise)")
    # plt.show()

    # plt.figure(figsize=(6,5))
    # sns.scatterplot(
    #     data=df_analysis, x="UMAP1", y="UMAP2",
    #     hue="cluster_label_dtw_kmeans", palette="tab20",
    #     s=10, alpha=0.8, edgecolor=None, legend=False
    # )
    # plt.title("UMAP colored by Kmeans clusters (−1 = noise)")
    # plt.show()

    # plt.figure(figsize=(6,5))
    # sns.scatterplot(
    #     data=df_analysis, x="DTW_UMAP1", y="DTW_UMAP2",
    #     hue="cluster_label_dtw_hdbscan", palette="tab20",
    #     s=10, alpha=0.8, edgecolor=None, legend=False
    # )
    # plt.title("UMAP colored by HDBSCAN clusters (−1 = noise)")
    # plt.show()

    # plt.figure(figsize=(6,5))
    # sns.scatterplot(
    #     data=df_analysis, x="DTW_UMAP1", y="DTW_UMAP2",
    #     hue="cluster_label_dtw_kmeans", palette="tab20",
    #     s=10, alpha=0.8, edgecolor=None, legend=False
    # )
    # plt.title("UMAP colored by Kmeans clusters (−1 = noise)")
    # plt.show()
    
    
    
    if "cluster_label_hdbscan" in df_analysis.columns:
        # HDBSCAN uses -1 for noise; use a palette that includes it
        plt.figure(figsize=(6,5))
        sns.scatterplot(
            data=df_analysis, x="PC1", y="PC2",
            hue="cluster_label_hdbscan", palette="tab20",
            s=10, alpha=0.8, edgecolor=None, legend=False
        )
        plt.title("PCA colored by HDBSCAN clusters (−1 = noise)")
        plt.show()


    # ---------- 4) Optional: PCA explained variance curve ----------
    if hasattr(pca_model, "explained_variance_ratio_"):
        plt.figure(figsize=(6,4))
        evr = pca_model.explained_variance_ratio_
        plt.plot(np.arange(1, len(evr)+1), np.cumsum(evr)*100, marker="o")
        plt.xlabel("Number of PCA components")
        plt.ylabel("Cumulative variance explained (%)")
        plt.title("PCA explained variance")
        plt.grid(alpha=0.2)
        plt.show()

    # ---------- 5) Optional: cluster size bars ----------
    if "cluster_label_kmeans" in df_analysis.columns:
        plt.figure(figsize=(6,3))
        k_counts = df_analysis["cluster_label_kmeans"].value_counts().sort_index()
        sns.barplot(x=k_counts.index, y=k_counts.values, color="tab:blue")
        plt.title("KMeans cluster sizes")
        plt.xlabel("Cluster"); plt.ylabel("Count")
        plt.show()

    if "cluster_label_hdbscan" in df_analysis.columns:
        plt.figure(figsize=(6,3))
        h_counts = df_analysis["cluster_label_hdbscan"].value_counts().sort_index()
        sns.barplot(x=h_counts.index.astype(str), y=h_counts.values, color="tab:green")
        plt.title("HDBSCAN cluster sizes (−1 = noise)")
        plt.xlabel("Cluster"); plt.ylabel("Count")
        plt.show()
        


        
    
    #############
    ### DTW
    #############
    # dtw_features=[
    #     # "elongation",
    #     # "sphericity",
    #     "percentage_dead_mask",
    #     # "nr_dead_mask_pixels",
    #     "organoid_contact_pixels",
    #     "tcell_contact_pixels",
    #     # "displacement",
    #     "mean_square_displacement",
    #     "speed",
    #     # # "dead",
    #     # "active_tcell_contact",
    #     # "position_t"
    # ]
    # non_binary = [c for c in dtw_features if "contact" not in c]
    # dtw_result = compute_dtw_window_clusters(
    #     df_tracks=df_tracks,                     # from your code
    #     df_windows=df_windows_descriptive,       
    #     features=dtw_features,
    #     non_binary_features=non_binary,
    #     umap_n_neighbors=50,
    #     umap_min_dist=0.1,
    #     random_state=seed,
    #     sample_frac=None,        # set e.g. 0.2 if it’s too big
    #     max_windows=None,         # or set e.g. 4000 to cap
    #     clusterer=hdbscan_clusterer,
    #     out_col_name="cluster_label_dtw_hdbscan"
    # )
    # join_keys = [k for k in ["sample_name","TrackID","sub_TrackID","window_start_position_t","window_end_position_t"]
    #              if k in dtw_result.columns and k in df_analysis.columns]
    # df_analysis = df_analysis.merge(
    #     dtw_result[join_keys + ["DTW_UMAP1","DTW_UMAP2","cluster_label_dtw_hdbscan"]],
    #     on=join_keys,
    #     how="left"
    # )