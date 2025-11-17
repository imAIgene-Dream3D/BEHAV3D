from behav3d.preprocessing.segmentation.unet_segmentation import run_behav3d_unet_segmentation

from behav3d.preprocessing.tracking.laptracking import run_tcell_laptracking
from behav3d.utils.tracking import convert_all_tracked_images_to_csv

from behav3d.analysis.feature_extraction import run_feature_extraction
from behav3d.analysis.tcell_analysis import filter_tcell_tracks, run_tcell_analysis
from behav3d.analysis import summarize_track_features
from behav3d.analysis.organoid_analysis import filter_organoid_tracks, run_organoid_analysis
from behav3d.analysis.backprojection import backproject_mean_features_behav3d

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

seed = 123
random.seed(seed)
np.random.seed(seed)

# ssd_dir = r"/Volumes/T7_Sam/"
ssd_dir = r"F:/"
ssd_dir = Path(ssd_dir)
output_dir = Path(ssd_dir, r"BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE")
metadata_csv_path = Path(ssd_dir, r"BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE/metadata.csv")

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
window_size=100
chosen_intervals = 50

# window_size=None
# chosen_intervals = None

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
df_windows_descriptive = create_windowed_track_dataset(
    df_tracks=df_positions,
    columns_to_summarize=features,
    window_size = window_size,
    step_size = chosen_intervals,
    time_col = "position_t",
    id_cols = ["sample_name", "TrackID"],
)

non_feature_cols = [
    "sample_name",
    "TrackID",
    "sub_TrackID",
    "window_start_position_t",
    "window_end_position_t",
    "window_length_frames",
]

# Drop anything that ends with "_signal_type" or other metadata
feature_cols = [
    col for col in df_windows_descriptive.columns
    if (col not in non_feature_cols)
    and (not col.endswith("_signal_type"))
]

### 1) Reduce redundancy of similar features
X_reduced, dropped, report = drop_highly_correlated_features(
    df=df_windows_descriptive,
    feature_cols=feature_cols,
    threshold=0.95
)
print(f"Dropped {len(dropped)} features.")
# report  # see which were kept vs dropped and why

df_analysis = df_windows_descriptive.copy()

### 2) Scale features and run PCA
X_scaled = StandardScaler().fit_transform(X_reduced)
pca_model = PCA(n_components=0.95, random_state=seed)
X_pca = pca_model.fit_transform(X_scaled)
df_analysis["PC1"] = X_pca[:, 0]
df_analysis["PC2"] = X_pca[:, 1]


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
df_analysis["UMAP1"] = umap_embedding[:,0]
df_analysis["UMAP2"] = umap_embedding[:,1]

### 4) Run Leiden clustering
labels_leiden = run_leiden_clustering(
    X=X_pca,                # or umap_embedding, or X_scaled
    n_neighbors=80,
    metric="cosine",
    resolution=0.35,
    random_state=seed,
)
df_analysis["cluster_label_leiden"] = labels_leiden
df_analysis["ClusterID"] = labels_leiden
df_analysis["ClusterID"].value_counts()

###### ANALYSIS VALUES
plot_number_per_clusters(df_analysis)
plot_per_cluster_proportions(df_analysis)

# ---------- 1) PCA & UMAP quick looks ----------
plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
sns.scatterplot(data=df_analysis, x="PC1", y="PC2", s=8, alpha=0.6, edgecolor=None)
plt.title("PCA (PC1 vs PC2)")

plt.subplot(1,2,2)
sns.scatterplot(data=df_analysis, x="UMAP1", y="UMAP2", s=8, alpha=0.6, edgecolor=None)
plt.title("UMAP (2D)")
plt.tight_layout()
plt.show()

plot_feature_cluster_heatmap(
    df_analysis,
    feature_cols=feature_cols,
    cluster_col="ClusterID",
    figsize=(8.27, 11.69),
)




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
    
plot_umap_feature_grid(df_analysis, feature_cols=feature_cols)

plot_feature_cluster_heatmap(
    df_analysis,
    feature_cols=feature_cols,
    cluster_col="ClusterID",
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

plt.figure(figsize=(6,5))
sns.scatterplot(
    data=df_analysis, x="UMAP1", y="UMAP2",
    hue="cluster_label_leiden", palette="tab20",
    s=10, alpha=0.8, edgecolor=None
)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.title("UMAP colored by Leiden clusters (−1 = noise)")


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