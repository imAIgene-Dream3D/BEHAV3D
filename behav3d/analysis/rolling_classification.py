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
import tqdm

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

from pathlib import Path
from behav3d.utils import format_time
from behav3d.utils.filtering import plot_filter_count

import yaml
import time
import seaborn as sns
import hdbscan

output_dir = r"/Volumes/T7_sam/BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE"
metadata_csv_path = r"/Volumes/T7_sam/BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE/metadata.csv"
metadata = load_behav3d_metadata(metadata_csv_path)


# from tsfresh import extract_features

analysis_outdir = Path(output_dir, "analysis", "tcell")
feature_outdir = Path(analysis_outdir, "track_features")
df_tracks_path = Path(feature_outdir, f"BEHAV3D_tcell_combined_track_features_filtered.csv")
# df_tracks_path = Path(feature_outdir, f"BEHAV3D_tcell_combined_track_features.csv")
df_tracks_orig = pd.read_csv(df_tracks_path)
df_tracks_orig=df_tracks_orig.sort_values(by=["sample_name", "TrackID", "position_t"])

# --- Create descriptive features per value ---
window_size=100
chosen_intervals = 50

window_size=None
chosen_intervals = None


groupby=["sample_name", "TrackID"]
features=[
    # "elongation",
    # "sphericity",
    "percentage_dead_mask",
    "nr_dead_mask_pixels",
    "organoid_contact_pixels",
    "tcell_contact_pixels",
    "displacement",
    "mean_square_displacement",
    "speed",
    # # "dead",
    # "active_tcell_contact",
    # "position_t"
]


df_tracks = df_tracks_orig[features+["sample_name", "TrackID", "position_t"]]
df_windows_descriptive = create_windowed_track_dataset(
    df_tracks=df_tracks,
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

# 1) Reduce redundancy
X_reduced, dropped, report = drop_highly_correlated_features(
    df=df_windows_descriptive,
    feature_cols=feature_cols,
    threshold=0.95
)
print(f"Dropped {len(dropped)} features.")
report  # see which were kept vs dropped and why

df_analysis = df_windows_descriptive.copy()
# --- PCA ---
X_scaled = StandardScaler().fit_transform(X_reduced)

pca_model = PCA(n_components=0.95, random_state=123)
X_pca = pca_model.fit_transform(X_scaled)

# df_pca = pd.DataFrame(X_pca[:, :2], columns=["PC1", "PC2"], index=df_windows_descriptive.index)
# df_analysis = pd.concat([df_windows_descriptive, df_pca], axis=1)
df_analysis["PC1"] = X_pca[:, 0]
df_analysis["PC2"] = X_pca[:, 1]
# --- Clustering ---
# Option A: KMeans
n_clusters = 4  # you can tune this
kmeans = KMeans(n_clusters=n_clusters, random_state=123, n_init="auto")
cluster_labels = kmeans.fit_predict(X_pca)
df_analysis["cluster_label_kmeans"] = cluster_labels


# Option B: (alternative) HDBSCAN (for variable density clusters)
clusterer = hdbscan.HDBSCAN(min_cluster_size=50, metric='euclidean')
cluster_labels = clusterer.fit_predict(X_pca)
df_analysis["cluster_label_hdbscan"] = cluster_labels

# --- UMAP ---
umap_model = umap.UMAP(
            n_components=2, 
            n_neighbors=100, 
            min_dist=0.1, 
            # init="random", 
            metric = "cosine",
            random_state=123,
            )

umap_embedding = umap_model.fit_transform(X_pca)
df_analysis["UMAP1"] = umap_embedding[:,0]
df_analysis["UMAP2"] = umap_embedding[:,1]
    
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

# ---------- 2) PCA colored by clusters ----------
if "cluster_label_kmeans" in df_analysis.columns:
    plt.figure(figsize=(6,5))
    sns.scatterplot(
        data=df_analysis, x="PC1", y="PC2",
        hue="cluster_label_kmeans", palette="tab20",
        s=10, alpha=0.8, edgecolor=None, legend=False
    )
    plt.title("PCA colored by KMeans clusters")
    plt.show()

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

# ---------- 3) UMAP colored by clusters ----------
if "cluster_label_kmeans" in df_analysis.columns:
    plt.figure(figsize=(6,5))
    sns.scatterplot(
        data=df_analysis, x="UMAP1", y="UMAP2",
        hue="cluster_label_kmeans", palette="tab20",
        s=10, alpha=0.8, edgecolor=None, legend=False
    )
    plt.title("UMAP colored by KMeans clusters")
    plt.show()

if "cluster_label_hdbscan" in df_analysis.columns:
    plt.figure(figsize=(6,5))
    sns.scatterplot(
        data=df_analysis, x="UMAP1", y="UMAP2",
        hue="cluster_label_hdbscan", palette="tab20",
        s=10, alpha=0.8, edgecolor=None, legend=False
    )
    plt.title("UMAP colored by HDBSCAN clusters (−1 = noise)")
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
    
# df_umap = pd.merge(df_tracks, umap_embedding, how="left")
# df_umap = pd.merge(df_rolling, umap_embedding, how="left")
for feature in feature_cols:
    sns.scatterplot(
        data=df_analysis,
        x="UMAP1",
        y="UMAP2",
        hue=feature,
        s=5,
        alpha=0.5,
        palette="viridis",
        # legend=False
    )
    plt.show()
        
    
from pandas.api.types import is_numeric_dtype
def plot_umap_feature_grid(
    df: pd.DataFrame,
    feature_cols: list[str],
    x_col: str = "UMAP1",
    y_col: str = "UMAP2",
    ncols: int = 4,
    max_plots: int | None = None,
    point_size: int = 5,
    alpha: float = 0.5,
    numeric_cmap: str = "viridis",
    categorical_palette: str = "tab20",
    add_colorbar: bool = True,
    page: int = 0,   # for pagination: 0-based page index
    ):
    """
    Creates a multi-row, multi-column grid of UMAP scatterplots colored by each feature in feature_cols.
    Filters out non-scalar or missing features automatically. Supports pagination via `page`.
    """
    # Filter valid features (exist, scalar, not all NaN)
    valid = []
    for c in feature_cols:
        if c in df.columns and _is_scalar_series(df[c]) and df[c].notna().any():
            valid.append(c)
    if max_plots is not None:
        valid = valid[:max_plots]

    if len(valid) == 0:
        raise ValueError("No valid features to plot.")

    n = len(valid)
    nrows = math.ceil(n / ncols)

    # Pagination support: choose a slice of features per page
    per_page = nrows * ncols
    start = page * per_page
    end = min(start + per_page, len(valid))
    feats = valid[start:end]
    if len(feats) == 0:
        raise ValueError(f"No features to plot on page {page} (only {math.ceil(len(valid)/per_page)} page(s) available).")

    # Axes limits shared across panels
    x_min, x_max = df[x_col].min(), df[x_col].max()
    y_min, y_max = df[y_col].min(), df[y_col].max()

    # Build grid
    fig, axes = plt.subplots(
        nrows=math.ceil(len(feats)/ncols),
        ncols=ncols,
        figsize=(4*ncols, 3.5*math.ceil(len(feats)/ncols)),
        squeeze=False,
        constrained_layout=True
    )

    for i, feat in enumerate(feats):
        r, c = divmod(i, ncols)
        ax = axes[r, c]

        s = df[feat]
        # Numeric vs categorical handling
        if is_numeric_dtype(s):
            # Numeric: use matplotlib scatter for easy colorbar handling
            sc = ax.scatter(
                df[x_col], df[y_col],
                s=point_size, alpha=alpha,
                c=s, cmap=numeric_cmap, edgecolors="none"
            )
            if add_colorbar:
                cb = plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
                cb.ax.tick_params(labelsize=8)
        else:
            # Categorical: enforce category dtype and use seaborn palette
            s_cat = s.astype("category")
            tmp = df.copy()
            tmp[feat] = s_cat
            sns.scatterplot(
                data=tmp, x=x_col, y=y_col, hue=feat,
                palette=categorical_palette, s=point_size, alpha=alpha,
                legend=False, ax=ax
            )

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_title(feat, fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("")

    # Hide any unused axes (if grid not full)
    total_cells = axes.size
    for j in range(len(feats), total_cells):
        r, c = divmod(j, ncols)
        axes[r, c].axis("off")

    # Add a common title
    fig.suptitle("UMAP colored by features", fontsize=14)
    plt.show()     
      
      
        
def dynamic_time_warping():
    """
    DYNAMIC TIME WARPING + UMAP
    """
    ### Z-scale certain features
    scaler = StandardScaler()
    # scaler = RobustScaler()
    ## select all items in list not *contact*

    non_binary_features = [x for x in features if "contact" not in x]

    df_tracks[non_binary_features]= pd.DataFrame(
        scaler.fit_transform(df_tracks[non_binary_features]), 
        columns=df_tracks[non_binary_features].columns
        )


    chosen_intervals = 75
    dtw_input_tracks=[]
    dtw_rownames=[]
    unique_tracks = df_tracks.groupby(['TrackID', 'sample_name'])
    for (TrackID, sample_name), group in unique_tracks:
        track_features = group[features].to_numpy().astype(np.double)
        
        n_frames = track_features.shape[0]

        if n_frames < window_size:
            continue  # skip if the track is too short

        for start in range(0, n_frames - window_size + 1, chosen_intervals):
            window = track_features[start:start + window_size]
            dtw_input_tracks.append(window)
            position_t = int(group["position_t"].to_numpy()[start + window_size -1])
            dtw_rownames.append(f"{TrackID}--{sample_name}--{position_t}")

    ### Select random samples from data
    # sampled_indices = random.sample(range(len(dtw_input_tracks)), 20000)

    # dtw_input_tracks = [dtw_input_tracks[i] for i in sampled_indices]
    # dtw_rownames = [dtw_rownames[i] for i in sampled_indices]


    dtw_distance_matrix = dtw_ndim.distance_matrix_fast(dtw_input_tracks)
    dtw_distance_matrix = pd.DataFrame(dtw_distance_matrix, index=dtw_rownames, columns=dtw_rownames)

    umap_model = umap.UMAP(
            n_components=2, 
            n_neighbors=15, 
            min_dist=0.1, 
            init="random", 
            random_state=123,
            metric="precomputed", 
            )

    umap_embedding = umap_model.fit_transform(dtw_distance_matrix.values)
    umap_embedding = pd.DataFrame(umap_embedding, columns=['UMAP1', 'UMAP2'])
    umap_embedding[['TrackID', 'sample_name', 'position_t']] = pd.DataFrame(
        [string.split('--') for string in dtw_distance_matrix.index]
        )
    umap_embedding["TrackID"] = umap_embedding["TrackID"].astype(np.int64)
    umap_embedding["position_t"] = umap_embedding["position_t"].astype(np.float64)
    umap_embedding = umap_embedding.sort_values(by=["sample_name", "TrackID", "position_t"])

    keys = ["sample_name", "TrackID", "position_t"]
    rolling_for_sampled_windows = df_rolling.merge(
        umap_embedding[keys].drop_duplicates(),
        on=keys,
        how="inner"
    )

    df_umap = umap_embedding.merge(
        rolling_for_sampled_windows[keys + features],
        on=keys,
        how="left"
    )

    # df_umap = pd.merge(df_tracks, umap_embedding, how="left")
    # df_umap = pd.merge(df_rolling, umap_embedding, how="left")
    for feature in features:
        sns.scatterplot(
            data=df_umap,
            x="UMAP1",
            y="UMAP2",
            hue=feature,
            s=20,
            alpha=0.5,
            palette="viridis",
            # legend=False
        )
        plt.show()

"""
#############################

TSFRESH

#############################
"""


from tsfresh import select_features
from tsfresh.utilities.dataframe_functions import impute

id_columns = ["sample_name", "TrackID"]
features=[
    # "elongation",
    # "sphericity",
    "percentage_dead_mask",
    # "nr_dead_mask_pixels",
    "organoid_contact",
    "tcell_contact",
    # "displacement",
    # "mean_square_displacement",
    "speed",
    # # "dead",
    # "active_tcell_contact",
    # "position_t"
]

df_tracks_tsfresh = df_tracks.copy()

min_track_length=100
max_track_length=100

if min_track_length is not None:
    df_tracks_tsfresh=df_tracks_tsfresh.groupby(id_columns).filter(lambda group: len(group) >= min_track_length).reset_index(drop=True)
if max_track_length is not None:
    df_tracks_tsfresh=df_tracks_tsfresh.groupby(id_columns).apply(lambda group: group.iloc[:max_track_length]).reset_index(drop=True)
    
    
df_tracks_tsfresh["composite_id"] = (
    df_tracks_tsfresh["sample_name"] + "--" + df_tracks_tsfresh["TrackID"].astype(str)
)

df_tracks_tsfresh = df_tracks_tsfresh[features+["composite_id","position_t"]]
df_tracks_tsfresh["organoid_contact"] = df_tracks_tsfresh["organoid_contact"].astype(np.int16)
df_tracks_tsfresh["tcell_contact"] = df_tracks_tsfresh["tcell_contact"].astype(np.int16)

### Z-scale certain features
scaler = StandardScaler()
# scaler = RobustScaler()
df_tracks_tsfresh[features] = pd.DataFrame(scaler.fit_transform(pd.DataFrame(df_tracks_tsfresh[features])), columns=features)


extracted_features = extract_features(
    df_tracks_tsfresh, 
    column_id="composite_id", 
    column_sort="position_t",
    impute_function=impute,
    n_jobs=0
    )

extracted_features_scaled = extracted_features.copy()
scaler = StandardScaler()
extracted_features_scaled = scaler.fit_transform(extracted_features)
extracted_features_scaled = pd.DataFrame(extracted_features_scaled, 
                                         index=extracted_features.index, 
                                         columns=extracted_features.columns)
# Start with your extracted tsfresh features
features = extracted_features_scaled.copy()  # just to be safe

# 1. Drop columns that are all NaN
features_clean = features.dropna(axis=1, how='all')

# 2. Drop constant columns (no variance)
features_clean = features_clean.loc[:, features_clean.nunique() > 1]

# 3. Drop highly correlated columns (e.g., correlation > 0.95)
corr_matrix = features_clean.corr().abs()
upper_triangle = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
to_drop = [col for col in upper_triangle.columns if any(upper_triangle[col] > 0.95)]
features_clean = features_clean.drop(columns=to_drop)

# # 4. (Optional) Standardize features (for later clustering)
# scaler = StandardScaler()
# features_scaled = scaler.fit_transform(features_clean)

# # 5. Keep as DataFrame with original index and column names
# features_clean_scaled = pd.DataFrame(features_scaled, 
#                                      index=features_clean.index, 
#                                      columns=features_clean.columns)

# extracted_features.index.name = 'composite_id'
# features_clean = features_clean.reset_index()
# df_extracted_features.index.name = 'index'
# summarized_features = df_tracks_tsfresh.groupby(["composite_id"])[features].mean().reset_index()
# summarized_features = extracted_features.groupby(["composite_id"]).mean().reset_index()
df_extracted_features = features_clean.copy()
df_extracted_features = df_extracted_features.reset_index()

df_extracted_features[["sample_name", "TrackID"]] = df_extracted_features["composite_id"].str.split("--", expand=True)
df_extracted_features["sample_name"] = df_extracted_features["sample_name"].astype(str)
df_extracted_features["TrackID"] = df_extracted_features["TrackID"].astype(int)

umap_features = df_extracted_features.drop(columns=["composite_id", "sample_name", "TrackID"])
pca = PCA(n_components=20)
X_pca = pca.fit_transform(umap_features)

umap_model = umap.UMAP(
        n_components=2, 
        n_neighbors=30, 
        min_dist=0.1, 
        init="random", 
        random_state=123,
        # metric="precomputed", 
        )

# umap_embedding = umap_model.fit_transform(umap_features)
umap_embedding = umap_model.fit_transform(X_pca)
# df_extracted_features = df_extracted_features.reset_index().rename(columns={"index": "composite_id"})
df_extracted_features[["sample_name","TrackID"]] = df_extracted_features["composite_id"].str.split("--", expand=True)
df_extracted_features["sample_name"] = df_extracted_features["sample_name"].astype(str)
df_extracted_features["TrackID"] = df_extracted_features["TrackID"].astype(int)

umap_embedding = pd.DataFrame(umap_embedding, columns=['UMAP1', 'UMAP2'])

df_umap = umap_embedding
df_umap["sample_name"] = df_extracted_features["sample_name"]
df_umap["TrackID"] = df_extracted_features["TrackID"]
df_umap = pd.merge(df_umap, df_extracted_features, how="left", on=["TrackID", "sample_name"])

sns.scatterplot(
            data=df_umap,
            x="UMAP1",
            y="UMAP2",
            # hue=feature,
            s=40,
            alpha=0.95,
            palette="viridis",
            # legend=False
        )

summarized_features = df_tracks_tsfresh.groupby(["composite_id"])[features].mean().reset_index()

df_umap_plot = pd.merge(df_umap, summarized_features, how="left", on=["composite_id"])
for feature in features:
    sns.scatterplot(
            data=df_umap_plot,
            x="UMAP1",
            y="UMAP2",
            hue=feature,
            s=40,
            alpha=0.95,
            palette="viridis",
            # legend=False
        )
    plt.show()
