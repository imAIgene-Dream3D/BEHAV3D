from behav3d.preprocessing.segmentation.unet_segmentation import run_behav3d_unet_segmentation

from behav3d.preprocessing.tracking.laptracking import run_tcell_laptracking
from behav3d.utils.tracking import convert_all_tracked_images_to_csv

from behav3d.analysis.feature_extraction import run_feature_extraction
from behav3d.analysis.tcell_analysis import filter_tcell_tracks, summarize_track_features, run_tcell_analysis
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
from behav3d.utils.analysis import plot_filter_count

import yaml
import time
import seaborn as sns

output_dir = r"/Volumes/T7_sam/BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE"
metadata_csv_path = r"/Volumes/T7_sam/BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE/metadata.csv"
metadata = load_behav3d_metadata(metadata_csv_path)


from tsfresh import extract_features



window_size=50
groupby=["sample_name", "TrackID"]
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

analysis_outdir = Path(output_dir, "analysis", "tcell")
feature_outdir = Path(analysis_outdir, "track_features")
df_tracks_path = Path(feature_outdir, f"BEHAV3D_tcell_combined_track_features_filtered.csv")


df_tracks = pd.read_csv(df_tracks_path)
df_tracks=df_tracks.sort_values(by=["sample_name", "TrackID", "position_t"])

df_tracks = df_tracks[features+["sample_name", "TrackID", "position_t"]]

df_rolling = (
        df_tracks
        .groupby(groupby)[features]
        .rolling(window=window_size, min_periods=window_size)
        .mean()
        # .reset_index(drop=True)
        .reset_index()
        .drop(columns='level_2')
        # This gives you sample_name, TrackID, and the original index
    )
### Z-scale certain features
scaler = StandardScaler()
# scaler = RobustScaler()
## select all items in list not *contact*

non_binary_features = [x for x in features if "contact" not in x]

df_tracks[non_binary_features]= pd.DataFrame(
    scaler.fit_transform(df_tracks[non_binary_features]), 
    columns=df_tracks[non_binary_features].columns
    )


chosen_intervals = 10
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
df_umap = pd.merge(df_tracks, umap_embedding, how="left")

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
