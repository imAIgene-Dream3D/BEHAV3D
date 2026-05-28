"""
This script performs Dynamic Time Warpign to calculate distance between tracks.
It then fits this to a UMAP and perform K means clustering.
It then overlays the features back over the UMAP and creates
a heatmap with summarized feature values per cluster

-------------------------------------
--------------- INPUT ---------------
-------------------------------------

- BEHAV3D track features .csv
- umap_minimal_distance
- umap_n_neighbors
- nr_of_clusters

-------------------------------------
--------------- OUTPUT --------------
-------------------------------------

# Features of a track at each timepoint per sample in the metadata csv (.csv)
- See "FEATURES TRACKS"

# Combined summarized features for each track for all samples in metadata csv (.csv)
- 
"""
import argparse
import re
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
from behav3d.core.utils import format_time, expand_column_patterns
from behav3d.analysis.filtering import (
    round_legend_ticks
)
import yaml
import time
import seaborn as sns
# df_tracks=df_tracks[df_tracks["relative_time"]<=30]


def _handle_nan_in_distance_matrix(distance_matrix, context="distance matrix"):
    """
    Handle NaN values in a distance matrix by replacing with maximum finite value.
    
    Parameters
    ----------
    distance_matrix : pd.DataFrame or np.ndarray
        The distance matrix to check and fix
    context : str
        Description of where this matrix comes from (for logging)
        
    Returns
    -------
    pd.DataFrame or np.ndarray
        The distance matrix with NaN values replaced
        
    Notes
    -----
    If all values are NaN, falls back to 0 and logs a serious data quality warning.
    """
    is_dataframe = isinstance(distance_matrix, pd.DataFrame)
    
    if is_dataframe:
        nan_count = distance_matrix.isna().sum().sum()
        if nan_count == 0:
            return distance_matrix
        
        print(f"  ⚠️ Warning: {context} contains {nan_count} NaN values")
        print(f"  → Replacing with maximum finite value")
        
        max_val = distance_matrix.max().max()
        if pd.isna(max_val):
            print(f"  ⚠️ Warning: All values in {context} are NaN.")
            print("     → Falling back to 0 for all distances. Please check input data and DTW configuration.")
            return distance_matrix.fillna(0)
        return distance_matrix.fillna(max_val)
    else:
        # numpy array
        if not np.isnan(distance_matrix).any():
            return distance_matrix
        
        nan_count = np.isnan(distance_matrix).sum()
        print(f"  ⚠️ Warning: {context} contains {nan_count} NaN values")
        print(f"  → Replacing with maximum finite value")
        
        max_val = np.nanmax(distance_matrix)
        if np.isnan(max_val):
            print(f"  ⚠️ Warning: All values in {context} are NaN.")
            print("     → Falling back to 0 for all distances. Please check input data and DTW configuration.")
            return np.nan_to_num(distance_matrix, nan=0)
        return np.nan_to_num(distance_matrix, nan=max_val)
    
def run_tcell_analysis(
    config=None,
    output_dir=None,
    df_tracks_path=None,
    df_tracks_summarized_path=None,
    cell_type="tcell",
    columns_to_use=[
        "mean_square_displacement",
        "speed",
        "mean_dead_dye",
        "tcell_contact",
        "organoid_contact"
    ],
    columns_to_normalize=[
        "mean_square_displacement",
        "speed",
        "mean_dead_dye"
    ],
    # dtw_features=[
    #     "z_mean_square_displacement", 
    #     "z_speed", 
    #     "z_mean_dead_dye", 
    #     "tcell_contact", 
    #     "organoid_contact"
    #     ],
    umap_minimal_distance=None,
    umap_n_neighbors=None,
    nr_of_clusters=None,
    cluster_percentage_group_by=None,
    plot_results=True,
    seed=42,
    output_subdir_name="results",
    feature_scaling_preset=None,
    min_track_length=None,
    max_track_length=None,
    ):
    print(f"--------------- Performing {cell_type} behavioral analysis ---------------")
    start_time = time.time()
    # assert config is not None or all(
    #     [output_dir, umap_minimal_distance, umap_n_neighbors, nr_of_clusters]
    # ), "Either 'config' or 'output_dir, umap_minimal_distance, umap_n_neighbors, nr_of_clusters' parameters must be supplied"
        
    if output_dir is None:
        output_dir = config['output_dir']
    if feature_scaling_preset is None and config is not None:
        feature_scaling_preset = config.get("feature_scaling_preset", None)
        
    analysis_outdir = Path(output_dir, "analysis", cell_type)
    feature_outdir = Path(analysis_outdir, "track_features")
    
    if not analysis_outdir.exists():
        analysis_outdir.mkdir(parents=True)
    if not feature_outdir.exists():
        feature_outdir.mkdir(parents=True)    
    
    if df_tracks_path is None:
        df_tracks_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features_filtered.csv")
    if df_tracks_summarized_path is None:
        df_tracks_summarized_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features_summarized.csv")
    
    df_tracks = pd.read_csv(df_tracks_path, low_memory=False)
    if Path(df_tracks_summarized_path).exists():
        df_tracks_summarized = pd.read_csv(df_tracks_summarized_path, low_memory=False)
    else:
        df_tracks_summarized = summarize_feature_dtw_tracks(df_tracks)
    df_tracks=df_tracks.sort_values(by=["sample_name", "TrackID", "relative_time"])

    if min_track_length is not None or max_track_length is not None:
        df_tracks = filter_feature_dtw_track_lengths(
            df_tracks,
            min_track_length=min_track_length,
            max_track_length=max_track_length,
        )
        df_tracks_summarized = summarize_feature_dtw_tracks(df_tracks)
    
    if df_tracks_summarized["track_length"].nunique() != 1:
        print("Warning: The track lengths are not cut to similar length, this might influence dynamic time warping")
        print("Set 'min_track_length' and 'max_track_length' to the same value to create equal tracks")

    if feature_scaling_preset is None:
        non_wildcard_cols_to_normalize = []
        for col in columns_to_normalize:
            non_wildcard_cols_to_normalize += expand_column_patterns(col, df_tracks.columns)
        # print(f"{get_current_time()} - Perform z-normalization on selected feature columns")
        df_tracks=normalize_track_features(
            df_tracks, 
            columns=non_wildcard_cols_to_normalize
        )
        
        dtw_features = [x if x not in columns_to_normalize else f"z_{x}" for x in columns_to_use]
        plot_feature_cols = _resolve_selected_plot_feature_columns(
            df_tracks_summarized=df_tracks_summarized,
            df_tracks=df_tracks,
            selected_features=columns_to_use,
        )
    elif feature_scaling_preset == "original_behav3d":
        df_tracks, dtw_features = apply_original_behav3d_feature_scaling(
            df_tracks,
            cell_type=cell_type,
        )
        plot_feature_cols = _resolve_selected_plot_feature_columns(
            df_tracks_summarized=df_tracks_summarized,
            df_tracks=df_tracks,
            selected_features=[
                "mean_square_displacement",
                "speed",
                "mean_dead_dye",
                f"{cell_type}_contact",
                "organoid_contact",
            ],
        )
    else:
        raise ValueError(
            "Invalid feature_scaling_preset. Supported values are None and 'original_behav3d'."
        )
    
    dtw_distance_matrix = calculate_dtw(
        df_tracks, 
        features=dtw_features
        )
    
    umap_embedding = fit_umap(
        dtw_distance_matrix=dtw_distance_matrix,
        umap_n_neighbors=umap_n_neighbors,
        umap_minimal_distance=umap_minimal_distance,
        random_state=seed
    )
    
    df_clusters = cluster_umap(
        umap_embedding=umap_embedding,
        output_dir = output_dir,
        nr_of_clusters=nr_of_clusters,
        df_tracks=df_tracks,
        df_tracks_summarized=df_tracks_summarized,
        cell_type=cell_type,
        cluster_percentage_group_by=cluster_percentage_group_by,
        plot_results=plot_results,
        random_state=seed,
        output_subdir_name=output_subdir_name,
        plot_feature_cols=plot_feature_cols,
    )
    end_time = time.time()
    h,m,s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    return(df_clusters)


def filter_feature_dtw_track_lengths(
    df_tracks,
    min_track_length=None,
    max_track_length=None,
    group_cols=("sample_name", "TrackID"),
):
    """Filter tracks by number of timepoints and trim each track from its start."""
    group_cols = [col for col in group_cols if col in df_tracks.columns]
    if len(group_cols) == 0:
        raise ValueError("Track length filtering requires sample_name and/or TrackID columns.")

    time_cols = [col for col in ["relative_time", "position_t", "time"] if col in df_tracks.columns]
    sort_cols = group_cols + time_cols
    df_filtered = df_tracks.sort_values(sort_cols).reset_index(drop=True).copy()
    before_tracks = df_filtered[group_cols].drop_duplicates().shape[0]

    if min_track_length is not None:
        min_track_length = int(min_track_length)
        track_lengths = df_filtered.groupby(group_cols, observed=True)[group_cols[0]].transform("size")
        df_filtered = df_filtered.loc[track_lengths >= min_track_length].reset_index(drop=True)

    if max_track_length is not None:
        max_track_length = int(max_track_length)
        if max_track_length < 1:
            raise ValueError("max_track_length must be at least 1.")
        df_filtered = (
            df_filtered
            .groupby(group_cols, group_keys=False, observed=True)
            .head(max_track_length)
            .reset_index(drop=True)
        )

    after_tracks = df_filtered[group_cols].drop_duplicates().shape[0]
    print(
        "- Applied Feature DTW track length filtering: "
        f"{before_tracks} → {after_tracks} tracks"
    )
    if len(df_filtered) == 0:
        raise ValueError("No tracks remain after Feature DTW track length filtering.")
    return df_filtered


def summarize_feature_dtw_tracks(df_tracks):
    """Summarize the in-memory tracks used for Feature DTW plotting/clustering."""
    group_cols = ["sample_name", "TrackID"]
    missing = [col for col in group_cols if col not in df_tracks.columns]
    if missing:
        raise ValueError("Cannot summarize Feature DTW tracks; missing columns: " + ", ".join(missing))

    df_tracks = df_tracks.copy()
    for col in df_tracks.columns:
        if col.startswith("touching_"):
            df_tracks[col] = df_tracks[col].astype(str)

    grouped = df_tracks.groupby(group_cols, observed=True)
    df_summarized = grouped.size().reset_index(name="track_length")

    num_bool_cols = df_tracks.select_dtypes(include=[np.number, bool]).columns.drop(
        ["TrackID"], errors="ignore"
    )
    if len(num_bool_cols) > 0:
        means = grouped[num_bool_cols].mean().reset_index()
        means = means.rename(columns={col: f"mean_{col}" for col in num_bool_cols})
        df_summarized = pd.merge(df_summarized, means, on=group_cols, how="left")

    if "dead" in df_tracks.columns:
        df_summarized["dies"] = grouped["dead"].any().reset_index()["dead"]
    if "displacement_from_origin" in df_tracks.columns:
        df_summarized["displacement_from_origin"] = (
            grouped["displacement_from_origin"].max().reset_index()["displacement_from_origin"]
        )
    if "cumulative_displacement" in df_tracks.columns:
        df_summarized["cumulative_displacement"] = (
            grouped["cumulative_displacement"].last().reset_index()["cumulative_displacement"]
        )

    contact_cols = [
        col for col in df_tracks.columns
        if col.endswith("_contact") and not col.startswith("active_")
    ]
    for contact_col in contact_cols:
        active_col = f"active_{contact_col}"
        if active_col in df_tracks.columns:
            def _active_contact_mean(group, contact_col=contact_col, active_col=active_col):
                mask = group[contact_col].fillna(False).astype(bool)
                return group.loc[mask, active_col].mean() if mask.any() else 0

            result = grouped.apply(_active_contact_mean, include_groups=False)
            if isinstance(result, pd.DataFrame):
                result = result.iloc[:, 0]
            df_summarized[active_col] = result.values

    metadata_cols = ["TrackID", "sample_name"]
    for col in ["well", "exp_nr"]:
        if col in df_tracks.columns:
            metadata_cols.append(col)
    metadata_cols.extend([col for col in df_tracks.columns if col.endswith("_line_condition")])
    metadata_cols = list(dict.fromkeys(metadata_cols))

    df_trackinfo = df_tracks[metadata_cols].drop_duplicates()
    return pd.merge(df_trackinfo, df_summarized, how="left")


def _resolve_selected_plot_feature_columns(
    df_tracks_summarized,
    df_tracks,
    selected_features,
    summary_match_mode="all",
):
    """Resolve user-selected feature names to columns available for result PDFs."""
    if selected_features is None:
        selected_features = []
    elif isinstance(selected_features, str):
        selected_features = [selected_features]
    else:
        selected_features = list(selected_features)
    summary_cols = [str(c) for c in list(df_tracks_summarized.columns)]
    summary_set = set(summary_cols)
    track_cols = set(str(c) for c in list(df_tracks.columns))
    out = []
    seen = set()

    def _is_allowed_plot_column(col):
        return "_on_distance" not in str(col)

    def _add(col):
        col = str(col)
        if col in summary_set and col not in seen and _is_allowed_plot_column(col):
            out.append(col)
            seen.add(col)

    def _summary_matches_for_feature(raw_feature):
        token = str(raw_feature).strip()
        if token == "":
            return []
        escaped = re.escape(token)
        pattern = re.compile(rf"(^|_){escaped}($|_)")
        matched = [
            col for col in summary_cols
            if _is_allowed_plot_column(col) and bool(pattern.search(str(col)))
        ]
        if summary_match_mode == "all":
            return matched
        return matched[:1]

    for feature in selected_features:
        feature = str(feature)
        expanded = expand_column_patterns(feature, df_tracks_summarized.columns)
        exact_summary_matches = [
            str(col) for col in expanded
            if str(col) in summary_set and _is_allowed_plot_column(col)
        ]
        for col in exact_summary_matches:
            _add(col)

        base_feature = feature[2:] if feature.startswith("z_") else feature
        base_feature = str(base_feature).strip()
        if base_feature in summary_set:
            _add(base_feature)
            continue

        mean_feature = f"mean_{base_feature}"
        if mean_feature in summary_set:
            _add(mean_feature)
            continue

        if len(exact_summary_matches) == 0 and base_feature not in summary_set:
            for col in _summary_matches_for_feature(base_feature):
                _add(col)
            if base_feature in track_cols:
                for col in _summary_matches_for_feature(base_feature):
                    _add(col)
    return out




def calculate_dtw(
    df_tracks, 
    features=[
        "z_mean_square_displacement", 
        "z_speed", 
        "z_mean_dead_dye", 
        "tcell_contact", 
        "organoid_contact"
        ]
    ):
    
    print("- Calculating the dynamic time warping distance matrix")
    df_tracks=df_tracks.sort_values(by=["sample_name", "TrackID", "relative_time"])
    
    nr_tracks=len(df_tracks[["sample_name", "TrackID"]].drop_duplicates())
    nr_timepoints=len(df_tracks["relative_time"].unique())
    nr_features=len(features)
    
    # Check for NaN in selected features and warn user
    nan_counts = df_tracks[features].isna().sum()
    features_with_nan = nan_counts[nan_counts > 0]
    if len(features_with_nan) > 0:
        print(f"  ⚠️ Warning: Found NaN values in selected features:")
        for feat, count in features_with_nan.items():
            print(f"     - {feat}: {count} NaN values")
        print(f"  → Filling NaN with 0 for DTW calculation")
    
    # Fill NaN values with 0 for DTW calculation
    df_tracks_filled = df_tracks.copy()
    df_tracks_filled[features] = df_tracks_filled[features].fillna(0)

    dtw_input_tracks=[]
    dtw_rownames=[]
    unique_tracks = df_tracks_filled.groupby(['TrackID', 'sample_name'])
    for (TrackID, sample_name), group in unique_tracks:
        track_features = group[features].to_numpy().astype(np.double)
        dtw_rownames.append(f"{TrackID}--{sample_name}")
        dtw_input_tracks.append(track_features)
    
    dtw_distance_matrix = dtw_ndim.distance_matrix_fast(dtw_input_tracks)
    dtw_distance_matrix = pd.DataFrame(dtw_distance_matrix, index=dtw_rownames, columns=dtw_rownames)
    
    # Check and handle NaN values in the distance matrix
    dtw_distance_matrix = _handle_nan_in_distance_matrix(
        dtw_distance_matrix, context="DTW distance matrix"
    )
    
    return(dtw_distance_matrix)

def fit_umap(
    dtw_distance_matrix,
    config=None,
    umap_n_neighbors=None,
    umap_minimal_distance=None,
    random_state=None
    ):
    
    print("- Fitting the dynamic time warping to a UMAP")
    assert config is not None or all(
            [umap_n_neighbors, umap_minimal_distance]
        ), "Either 'config' or 'umap_n_neighbors and umap_minimal_distance' must be supplied"
            
    if umap_n_neighbors is None or umap_minimal_distance is None:
        assert config is not None, "Provide config or both umap_n_neighbors and umap_minimal_distance"
        umap_n_neighbors = config['umap_n_neighbors']
        umap_minimal_distance = config["umap_minimal_distance"]
    
    # Final safety check for NaN values in the distance matrix
    dtw_matrix_values = dtw_distance_matrix.values.copy()
    dtw_matrix_values = _handle_nan_in_distance_matrix(
        dtw_matrix_values, context="distance matrix (pre-UMAP)"
    )
        
    umap_model = umap.UMAP(
        n_components=2, 
        n_neighbors=umap_n_neighbors, 
        min_dist=umap_minimal_distance, 
        init="random", 
        random_state=random_state,
        metric="precomputed", 
        )
    umap_embedding = umap_model.fit_transform(dtw_matrix_values)
    umap_embedding = pd.DataFrame(umap_embedding, columns=['UMAP1', 'UMAP2'])
    umap_embedding[['TrackID', 'sample_name']] = pd.DataFrame(
        [string.split('--') for string in dtw_distance_matrix.index]
        )
    umap_embedding["TrackID"] = umap_embedding["TrackID"].astype(np.int64)
    return(umap_embedding)

def cluster_umap(
    umap_embedding,
    config=None,
    nr_of_clusters=None,
    df_tracks=None,
    df_tracks_summarized=None,
    random_state=None,
    output_dir = None,
    cell_type="tcell",
    cluster_percentage_group_by=None,
    plot_results=True,
    output_subdir_name="results",
    plot_feature_cols=None,
    ):
    
    assert config is not None or all(
        [output_dir, nr_of_clusters]
    ), "Either 'config' or 'output_dir, nr_of_clusters' parameters must be supplied"
      
    print("- Performing clustering on the UMAP data")
    if all([output_dir, nr_of_clusters]) is None:
        output_dir = Path(config['output_dir'])
        nr_of_clusters=config["nr_of_clusters"]
    
    analysis_outdir = Path(output_dir, "analysis", cell_type)
    feature_outdir = Path(analysis_outdir, "track_features")
    results_outdir = Path(analysis_outdir, str(output_subdir_name).strip() or "results")
    if not analysis_outdir.exists():
        analysis_outdir.mkdir(parents=True)
    if not feature_outdir.exists():
        feature_outdir.mkdir(parents=True)
    if not results_outdir.exists():
        results_outdir.mkdir(parents=True)
    
    if df_tracks is None:
        df_tracks_path = Path(feature_outdir, f"BEHAV3D_combined_track_features_filtered.csv")
        df_tracks = pd.read_csv(df_tracks_path)
    if df_tracks_summarized is None:
        df_tracks_summarized_path = Path(feature_outdir, f"BEHAV3D_combined_track_features_summarized.csv")
        df_tracks_summarized = pd.read_csv(df_tracks_summarized_path)
      
    # df_tracks=df_tracks.sort_values(by=["sample_name", "TrackID", "relative_time"])
    # TrackIDs = df_tracks[["sample_name", "TrackID"]].drop_duplicates().reset_index(drop=True)
    # df_umap = pd.concat([TrackIDs, umap_embedding], axis=1)
    # df_trackinfo = df_tracks[['TrackID', 'sample_name','well', 'exp_nr', 'organoid_line', 'tcell_line']].drop_duplicates()
    df_umap = pd.merge(df_tracks_summarized, umap_embedding, how="left")
    
    # Perform clustering
    scaler = StandardScaler()
    # umap_scaled = scaler.fit_transform(umap_embedding)  # Standardize UMAP coordinates
    kmeans = KMeans(n_clusters=nr_of_clusters, n_init=100, random_state=random_state)
    
    ### TODO Do the clustering based straight on the DTW distances
    ### Miguel suggested clusters are more meaningful before the umap embedding
    
    umap_embedding.columns.difference(df_tracks_summarized.columns)
    
    df_umap["ClusterID"] = kmeans.fit_predict(
        df_umap[umap_embedding.columns.difference(df_tracks_summarized.columns)]
        )
    
    # df_umap["ClusterID"] = kmeans.fit_predict(
    #     umap_embedding.drop(columns=["TrackID","sample_name"])
    #     )
    # df_umap["cluster2"] = kmeans.fit_predict(umap_embedding)
    
    # Set cluster index to start from 1 for backprojection purposes
    df_umap["ClusterID"]=df_umap["ClusterID"]+1
    df_umap["ClusterID"]=df_umap["ClusterID"].astype('category')
    
    df_umap_out_path = Path(results_outdir, f"BEHAV3D_{cell_type}_UMAP_clusters.csv")
    print(f"- Writing clustered tracks to {df_umap_out_path}")
    df_umap.to_csv(df_umap_out_path, sep=",", index=False)

    # Dynamically detect *_line_condition columns instead of hardcoding
    line_cols = [c for c in df_tracks_summarized.columns if c.endswith('_line_condition')]
    sample_cols = line_cols + ["sample_name"]
    if "well" in df_tracks_summarized.columns:
        sample_cols.append("well")
    if "exp_nr" in df_tracks_summarized.columns:
        sample_cols.append("exp_nr")
    if cluster_percentage_group_by is None:
        cluster_percentage_context_cols = []
    elif isinstance(cluster_percentage_group_by, str):
        cluster_percentage_context_cols = [cluster_percentage_group_by]
    else:
        cluster_percentage_context_cols = list(cluster_percentage_group_by)
    metadata_cols = [
        c for c in sample_cols + cluster_percentage_context_cols
        if c in df_umap.columns
    ]
    if plot_feature_cols is None:
        info_cols = list(df_umap.drop(columns=["TrackID", "UMAP1", "UMAP2", "ClusterID"]).columns)
    else:
        selected_plot_cols = [c for c in list(plot_feature_cols) if c in df_umap.columns]
        if len(selected_plot_cols) == 0:
            print(
                "  ⚠️ Warning: None of the selected DTW features could be resolved to summarized columns for plotting. "
                "Output PDFs will show metadata/context columns only."
            )
        info_cols = list(dict.fromkeys(metadata_cols + selected_plot_cols))
    
    cluster_UMAP_path = Path(results_outdir, f"BEHAV3D_{cell_type}_UMAP_clusters.pdf")
    print(f"- Plotting clustered UMAP plots with displayed Track features to {cluster_UMAP_path}")
    plot_feature_umap(
        df_umap=df_umap,
        info_cols=info_cols,
        sample_cols=sample_cols,
        outpath=cluster_UMAP_path,
        rows_per_page = 4,
        nr_cols = 2,
        rows_first_img = 2,
        figsize = (8.27, 11.69),
        plot_results=plot_results
    )
    
    ### Producing a heatmap of the summarized features again summarized over all tracks
    ### Belonging to that cluster
    cluster_features_heatmap_path = Path(results_outdir, f"BEHAV3D_{cell_type}_UMAP_cluster_feature_heatmap.pdf")
    print(f"- Plotting heatmaps with summarized cluster features to {cluster_features_heatmap_path}")
    plot_clustering_feature_heatmap(
        df_umap,
        info_cols,
        sample_cols,
        cluster_features_heatmap_path,
        rows_per_page = 7,
        nr_cols = 2,
        figsize = (8.27, 11.69),
        plot_results=plot_results
    )

    # Dynamically detect ALL line_condition columns (for grid plotting across all cell types)
    all_line_cols = [c for c in df_umap.columns if c.endswith('_line_condition')]
    
    # Prioritize THIS cell type's line_condition to be first (for rows in grid)
    # Look for line columns that match this cell_type
    this_cell_line_cols = [c for c in all_line_cols if f"_{cell_type}_" in c]
    other_line_cols = [c for c in all_line_cols if c not in this_cell_line_cols]
    
    # Arrange: this cell type's conditions first, then others
    line_cols = this_cell_line_cols + other_line_cols
    
    # Include cluster_percentage_group_by columns (e.g., exp_nr) if they exist
    if cluster_percentage_group_by is None:
        cluster_percentage_groups = []
    elif isinstance(cluster_percentage_group_by, str):
        cluster_percentage_groups = [cluster_percentage_group_by]
    else:
        cluster_percentage_groups = list(cluster_percentage_group_by)
    extra_group_cols = []
    for col in cluster_percentage_groups:
        if col in df_umap.columns and col not in line_cols:
            extra_group_cols.append(col)
    
    group_cols_sample = line_cols + extra_group_cols + ["sample_name", "ClusterID"]
    group_cols_sample_only = line_cols + extra_group_cols + ["sample_name"]
    
    df_clust_perc = (
        df_umap
        .groupby(group_cols_sample, observed=True)
        .size()
        .reset_index(name="count")
    )
    
    sample_totals = (
        df_clust_perc
        .groupby(group_cols_sample_only, observed=True)["count"]
        .sum()
        .reset_index(name="sample_total")
    )
    
    combo_totals = (
        df_clust_perc
        .groupby(line_cols, observed=True)["count"]
        .sum()
        .reset_index(name="combo_total")
    ) if line_cols else None
    
    df_clust_perc = df_clust_perc.merge(sample_totals, how="left", on=group_cols_sample_only)
    if combo_totals is not None:
        df_clust_perc = df_clust_perc.merge(combo_totals, how="left", on=line_cols)
    df_clust_perc["percentage"] = (df_clust_perc["count"] / df_clust_perc["sample_total"])
    
    cluster_percentage_plot_prefix = Path(results_outdir, f"BEHAV3D_{cell_type}_UMAP_cluster_percentages")
    print(f"- Plotting percentage plots of each cluster per sample to {cluster_percentage_plot_prefix}_*.pdf")

    plot_cluster_percentage_bars(
        df_clust_perc,
        cluster_percentage_plot_prefix,
        group_by_columns=cluster_percentage_group_by,
        plot_results=plot_results
    )
        
    df_clust_perc = df_clust_perc.reset_index(drop=True)
    df_clust_perc_out_path = Path(results_outdir, f"BEHAV3D_{cell_type}_UMAP_cluster_percentages.csv")
    print(f"- Writing cluster percentages to {df_clust_perc_out_path}")
    df_clust_perc.to_csv(df_clust_perc_out_path, sep=",", index=False)
    
    df_clust_tracks_out_path = Path(results_outdir, f"BEHAV3D_{cell_type}_combined_track_features_clustered.csv")
    print(f"- Writing clustered track info to {df_clust_tracks_out_path}")
    df_tracks = pd.merge(df_tracks, df_umap[["TrackID", "ClusterID"]], on='TrackID', how='left')
    
    df_tracks.to_csv(df_clust_tracks_out_path, sep=",", index=False)
    return()

def _minmax_scale(series):
    series = pd.to_numeric(series, errors="coerce").astype("float64")
    min_value = series.min(skipna=True)
    max_value = series.max(skipna=True)
    if pd.isna(min_value) or pd.isna(max_value) or max_value == min_value:
        return pd.Series(0.0, index=series.index)
    return (series - min_value) / (max_value - min_value)

def _original_behav3d_quantile_scale(values):
    values = pd.to_numeric(values, errors="coerce").astype("float64")
    std = values.std()
    if pd.isna(std) or std == 0:
        z_values = pd.Series(0.0, index=values.index)
    else:
        z_values = (values - values.mean()) / std

    q75 = z_values.quantile(0.75)
    min_z = z_values.min(skipna=True)
    if pd.isna(q75) or pd.isna(min_z):
        q_values = pd.Series(0.0, index=values.index)
    else:
        q_values = z_values.where(z_values > q75, min_z)

    q_values = _minmax_scale(q_values)
    max_quantile = q_values.quantile(0.9999999)
    if pd.notna(max_quantile) and max_quantile != 0:
        q_values = q_values / max_quantile

    return q_values

def apply_original_behav3d_feature_scaling(df_tracks, cell_type="tcell"):
    """
    Add the original BEHAV3D/R-style scaled feature columns used for feature DTW.
    """
    contact_col = f"{cell_type}_contact"
    required_cols = [
        "mean_square_displacement",
        "speed",
        "mean_dead_dye",
        "organoid_contact",
        contact_col,
        "exp_nr",
    ]
    missing_cols = [col for col in required_cols if col not in df_tracks.columns]
    if missing_cols:
        raise ValueError(
            "feature_scaling_preset='original_behav3d' requires missing columns: "
            + ", ".join(missing_cols)
        )

    print("- Applying original BEHAV3D feature scaling preset")
    df_scaled = df_tracks.copy()
    quantile_features = {
        "mean_square_displacement": "q_mean_square_displacement",
        "speed": "q_speed",
        "mean_dead_dye": "q_mean_dead_dye",
    }

    group_key = "exp_nr"
    for source_col, out_col in quantile_features.items():
        df_scaled[out_col] = (
            df_scaled
            .groupby(group_key, group_keys=False)[source_col]
            .transform(_original_behav3d_quantile_scale)
        )

    contact_features = {
        "organoid_contact": "s_organoid_contact",
        contact_col: f"s_{contact_col}",
    }
    for source_col, out_col in contact_features.items():
        df_scaled[out_col] = (
            df_scaled
            .groupby(group_key, group_keys=False)[source_col]
            .transform(_minmax_scale)
        )

    dtw_features = [
        "q_mean_square_displacement",
        "q_speed",
        "q_mean_dead_dye",
        "s_organoid_contact",
        f"s_{contact_col}",
    ]
    return df_scaled, dtw_features

def normalize_track_features(
    df_tracks,
    columns = [
        "mean_square_displacement",
        "speed",
        "mean_dead_dye"
    ]
    ):
    
    for column_name in columns:
        # Skip if column doesn't exist
        if column_name not in df_tracks.columns:
            continue
        
        # Skip non-numeric columns (e.g., orientation_vector which is a string representation)
        if df_tracks[column_name].dtype == 'object':
            # Try to convert to numeric, skip if it fails
            try:
                df_tracks[column_name] = pd.to_numeric(df_tracks[column_name], errors='coerce')
            except:
                print(f"Skipping normalization of non-numeric column: {column_name}")
                continue
        
        # Check if column has valid numeric data and variance
        if df_tracks[column_name].std() > 0:
            df_tracks.loc[:, f'z_{column_name}'] = df_tracks[column_name].transform(
                lambda x: (x - x.mean()) / x.std() if x.std() > 0 else 0
            )
        else:
            print(f"Skipping normalization of zero-variance column: {column_name}")
            df_tracks[f'z_{column_name}'] = 0
    
    return (df_tracks)

def plot_cluster_percentage_bars(
    df_clust_perc: pd.DataFrame,
    outprefix,
    group_by_columns=None,
    plot_results=True,
    grid_figsize=(20, 10),
    ncols_samples=3,
    sample_panel_width=6.0,
    sample_row_height=2.8,
    title_fontsize=24,
    label_fontsize=12,
    info_fontsize=10
):
    """
    Produce a PDF with cluster percentage visualizations.

    Page(s) 1+: 'Combined' view (samples pooled). Grid layout:
                rows = this cell_type's line_condition values
                columns = values of user-selected grouping column
                One page per selected grouping column.
                If group_by_columns is empty/None, skip these pages.

    Last Page: 'Per-sample' view. Panels laid out in 3 columns (configurable).
               Each panel is a single horizontal stacked bar of ClusterID percentages
               for one sample_name.
    """
    outprefix = Path(outprefix)
    outpdf = outprefix.with_suffix(".pdf")

    # Ensure group_by_columns is always a list
    if group_by_columns is None:
        group_by_columns = []
    elif isinstance(group_by_columns, str):
        group_by_columns = [group_by_columns]
    
    # Filter to columns that actually exist in the data
    group_by_columns = [c for c in group_by_columns if c in df_clust_perc.columns]

    # Dynamically detect all *_line_condition columns
    # The first one is this cell type's column (cluster_umap orders them: this cell type first)
    all_line_cols = [c for c in df_clust_perc.columns if c.endswith('_line_condition')]
    this_cell_line_col = all_line_cols[0] if all_line_cols else None
    
    # Exclude the row column from group_by_columns (can't be both rows and columns)
    if this_cell_line_col:
        group_by_columns = [c for c in group_by_columns if c != this_cell_line_col]
    
    # ------------- Prepare pooled data for grid pages -------------
    # Include group_by_columns in the pooled grouping so we can filter by them
    pooled_group_cols = all_line_cols + group_by_columns + ['ClusterID']
    # Remove duplicates while preserving order
    pooled_group_cols = list(dict.fromkeys(pooled_group_cols))
    
    pooled = (
        df_clust_perc
        .groupby(pooled_group_cols, as_index=False, observed=True)['count']
        .sum()
    )
    # Total count per line_condition combo (across ClusterID)
    pooled_total_cols = all_line_cols + group_by_columns
    pooled_total_cols = list(dict.fromkeys(pooled_total_cols))
    
    if pooled_total_cols:
        pooled_combo_totals = (
            pooled.groupby(pooled_total_cols, observed=True)['count']
            .sum()
            .reset_index(name='combo_total_pooled')
        )
        pooled = pooled.merge(pooled_combo_totals, on=pooled_total_cols, how='left')
    else:
        pooled['combo_total_pooled'] = pooled['count'].sum()
    pooled['percentage_pooled'] = pooled['count'] / pooled['combo_total_pooled']

    # ------------- Prepare per-sample data -------------
    per_sample = (
        df_clust_perc
        .groupby(['sample_name', 'ClusterID'], as_index=False, observed=True)['count']
        .sum()
    )
    per_sample_totals = (
        per_sample.groupby('sample_name', observed=True)['count']
        .sum()
        .reset_index(name='sample_total_all')
    )
    per_sample = per_sample.merge(per_sample_totals, on='sample_name', how='left')
    per_sample['percentage_overall'] = per_sample['count'] / per_sample['sample_total_all']
    samples = per_sample['sample_name'].drop_duplicates().tolist()

    # Common ClusterID ordering for consistent legends across pages
    cluster_order = (
        df_clust_perc['ClusterID']
        .drop_duplicates()
        .sort_values()
        .tolist()
    )
    
    # ------------- Build list of grid pages to generate -------------
    grids_to_plot = []
    
    if group_by_columns and this_cell_line_col:
        # One grid page per selected column
        line_values_rows = df_clust_perc[this_cell_line_col].drop_duplicates().sort_values().tolist()
        nrows = len(line_values_rows)
        
        for col_for_columns in group_by_columns:
            if col_for_columns in df_clust_perc.columns:
                col_values = df_clust_perc[col_for_columns].drop_duplicates().sort_values().tolist()
                grids_to_plot.append({
                    'row_col': this_cell_line_col,
                    'row_values': line_values_rows,
                    'col_col': col_for_columns,
                    'col_values': col_values,
                    'nrows': nrows,
                    'ncols': len(col_values)
                })
    
    
    with PdfPages(outpdf) as pdf:
        for grid_config in grids_to_plot:
            row_col = grid_config['row_col']
            row_values = grid_config['row_values']
            col_col = grid_config['col_col']
            col_values = grid_config['col_values']
            nrows = grid_config['nrows']
            ncols = grid_config['ncols']
            
            # ------------------------------- GRID PAGE -------------------------------
            fig1, axes1 = plt.subplots(
                nrows,
                ncols,
                figsize=grid_figsize,
                sharex=True,
                sharey=True,
                squeeze=False
            )

            legend_handles_1, legend_labels_1 = None, None

            for i, row_val in enumerate(row_values):
                for j, col_val in enumerate(col_values):
                    ax = axes1[i, j]

                    # Build subset filter
                    subset = pooled[
                        (pooled[row_col] == row_val) &
                        (pooled[col_col] == col_val)
                    ]
                    
                    # Set column headers (top row)
                    if i == 0:
                        ax.set_title(f'{col_val}', fontsize=title_fontsize)
                    
                    # Set row labels (left side)
                    if j == 0:
                        ax.text(
                            -0.05, 0.5, f'{row_val}',
                            ha='right', va='center',
                            rotation=0, transform=ax.transAxes,
                            fontsize=title_fontsize
                        )

                    if subset.empty:
                        ax.axis('off')
                        continue

                    # Ensure one row per ClusterID
                    subset_unique = subset.drop_duplicates(subset=['ClusterID'], keep='first')
                    
                    row = (subset_unique
                           .set_index('ClusterID')
                           .reindex(cluster_order)
                           ['percentage_pooled']
                           .fillna(0.0))
                    pivot = row.to_frame().T
                    pivot.index = ['combined']

                    bar_ax = pivot.plot(kind='barh', stacked=True, ax=ax, legend=False)

                    if legend_handles_1 is None:
                        legend_handles_1, legend_labels_1 = bar_ax.get_legend_handles_labels()

                    num_cells = int(subset['combo_total_pooled'].iloc[0]) if not subset.empty else 0
                    ax.text(
                        0.5, 0.08, f'# Cells (pooled): {num_cells}',
                        ha='center', va='center',
                        transform=ax.transAxes, fontsize=info_fontsize
                    )

                    ax.axis('off')
                    
            row_col_name = row_col.replace('_line_condition', '').replace('_', ' ').title()
            col_col_name = col_col.replace('_line_condition', '').replace('_', ' ').title()
            title_text = f'Cluster percentages: {row_col_name} × {col_col_name}'
            fig1.suptitle(title_text, fontsize=title_fontsize, y=0.995)

            if legend_handles_1 and legend_labels_1:
                fig1.legend(
                    legend_handles_1, legend_labels_1,
                    fontsize=label_fontsize,
                    title='ClusterID', title_fontsize=label_fontsize,
                    bbox_to_anchor=(0.92, 0.5), loc='center left'
                )
                fig1.tight_layout(rect=[0, 0, 0.88, 0.96])
            else:
                fig1.tight_layout(rect=[0, 0, 1, 0.96])

            if plot_results:
                plt.show()
            pdf.savefig(fig1, bbox_inches='tight')
            plt.close(fig1)

        # ------------------------------- PER-SAMPLE PAGE -------------------------------
        n_panels = len(samples)
        ncols = ncols_samples
        nrows = int(np.ceil(n_panels / ncols))
        fig_w = sample_panel_width * ncols
        fig_h = sample_row_height * nrows

        fig2, axes2 = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h), squeeze=False)

        legend_handles_2, legend_labels_2 = None, None

        for k, sample in enumerate(samples):
            i, j = divmod(k, ncols)
            ax = axes2[i, j]

            sub = per_sample[per_sample['sample_name'] == sample]
            if sub.empty:
                ax.axis('off')
                continue

            # Ensure one row per ClusterID
            sub_unique = sub.drop_duplicates(subset=['ClusterID'], keep='first')
            
            row = (sub_unique
                   .set_index('ClusterID')
                   .reindex(cluster_order)
                   ['percentage_overall']
                   .fillna(0.0))
            pivot = row.to_frame().T
            pivot.index = [sample]

            bar_ax = pivot.plot(kind='barh', stacked=True, ax=ax, legend=False)

            if legend_handles_2 is None:
                legend_handles_2, legend_labels_2 = bar_ax.get_legend_handles_labels()

            num_cells = int(sub['sample_total_all'].iloc[0]) if not sub.empty else 0
            ax.set_title(sample, fontsize=label_fontsize)
            ax.text(
                0.5, 0.08, f'# Cells: {num_cells}',
                ha='center', va='center',
                transform=ax.transAxes, fontsize=info_fontsize
            )

            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_yticks([])
            ax.set_xticks([])
            ax.spines[['top','right','bottom','left']].set_visible(False)

        # Hide unused axes
        for k in range(n_panels, nrows * ncols):
            i, j = divmod(k, ncols)
            axes2[i, j].axis('off')

        fig2.suptitle('Per-sample cluster composition', fontsize=title_fontsize, y=0.995)

        if legend_handles_2 and legend_labels_2:
            fig2.legend(
                legend_handles_2, legend_labels_2,
                fontsize=label_fontsize,
                title='ClusterID', title_fontsize=label_fontsize,
                bbox_to_anchor=(0.92, 0.5), loc='center left'
            )
            fig2.tight_layout(rect=[0, 0, 0.88, 0.96])
        else:
            fig2.tight_layout(rect=[0, 0, 1, 0.96])

        if plot_results:
            plt.show()
        pdf.savefig(fig2, bbox_inches='tight')
        plt.close(fig2)
                      
def plot_feature_umap(
    df_umap,
    info_cols,
    sample_cols,
    outpath,
    rows_per_page = 4,
    nr_cols = 2,
    rows_first_img = 2,
    figsize = (8.27, 11.69),
    plot_results=True
    ):
    n_plots = len(info_cols)
    n_rows = (n_plots // nr_cols) + (1 if n_plots % nr_cols != 0 else 0) + rows_first_img
    nr_pages = (n_rows // rows_per_page) + (1 if n_rows % rows_per_page != 0 else 0)

    # Create PDF file
    
    with PdfPages(outpath) as pdf:
        plot_idx = 0  # Track which plot we are adding

        for page in range(nr_pages):
            fig = plt.figure(figsize=figsize)
            gs = GridSpec(rows_per_page, nr_cols, figure=fig, wspace=0.3)

            # First image on the first page
            if page == 0:
                ax = fig.add_subplot(gs[:rows_first_img, :])
                sns.scatterplot(
                    data=df_umap,
                    x="UMAP1",
                    y="UMAP2",
                    hue="ClusterID",
                    s=20,
                    alpha=0.5,
                    palette="Set1",
                    ax=ax,
                )
                ax.legend(
                    loc="upper left", 
                    prop={'size': 8}, 
                    bbox_to_anchor=(1, 1), 
                    borderpad=0.3, 
                    labelspacing=0.4,
                    columnspacing=0.1,
                    frameon=False
                    )
                # legend.set_title("ClusterID", prop={'size': 9})
                ax.set_title("ClusterID", fontsize=10, loc='center')
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                ax.set_xlabel("")
                ax.set_ylabel("")
                # plot_idx += 1  # Increment for the next plots
                if plot_results:
                    plt.show()
            # Remaining plots
            remaining_axes = [
                fig.add_subplot(gs[i, j])
                for i in range(rows_first_img if page == 0 else 0, rows_per_page)
                for j in range(nr_cols)
            ]

            for ax in remaining_axes:
                if plot_idx >= len(info_cols):
                    ax.remove()  # Remove empty axes
                    continue
                colorcol = info_cols[plot_idx]
                if colorcol in sample_cols or df_umap.dtypes[colorcol]==bool:
                    sns.scatterplot(
                        data=df_umap,
                        x="UMAP1",
                        y="UMAP2",
                        s=20,
                        hue=colorcol,
                        alpha=0.5,
                        palette="Set2",
                        ax=ax
                    )
                else:
                    sns.scatterplot(
                        data=df_umap,
                        x="UMAP1",
                        y="UMAP2",
                        s=20,
                        hue=colorcol,
                        palette="viridis",
                        alpha=0.5,
                        ax=ax
                    )

                # Reduce legend size & move outside plot
                ax.legend(
                    loc="upper left", 
                    prop={'size': 8}, 
                    bbox_to_anchor=(1, 1), 
                    borderpad=0.3, 
                    labelspacing=0.4,
                    columnspacing=0.1,
                    frameon=False
                    )
                # legend.set_title(colorcol, prop={'size': 9})
                ax.set_title(colorcol, fontsize=10, loc='center')
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                ax.set_xlabel("")
                ax.set_ylabel("")

                plot_idx += 1  # Move to the next plot

            # Save and close the figure for this page
            # fig.tight_layout()
            fig.subplots_adjust(left=0.05, right=0.85, top=0.95, bottom=0.05)
            # fig.tight_layout(pad=2.0)
            # fig.set_constrained_layout(True)
            # plt.show()
            pdf.savefig(fig, dpi=600)
            plt.close(fig)

def plot_clustering_feature_heatmap(
    df_umap,
    info_cols,
    sample_cols,
    outpath,
    rows_per_page = 7,
    nr_cols = 2,
    figsize = (8.27, 11.69),
    plot_results=True,
    show_points=False,
    point_alpha=0.5,
    point_size=8,
    mean_marker_size=60,
):
    """
    Produce a PDF with:
      • Page 1: full-page min–max scaled heatmap of cluster means.
      • Subsequent pages: per-feature violin plots tiled across pages.
    """

    info_cols   = list(info_cols) if info_cols is not None else []
    sample_cols = list(sample_cols) if sample_cols is not None else []

    # ---- Helper ----
    def _round_legend_ticks(max_val):
        try:
            return round_legend_ticks(max_val)
        except Exception:
            if not np.isfinite(max_val) or max_val <= 0:
                return 1.0
            magnitude = 10.0 ** np.floor(np.log10(max_val))
            return float(np.ceil(max_val / magnitude) * magnitude)

    # ---- Cluster means ----
    _cols_for_means = list(dict.fromkeys(list(info_cols) + ["ClusterID"]))
    df_for_means = (
        df_umap[_cols_for_means]
        .drop(columns=sample_cols, errors="ignore")
    )
    cluster_means = (
        df_for_means
        .groupby("ClusterID", observed=False)
        .mean(numeric_only=True)
        .reset_index()
    )

    # ---- Min-max scaling ----
    cluster_means_scaled = cluster_means.copy()
    scale_columns = [c for c in cluster_means.columns if c != "ClusterID"]

    X = cluster_means_scaled[scale_columns].apply(pd.to_numeric, errors="coerce")
    X = X.replace([np.inf, -np.inf], np.nan)
    all_nan_cols = X.columns[X.isna().all()].tolist()
    if all_nan_cols:
        X = X.drop(columns=all_nan_cols)
        scale_columns = [c for c in scale_columns if c not in all_nan_cols]

    if len(scale_columns) > 0:
        X_filled = X.copy()
        med = X_filled.median(numeric_only=True)
        X_filled = X_filled.fillna(med)
        cluster_means_scaled[scale_columns] = MinMaxScaler().fit_transform(X_filled[scale_columns])
        df_heatmap_scaled = cluster_means_scaled.melt(id_vars="ClusterID", var_name="var", value_name="AU")
        overall_heatmap_data = df_heatmap_scaled.pivot(index="var", columns="ClusterID", values="AU")
    else:
        overall_heatmap_data = pd.DataFrame()

    # ---- Prepare violin plot data ----
    value_cols = [c for c in info_cols if c not in set(sample_cols)]
    df_values = df_umap[["ClusterID"] + value_cols].copy()
    for c in value_cols:
        df_values[c] = pd.to_numeric(df_values[c], errors="coerce")
    df_values.replace([np.inf, -np.inf], np.nan, inplace=True)
    df_long = df_values.melt(id_vars="ClusterID", var_name="var", value_name="value")

    cluster_order = sorted(df_values["ClusterID"].dropna().unique().tolist())
    feat_names = [c for c in value_cols if c in df_long["var"].unique()]
    n_plots = len(feat_names)
    rows_for_plots = (n_plots + nr_cols - 1) // nr_cols
    nr_pages = max(1, (rows_for_plots + rows_per_page - 1) // rows_per_page)

    with PdfPages(outpath) as pdf:
        # ---- Page 1: Full-page heatmap ----
        fig, ax = plt.subplots(figsize=figsize)
        if not overall_heatmap_data.empty:
            try:
                heatmap = sns.heatmap(
                    overall_heatmap_data,
                    ax=ax,
                    cmap="viridis",
                    cbar=True,
                    yticklabels=True
                )
                ax.set_title("Min–Max Scaled Cluster Means", fontsize=14, pad=14)
                ax.set_xlabel("ClusterID", fontsize=10)
                ax.set_ylabel("", fontsize=10)
                ax.tick_params(axis="y", labelsize=6)
                ax.tick_params(axis="x", labelsize=8, rotation=0)
                cbar = heatmap.collections[0].colorbar
                cbar.ax.tick_params(labelsize=8)
                fig.tight_layout(pad=2.0)
            except Exception:
                ax.text(0.5, 0.5, "Overview heatmap unavailable", ha="center", va="center")
                ax.axis("off")
        else:
            ax.text(0.5, 0.5, "No features available for overview scaling", ha="center", va="center")
            ax.axis("off")

        pdf.savefig(fig, dpi=600)
        plt.close(fig)

        # ---- Remaining pages: Violin plots ----
        plot_idx = 0
        for page in range(nr_pages):
            fig = plt.figure(figsize=figsize)
            gs = GridSpec(rows_per_page, nr_cols, figure=fig, hspace=1.5, wspace=0.3)
            remaining_axes = [
                fig.add_subplot(gs[r, c])
                for r in range(rows_per_page)
                for c in range(nr_cols)
            ]

            for ax in remaining_axes:
                if plot_idx >= n_plots:
                    ax.remove()
                    continue

                feat = feat_names[plot_idx]
                sub = df_long.loc[df_long["var"] == feat, ["ClusterID", "value"]].dropna(subset=["ClusterID", "value"])
                if sub.empty:
                    ax.text(0.5, 0.5, f"{feat}\n(no finite data)", ha="center", va="center")
                    ax.axis("off"); plot_idx += 1; continue

                try:
                    sns.violinplot(data=sub, x="ClusterID", y="value", order=cluster_order, inner=None, ax=ax, cut=0)
                except Exception:
                    ax.text(0.5, 0.5, f"{feat}\n(plot unavailable)", ha="center", va="center")
                    ax.axis("off"); plot_idx += 1; continue

                if show_points:
                    sns.stripplot(
                        data=sub,
                        x="ClusterID",
                        y="value",
                        order=cluster_order,
                        ax=ax,
                        dodge=False,
                        jitter=0.2,
                        alpha=point_alpha,
                        size=point_size
                    )

                means = sub.groupby("ClusterID", observed=False)["value"].mean().reindex(cluster_order)
                ax.scatter(
                    np.arange(len(cluster_order)),
                    means.values,
                    s=mean_marker_size,
                    edgecolor="black",
                    linewidths=0.8,
                    zorder=3
                )

                ax.set_title(feat, fontsize=9)
                ax.set_xlabel("ClusterID", fontsize=8)
                ax.set_ylabel("Value", fontsize=8)
                ax.tick_params(axis="x", rotation=0, labelsize=7)
                ax.tick_params(axis="y", labelsize=7)
                plot_idx += 1

            fig.subplots_adjust(left=0.20, right=0.98, top=0.95, bottom=0.08)
            pdf.savefig(fig, dpi=600)
            plt.close(fig)

    if plot_results:
        print(f"Saved PDF to: {outpath}")
        
# def plot_clustering_feature_heatmap(
#     df_umap,
#     info_cols,
#     sample_cols,
#     outpath,
#     rows_per_page = 7,
#     nr_cols = 2,
#     rows_first_img = 4,
#     figsize = (8.27, 11.69),
#     plot_results=True
# ):
#     """
#     Produce a PDF with:
#       • an overview (min-max scaled) heatmap of cluster means (top of first page)
#       • per-feature heatmaps (original value scaling), tiled below & across pages.

#     Handles non-finite values robustly and avoids MinMaxScaler errors by:
#       • coercing to numeric,
#       • replacing ±inf with NaN,
#       • dropping all-NaN columns,
#       • filling remaining NaNs with column median before scaling.
#     """
#     import numpy as np
#     import pandas as pd
#     import seaborn as sns
#     import matplotlib.pyplot as plt
#     from matplotlib.backends.backend_pdf import PdfPages
#     from matplotlib.gridspec import GridSpec
#     from sklearn.preprocessing import MinMaxScaler

#     # Accept iterables for cols
#     info_cols   = list(info_cols) if info_cols is not None else []
#     sample_cols = list(sample_cols) if sample_cols is not None else []

#     # Helper: tolerant tick rounding if round_legend_ticks doesn't exist
#     def _round_legend_ticks(max_val):
#         try:
#             return round_legend_ticks(max_val)  # provided elsewhere in your codebase
#         except Exception:
#             # Simple fallback: round up to 1-2 significant digits
#             if not np.isfinite(max_val) or max_val <= 0:
#                 return 1.0
#             magnitude = 10.0 ** np.floor(np.log10(max_val))
#             return float(np.ceil(max_val / magnitude) * magnitude)

#     # ---- Compute cluster means (numeric only) ----
#     # Drop sample_cols silently if absent, and be explicit about observed
#     df_for_means = (
#         df_umap[list(info_cols) + ["ClusterID"]]
#         .drop(columns=sample_cols, errors="ignore")
#     )

#     # Group & average (numeric columns only). observed=False to silence FutureWarning.
#     cluster_means = (
#         df_for_means
#         .groupby("ClusterID", observed=False)
#         .mean(numeric_only=True)
#         .reset_index()
#     )

#     # Long form for per-feature heatmaps (original values)
#     df_heatmap = cluster_means.melt(id_vars="ClusterID", var_name="var", value_name="value")

#     # ---- Overview heatmap: min-max scale each feature across clusters ----
#     cluster_means_scaled = cluster_means.copy()
#     scale_columns = [c for c in cluster_means.columns if c != "ClusterID"]

#     # Coerce to numeric, clean infinities, drop all-NaN columns, fill NaNs with median
#     X = cluster_means_scaled[scale_columns].apply(pd.to_numeric, errors="coerce")
#     X = X.replace([np.inf, -np.inf], np.nan)

#     # Drop columns that are entirely NaN
#     all_nan_cols = X.columns[X.isna().all()].tolist()
#     if all_nan_cols:
#         # Remove from both X and list of columns to scale
#         X = X.drop(columns=all_nan_cols)
#         scale_columns = [c for c in scale_columns if c not in all_nan_cols]

#     # Only scale if anything remains
#     if len(scale_columns) > 0:
#         # Fill remaining NaNs per-column with median
#         X_filled = X.copy()
#         # Use numeric_only to silence pandas warnings on empty slices
#         med = X_filled.median(numeric_only=True)
#         X_filled = X_filled.fillna(med)
#         # Handle constant columns to avoid 0/0 in MinMaxScaler (it handles it, but keep clean)
#         if len(scale_columns) > 0:
#             cluster_means_scaled[scale_columns] = MinMaxScaler().fit_transform(X_filled[scale_columns])

#         df_heatmap_scaled = cluster_means_scaled.melt(id_vars="ClusterID", var_name="var", value_name="AU")
#         overall_heatmap_data = df_heatmap_scaled.pivot(index="var", columns="ClusterID", values="AU")
#     else:
#         # No valid numeric columns to scale -> create empty frame so we skip overview plot
#         overall_heatmap_data = pd.DataFrame()

#     # ---- Pagination math for per-feature plots ----
#     feat_names = df_heatmap["var"].unique().tolist()
#     n_plots = len(feat_names)

#     # rows needed for all per-feature plots (not counting the overview block)
#     rows_for_plots = (n_plots + nr_cols - 1) // nr_cols
#     # include the reserved rows for the overview on the first page
#     total_rows = rows_first_img + rows_for_plots
#     nr_pages = max(1, (total_rows + rows_per_page - 1) // rows_per_page)

#     with PdfPages(outpath) as pdf:
#         plot_idx = 0  # index into feat_names

#         for page in range(nr_pages):
#             fig = plt.figure(figsize=figsize)
#             gs = GridSpec(rows_per_page, nr_cols, figure=fig, hspace=1.5, wspace=0.3)

#             # ---- First page overview heatmap (if we have anything to show) ----
#             start_row = 0
#             if page == 0 and not overall_heatmap_data.empty:
#                 ax = fig.add_subplot(gs[:rows_first_img, :])
#                 # Guard against all-NaN rows/cols
#                 try:
#                     heatmap = sns.heatmap(overall_heatmap_data, ax=ax, cmap="viridis", cbar=True)
#                     ax.set_title("Min–Max scaled heatmap", fontsize=16)
#                     ax.set_xlabel("ClusterID")
#                     ax.set_ylabel("")
#                     ax.tick_params(axis="y", labelsize=8)
#                     cbar = heatmap.collections[0].colorbar
#                     cbar.ax.tick_params(labelsize=12)
#                 except Exception as _:
#                     ax.text(0.5, 0.5, "Overview heatmap unavailable", ha="center", va="center")
#                     ax.axis("off")
#                 start_row = rows_first_img
#             elif page == 0:
#                 # reserve rows but put a note if no overview
#                 ax = fig.add_subplot(gs[:rows_first_img, :])
#                 ax.text(0.5, 0.5, "No features available for overview scaling", ha="center", va="center")
#                 ax.axis("off")
#                 start_row = rows_first_img

#             # ---- Remaining per-feature heatmaps on this page ----
#             remaining_axes = [
#                 fig.add_subplot(gs[r, c])
#                 for r in range(start_row, rows_per_page)
#                 for c in range(nr_cols)
#             ]

#             for ax in remaining_axes:
#                 if plot_idx >= n_plots:
#                     ax.remove()
#                     continue

#                 col = feat_names[plot_idx]
#                 col_df = df_heatmap[df_heatmap["var"] == col]
#                 col_pivot = col_df.pivot(index="var", columns="ClusterID", values="value")

#                 # Determine color scale: vmin at 0 for readability (as before), vmax rounded
#                 vmax_raw = pd.to_numeric(col_df["value"], errors="coerce").replace([np.inf, -np.inf], np.nan).max()
#                 if not np.isfinite(vmax_raw) or vmax_raw <= 0:
#                     vmin, vmax = 0.0, 1.0
#                 else:
#                     vmin, vmax = 0.0, float(_round_legend_ticks(float(vmax_raw)))

#                 try:
#                     heatmap = sns.heatmap(col_pivot, ax=ax, cmap="viridis", cbar=True, vmin=vmin, vmax=vmax)
#                 except Exception:
#                     # If data is all-NaN or bad, show a placeholder
#                     ax.text(0.5, 0.5, f"{col}\n(no finite data)", ha="center", va="center")
#                     ax.axis("off")
#                     plot_idx += 1
#                     continue

#                 ax.set_title(col)
#                 ax.set_xlabel("ClusterID")
#                 ax.set_ylabel("")
#                 ax.set_yticks([])  # keep compact as in your original
#                 cbar = heatmap.collections[0].colorbar
#                 cbar.ax.tick_params(labelsize=8)
#                 try:
#                     cbar.set_ticks([vmin, (vmin + vmax) / 2.0, vmax])
#                 except Exception:
#                     pass

#                 plot_idx += 1

#             # Layout and save page
#             fig.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.05)
#             pdf.savefig(fig, dpi=600)
#             plt.close(fig)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Input parameters for automatic data transfer.')
    parser.add_argument('-c', '--config', type=str, help='path to a config.yml file that stores all required paths', required=False)
    args = parser.parse_args()
    with open(args.config, "r") as parameters:
        config=yaml.load(parameters, Loader=yaml.SafeLoader)
    # metadata = pd.read_csv(config["metadata_csv_path"])
    run_tcell_analysis(config)
