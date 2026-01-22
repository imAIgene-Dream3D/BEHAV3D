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
    plot_filter_count, 
    filter_by_full_duration, 
    filter_minimal_track_length, 
    trim_to_maximal_track_length,
    plot_dead_dye_distribution,
    plot_touching_nontouching_distribution,
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
    seed=42
    ):
    print(f"--------------- Performing {cell_type} behavioral analysis ---------------")
    start_time = time.time()
    # assert config is not None or all(
    #     [output_dir, umap_minimal_distance, umap_n_neighbors, nr_of_clusters]
    # ), "Either 'config' or 'output_dir, umap_minimal_distance, umap_n_neighbors, nr_of_clusters' parameters must be supplied"
        
    if output_dir is None:
        output_dir = config['output_dir']
        
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
    df_tracks_summarized = pd.read_csv(df_tracks_summarized_path, low_memory=False)
    df_tracks=df_tracks.sort_values(by=["sample_name", "TrackID", "relative_time"])
    
    if df_tracks_summarized["track_length"].nunique() != 1:
        print("Warning: The track lengths are not cut to similar length, this might influence dynamic time warping")
        print("Set 'min_track_length' and 'max_track_length' to the same value to create equal tracks")

    non_wildcard_cols_to_normalize = []
    for col in columns_to_normalize:
        non_wildcard_cols_to_normalize += expand_column_patterns(col, df_tracks.columns)
    # print(f"{get_current_time()} - Perform z-normalization on selected feature columns")
    df_tracks=normalize_track_features(
        df_tracks, 
        columns=non_wildcard_cols_to_normalize
    )
    
    dtw_features = [x if x not in columns_to_normalize else f"z_{x}" for x in columns_to_use]
    
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
        random_state=seed
    )
    end_time = time.time()
    h,m,s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    return(df_clusters)

def filter_cell_tracks(
    metadata,
    config=None,
    output_dir=None,
    exp_duration=None,
    min_track_length=None,
    max_track_length=None,
    filter_t0_dead=True,
    cell_type="tcell",
    time_type="frames", #can also be "relative_time"
    plot_results=True,
    df_input_path=None  # Optional: path to input CSV (e.g., advanced features CSV)
    ):
    """
    This code filters tracks based on supplied parameters in the config.yml
    
    Filtering is based on:
    - Maximum experiment length (exp_duration)    
    - Minimum track length (min_track_length)
    - Tracks starting at timepoint 1 with a dead dye mean over the dead_dye_threshold (dead_dye_threshold)
    
    Additonally, all tracks are cut down to:
    - Maximum track length (max_track_length)
    
    Parameters
    ----------
    df_input_path : Path or str, optional
        Path to input CSV file. If provided, reads from this file instead of the default
        combined_track_features.csv. Useful for reading advanced features CSV with active killing data.
    
    Output:
    - A .csv file containing filtered tracks from all experiments
    """
    
    # assert config is not None or all(
    #     [output_dir, exp_duration, min_track_length, max_track_length]
    # ), "Either 'config' or 'output_dir, tcell_exp_duration, tcell_min_track_length and tcell_max_track_length' parameters must be supplied"
    
    start_time = time.time()
    
    print(f"--------------- Filtering tracks ---------------")
    print(f"Basing filtering on {time_type}")
    # if not all([output_dir, exp_duration, min_track_length, max_track_length]):
    #     output_dir = config['output_dir']
    #     exp_duration = config[f'{cell_type}_exp_duration']
    #     min_track_length = config[f'{cell_type}_min_track_length']
    #     max_track_length = config[f'{cell_type}_max_track_length']
    
    analysis_outdir = Path(output_dir, "analysis", cell_type)
    feature_outdir = Path(analysis_outdir, "track_features")
    qc_outdir = Path(analysis_outdir, "quality_control")

    if not analysis_outdir.exists():
        analysis_outdir.mkdir(parents=True)
    if not feature_outdir.exists():
        feature_outdir.mkdir(parents=True)
    if not qc_outdir.exists():
        qc_outdir.mkdir(parents=True)
    
    # Use provided input path or default to combined_track_features.csv
    if df_input_path is not None:
        df_all_tracks_path = Path(df_input_path)
        print(f"  Using input file: {df_all_tracks_path.name}")
    else:
        df_all_tracks_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features.csv")
    
    if not df_all_tracks_path.exists():
        raise FileNotFoundError(
            f"Track features file not found: {df_all_tracks_path}\n"
            f"Please run Feature Extraction for {cell_type} first."
        )
    
    df_all_tracks = pd.read_csv(df_all_tracks_path)
    
    df_all_tracks['sample_name'] = df_all_tracks['sample_name'].astype(str)
    metadata['sample_name'] = metadata['sample_name'].astype(str)

    # Dynamically detect *_line_condition columns from metadata
    line_condition_cols = [c for c in metadata.columns if c.endswith('_line_condition')]
    group_cols = ['TrackID', 'sample_name'] + line_condition_cols + ['exp_nr', 'well']

    cols_present = [c for c in group_cols if c in metadata.columns]
    metadata_info = metadata.loc[:, cols_present].copy()
    df_all_tracks_filt = pd.merge(df_all_tracks, metadata_info, how="left", on="sample_name")
    
    # Function to count the number of unique tracks in the DataFrame
    def count_tracks(df_all_tracks, col_name="nr_tracks", df_track_counts=None):
        # Dynamically detect which grouping columns are present (including *_line_condition columns)
        line_cols_in_tracks = [c for c in df_all_tracks.columns if c.endswith('_line_condition')]
        potential_group_cols = ['sample_name'] + line_cols_in_tracks + ['exp_nr', 'well']
        group_cols_for_count = [c for c in potential_group_cols if c in df_all_tracks.columns]
        
        nr_tracks=df_all_tracks.groupby(group_cols_for_count).agg(
            nr_tracks=pd.NamedAgg(column='TrackID', aggfunc='nunique')
        ).reset_index()
        nr_tracks=nr_tracks.rename(columns={"nr_tracks":col_name})
        if df_track_counts is None:
            return(nr_tracks)
        else:
            return(pd.merge(df_track_counts, nr_tracks, how="left"))
    
    # Counting the nr of tracks before filtering
    df_track_counts=count_tracks(df_all_tracks_filt, col_name="nr_tracks_before_filtering")
    
    # Filtering the tracks based on the total experimental duration
    # Any timepoint after this will be filtered out 
    if time_type=="real_time" or time_type=="hours":
        time_column="time"  # 'time' column contains actual hours
        # Don't subtract 1 for real time units (hours)
    else:  # frames
        time_column="position_t"
        # Only subtract 1 if values are not None (frames are 0-indexed)
        if exp_duration is not None:
            exp_duration = exp_duration - 1
        if min_track_length is not None:
            min_track_length = min_track_length - 1
        if max_track_length is not None:
            max_track_length = max_track_length - 1
    
    df_all_tracks_filt = df_all_tracks_filt.sort_values(group_cols + [time_column]).reset_index(drop=True)

    df_all_tracks_filt=filter_by_full_duration(
        df=df_all_tracks_filt,
        time_column=time_column,
        exp_duration=exp_duration
        )
    
    df_track_counts=count_tracks(
        df_all_tracks_filt, 
        col_name="nr_tracks_exp_duration", 
        df_track_counts=df_track_counts
        )

    # Filtering out tracks under specific track length
    df_all_tracks_filt = filter_minimal_track_length(
        df=df_all_tracks_filt,
        min_track_length=min_track_length,
        time_column=time_column,
        group_cols=group_cols
    )

    # Cutting down tracks to maximal track length    
    df_all_tracks_filt = trim_to_maximal_track_length(
        df=df_all_tracks_filt,
        max_track_length=max_track_length,
        time_column=time_column,
        group_cols=group_cols
    )
    
    df_track_counts = count_tracks(
            df_all_tracks_filt, 
            col_name="nr_tracks_min_track_length", 
            df_track_counts=df_track_counts
        )
    
    # Plot the number of cells having contact with other cell types
    # Dynamically detect all *_contact columns
    contact_columns = [col for col in df_all_tracks_filt.columns if col.endswith('_contact') and not col.startswith('active_')]
    for contact_col in contact_columns:
        # Extract the target cell type from column name (e.g., "organoid1_contact" -> "organoid1")
        target_type = contact_col.replace('_contact', '')
        plot_touching_outpath = Path(qc_outdir, f"BEHAV3D_{target_type}_touching_distribution.pdf")
        print(f"- Plotting {target_type} touching distribution to {plot_touching_outpath}")
        plot_touching_nontouching_distribution(
            df_all_tracks_filt, 
            outpath=plot_touching_outpath,
            contact_column=contact_col,
            nr_cols=3,
            rows_per_page=3,
            )
    
    filter_cols=["nr_tracks_before_filtering", "nr_tracks_exp_duration", "nr_tracks_min_track_length"]
    if "mean_dead_dye" in df_all_tracks_filt.columns:
        # Plot the distribution of dead dye intensity of all timepoints and at timepoint 1
        # Can be used to aid in the choice of dead_dye_threshold
        plot_dead_dye_distr_outpath = Path(qc_outdir, f"BEHAV3D_dead_dye_distribution.pdf")
        print(f"- Plotting dead dye distribution at all timepoints to {plot_dead_dye_distr_outpath}")
        plot_dead_dye_distribution(
            df_all_tracks_filt,
            outpath=plot_dead_dye_distr_outpath,
            nr_cols=2,
            rows_per_page=2
            )
    
        plot_dead_dye_distr_t0_outpath = Path(qc_outdir, f"BEHAV3D_dead_dye_distribution_t0.pdf")
        print(f"- Plotting dead dye distribution at timepoint 1 to {plot_dead_dye_distr_outpath}")
        plot_dead_dye_distribution(
            df_all_tracks_filt[df_all_tracks_filt["relative_time"]==1],
            outpath=plot_dead_dye_distr_t0_outpath,
            nr_cols=2,
            rows_per_page=2
            )
    
        # Filter out all T cells that are dead based on the threshold at the first timepoint of a track
        if filter_t0_dead:
            assert 'dead' in df_all_tracks_filt.columns, "The column 'dead' is not present in the DataFrame, but filter_t0_dead is supplied"
            dead_t0 = df_all_tracks_filt[
                (df_all_tracks_filt["relative_time"]==1) & 
                (df_all_tracks_filt["dead"])
                ][["TrackID","sample_name"]]
            df_all_tracks_filt=df_all_tracks_filt[~df_all_tracks_filt.set_index(['TrackID', 'sample_name']).index.isin(dead_t0.set_index(['TrackID', 'sample_name']).index)]
            filter_cols.append("nr_tracks_dead_t1")   
        df_track_counts=count_tracks(df_all_tracks_filt, col_name="nr_tracks_dead_t1", df_track_counts=df_track_counts)

    plot_filter_count_outpath = Path(qc_outdir, f"BEHAV3D_filter_counts.pdf")
    print(f"- Plotting track counts after filtering steps to {plot_filter_count_outpath}")
    plot_filter_count(
        df_track_counts,
        outpath=plot_filter_count_outpath,
        nr_cols=3,
        rows_per_page = 3,
        filter_cols=filter_cols,
        plot_results=plot_results
    )
    
    # Write the filtered tracks to a .csv
    filt_tracks_out_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features_filtered.csv")
    print(f"- Writing filtered tracks to {filt_tracks_out_path}")
    df_all_tracks_filt = df_all_tracks_filt.sort_values(
        by=["sample_name", "TrackID", "position_t"], 
        ascending=[True, True, True]  # example: last column sorted descending
    )    
    new_order = ["sample_name"] + [c for c in df_all_tracks_filt.columns.tolist() if c != "sample_name"]
    df_all_tracks_filt = df_all_tracks_filt[new_order]
    
    df_all_tracks_filt.to_csv(filt_tracks_out_path, sep=",", index=False)
    end_time = time.time()
    h,m,s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    return(df_all_tracks_filt)

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

def rolling_classification(
    df_tracks,
    window_size=20,
    groupby=["sample_name", "TrackID"],
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
    ):
    df_tracks = df_tracks.sort_values(by=groupby + ["position_t"])

    # df_rolling = df_tracks.groupby(groupby)[features].rolling(window=window_size, min_periods=window_size).mean().reset_index().drop(columns='level_2')
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
    df_rolling["position_t"] = df_tracks["position_t"].values
    # df_rolling = pd.concat([df_tracks[['sample_name', 'TrackID', 'position_t']], df_rolling], axis=1)

    
    # df_rolling = df_tracks[features].rolling(window=window_size, min_periods=window_size).mean()
    df_rolling = df_rolling.dropna()
    
    df_rolling_features = df_rolling[features]
    
    scaler = StandardScaler()
    scaler = RobustScaler()
    df_rolling_features = pd.DataFrame(scaler.fit_transform(df_rolling_features), columns=df_rolling_features.columns)

    n_components = 2
    # Step 3: Dimensionality reductions
    pca_model = PCA(n_components=n_components)
    pca_result = pca_model.fit_transform(df_rolling_features)

    tsne_model = TSNE(n_components=n_components, random_state=123, perplexity=30)
    tsne_result = tsne_model.fit_transform(df_rolling_features)

    umap_model = umap.UMAP(n_components=n_components, n_neighbors=15, min_dist=0.1, init="random", random_state=123)
    umap_result = umap_model.fit_transform(df_rolling_features)
    
    # Step 4: Plot all three embeddings side by side
    fig, axs = plt.subplots(1, 3, figsize=(18, 5))
    embeddings = [pca_result, umap_result]
    titles = ["PCA", "UMAP"]
    
    if n_components == 2:
        for ax, emb, title in zip(axs, embeddings, titles):
            scatter = ax.scatter(emb[:, 0], emb[:, 1], c=np.arange(len(emb)), cmap="viridis", s=2, alpha=0.7)
            ax.set_title(title)
            ax.set_xlabel(f'{title} 1')
            ax.set_ylabel(f'{title} 2')
        
        plt.tight_layout()
        plt.show()
    else:
        for i, (emb, title) in enumerate(zip(embeddings, titles), 1):
            ax = fig.add_subplot(1, 3, i, projection='3d')
            scatter = ax.scatter(emb[:, 0], emb[:, 1], emb[:, 2], c=np.arange(len(emb)), cmap="viridis", s=2, alpha=0.7)
            ax.set_title(title)
            ax.set_xlabel(f'{title} 1')
            ax.set_ylabel(f'{title} 2')
            ax.set_zlabel(f'{title} 3')

        plt.tight_layout()
        plt.show()

    n_components=2
    n_neighbors=15
    min_dist=0.1
    
    umap_model = umap.UMAP(
        n_components=n_components, 
        n_neighbors=n_neighbors, 
        min_dist=min_dist, 
        init="random", 
        random_state=123,
        metric="euclidean", 
        )
    
    pca = PCA(n_components=min(10, len(features)))
    # pca = PCA(n_components=0.95)
    X_pca = pca.fit_transform(df_rolling_features[features])
    umap_result = umap_model.fit_transform(X_pca)
    
    umap_result = umap_model.fit_transform(df_rolling_features)
    
    fig = plt.figure()
    if n_components == 1:
        ax = fig.add_subplot(111)
        ax.scatter(umap_result[:,0], range(len(umap_result)), c=np.arange(len(umap_result)))
    if n_components == 2:
        ax = fig.add_subplot(111)
        ax.scatter(umap_result[:,0], umap_result[:,1], c=np.arange(len(umap_result)), s=1, alpha=0.5)
    if n_components == 3:
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(umap_result[:,0], umap_result[:,1], umap_result[:,2], c=np.arange(len(umap_result)), s=1, alpha=0.5)

    clusterer = HDBSCAN(
        min_cluster_size=100,         # Minimum size of a cluster
        # min_samples=20,               # Controls outlier sensitivity (higher = stricter)
        # cluster_selection_epsilon=0.0,# Optional: smooths cluster boundaries
        # cluster_selection_method='eom',  # Default is good; use 'leaf' if you want more granularity
        )
    cluster_labels = clusterer.fit_predict(umap_result)

    ### TODO Do the clustering based straight on the DTW distances
    ### Miguel suggested clusters are more meaningful before the umap embedding
    df_rolling_features['cluster'] = cluster_labels

    # Plot UMAP colored by HDBSCAN clusters
    plt.figure(figsize=(10, 8))
    unique_labels = np.unique(cluster_labels)

    # Assign a color for each cluster (noise will be colored gray)
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_labels)))
    for i, label in enumerate(unique_labels):
        mask = cluster_labels == label
        plt.scatter(
            umap_result[mask, 0],
            umap_result[mask, 1],
            s=5,
            color='gray' if label == -1 else colors[i],
            label=f'Noise' if label == -1 else f'Cluster {label}',
            alpha=0.8
        )

    plt.title('UMAP Projection with HDBSCAN Clusters')
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    plt.legend(markerscale=4, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()
    
    # Plot UMAP colored by each selected feature
    for feature in df_rolling_features.columns:
        plt.figure(figsize=(8, 6))
        plt.scatter(
            umap_result[:, 0], 
            umap_result[:, 1], 
            c=df_rolling_features[feature], 
            # cmap='viridis', 
            s=5, 
            alpha=0.8
        )
        plt.colorbar(label=feature)
        plt.title(f'UMAP colored by {feature}')
        plt.xlabel('UMAP 1')
        plt.ylabel('UMAP 2')
        plt.tight_layout()
        plt.show()
    
    labels = KMeans(n_clusters=5).fit_predict(df_rolling_features)  # on raw 2 features
    plt.scatter(df_rolling_features.iloc[:,0], df_rolling_features.iloc[:,1], c=labels, cmap="Spectral", s=2)
    
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
    plot_results=True
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
    results_outdir = Path(analysis_outdir, "results")
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
    info_cols = df_umap.drop(columns=["TrackID", "UMAP1", "UMAP2", "ClusterID"]).columns
    
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
    extra_group_cols = []
    if cluster_percentage_group_by:
        for col in cluster_percentage_group_by:
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
    df_for_means = (
        df_umap[list(info_cols) + ["ClusterID"]]
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