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
from behav3d.utils import format_time, expand_column_patterns
from behav3d.utils.filtering import plot_filter_count, filter_by_full_duration, filter_minimal_track_length, trim_to_maximal_track_length
import yaml
import time
import seaborn as sns
# df_tracks=df_tracks[df_tracks["relative_time"]<=30]

def run_tcell_analysis(
    config=None,
    output_dir=None,
    df_tracks_path=None,
    df_tracks_summarized_path=None,
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
    plot_results=True,
    seed=42
    ):
    print(f"--------------- Performing T-cell behavioral analysis ---------------")
    start_time = time.time()
    # assert config is not None or all(
    #     [output_dir, umap_minimal_distance, umap_n_neighbors, nr_of_clusters]
    # ), "Either 'config' or 'output_dir, umap_minimal_distance, umap_n_neighbors, nr_of_clusters' parameters must be supplied"
        
    if output_dir is None:
        output_dir = config['output_dir']
        
    analysis_outdir = Path(output_dir, "analysis", "tcell")
    feature_outdir = Path(analysis_outdir, "track_features")
    
    if not analysis_outdir.exists():
        analysis_outdir.mkdir(parents=True)
    if not feature_outdir.exists():
        feature_outdir.mkdir(parents=True)    
    
    if df_tracks_path is None:
        df_tracks_path = Path(feature_outdir, f"BEHAV3D_tcell_combined_track_features_filtered.csv")
    if df_tracks_summarized_path is None:
        df_tracks_summarized_path = Path(feature_outdir, f"BEHAV3D_tcell_combined_track_features_summarized.csv")
    
    df_tracks = pd.read_csv(df_tracks_path, low_memory=False)
    df_tracks_summarized = pd.read_csv(df_tracks_summarized_path, low_memory=False)
    df_tracks=df_tracks.sort_values(by=["sample_name", "TrackID", "relative_time"])
    
    if df_tracks_summarized["track_length"].nunique() != 1:
        print("Warning: The track lengths are not cut to similar length, this might influence dynamic time warping")
        print("Set 'tcell_min_track_length' and 'tcell_max_track_length' to the same value to create equal tracks")

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
        plot_results=plot_results,
        random_state=seed
    )
    end_time = time.time()
    h,m,s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    return(df_clusters)

def filter_tcell_tracks(
    metadata,
    config=None,
    output_dir=None,
    exp_duration=None,
    min_track_length=None,
    max_track_length=None,
    filter_t0_dead=True,
    cell_type="tcell",
    time_type="frames", #can also be "relative_time"
    plot_results=True
    ):
    """
    This code filters tracks based on supplied parameters in the config.yml
    
    Filtering is based on:
    - Maximum experiment length (tcell_exp_duration)    
    - Minimum track length (tcell_min_track_length)
    - Tracks starting at timepoint 1 with a dead dye mean over the dead_dye_threshold (dead_dye_threshold)
    
    Additonally, all tracks are cut down to:
    - Maximum track length (tcell_max_track_length)
    
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
    
    df_all_tracks_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features.csv")
    df_all_tracks = pd.read_csv(df_all_tracks_path)
    
    df_all_tracks['sample_name'] = df_all_tracks['sample_name'].astype(str)
    metadata['sample_name'] = metadata['sample_name'].astype(str)

    group_cols = ['TrackID', 'sample_name', 'organoid_line', 'tcell_line', 'exp_nr', 'well']

    cols_present = [c for c in group_cols if c in metadata.columns]
    metadata_info = metadata.loc[:, cols_present].copy()
    df_all_tracks_filt = pd.merge(df_all_tracks, metadata_info, how="left", on="sample_name")
    
    # Function to count the number of unique tracks in the DataFrame
    def count_tracks(df_all_tracks, col_name="nr_tracks", df_track_counts=None):
        nr_tracks=df_all_tracks.groupby([
            'sample_name', 'organoid_line', 'tcell_line', 'exp_nr', 'well']
            ).agg(nr_tracks=pd.NamedAgg(column='TrackID', aggfunc='nunique')).reset_index()
        nr_tracks=nr_tracks.rename(columns={"nr_tracks":col_name})
        if df_track_counts is None:
            return(nr_tracks)
        else:
            return(pd.merge(df_track_counts, nr_tracks, how="left"))
    
    # Counting the nr of tracks before filtering
    df_track_counts=count_tracks(df_all_tracks_filt, col_name="nr_tracks_before_filtering")
    
    # Filtering the tracks based on the total experimental duration
    # Any timepoint after this will be filtered out 
    if time_type=="real_time":
        time_column="time"
    else:
        time_column="position_t"
        exp_duration = exp_duration-1
        min_track_length = min_track_length-1
        max_track_length = max_track_length-1
    
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
    
    # Plot the number of cells having contact with another T cell and Organoid for analysis
    # of the set contact_threshold
    if "tcell_contact" in df_all_tracks_filt.columns:
        plot_tcell_touching_outpath = Path(qc_outdir, f"BEHAV3D_tcell_touching_distribution.pdf")
        print(f"- Plotting tcell touching distribution to {plot_tcell_touching_outpath}")
        plot_touching_nontouching_distribution(
            df_all_tracks_filt, 
            outpath=plot_tcell_touching_outpath,
            contact_column="tcell_contact",
            nr_cols=3,
            rows_per_page = 3,
            )
    
    if "organoid_contact" in df_all_tracks_filt.columns:
        plot_organoid_touching_outpath = Path(qc_outdir, f"BEHAV3D_organoid_touching_distribution.pdf")
        print(f"- Plotting organoid touching distribution to {plot_organoid_touching_outpath}")
        plot_touching_nontouching_distribution(
            df_all_tracks_filt, 
            outpath=plot_organoid_touching_outpath,
            contact_column="organoid_contact",
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

    # dtw_input_tracks = np.empty((nr_tracks, nr_timepoints, nr_features),dtype=np.double)
    
    dtw_input_tracks=[]
    dtw_rownames=[]
    unique_tracks = df_tracks.groupby(['TrackID', 'sample_name'])
    for (TrackID, sample_name), group in unique_tracks:
        track_features = group[features].to_numpy().astype(np.double)
        dtw_rownames.append(f"{TrackID}--{sample_name}")
        # dtw_track = {
        #     "TrackID": TrackID,
        #     "sample_name": sample_name,
        #     "dtw_input": track_features
        # }
        # dtw_track = pd.DataFrame(dtw_track)
        dtw_input_tracks.append(track_features)
    
    dtw_distance_matrix = dtw_ndim.distance_matrix_fast(dtw_input_tracks)
    dtw_distance_matrix = pd.DataFrame(dtw_distance_matrix, index=dtw_rownames, columns=dtw_rownames)
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
        
    umap_model = umap.UMAP(
        n_components=2, 
        n_neighbors=umap_n_neighbors, 
        min_dist=umap_minimal_distance, 
        init="random", 
        random_state=random_state,
        metric="precomputed", 
        )
    umap_embedding = umap_model.fit_transform(dtw_distance_matrix.values)
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
    plot_results=True
    ):
    
    assert config is not None or all(
        [output_dir, nr_of_clusters]
    ), "Either 'config' or 'output_dir, nr_of_clusters' parameters must be supplied"
      
    print("- Performing clustering on the UMAP data")
    if all([output_dir, nr_of_clusters]) is None:
        output_dir = Path(config['output_dir'])
        nr_of_clusters=config["nr_of_clusters"]
    
    tcell_outdir = Path(output_dir, "analysis", "tcell")
    feature_outdir = Path(tcell_outdir, "track_features")
    results_outdir = Path(tcell_outdir, "results")
    if not tcell_outdir.exists():
        tcell_outdir.mkdir(parents=True)
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
    
    df_umap_out_path = Path(results_outdir, f"BEHAV3D_tcell_UMAP_clusters.csv")
    print(f"- Writing clustered tracks to {df_umap_out_path}")
    df_umap.to_csv(df_umap_out_path, sep=",", index=False)

    sample_cols = ["organoid_line", "tcell_line"]
    info_cols = df_umap.drop(columns=["TrackID", "sample_name", "well", "exp_nr", "UMAP1", "UMAP2", "ClusterID"]).columns
    
    cluster_UMAP_path = Path(results_outdir, f"BEHAV3D_tcell_UMAP_clusters.pdf")
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
    cluster_features_heatmap_path = Path(results_outdir, f"BEHAV3D_tcell_UMAP_cluster_feature_heatmap.pdf")
    print(f"- Plotting heatmaps with summarized cluster features to {cluster_features_heatmap_path}")
    plot_clustering_feature_heatmap(
        df_umap,
        info_cols,
        sample_cols,
        cluster_features_heatmap_path,
        rows_per_page = 7,
        nr_cols = 2,
        rows_first_img = 4,
        figsize = (8.27, 11.69),
        plot_results=plot_results
    )

    df_clust_perc = df_umap.groupby(["organoid_line", "tcell_line", "ClusterID"]).size().reset_index(name='count')
    total_counts = df_clust_perc.groupby(['organoid_line', 'tcell_line'])['count'].sum().reset_index(name='total_count')
    
    df_clust_perc = pd.merge(df_clust_perc, total_counts)
    df_clust_perc["percentage"] = (df_clust_perc['count'] / df_clust_perc['total_count'])
    
    cluster_percentage_plot_path = Path(results_outdir, f"BEHAV3D_tcell_UMAP_cluster_percentages.pdf")
    print(f"- Plotting percentage plots of each cluster per combination of T-cell and organoid line to {cluster_percentage_plot_path}")
    plot_cluster_percentage_bars(
        df_clust_perc,
        cluster_percentage_plot_path,
        plot_results=plot_results
    )
        
    df_clust_perc = df_clust_perc.reset_index(drop=True)
    df_clust_perc_out_path = Path(results_outdir, f"BEHAV3D_tcell_UMAP_cluster_percentages.csv")
    print(f"- Writing cluster percentages to {df_clust_perc_out_path}")
    df_clust_perc.to_csv(df_clust_perc_out_path, sep=",", index=False)
    
    df_clust_tracks_out_path = Path(results_outdir, f"BEHAV3D_tcell_combined_track_features_clustered.csv")
    print(f"- Writing clustered track info to {df_clust_perc_out_path}")
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
        df_tracks.loc[:, f'z_{column_name}'] = df_tracks[column_name].transform(lambda x: (x - x.mean()) / x.std())
    return (df_tracks)

def plot_cluster_percentage_bars(
    df_clust_perc,
    outpath,
    plot_results=True
    ):
    with PdfPages(outpath) as pdf:
        tcell_lines = df_clust_perc['tcell_line'].unique()
        organoid_lines = df_clust_perc['organoid_line'].unique()
        
        fig, axes = plt.subplots(
            len(tcell_lines), 
            len(organoid_lines), 
            figsize=(20, 10), sharex=True, sharey=True)

        axes = np.atleast_2d(axes)
        # Plot horizontal stacked bar charts
        for i, tcell_line in enumerate(tcell_lines):
            for j, organoid_line in enumerate(organoid_lines):
                ax = axes[i, j]
                subset = df_clust_perc[(df_clust_perc['tcell_line'] == tcell_line) & (df_clust_perc['organoid_line'] == organoid_line)]
                if i == 0:
                    ax.set_title(f'{organoid_line}', fontsize=30)
                if j == 0:
                    ax.set_ylabel(f'{tcell_line}', fontsize=30) 
                if subset.empty:
                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)
                    ax.spines['left'].set_visible(False)
                    ax.spines['bottom'].set_visible(False)
                    continue
                subset_pivot = subset.pivot(index='tcell_line', columns='ClusterID', values='percentage').fillna(0)
                subset_pivot.plot(kind='barh', stacked=True, ax=ax, legend=False)
                if i == 0:
                    ax.set_title(f'{organoid_line}', fontsize=30)
                if j == 0:
                    print(tcell_line)
                    ax.set_ylabel(f'{tcell_line}', fontsize=30) 
    
                num_cells = subset['count'].sum()
                ax.text(
                    0.5, 
                    0.1, 
                    f'# Cells: {num_cells}', 
                    ha='center', 
                    va='center', 
                    transform=ax.transAxes, 
                    fontsize=20
                    )
                
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['left'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                
                # ax.set_xlabel('Percentage')

         # Add legend
        handles, labels = ax.get_legend_handles_labels()
        fig.legend(
            handles, 
            labels, 
            fontsize=30,
            title='ClusterID', 
            title_fontsize=30,
            bbox_to_anchor=(0.9, 0.5), 
            loc='center left')
        fig.tight_layout(rect=[0, 0, 0.85, 1])
        if plot_results:
            plt.show()
        pdf.savefig(fig, bbox_inches='tight', dpi=600)
        plt.close(fig)
            
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
                    s=40,
                    alpha=0.95,
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
                        s=40,
                        hue=colorcol,
                        alpha=0.8,
                        palette="Set2",
                        ax=ax
                    )
                else:
                    sns.scatterplot(
                        data=df_umap,
                        x="UMAP1",
                        y="UMAP2",
                        s=40,
                        hue=colorcol,
                        palette="viridis",
                        alpha=0.6,
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
    rows_first_img = 4,
    figsize = (8.27, 11.69),
    plot_results=True,
    show_points=False,        # overlay individual samples
    point_alpha=0.5,         # transparency for individual samples
    point_size=8,            # size for individual points
    mean_marker_size=60,     # size for mean markers
):
    """
    Produce a PDF with:
      • an overview (min-max scaled) heatmap of cluster means (top of first page)
      • per-feature VIOLIN plots (original value scaling) tiled below & across pages,
        showing the distribution of every sample in each cluster and the cluster mean.

    Robust to non-finite values by:
      • coercing to numeric,
      • replacing ±inf with NaN,
      • dropping all-NaN columns for the overview scaler,
      • filling remaining NaNs with column median before scaling (overview only).

    Parameters added:
      • show_points (bool): overlay each sample as jittered points on top of violins.
      • point_alpha (float), point_size (int): styling for sample points.
      • mean_marker_size (int): styling for the mean marker per cluster.
    """
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.gridspec import GridSpec
    from sklearn.preprocessing import MinMaxScaler

    # Accept iterables for cols
    info_cols   = list(info_cols) if info_cols is not None else []
    sample_cols = list(sample_cols) if sample_cols is not None else []

    # Helper: tolerant tick rounding if round_legend_ticks doesn't exist
    def _round_legend_ticks(max_val):
        try:
            return round_legend_ticks(max_val)  # provided elsewhere in your codebase
        except Exception:
            # Simple fallback: round up to 1-2 significant digits
            if not np.isfinite(max_val) or max_val <= 0:
                return 1.0
            magnitude = 10.0 ** np.floor(np.log10(max_val))
            return float(np.ceil(max_val / magnitude) * magnitude)

    # ---- Compute cluster means (numeric only) for overview ----
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

    # ---- Overview heatmap: min-max scale each feature across clusters ----
    cluster_means_scaled = cluster_means.copy()
    scale_columns = [c for c in cluster_means.columns if c != "ClusterID"]

    # Coerce to numeric, clean infinities, drop all-NaN columns, fill NaNs with median
    X = cluster_means_scaled[scale_columns].apply(pd.to_numeric, errors="coerce")
    X = X.replace([np.inf, -np.inf], np.nan)

    # Drop columns that are entirely NaN
    all_nan_cols = X.columns[X.isna().all()].tolist()
    if all_nan_cols:
        X = X.drop(columns=all_nan_cols)
        scale_columns = [c for c in scale_columns if c not in all_nan_cols]

    # Build overview pivot if we have anything left to scale
    if len(scale_columns) > 0:
        X_filled = X.copy()
        med = X_filled.median(numeric_only=True)
        X_filled = X_filled.fillna(med)
        cluster_means_scaled[scale_columns] = MinMaxScaler().fit_transform(X_filled[scale_columns])

        df_heatmap_scaled = cluster_means_scaled.melt(id_vars="ClusterID", var_name="var", value_name="AU")
        overall_heatmap_data = df_heatmap_scaled.pivot(index="var", columns="ClusterID", values="AU")
    else:
        overall_heatmap_data = pd.DataFrame()

    # ---- Long-form per-sample data for violin plots ----
    # Take original per-sample values, coerce numeric, keep ClusterID
    value_cols = [c for c in info_cols if c not in set(sample_cols)]
    df_values = df_umap[["ClusterID"] + value_cols].copy()

    # Clean values
    for c in value_cols:
        df_values[c] = pd.to_numeric(df_values[c], errors="coerce")
    df_values.replace([np.inf, -np.inf], np.nan, inplace=True)

    # Melt to long form: one row per sample-feature
    df_long = df_values.melt(id_vars="ClusterID", var_name="var", value_name="value")
    # Keep cluster ordering stable and explicit
    cluster_order = sorted(df_values["ClusterID"].dropna().unique().tolist())

    # Features list (preserve input order from info_cols but only those present)
    feat_names = [c for c in value_cols if c in df_long["var"].unique()]
    n_plots = len(feat_names)

    # ---- Pagination math for per-feature plots ----
    rows_for_plots = (n_plots + nr_cols - 1) // nr_cols
    total_rows     = rows_first_img + rows_for_plots
    nr_pages       = max(1, (total_rows + rows_per_page - 1) // rows_per_page)

    with PdfPages(outpath) as pdf:
        plot_idx = 0  # index into feat_names

        for page in range(nr_pages):
            fig = plt.figure(figsize=figsize)
            gs = GridSpec(rows_per_page, nr_cols, figure=fig, hspace=1.5, wspace=0.3)

            # ---- First page overview heatmap (if available) ----
            start_row = 0
            if page == 0 and not overall_heatmap_data.empty:
                ax = fig.add_subplot(gs[:rows_first_img, :])
                try:
                    heatmap = sns.heatmap(overall_heatmap_data, ax=ax, cmap="viridis", cbar=True)
                    ax.set_title("Min–Max scaled heatmap", fontsize=16)
                    ax.set_xlabel("ClusterID")
                    ax.set_ylabel("")
                    ax.tick_params(axis="y", labelsize=8)
                    cbar = heatmap.collections[0].colorbar
                    cbar.ax.tick_params(labelsize=12)
                except Exception:
                    ax.text(0.5, 0.5, "Overview heatmap unavailable", ha="center", va="center")
                    ax.axis("off")
                start_row = rows_first_img
            elif page == 0:
                ax = fig.add_subplot(gs[:rows_first_img, :])
                ax.text(0.5, 0.5, "No features available for overview scaling", ha="center", va="center")
                ax.axis("off")
                start_row = rows_first_img

            # ---- Remaining per-feature VIOLIN plots on this page ----
            remaining_axes = [
                fig.add_subplot(gs[r, c])
                for r in range(start_row, rows_per_page)
                for c in range(nr_cols)
            ]

            for ax in remaining_axes:
                if plot_idx >= n_plots:
                    ax.remove()
                    continue

                feat = feat_names[plot_idx]
                sub = df_long.loc[df_long["var"] == feat, ["ClusterID", "value"]].copy()

                # Drop rows with no cluster or value
                sub = sub.dropna(subset=["ClusterID", "value"])
                if sub.empty:
                    ax.text(0.5, 0.5, f"{feat}\n(no finite data)", ha="center", va="center")
                    ax.axis("off")
                    plot_idx += 1
                    continue

                # Violin plot per cluster (distribution)
                try:
                    sns.violinplot(
                        data=sub,
                        x="ClusterID",
                        y="value",
                        order=cluster_order,
                        inner=None,
                        ax=ax,
                        cut=0
                    )
                except Exception:
                    # Fallback if seaborn struggles for some reason
                    ax.text(0.5, 0.5, f"{feat}\n(plot unavailable)", ha="center", va="center")
                    ax.axis("off")
                    plot_idx += 1
                    continue

                # Overlay all points (each sample), jittered
                if show_points:
                    try:
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
                    except Exception:
                        pass

                # Overlay mean marker per cluster
                try:
                    means = sub.groupby("ClusterID", observed=False)["value"].mean().reindex(cluster_order)
                    # X positions correspond to categorical ticks [0..n-1]
                    x_pos = np.arange(len(cluster_order))
                    ax.scatter(
                        x_pos,
                        means.values,
                        s=mean_marker_size,
                        edgecolor="black",
                        linewidths=0.8,
                        zorder=3
                    )
                except Exception:
                    pass

                ax.set_title(feat)
                ax.set_xlabel("ClusterID")
                ax.set_ylabel("Value")
                ax.tick_params(axis="x", rotation=0)

                plot_idx += 1

            # Layout and save page
            fig.subplots_adjust(left=0.12, right=0.98, top=0.95, bottom=0.08)
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

def plot_dead_dye_distribution(
    df_tracks,
    outpath,
    nr_cols=3,
    rows_per_page = 3,
    figsize=(8.27, 11.69),
    plot_results=True
    ):
    """
    Create a violin plot with an underlying scatterplot that provides an
    overview of the mean dead dye intensity per segment at the first timepoint
    in each experiment
    """
    sample_names = df_tracks["sample_name"].unique()
    n_plots = len(sample_names)
    n_rows = (n_plots // nr_cols) + (1 if n_plots % nr_cols != 0 else 0)
    nr_pages = (n_rows // rows_per_page) + (1 if n_rows % rows_per_page != 0 else 0)
    
    with PdfPages(outpath) as pdf:
        # df_time1 = df_tracks[df_tracks["relative_time"]==1]
        plot_idx = 0  # Track which plot we are adding
        for page in range(nr_pages):
            fig = plt.figure(figsize=figsize)
            gs = GridSpec(rows_per_page, nr_cols, figure=fig, wspace=0.3, hspace=0.1)
            remaining_axes = [
                fig.add_subplot(gs[i, j]) 
                for i in range(rows_per_page)
                for j in range(nr_cols)
                ]

            for ax in remaining_axes:
                if plot_idx >= n_plots:
                    ax.remove()  # Remove empty axes
                    continue
                
                sample = sample_names[plot_idx]
                df_subset = df_tracks[(df_tracks["sample_name"] == sample)]

                # Violin plot
                sns.violinplot(
                    data=df_subset, 
                    # x='sample_name', 
                    y='mean_dead_dye', 
                    dodge=False, 
                    inner=None,
                    ax=ax
                    )

                # Jitter plot (scatter over the violin)
                sns.stripplot(
                    data=df_subset, 
                    # x='sample_name', 
                    y='mean_dead_dye', 
                    color='black', 
                    alpha=0.5, 
                    jitter=True,
                    ax=ax
                    )
                ax.set_title(sample, fontsize=10, loc='center')
                ax.set_xticks([])
                ax.set_xticklabels([])
                ax.set_xlabel("")
                ax.set_ylabel("Mean Dead Dye")
                ax.grid(True, linestyle="--", alpha=0.7)
                sns.despine()
                plot_idx += 1
                
            # plt.show(fig)
            fig.subplots_adjust(left=0.05, right=0.85, top=0.95, bottom=0.05)
            pdf.savefig(fig, bbox_inches='tight', dpi=600)
            plt.close(fig)
        
def plot_touching_nontouching_distribution(
    df_tracks,
    outpath,
    contact_column='organoid_contact',
    nr_cols=3,
    rows_per_page = 3,
    figsize=(8.27, 11.69),
    plot_results=True
    ):
    """
    Create a barplot that provides an overview of how many cells make contact with
    other organoids/cells
    """
    
    sample_names = df_tracks["sample_name"].unique()
    n_plots = len(sample_names)
    n_rows = (n_plots // nr_cols) + (1 if n_plots % nr_cols != 0 else 0)
    nr_pages = (n_rows // rows_per_page) + (1 if n_rows % rows_per_page != 0 else 0)
        
    with PdfPages(outpath) as pdf:
        # organoid_lines = df_tracks["organoid_line"].unique()
        plot_idx = 0  # Track which plot we are adding
        for page in range(nr_pages):
            fig = plt.figure(figsize=figsize)
            gs = GridSpec(rows_per_page, nr_cols, figure=fig, wspace=0.5, hspace=0.3)

            remaining_axes = [
                fig.add_subplot(gs[i, j]) 
                for i in range(rows_per_page)
                for j in range(nr_cols)
                ]

            for ax in remaining_axes:
                if plot_idx >= n_plots:
                    ax.remove()  # Remove empty axes
                    continue
                
                sample = sample_names[plot_idx]
                df_subset = df_tracks[(df_tracks["sample_name"] == sample)]

                sns.histplot(
                    df_subset, 
                    x=contact_column, 
                    multiple="dodge", 
                    shrink=0.8, 
                    discrete=True, 
                    ax=ax
                    )
                ax.set_title(f"{sample}", fontsize=10)
                plot_idx += 1
            # Add a global title
            if contact_column == "organoid_contact":
                fig.suptitle("Touching vs. Non-touching Organoids", fontsize=14, fontweight="bold")
            elif contact_column == "tcell_contact":
                fig.suptitle("Touching vs. Non-touching T cells", fontsize=14, fontweight="bold")

            # Adjust layout and save
            # plt.tight_layout(rect=[0, 0, 1, 0.96])
            fig.subplots_adjust(left=0.05, right=0.85, top=0.95, bottom=0.05)
            # plt.show(fig)
            pdf.savefig(fig, bbox_inches='tight', dpi=600)
            plt.close(fig)
            # plt.savefig("output.pdf", bbox_inches='tight', dpi=600)

def round_legend_ticks(value):
    if value <= 1.0:
        # return np.ceil(value * 10) / 10
        return 1.0
    elif value <= 100:
        return np.ceil(value / 10) * 10
    elif value <= 10000:
        return np.ceil(value / 500) * 500
    return value

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Input parameters for automatic data transfer.')
    parser.add_argument('-c', '--config', type=str, help='path to a config.yml file that stores all required paths', required=False)
    args = parser.parse_args()
    with open(args.config, "r") as parameters:
        config=yaml.load(parameters, Loader=yaml.SafeLoader)
    # metadata = pd.read_csv(config["metadata_csv_path"])
    run_tcell_analysis(config)