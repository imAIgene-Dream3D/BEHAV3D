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
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler, MinMaxScaler

from pathlib import Path
from behav3d import format_time
import yaml
import time
import seaborn as sns
# df_tracks=df_tracks[df_tracks["relative_time"]<=30]

def run_tcell_analysis(
    config=None,
    output_dir=None,
    df_tracks_path=None,
    df_tracks_summarized_path=None,
    dtw_features=[
        "z_mean_square_displacement", 
        "z_speed", 
        "z_mean_dead_dye", 
        "tcell_contact", 
        "organoid_contact"
        ],
    umap_minimal_distance=None,
    umap_n_neighbors=None,
    nr_of_clusters=None
    ):
    print(f"--------------- Performing T-cell behavioral analysis ---------------")
    start_time = time.time()
    assert config is not None or all(
        [output_dir, umap_minimal_distance, umap_n_neighbors, nr_of_clusters]
    ), "Either 'config' or 'output_dir, umap_minimal_distance, umap_n_neighbors, nr_of_clusters' parameters must be supplied"
        
    if output_dir is None:
        output_dir = config['output_dir']
    analysis_outdir = Path(output_dir, "analysis", "tcells")
    feature_outdir = Path(analysis_outdir, "track_features")
    
    if not analysis_outdir.exists():
        analysis_outdir.mkdir(parents=True)
    if not feature_outdir.exists():
        feature_outdir.mkdir(parents=True)    
    
    if df_tracks_path is None:
        df_tracks_path = Path(feature_outdir, f"BEHAV3D_combined_track_features_filtered.csv")
    if df_tracks_summarized_path is None:
        df_tracks_summarized_path = Path(feature_outdir, f"BEHAV3D_combined_track_features_summarized.csv")
    
    df_tracks = pd.read_csv(df_tracks_path)
    df_tracks_summarized = pd.read_csv(df_tracks_summarized_path)
    df_tracks=df_tracks.sort_values(by=["sample_name", "TrackID", "relative_time"])
    
    if df_tracks_summarized["track_length"].nunique() != 1:
        print("Warning: The track lengths are not cut to similar length, this might influence dynamic time warping")
        print("Set 'tcell_min_track_length' and 'tcell_max_track_length' to the same value to create equal tracks")

    dtw_distance_matrix = calculate_dtw(
        df_tracks, 
        features=dtw_features
        )
    
    umap_embedding = fit_umap(
        dtw_distance_matrix=dtw_distance_matrix,
        umap_n_neighbors=umap_n_neighbors,
        umap_minimal_distance=umap_minimal_distance,
        random_state=None
    )
    
    df_clusters = cluster_umap(
        umap_embedding=umap_embedding,
        output_dir = output_dir,
        nr_of_clusters=nr_of_clusters,
        df_tracks=df_tracks,
        df_tracks_summarized=df_tracks_summarized,
        random_state=None
    )
    end_time = time.time()
    h,m,s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    return(df_clusters)

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

#  # dtw_input_tracks = np.empty((nr_tracks, nr_timepoints, nr_features),dtype=np.double)
#     colnames=['TrackID', 'sample_name']+features
#     dtw_input_tracks=pd.DataFrame(columns=colnames)
#     unique_tracks = df_tracks.groupby(['TrackID', 'sample_name'])
#     for (TrackID, sample_name), group in unique_tracks:
#         track_features = group[features].to_numpy().astype(np.double)
#         # dtw_input_tracks=pd.concat([dtw_input_tracks, group[colnames]])
#         # track_features = group[['TrackID', 'sample_name']+features]
#         # track_features = group[features].to_numpy().astype(np.double)
#         # dtw_input_tracks.append(track_features)
#     dtw_input = dtw_input_tracks.drop(columns=['TrackID', 'sample_name']).to_numpy().astype(np.double)
#     dtw_distance_matrix = dtw_ndim.distance_matrix_fast(dtw_input)
#     dtw_distance_matrix = pd.concat([pd.DataFrame(dtw_distance_matrix)
#     return(dtw_distance_matrix)

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
            
    if all([umap_n_neighbors, umap_minimal_distance]) is None:
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
    output_dir = None
    ):
    
    assert config is not None or all(
        [output_dir, nr_of_clusters]
    ), "Either 'config' or 'output_dir, nr_of_clusters' parameters must be supplied"
      
    print("- Performing clustering on the UMAP data")
    if all([output_dir, nr_of_clusters]) is None:
        output_dir = Path(config['output_dir'])
        nr_of_clusters=config["nr_of_clusters"]
    
    tcell_outdir = Path(output_dir, "analysis", "tcells")
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
    df_umap["ClusterID"] = kmeans.fit_predict(umap_embedding.drop(columns=["TrackID","sample_name"]))
    # df_umap["cluster2"] = kmeans.fit_predict(umap_embedding)
    
    # Set cluster index to start from 1 for backprojection purposes
    df_umap["ClusterID"]=df_umap["ClusterID"]+1
    df_umap["ClusterID"]=df_umap["ClusterID"].astype('category')
    
    df_umap_out_path = Path(results_outdir, f"BEHAV3D_UMAP_clusters.csv")
    print(f"- Writing clustered tracks to {df_umap_out_path}")
    df_umap.to_csv(df_umap_out_path, sep=",", index=False)

    print("- Producing clustered UMAP plots with displayed Track features")
    sample_cols = ["organoid_line", "tcell_line"]
    info_cols = df_umap.drop(columns=["TrackID", "sample_name", "well", "exp_nr", "UMAP1", "UMAP2", "ClusterID"]).columns
    
    cluster_UMAP_path = Path(results_outdir, f"BEHAV3D_UMAP_clusters.pdf")
    create_umap_plot(
        df_umap=df_umap,
        info_cols=info_cols,
        sample_cols=sample_cols,
        outpath=cluster_UMAP_path,
        rows_per_page = 4,
        nr_cols = 2,
        rows_first_img = 2,
        figsize = (8.27, 11.69)
    )
    
    ### Producing a heatmap of the summarized features again summarized over all tracks
    ### Belonging to that cluster
    print("- Producing heatmaps with summarized cluster features")
    cluster_features_heatmap_path = Path(results_outdir, f"BEHAV3D_UMAP_cluster_feature_heatmap.pdf")
    create_heatmap_plot(
        df_umap,
        info_cols,
        sample_cols,
        cluster_features_heatmap_path,
        rows_per_page = 7,
        nr_cols = 2,
        rows_first_img = 4,
        figsize = (8.27, 11.69)
    )

    print("- Producing percentage plots of each cluster per combination of T-cell and organoid line")
    df_clust_perc = df_umap.groupby(["organoid_line", "tcell_line", "ClusterID"]).size().reset_index(name='count')
    total_counts = df_clust_perc.groupby(['organoid_line', 'tcell_line'])['count'].sum().reset_index(name='total_count')
    
    df_clust_perc = pd.merge(df_clust_perc, total_counts)
    df_clust_perc["percentage"] = (df_clust_perc['count'] / df_clust_perc['total_count'])
    
    cluster_percentage_plot_path = Path(results_outdir, f"BEHAV3D_UMAP_cluster_percentages.pdf")
    create_percentage_plot(
        df_clust_perc,
        cluster_percentage_plot_path
    )
        
    df_clust_perc = df_clust_perc.reset_index(drop=True)
    df_clust_perc_out_path = Path(results_outdir, f"BEHAV3D_UMAP_cluster_percentages.csv")
    print(f"- Writing summarized tracks to {df_clust_perc_out_path}")
    df_clust_perc.to_csv(df_clust_perc_out_path, sep=",", index=False)
    
    return()

def create_percentage_plot(
    df_clust_perc,
    outpath
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
        plt.show()
        pdf.savefig(fig, bbox_inches='tight', dpi=600)
        plt.close(fig)
            
def create_umap_plot(
    df_umap,
    info_cols,
    sample_cols,
    outpath,
    rows_per_page = 4,
    nr_cols = 2,
    rows_first_img = 2,
    figsize = (8.27, 11.69)
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
            plt.show()
            pdf.savefig(fig, dpi=600)
            plt.close(fig)

def create_heatmap_plot(
    df_umap,
    info_cols,
    sample_cols,
    outpath,
    rows_per_page = 7,
    nr_cols = 2,
    rows_first_img = 4,
    figsize = (8.27, 11.69)
    ):

    # Producing heatmaps with summarized cluster features
    cluster_means = df_umap[list(info_cols)+["ClusterID"]].drop(columns=sample_cols).groupby('ClusterID').mean().reset_index()
    df_heatmap = cluster_means.melt(id_vars='ClusterID', var_name='var', value_name='value')

    # Plot an overall heatmap where every feature is scaled from 0 to 1
    cluster_means_scaled = cluster_means.copy()
    scale_columns = cluster_means_scaled.columns[cluster_means.columns != "ClusterID"]
    cluster_means_scaled[scale_columns] = MinMaxScaler().fit_transform(cluster_means_scaled[scale_columns])
    df_heatmap_scaled = cluster_means_scaled.melt(id_vars='ClusterID', var_name='var', value_name='AU')

    overall_heatmap_data = df_heatmap_scaled.pivot(index="var", columns="ClusterID", values="AU")
    
    # Plot the heatmap separated per feature, but scale in original values
    columns = df_heatmap["var"].unique()
    
    n_plots = len(columns)
    n_rows = (n_plots // nr_cols) + (1 if n_plots % nr_cols != 0 else 0) + rows_first_img
    nr_pages = (n_rows // rows_per_page) + (1 if n_rows % rows_per_page != 0 else 0)
    with PdfPages(outpath) as pdf:
        plot_idx = 0  # Track which plot we are adding

        for page in range(nr_pages):
            fig = plt.figure(figsize=figsize)
            gs = GridSpec(rows_per_page, nr_cols, figure=fig, hspace=1.5, wspace=0.3)

            # First image on the first page
            if page == 0:
                ax = fig.add_subplot(gs[:rows_first_img, :])
                heatmap = sns.heatmap(overall_heatmap_data, ax=ax, cmap="viridis", cbar=True)
                ax.set_title("Min-Max scaled heatmap", fontsize=16)
                ax.set_xlabel("ClusterID")
                ax.set_ylabel("")
                ax.tick_params(axis='y', labelsize=8)
                cbar = heatmap.collections[0].colorbar
                cbar.ax.tick_params(labelsize=12)
            # Remaining plots
            remaining_axes = [
                fig.add_subplot(gs[i, j])
                for i in range(rows_first_img if page == 0 else 0, rows_per_page)
                for j in range(nr_cols)
            ]

            for ax in remaining_axes:
                if plot_idx >= len(columns):
                    ax.remove()  # Remove empty axes
                    continue
                col = columns[plot_idx]
                col_heatmap = df_heatmap[df_heatmap["var"] == col].pivot(index="var", columns="ClusterID", values="value")
                vmin=0
                vmax = round_legend_ticks(df_heatmap[df_heatmap["var"] == col]["value"].max())
                
                heatmap = sns.heatmap(col_heatmap, ax=ax, cmap="viridis", cbar=True, vmin=0, vmax=vmax)
                ax.set_title(col)
                ax.set_xlabel("ClusterID")
                ax.set_ylabel("")
                # ax.set_xticks([])
                ax.set_yticks([])
                # ax.set_xticklabels([])
                ax.set_yticklabels([])
                cbar = heatmap.collections[0].colorbar
                cbar.ax.tick_params(labelsize=8)  # Adjust font size of colorbar ticks
                cbar.set_ticks([vmin, vmax/2, vmax])
                plot_idx += 1  # Move to the next plot

            # Save and close the figure for this page
            # fig.tight_layout()
            fig.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.05)
            # fig.tight_layout(pad=2.0)
            # fig.set_constrained_layout(True)
            plt.show()
            pdf.savefig(fig, dpi=600)
            plt.close(fig)

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