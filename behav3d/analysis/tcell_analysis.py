import argprase
from dtaidistance import dtw, dtw_ndim
import pandas as pd
import numpy as np
import umap
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import patchworklib as pw
import random
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from plotnine import *
from pathlib import Path
import yaml
# df_tracks=df_tracks[df_tracks["relative_time"]<=30]

def run_tcell_analysis(
    config,
    df_tracks_path=None,
    df_tracks_summarized_path =None 
    ):
    
    output_dir = config["output_dir"]
    if df_tracks_path is None:
        df_tracks_path = Path(output_dir, f"BEHAV3D_combined_track_features_filtered.csv")
    if df_tracks_summarized_path is None:
        df_tracks_summarized_path = Path(output_dir, f"BEHAV3D_combined_track_features_summarized.csv")
    
    df_tracks = pd.read_csv(df_tracks_path)
    df_tracks_summarized = pd.read_csv(df_tracks_summarized_path)
    df_tracks=df_tracks.sort_values(by=["sample_name", "TrackID", "relative_time"])
    
    dtw_distance_matrix = calculate_dtw(
        df_tracks
        )
    
    umap_embedding = fit_umap(
        dtw_distance_matrix=dtw_distance_matrix,
        config=config,
        random_state=None
    )
    
    df_clusters = cluster_umap(
        umap_embedding,
        
    )
    
def calculate_dtw(
    df_tracks, 
    features=[
        "mean_square_displacement", 
        "speed", 
        "dead_dye_mean", 
        "tcell_contact", 
        "organoid_contact"
        ]
    ):
    df_tracks=df_tracks.sort_values(by=["sample_name", "TrackID", "relative_time"])
    
    nr_tracks=len(df_tracks[["sample_name", "TrackID"]].drop_duplicates())
    nr_timepoints=len(df_tracks["relative_time"].unique())
    nr_features=len(features)

    dtw_input_tracks = np.empty((nr_tracks, nr_timepoints, nr_features),dtype=np.double)
    for i, feature in enumerate(features):
        pivot_df = df_tracks.pivot(
            index=['sample_name', 'TrackID'], 
            columns='relative_time', 
            values=feature
            )
        numpy_array = pivot_df.values.astype(np.double)
        dtw_input_tracks[:, :, i] = numpy_array
    dtw_input_tracks=dtw_input_tracks.astype(np.double)

    dtw_distance_matrix = dtw_ndim.distance_matrix_fast(dtw_input_tracks)
    return(dtw_distance_matrix)

def fit_umap(
    dtw_distance_matrix,
    config,
    random_state=None
    ):
    
    umap_model = umap.UMAP(
        n_components=2, 
        n_neighbors=config["umap_n_neighbors"], 
        min_dist=config["umap_minimal_distance"], 
        init="random", 
        random_state=random_state,
        metric="precomputed", 
        )
    umap_embedding = umap_model.fit_transform(dtw_distance_matrix)
    umap_embedding = pd.DataFrame(umap_embedding, columns=['UMAP1', 'UMAP2'])  
    return(umap_embedding)

def cluster_umap(
    umap_embedding,
    config,
    df_tracks=None,
    df_track_summarized=None,
    nr_clusters=6,
    random_state=None
    ):
    
    output_dir = config['output_dir']
    if df_tracks is None:
        df_tracks_path = Path(output_dir, f"BEHAV3D_combined_track_features_filtered.csv")
        df_tracks = pd.read_csv(df_tracks_path)
    if df_tracks_summarized is None:
        df_tracks_summarized_path = Path(output_dir, f"BEHAV3D_combined_track_features_summarized.csv")
        df_tracks_summarized = pd.read_csv(df_tracks_summarized_path)
      
    df_tracks=df_tracks.sort_values(by=["sample_name", "TrackID", "relative_time"])
    TrackIDs = df_tracks[["sample_name", "TrackID"]].drop_duplicates().reset_index(drop=True)
    df_umap = pd.concat([TrackIDs, umap_embedding], axis=1)
    # df_trackinfo = df_tracks[['TrackID', 'sample_name','well', 'exp_nr', 'organoid_line', 'tcell_line']].drop_duplicates()
    df_umap = pd.merge(df_track_summarized, df_umap, how="left")
    
    # Perform clustering
    scaler = StandardScaler()
    # umap_scaled = scaler.fit_transform(umap_embedding)  # Standardize UMAP coordinates
    kmeans = KMeans(n_clusters=config["nr_of_clusters"], n_init=100, random_state=random_state)
    df_umap["cluster"] = kmeans.fit_predict(umap_embedding)
    # df_umap["cluster2"] = kmeans.fit_predict(umap_embedding)
    df_umap["cluster"]=df_umap["cluster"].astype('category')
    
    df_umap_out_path = Path(output_dir, f"BEHAV3D_UMAP_clusters.csv")
    print(f"- Writing summarized tracks to {df_umap_out_path}")
    df_umap.to_csv(df_umap_out_path, sep=",", index=False)

    umap_plots = []
    
    info_cols = df_umap.drop(columns=["TrackID", "well", "exp_nr", "UMAP1", "UMAP2", "cluster"]).columns
    
    cluster_plot = (
            ggplot(df_umap, aes(x='UMAP1', y='UMAP2', color="cluster")) +
            geom_point(size=2, alpha=0.6) +
            labs(color="cluster") +
            labs(title="cluster") +
            labs(x="", y="") +
            theme_light(base_size=20) +
            theme_bw() +
            theme(axis_text_x=element_blank(), axis_text_y=element_blank(), aspect_ratio=1) +
            # scale_color_discrete() +
            coord_fixed()
            # facet_grid('organoid_line ~ sample_name')
        )
    
    for colorcol in info_cols:
        plot = (
            ggplot(df_umap, aes(x='UMAP1', y='UMAP2', color=colorcol)) +
            geom_point(size=2, alpha=0.6) + #show_legend=False
            labs(color=colorcol) +
            labs(title=colorcol) +
            labs(x="", y="") +
            theme_light(base_size=20) +
            theme_bw() +
            theme(axis_text_x=element_blank(), axis_text_y=element_blank(), aspect_ratio=1) +
            # scale_color_discrete() +
            coord_fixed()
            # facet_grid('organoid_line ~ sample_name')
        )
        umap_plots.append(pw.load_ggplot(plot, figsize=[4,4]))
    
    nr_cols = 2
    comb_plot=(pw.load_ggplot(cluster_plot, figsize=[4,4]))
    for idr in range(0, len(info_cols), nr_cols):
        rowplots = umap_plots[idr]
        for idc in range(1,nr_cols):
            if idc+idr+1 <= len(info_cols):
                rowplots |= umap_plots[idc+idr]
            else:
                rowplots |= pw.load_ggplot((ggplot()+geom_blank()+theme_void()), figsize=[4,4])
            comb_plot /= (rowplots)
    
    cluster_UMAP_path = Path(output_dir, f"BEHAV3D_UMAP_clusters.pdf")
    comb_plot.savefig(cluster_UMAP_path)
    
    cluster_means = df_umap.drop(
        columns=[
            "sample_name", 
            "organoid_line", 
            "tcell_line",
            "well",
            "exp_nr",
            "UMAP1",
            "UMAP2",
            "TrackID"
            ]).groupby('cluster').mean().reset_index()
    df_heatmap = cluster_means.melt(id_vars='cluster', var_name='var', value_name='value')
    heatmap_plot = (
        ggplot(df_heatmap, aes(x='cluster', y='var', fill='value')) +
        geom_tile() +
        scale_fill_gradient(low="blue", high="yellow") +  # Adjust the color scale as needed
        labs(x='Cluster', y='Mean Columns', fill='Value') +
        theme_minimal()
    )
    cluster_features_heatmap_path = Path(output_dir, f"BEHAV3D_UMAP_cluster_feature_heatmap.pdf")
    heatmap_plot.save(cluster_features_heatmap_path, width=8, height=8)
    
    df_clust_perc = df_umap.groupby(["organoid_line", "tcell_line", "cluster"]).size().reset_index(name='count')
    total_counts = df_clust_perc.groupby(['organoid_line', 'tcell_line'])['count'].sum().reset_index(name='total_count')
    total_counts["total_count"] = '# Cells: ' + total_counts["total_count"].astype("string")
    cluster_counts = pd.merge(df_clust_perc, total_counts)
    df_clust_perc["percentage"] = (df_clust_perc['count'] / df_clust_perc['count'].sum())

    total_counts_facet = cluster_counts.groupby(['organoid_line', 'tcell_line'])['total_count'].mean().reset_index()
    
    plot = (
        ggplot(df_clust_perc) + 
        geom_col(aes(x=0, y='percentage', fill='cluster')) + 
        # geom_text(aes(x=0, y=1,label='sum(count)'), va='bottom', size=8, position='identity') +
        theme_void() + 
        facet_grid('tcell_line ~ organoid_line', scales='free_y') +
        labs(x='', y='Cluster', title='Horizontal Stacked Bar Chart') +
        coord_flip() +
        theme(aspect_ratio = 0.2) +
        # geom_text(data=total_counts, aes(x=0, y=1,label='sum(count)'), va='bottom', size=8, position='identity') 
        geom_text(aes(label="total_count", x=-1.0, y=0.5), data=total_counts, va='bottom', size=8, position='identity')
        # annotate('text', x=0.5, y=105, label=lambda d: f"Total Count: {d['total_count'].iloc[0]}", size=12)
    )
    cluster_percentage_plot_path = Path(output_dir, f"BEHAV3D_UMAP_cluster_percentages.pdf")
    plot.save(cluster_percentage_plot_path, width=8, height=8)

    df_clust_perc = df_clust_perc.reset_index(drop=True)
    df_clust_perc_out_path = Path(output_dir, f"BEHAV3D_UMAP_cluster_percentages.csv")
    print(f"- Writing summarized tracks to {df_clust_perc_out_path}")
    df_clust_perc.to_csv(df_clust_perc_out_path, sep=",", index=False)
    
    return()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Input parameters for automatic data transfer.')
    parser.add_argument('-c', '--config', type=str, help='path to a config.yml file that stores all required paths', required=False)
    args = parser.parse_args()
    with open(args.config, "r") as parameters:
        config=yaml.load(parameters, Loader=yaml.SafeLoader)
    # metadata = pd.read_csv(config["metadata_csv_path"])
    run_tcell_analysis(config)