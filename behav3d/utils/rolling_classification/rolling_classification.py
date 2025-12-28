from behav3d.utils import load_behav3d_metadata, check_behav3d_metadata

import pandas as pd
from pathlib import Path

import pandas as pd
import numpy as np

import umap

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import random

from sklearn.cluster import KMeans, HDBSCAN
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from sklearn.decomposition import PCA
from sklearn.feature_selection import VarianceThreshold


from pathlib import Path

import seaborn as sns
# import hdbscan

from behav3d.utils.rolling_classification import *
from behav3d.utils.rolling_classification.features import *
from behav3d.utils.rolling_classification.clustering import *
from behav3d.utils.rolling_classification.plotting import *
from behav3d.utils.rolling_classification.videos import *
# %matplotlib inline

seed = 123
random.seed(seed)
np.random.seed(seed)

if __name__ == "__main__":
    # ssd_dir = r"/Volumes/T7_Sam/"
    ssd_dir = r"F:/"
    ssd_dir = Path(ssd_dir)
    
    output_dir = Path(ssd_dir, r"BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE")
    metadata_csv_path = Path(ssd_dir, r"BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE/metadata.csv")

    # output_dir = Path(ssd_dir, r"BHVD_BEHAV3D/BEHAV3D_python/runs/CombinedAnalysis_AmberMacrophage")
    # metadata_csv_path = Path(ssd_dir, r"BHVD_BEHAV3D/BEHAV3D_python/runs/CombinedAnalysis_AmberMacrophage/metadata.csv")


    # outfolder = r"C:/Users/Samde/Downloads"
    outfolder = Path(ssd_dir, r"BHVD_BEHAV3D\BEHAV3D_python\rolling_classification")
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

    window_size = 10
    # window_size = None

    groupby=["sample_name", "TrackID"]
    features=[
        "percentage_dead_mask",
        # "mean_dead_dye",
        # "nr_dead_mask_pixels",
        "organoid_contact_pixels",
        "tcell_contact_pixels",
        # "mean_square_displacement",
        "speed",
        # "directional_persistence",
        # "volume",
        # "extent",
        # "elongation",
        # "sphericity",
        # "solidity",
    ]

    # df_tracks = df_positions[["sample_name", "TrackID", "position_t"]+features]
    df_windows_descriptive = create_descriptive_track_dataset(
        df_tracks=df_positions,
        columns_to_summarize=features,
        window_size = window_size,
        step_size = 1,
        time_col = "position_t",
        id_cols = ["sample_name", "TrackID"],
    )

    df_windows_descriptive.to_csv(Path(outfolder,r"df_windows_descriptive.csv"), index=False)
    # df_windows_descriptive = pd.read_csv(Path(outfolder,r"df_windows_descriptive.csv"))
    
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
    # df_analysis, dropped, report = drop_highly_correlated_features(
    #     df=df_analysis,
    #     feature_cols=descriptive_feature_cols,
    #     threshold=0.95
    # )
    # print(f"Dropped {len(dropped)} highly correlated features.")
    # descriptive_feature_cols = [c for c in descriptive_feature_cols if c not in dropped]
    
    
    ### Remove features with no variance
    selector = VarianceThreshold(threshold=1e-4)
    selector.fit(df_analysis[descriptive_feature_cols])

    keep_mask = selector.get_support()
    kept_features = df_analysis[descriptive_feature_cols].columns[keep_mask].tolist()
    dropped_low_var = df_analysis[descriptive_feature_cols].columns[~keep_mask].tolist()
    print(f"Dropped {len(dropped_low_var)} low-variance features.")
    descriptive_feature_cols = [c for c in descriptive_feature_cols if c not in dropped_low_var]
    
    scaler = StandardScaler().fit(df_analysis[descriptive_feature_cols])
    df_analysis[descriptive_feature_cols] = scaler.transform(df_analysis[descriptive_feature_cols])
    
    adata_full = df_to_adata(df_analysis, descriptive_feature_cols, obs_cols=non_feature_cols)
    adata_full.uns["preprocessing"] = {
        "kept_features": list(descriptive_feature_cols),
        "scaler": {
            "mean": scaler.mean_.astype(float),
            "scale": scaler.scale_.astype(float),
        }}

    # df_analysis[descriptive_feature_cols] = scaler.fit_transform(df_analysis[descriptive_feature_cols])
    # --- 3) Subset windows (custom) ---
    chosen_intervals = 30
    # chosen_intervals = window_size
    df_leiden_subset = subset_windowed_tracks(
        df_windowed=df_analysis,
        step_size=chosen_intervals,
    )
    
    
    adata_sub  = df_to_adata(df_leiden_subset, descriptive_feature_cols, obs_cols=non_feature_cols)
    adata_sub.uns["preprocessing"] = {
        "kept_features": list(descriptive_feature_cols),
        "scaler": {
            "mean": scaler.mean_.astype(float),
            "scale": scaler.scale_.astype(float),
        }}
    # sc.pp.scale(adata_full, zero_center=True, max_value=None)
    # scaled_subset_X = adata_full[adata_sub.obs_names, adata_sub.var_names].X
    # adata_sub.X = scaled_subset_X.copy()
    
    """
    LEIDEN CLUSTERING
    """
    adata_sub = run_pca(
        adata_sub,
        pca_var=0.95,
        ncomps=len(descriptive_feature_cols),
        svd_solver='full', 
        random_state=seed
    )
    # pca_var=0.95
    # # sc.pp.pca(adata_sub, n_comps=None, )
    # sc.pp.pca(
    #     adata_sub, 
    #     zero_center=True,
    #     n_comps=len(descriptive_feature_cols), 
    #     svd_solver='full', 
    #     random_state=seed,
    #     )

    # var_ratio = adata_sub.uns["pca"]["variance_ratio"]
    # n_pcs = int(np.searchsorted(np.cumsum(var_ratio), pca_var) + 1)
    # adata_sub.obsm["X_pca"] = adata_sub.obsm["X_pca"][:, :n_pcs]
    
    # adata_sub = run_leiden_clustering(
    #     adata_sub, 
    #     resolution="auto", 
    #     stability_resolutions=(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0), # Gives 0.5 as result
    #     key_added="ClusterID",
    #     random_state=seed
    #     )

    adata_sub = run_leiden_clustering(
        adata_sub, 
        n_neighbors=50,
        resolution=0.25, 
        metric="euclidean",
        method="umap",
        use_rep="X_pca",
        key_added="ClusterID",
        random_state=seed
        )
    
    adata_sub = merge_small_clusters(
        adata_sub,
        key="ClusterID",
        min_size=500,
        use_rep="X_pca",  # or "X" or "X_umap" depending on what you trust
    )

    sc.tl.umap(
        adata_sub,
        min_dist=0.1,
        random_state=seed,
    )
    
    sc.pl.umap(adata_sub, color="ClusterID", size=2, alpha=0.5)
    
    feature_dict = {
        f: [c for c in descriptive_feature_cols if c.startswith(f + "_")]
        for f in features
    }
    known_prefixes = tuple(f + "_" for f in features)
    feature_dict["other"] = [c for c in descriptive_feature_cols if not c.startswith(known_prefixes)]
    
    sc.tl.dendrogram(adata_sub, groupby="ClusterID")
    sc.pl.heatmap(
        adata_sub,
        var_names=feature_dict,
        var_group_labels=list(feature_dict.keys()),
        var_group_rotation=0,
        groupby="ClusterID",
        standard_scale="var",
        figsize=(25, 20),
        swap_axes=True,
        dendrogram=True,
        show_gene_labels=True,
    )
    
    sc.pl.matrixplot(
        adata_sub,
        var_names=feature_dict,    
        groupby="ClusterID",
        standard_scale="var",
        figsize=(20, 20),
        swap_axes=True,
        dendrogram=True
    )
    
    sc.tl.ingest(
        adata_full,
        adata_sub,
        obs="ClusterID",
        embedding_method="umap"  # also transfers UMAP coords into adata_full.obsm["X_umap"]
    )
    
    adata_full = filter_short_state_runs(
        adata_full,
        cluster_key="ClusterID",
        id_cols = ["sample_name", "TrackID"],
        time_key = "position_t",
        length_removed = 1,
        new_key = "ClusterID_filt1"
    )
    
    adata_full = filter_short_state_runs(
        adata_full,
        cluster_key="ClusterID",
        id_cols = ["sample_name", "TrackID"],
        time_key = "position_t",
        length_removed = 2,
        new_key = "ClusterID_filt2"
    )
    
    adata_full = filter_short_state_runs(
        adata_full,
        cluster_key="ClusterID",
        id_cols = ["sample_name", "TrackID"],
        time_key = "position_t",
        length_removed = 3,
        new_key = "ClusterID_filt3"
    )
    
    adata_full = filter_short_state_runs(
        adata_full,
        cluster_key="ClusterID",
        id_cols = ["sample_name", "TrackID"],
        time_key = "position_t",
        length_removed = 5,
        new_key = "ClusterID_filt5"
    )
    
    adata_full.write(Path(outfolder, "adata_full.h5ad"), compression="gzip") 
    adata_sub.write(Path(outfolder, "adata_sub.h5ad")) 
    
    # write results back
    df_analysis = adata_add_back_to_df(
        df_analysis, adata_full,
        cols_from_obs=["ClusterID"],
    )
    
    df_leiden_subset = adata_add_back_to_df(
        df_leiden_subset, adata_sub,
        cols_from_obs=["ClusterID"],
    )
    
    ###### ANALYSIS VALUES
    plot_number_per_clusters(df_leiden_subset, cluster_col="ClusterID")
    plot_per_cluster_proportions(df_leiden_subset)
    

    compute_cluster_transition_matrix(
        adata=adata_full,
        cluster_key="ClusterID",
        id_cols = ["sample_name", "TrackID"],
        plot=True
    )
    
    _,transmat = compute_cluster_transition_matrix(
        adata=adata_full,
        cluster_key="ClusterID",
        id_cols = ["sample_name", "TrackID"],
        only_transitions=False,
        plot=True
    )
    
    transmat_nr,transmat = compute_cluster_transition_matrix(
        adata=adata_full,
        cluster_key="ClusterID",
        id_cols = ["sample_name", "TrackID"],
        only_transitions=True,
        plot=True
    )
    
    plot_exemplar_track_bars(
        adata_full,
        n_tracks=50,
        track_key="TrackID",
        time_key="position_t",
        state_key="ClusterID",
        seed=seed,
        cmap_name="tab20",
        ax=None,
    )
    
    """
    ################################################
    HMM STATE CLASSIFICATION (on descriptive data)
    ################################################
    """
    
    
    df_hmm_descriptive = df_analysis.copy()
    # scaled_X = adata_full[:, descriptive_feature_cols].X.toarray()
    # df_analysis_hmm.loc[:, descriptive_feature_cols] = scaled_X

    df_hmm_descriptive, sampled_keys = subset_full_tracks(
        df=df_hmm_descriptive,
        fraction=0.5,
        random_state=seed,
        id_cols=["sample_name", "TrackID"],
        return_selected_keys=True
    )
    
    # df_hmm_descriptive, hmm_model, selection_df = run_hmm_state_classification(
    #     df_features=df_hmm_descriptive,
    #     feature_cols=descriptive_feature_cols,
    #     k_min = 3,
    #     k_max = 15,
    #     n_states="auto", #10 came out of it
    #     out_col_name="hmm_state"
    # )
    
    
    df_hmm_descriptive, hmm_model_descriptive, selection_df = run_sticky_hmm_state_classification(
        df_features=df_hmm_descriptive,
        feature_cols=descriptive_feature_cols,
        n_states=8, #6 came out of it
        covariance_type="diag",
        stickiness_kappa=4,
        out_col_name="hmm_state",
        random_state=seed  
    )
    
    # Apply HMM to full dataset
    df_analysis, _, _ = run_sticky_hmm_state_classification(
        df_features=df_analysis,
        feature_cols=descriptive_feature_cols,
        model=hmm_model_descriptive,
        out_col_name="hmm_state"
    )
     
    df_hmm_descriptive["hmm_state"] =  df_hmm_descriptive["hmm_state"].astype(str).astype("category")
    df_analysis["hmm_state"] =  df_analysis["hmm_state"].astype(str).astype("category")
    
    merge_pandas_cols_into_obs_anndata(
        cols=["hmm_state"],
        adata=adata_sub,
        df_analysis=df_analysis,
        on=["sample_name", "TrackID", "position_t"],    
    )
   
    merge_pandas_cols_into_obs_anndata(
        cols=["hmm_state"],
        adata=adata_full,
        df_analysis=df_analysis,
        on=["sample_name", "TrackID", "position_t"],    
    )
    
    sc.pl.umap(adata_sub, color="hmm_state", size=2)

    adata_full = filter_short_state_runs(
        adata_full,
        cluster_key="hmm_state",
        id_cols = ["sample_name", "TrackID"],
        time_key = "position_t",
        length_removed = 1,
        new_key = "hmm_state_filt"
    )
 
    compute_cluster_transition_matrix(
        adata=adata_full,
        cluster_key="hmm_state",
        id_cols = ["sample_name", "TrackID"],
        plot=True
    )
    
    compute_cluster_transition_matrix(
        adata=adata_full,
        cluster_key="hmm_state",
        id_cols = ["sample_name", "TrackID"],
        only_transitions=True,
        plot=True
    )
    
    _,mat = compute_cluster_transition_matrix(
        adata=adata_full,
        cluster_key="hmm_state_filt",
        id_cols = ["sample_name", "TrackID"],
        only_transitions=True,
        plot=True
    )
    
    plot_hmm_transition_matrix(hmm_model_descriptive)

    #### PROJECT ON DESCRIPTIVE DATA
    sc.tl.dendrogram(adata_sub, groupby="hmm_state")
    sc.pl.heatmap(
        adata_sub,
        var_names=descriptive_feature_cols,
        # var_group_labels=list(feature_dict.keys()),
        groupby="hmm_state",
        standard_scale="var",
        figsize=(40, 20),
        swap_axes=True,
        dendrogram=True,
        show_gene_labels=True
    )
    
    sc.pl.matrixplot(
        adata_sub,
        var_names=descriptive_feature_cols,                 # dict: {group_label: [genes...]}
        # var_group_labels=list(feature_dict.keys()),
        groupby="hmm_state",
        standard_scale="var",
        figsize=(20, 20),
        swap_axes=True,
        dendrogram=True
    )
    
    """
    ################################################
    sticky HMM STATE CLASSIFICATION (on raw timepoint data)
    ################################################
    """
    
        
    df_tracks = df_positions[["sample_name", "TrackID", "position_t"]+features].copy()
    nonbinary_cols = select_nonbinary_columnnames(
        df=df_tracks, 
        cols= features
        )
    df_tracks[nonbinary_cols]=df_tracks[nonbinary_cols].astype(float)
    
    df_tracks.loc[:, nonbinary_cols] = StandardScaler().fit_transform(df_tracks.loc[:, nonbinary_cols])
    
    df_sticky_hmm_raw = subset_full_tracks(
        df=df_tracks,
        id_cols=["sample_name", "TrackID"],
        sampled_keys=sampled_keys
    )
    
    # df_sticky_hmm_raw, hsmm_model_raw, selection_df = run_sticky_hmm_state_classification(
    #     df_features=df_sticky_hmm_raw,
    #     feature_cols=features,
    #     k_min = 3,
    #     k_max = 15,
    #     n_states="auto", #6 came out of it
    #     # n_states=10,
    #     out_col_name="hsmm_state_instant",
    #     random_state=seed  
    # )
     
    df_sticky_hmm_raw, sticky_hmm_model_raw, selection_df = run_sticky_hmm_state_classification(
        df_features=df_sticky_hmm_raw,
        feature_cols=features,
        n_states=9,
        out_col_name="hmm_state_raw",
        covariance_type="diag",
        stickiness_kappa=8,
        random_state=seed  
    )
    
    # Apply HMM to full dataset
    df_tracks, _, _ = run_sticky_hmm_state_classification(
        df_features=df_tracks,
        feature_cols=features,
        model=sticky_hmm_model_raw,
        out_col_name="hmm_state_raw"
    )
    
    df_sticky_hmm_raw["hmm_state_raw"] =  df_sticky_hmm_raw["hmm_state_raw"].astype(str).astype("category")
    df_tracks["hmm_state_raw"] =  df_tracks["hmm_state_raw"].astype(str).astype("category")
    
    merge_pandas_cols_into_obs_anndata(
        cols=["hmm_state_raw"],
        adata=adata_sub,
        df_analysis=df_tracks,
        on=["sample_name", "TrackID", "position_t"],    
    )
   
    merge_pandas_cols_into_obs_anndata(
        cols=["hmm_state_raw"],
        adata=adata_full,
        df_analysis=df_tracks,
        on=["sample_name", "TrackID", "position_t"],    
    )
    
    sc.pl.umap(adata_sub, color="hmm_state_raw", size=6)

    adata_sub_hmm_raw = df_to_adata(
        df_sticky_hmm_raw, features, obs_cols=["sample_name", "TrackID", "position_t", "hmm_state_raw"]
    )
    
    #### PROJECT ON DESCRIPTIVE DATA
    sc.tl.dendrogram(adata_sub, groupby="hmm_state")
    sc.pl.heatmap(
        adata_sub,
        var_names=descriptive_feature_cols,
        # var_group_labels=list(feature_dict.keys()),
        groupby="hmm_state_raw",
        standard_scale="var",
        figsize=(40, 20),
        swap_axes=True,
        dendrogram=True,
        show_gene_labels=True
    )
    
    sc.pl.matrixplot(
        adata_sub,
        var_names=descriptive_feature_cols,                 # dict: {group_label: [genes...]}
        # var_group_labels=list(feature_dict.keys()),
        groupby="hmm_state_raw",
        standard_scale="var",
        figsize=(20, 20),
        swap_axes=True,
        dendrogram=True
    )
      
    plot_number_per_clusters(df_hmm_raw, cluster_col="hmm_state")
    plot_hmm_transition_matrix(hmm_model_raw)
    plot_hmm_transition_matrix(sticky_hmm_model_raw)
    
    
    # Contingency table
    pd.crosstab(adata_full.obs['ClusterID'], adata_full.obs['hmm_state_raw'], normalize='index')
    
    
    adata_full.write("/Users/s.deblank-3/Downloads/adata_full.h5ad") 
    adata_sub.write("/Users/s.deblank-3/Downloads/adata_sub.h5ad") 
    
    adata_full.write(r"C:\Users\Samde/Downloads/adata_full.h5ad") 
    adata_sub.write(r"C:\Users\Samde/Downloads/adata_sub.h5ad")
    
    adata = sc.read_h5ad(Path(outfolder,"adata_full.h5ad"))
    adata_sub = sc.read_h5ad(Path(outfolder,"adata_sub.h5ad"))
    adata.write_zarr("my_data.zarr")
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    from behav3d.utils.fileio import load_image, load_zarr, save_as_zarr, append_to_zarr
    import scanpy as sc
    from pathlib import Path
    
    outfolder = r"C:/Users/Samde/Downloads"
    adata = sc.read_h5ad(Path(outfolder,"adata_cluster_from_1.h5ad"))
    
    adata.obs["ClusterID"]=adata.obs["ClusterID"].astype(np.int64)+1
    adata.obs["hmm_state_raw"]=adata.obs["hmm_state_raw"].astype(np.int64)+1
    adata.obs["hmm_state"]=adata.obs["hmm_state"].astype(np.int64)+1
    
    img = load_zarr(r"F:\BHVD_BEHAV3D\BEHAV3D_python\runs\ROCHE\images\ROCHE_JM1_Exp042-8_Img02_10T_HER2I\ROCHE_JM1_Exp042-8_Img02_10T_HER2I_tcell_tracked.zarr")
    raw = load_zarr(r"F:\BHVD_BEHAV3D\BEHAV3D_python\runs\ROCHE\images\ROCHE_JM1_Exp042-8_Img02_10T_HER2I\ROCHE_JM1_Exp042-8_Img02_10T_HER2I.zarr")
    
    hmm_descr_path = Path(outfolder, "hmm_descr_state_overlay.zarr")
    hmm_instant_path=Path(outfolder, "hmm_raw_state_overlay.zarr")
    leiden_path = Path(outfolder, "leiden_overlay.zarr")
    
    hmm_instant_img = relabel_from_adata(
        img,
        adata,
        obs_col="hmm_state",
        sample_name = "ROCHE_JM1_Exp042-8_Img02_10T_HER2I"
    )
    save_as_zarr(hmm_instant_img, hmm_instant_path)

    hmm_desc_img = relabel_from_adata(
        img,
        adata,
        obs_col="hmm_state_raw",
        sample_name = "ROCHE_JM1_Exp042-8_Img02_10T_HER2I"
    )
    save_as_zarr(hmm_desc_img, hmm_descr_path)

    leiden_img = relabel_from_adata(
        img,
        adata,
        obs_col="ClusterID",
        sample_name = "ROCHE_JM1_Exp042-8_Img02_10T_HER2I"
    )
    save_as_zarr(leiden_img, leiden_path)


    leiden_img = load_zarr(leiden_path)
    leiden_img = np.repeat(np.expand_dims(leiden_img, 1), axis=1, repeats=3)
    
    hmm_desc_img = load_zarr(hmm_descr_path)
    hmm_desc_img = np.repeat(np.expand_dims(hmm_desc_img, 1), axis=1, repeats=3)
    
    hmm_instant_img = load_zarr(hmm_instant_path)
    hmm_instant_img = np.repeat(np.expand_dims(hmm_instant_img, 1), axis=1, repeats=3)

    img = np.repeat(np.expand_dims(img, 1), axis=1, repeats=3)
    
    import napari
    viewer = napari.Viewer()
    viewer.add_image(raw)
    viewer.add_labels(img)
    viewer.add_labels(leiden_img, name="leiden")
    viewer.add_labels(hmm_desc_img, name="hmm_descr")
    viewer.add_labels(hmm_instant_img, name="hmm_instant")

    viewer.close()


    plot_cluster_feature_heatmap(
        df_analysis_subset,
        feature_cols=descriptive_feature_cols,
        cluster_col="ClusterID",
        # figsize=(8.27, 11.69),
    )
    
    plot_feature_cluster_heatmap(
        df_analysis_subset,
        feature_cols=descriptive_feature_cols,
        cluster_col="ClusterID",
        figsize=(8.27, 11.69),
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
        out_dir = outfolder,
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
        out_dir=outfolder,
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
        clusters=[1],  
        df_windows=df_analysis_subset,
        df_positions=df_positions,
        output_folder=output_dir,  # folder containing images/<sample>/<sample>.zarr
        out_dir=outfolder,
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
    
    plt.figure(figsize=(6,5))
    ax = sns.scatterplot(
    data=df_analysis_subset, x="UMAP1", y="UMAP2",
    hue="cluster_label_leiden", palette="tab20",
    s=5, alpha=0.3, edgecolor=None  
    )
    # bigger legend dots without changing plot dots
    leg = ax.legend(
        bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.,
        markerscale=6,  # <-- enlarge legend markers
        scatterpoints=1 # keep one dot per legend entry
    )
    plt.title("UMAP colored by Leiden clusters (−1 = noise)")
    plt.tight_layout()
    
    
    
    
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
    
    """
    ################################################
    HMM STATE CLASSIFICATION (on raw timepoint data)
    ################################################
    """
    
    
    df_tracks = df_positions[["sample_name", "TrackID", "position_t"]+features].copy()
    nonbinary_cols = select_nonbinary_columnnames(
        df=df_tracks, 
        cols= features
        )
    df_tracks[nonbinary_cols]=df_tracks[nonbinary_cols].astype(float)
    
    df_tracks.loc[:, nonbinary_cols] = StandardScaler().fit_transform(df_tracks.loc[:, nonbinary_cols])
    
    df_hmm_raw, sampled_keys = subset_full_tracks(
        df=df_tracks,
        id_cols=["sample_name", "TrackID"],
        return_selected_keys=True
        # sampled_keys=sampled_keys
    )
    
    # df_hmm_raw, hmm_model_raw, selection_df = run_hmm_state_classification(
    #     df_features=df_hmm_raw,
    #     feature_cols=features,
    #     k_min = 3,
    #     k_max = 15,
    #     n_states="auto", #6 came out of it
    #     # n_states=10,
    #     out_col_name="hmm_state",
    #     random_state=seed  
    # )
     
    df_hmm_raw, hmm_model_raw, selection_df = run_hmm_state_classification(
        df_features=df_hmm_raw,
        feature_cols=features,
        # k_min = 2,
        # k_max = 12,
        # n_states="auto" #6 came out of it
        n_states=10,
        out_col_name="hmm_state",
        random_state=seed  
    )
    
    # Apply HMM to full dataset
    df_tracks, _, _ = run_hmm_state_classification(
        df_features=df_tracks,
        feature_cols=features,
        model=hmm_model_raw,
        out_col_name="hmm_state"
    )
    
    df_hmm_raw["hmm_state"] =  df_hmm_raw["hmm_state"].astype(str).astype("category")
    df_tracks["hmm_state"] =  df_tracks["hmm_state"].astype(str).astype("category")
    
    merge_pandas_cols_into_obs_anndata(
        cols=["hmm_state"],
        adata=adata_sub,
        df_analysis=df_tracks,
        on=["sample_name", "TrackID", "position_t"],    
    )
   
    merge_pandas_cols_into_obs_anndata(
        cols=["hmm_state"],
        adata=adata_full,
        df_analysis=df_tracks,
        on=["sample_name", "TrackID", "position_t"],    
    )
    
    sc.pl.umap(adata_sub, color="hmm_state", size=6)

    adata_sub_hmm_raw = df_to_adata(
        df_hmm_raw, features, obs_cols=["sample_name", "TrackID", "position_t", "hmm_state"]
    )
    
    sc.tl.dendrogram(adata_sub_hmm_raw, groupby="hmm_state")
    sc.pl.heatmap(
        adata_sub_hmm_raw,
        var_names=features,
        # var_group_labels=list(feature_dict.keys()),
        groupby="hmm_state",
        standard_scale="var",
        figsize=(8.27, 11.69),
        swap_axes=True,
        dendrogram=True,
        show_gene_labels=True
    )
    
    sc.pl.matrixplot(
        adata_sub_hmm_raw,
        var_names=features,                 # dict: {group_label: [genes...]}
        # var_group_labels=list(feature_dict.keys()),
        groupby="hmm_state",
        standard_scale="var",
        figsize=(8.27, 11.69),
        swap_axes=True,
        dendrogram=True
    )
    
    #### PROJECT ON DESCRIPTIVE DATA
    sc.tl.dendrogram(adata_sub, groupby="hmm_state")
    sc.pl.heatmap(
        adata_sub,
        var_names=descriptive_feature_cols,
        # var_group_labels=list(feature_dict.keys()),
        groupby="hmm_state",
        standard_scale="var",
        figsize=(11.69, 8.27),
        swap_axes=True,
        dendrogram=True,
        show_gene_labels=True
    )
    
    sc.pl.matrixplot(
        adata_sub,
        var_names=descriptive_feature_cols,                 # dict: {group_label: [genes...]}
        # var_group_labels=list(feature_dict.keys()),
        groupby="hmm_state",
        standard_scale="var",
        figsize=(8.27, 11.69),
        swap_axes=True,
        dendrogram=True
    )
      
    plot_number_per_clusters(df_hmm_raw, cluster_col="hmm_state")
    plot_hmm_transition_matrix(hmm_model_raw)
    
    
    
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
        

    """
    HMM STATE CLASSIFICATION (on windowed data)
    """
    df_hmm_descriptive = df_analysis.copy()
    # scaled_X = adata_full[:, descriptive_feature_cols].X.toarray()
    # df_analysis_hmm.loc[:, descriptive_feature_cols] = scaled_X

    df_hmm_descriptive, sampled_keys = subset_full_tracks(
        df=df_hmm_descriptive,
        fraction=0.2,
        random_state=seed,
        id_cols=["sample_name", "TrackID"],
        return_selected_keys=True
    )
    
    # df_hmm_descriptive, hmm_model, selection_df = run_hmm_state_classification(
    #     df_features=df_hmm_descriptive,
    #     feature_cols=descriptive_feature_cols,
    #     k_min = 3,
    #     k_max = 15,
    #     n_states="auto", #10 came out of it
    #     out_col_name="hmm_state"
    # )
    
    
    df_hmm_descriptive, hmm_model_descriptive, selection_df = run_hmm_state_classification(
        df_features=df_hmm_descriptive,
        feature_cols=descriptive_feature_cols,
        n_states=7, #6 came out of it
        out_col_name="hmm_state",
        random_state=seed  
    )
    
    # Apply HMM to full dataset
    df_analysis, _, _ = run_hmm_state_classification(
        df_features=df_analysis,
        feature_cols=descriptive_feature_cols,
        model=hmm_model_descriptive,
        out_col_name="hmm_state"
    )
    
    df_hmm_descriptive["hmm_state"] =  df_hmm_descriptive["hmm_state"].astype(str).astype("category")
    df_analysis["hmm_state"] =  df_analysis["hmm_state"].astype(str).astype("category")
    
    merge_pandas_cols_into_obs_anndata(
        cols=["hmm_state"],
        adata=adata_sub,
        df_analysis=df_analysis,
        on=["sample_name", "TrackID", "position_t"],    
    )
   
    merge_pandas_cols_into_obs_anndata(
        cols=["hmm_state"],
        adata=adata_full,
        df_analysis=df_analysis,
        on=["sample_name", "TrackID", "position_t"],    
    )
    
     # hmm_lookup = df_analysis[["sample_name", "TrackID", "position_t"] + ["hmm_state"]].copy()
    # obs_df = adata_full.obs.reset_index().rename(columns={"index": "_obs_index"})
    # obs_df = obs_df.drop(columns=["hmm_state"], errors="ignore")
    # obs_merged = obs_df.merge(hmm_lookup, on=["sample_name", "TrackID", "position_t"], how="left")
    # adata_full.obs["hmm_state"] = (
    #     obs_merged
    #     .set_index("_obs_index")
    #     .loc[adata_full.obs.index, "hmm_state"]
    #     .astype(str)
    #     .astype("category")
    # )
    
    # hmm_lookup = df_analysis_hmm[["sample_name", "TrackID", "position_t"] + ["hmm_state"]].copy()

    sc.pl.umap(adata_sub, color="hmm_state", size=6)
    
    sc.tl.dendrogram(adata_sub, groupby="hmm_state")
    sc.pl.heatmap(
        adata_sub,
        var_names=feature_dict,
        var_group_labels=list(feature_dict.keys()),
        groupby="hmm_state",
        standard_scale="var",
        figsize=(8.27, 11.69),
        swap_axes=True,
        dendrogram=True,
        show_gene_labels=True
    )
    
    sc.pl.matrixplot(
        adata_sub,
        var_names=feature_dict,                 # dict: {group_label: [genes...]}
        # var_group_labels=list(feature_dict.keys()),
        groupby="hmm_state",
        standard_scale="var",
        figsize=(8.27, 11.69),
        swap_axes=True,
        dendrogram=True
    )
    
    plot_number_per_clusters(df_hmm_descriptive, cluster_col="hmm_state")
    plot_hmm_transition_matrix(hmm_model_descriptive)
        
    
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