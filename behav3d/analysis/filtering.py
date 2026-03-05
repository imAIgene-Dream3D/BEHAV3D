
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pathlib import Path
import time
from behav3d.core.utils import format_time
import numpy as np

def filter_by_full_duration(df: pd.DataFrame,
                            exp_duration=None,
                            time_column="position_t"
                            ) -> pd.DataFrame:
    """
    Cut rows where time_column > exp_duration.
    """
    if exp_duration is not None:
        df=df.drop(df[df[time_column]>exp_duration].index)
    return df


def filter_minimal_track_length(
    df,
    min_track_length=None,
    time_column="position_t",
    group_cols=["TrackID", "sample_name"]
    ):

    first_t = (
            df.groupby(group_cols)[time_column]
            .transform("first")
        )

    # --- 2) Min track length on the chosen axis (keep groups with sufficient span) ---
    if min_track_length is not None:
        # span = last - first on the chosen axis
        span = (
            df.groupby(group_cols)[time_column]
            .transform("last") - first_t
        )
        df = df[span >= float(min_track_length)].reset_index(drop=True)
    return df


def trim_to_maximal_track_length(
    df,
    max_track_length=None,
    time_column="position_t",
    group_cols=["TrackID", "sample_name"]
    ):

    if max_track_length is not None:
        last_t = df.groupby(group_cols)[time_column].transform("last")
        within_last_window = (last_t - df[time_column]) <= float(max_track_length)
        df = df[within_last_window].reset_index(drop=True)
    return df


def plot_filter_count(
    df_track_counts,
    outpath,
    nr_cols=3,
    rows_per_page = 3,
    figsize=(8.27, 11.69),
    filter_cols=["nr_tracks_before_filtering", "nr_tracks_exp_duration", "nr_tracks_min_track_length", "nr_tracks_dead_t1"],
    plot_results=False
    ):
    
    sample_names = df_track_counts["sample_name"].unique()
    n_plots = len(sample_names)
    n_rows = (n_plots // nr_cols) + (1 if n_plots % nr_cols != 0 else 0)
    nr_pages = (n_rows // rows_per_page) + (1 if n_rows % rows_per_page != 0 else 0)
    
    with PdfPages(outpath) as pdf:
        plot_idx = 0
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
                    ax.remove()
                    continue
                
                sample = sample_names[plot_idx]
                df_subset = df_track_counts[df_track_counts["sample_name"] == sample]
                
                sns.barplot(
                    x=filter_cols, 
                    y=df_subset[filter_cols].values.flatten(),
                    ax=ax
                )
                
                ax.set_title(sample, fontsize=10, loc='center')
                ax.set_xlabel("")
                ax.set_ylabel("Count")
                
                # Force y-axis to use integer formatting (track counts should always be whole numbers)
                ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
                
                # Use plt.setp for better compatibility with matplotlib
                plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
                
                plot_idx += 1
            
            fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
            # plt.show()
            pdf.savefig(fig, bbox_inches='tight', dpi=600)
            plt.close(fig)


def filter_and_truncate_tracks_anndata(
    adata,
    groupby_cols = ["sample_name", "TrackID"],  # <-- composite track identity
    time_col: str = "position_t",                                  # <-- time/frame column
    min_length: int = 10,                                   # x
    max_length: int = 200                                   # y
):
    missing = [k for k in groupby_cols if k not in adata.obs.columns]
    if missing:
        raise KeyError(f"Missing groupby_cols in adata.obs: {missing}")

    if time_col not in adata.obs.columns:
        raise KeyError(f"'{time_col}' not found in adata.obs. Set time_col to your time column.")

    obs = adata.obs.copy()

    # --- 1) Filter groups with too few positions ---
    # size per group (composite track key)
    group_sizes = obs.groupby(list(groupby_cols), observed=True).size()
    keep_groups = group_sizes[group_sizes >= min_length].index

    # fast membership test: join on multi-index
    obs["_keep"] = obs.set_index(list(groupby_cols)).index.isin(keep_groups)
    obs = obs.loc[obs["_keep"]].copy()

    # --- 2) Keep last max_length per group by time ordering ---
    # stable tie-breaker (in case time_col has ties)
    obs["_orig_idx"] = np.arange(len(obs))
    obs = obs.sort_values(list(groupby_cols) + [time_col, "_orig_idx"])

    # rank within group
    obs["_rank"] = obs.groupby(list(groupby_cols), observed=True).cumcount()
    sizes = obs.groupby(list(groupby_cols), observed=True)["_rank"].transform("max") + 1
    obs = obs.loc[obs["_rank"] >= (sizes - max_length)]

    # subset AnnData with rows in this exact order
    idx = obs.index
    adata_out = adata[idx].copy()

    # cleanup helper cols if they persisted into adata_out.obs
    for col in ["_keep", "_orig_idx", "_rank"]:
        if col in adata_out.obs.columns:
            adata_out.obs.drop(columns=[col], inplace=True)

    return adata_out


def filter_tracks(
    metadata,
    output_dir=None,
    exp_duration=None,
    min_track_length=None,
    max_track_length=None,
    filter_t0_dead=True,
    min_size=None,
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
    
    start_time = time.time()
    
    print(f"--------------- Filtering tracks ---------------")
    print(f"Basing filtering on {time_type}")

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

    # Filter out tracks whose volume at t=1 is below min_size
    if min_size is not None and "volume" in df_all_tracks_filt.columns:
        t1_small = df_all_tracks_filt[
            (df_all_tracks_filt["relative_time"] == 1) &
            (df_all_tracks_filt["volume"] < min_size)
        ][["TrackID", "sample_name"]]
        before = df_all_tracks_filt[["TrackID", "sample_name"]].drop_duplicates().shape[0]
        df_all_tracks_filt = df_all_tracks_filt[
            ~df_all_tracks_filt.set_index(["TrackID", "sample_name"]).index.isin(
                t1_small.set_index(["TrackID", "sample_name"]).index
            )
        ]
        after = df_all_tracks_filt[["TrackID", "sample_name"]].drop_duplicates().shape[0]
        print(f"  Filtered by min size ({min_size}): {before} → {after} tracks")
        df_track_counts = count_tracks(
            df_all_tracks_filt,
            col_name="nr_tracks_min_size",
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

def subset_selection_of_tracks(
    df: pd.DataFrame,
    fraction: float = 0.1,
    id_cols: list[str] = ["sample_name", "TrackID"],
    random_state: int = 0,
    return_selected_keys: bool = False,
    sampled_keys=None,
):
    """
    Randomly subsample a fraction of unique tracks defined by id_cols,
    then return df restricted to those full tracks.

    If `selected_keys` is provided (a DataFrame containing columns `id_cols`),
    no random sampling is done; instead, df is restricted to the tracks
    matching those keys. This is useful for reproducing an earlier selection.
    """
    track_keys = df[id_cols].drop_duplicates()

    if sampled_keys is None:
        sampled_keys = track_keys.sample(frac=fraction, random_state=random_state)
    else:
        sampled_keys = sampled_keys[id_cols].drop_duplicates()

    mask = (
        df[id_cols]
        .merge(sampled_keys, on=id_cols, how="left", indicator=True)["_merge"]
        .eq("both")
    ).to_numpy()
    df_sub = df.loc[mask]

    if return_selected_keys:
        return df_sub, sampled_keys
    else:
        return df_sub

def subset_timepoints_from_tracks(
    df_windowed, 
    step_size, 
    time_col="position_t",
    id_cols=("sample_name", "TrackID"), 
    keep="first"
    ):
    
    df_windowed = df_windowed.dropna()
    step_size = max(int(step_size), 1)
    if step_size == 1:
        return df_windowed

    # ensure stable order within each track
    dfw = df_windowed.sort_values(list(id_cols) + [time_col], kind="mergesort").copy()

    if keep == "first":
        rank = dfw.groupby(list(id_cols), sort=False).cumcount()
        mask = (rank % step_size) == 0
    elif keep == "last":
        rank_from_end = dfw.groupby(list(id_cols), sort=False).cumcount(ascending=False)
        mask = (rank_from_end % step_size) == 0
    else:
        raise ValueError("keep must be 'first' or 'last'")

    return dfw.loc[mask]
    
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
                df_subset = df_tracks[(df_tracks["sample_name"] == sample)].copy()

                # Convert contact to binary: ANY value > 0 means touching (1), else not touching (0)
                # This handles both boolean (True/False) and numeric contact columns
                df_subset['contact_binary'] = (df_subset[contact_column].astype(float) > 0).astype(int)
                
                # Count touching vs not touching
                touching_counts = df_subset['contact_binary'].value_counts().reindex([0, 1], fill_value=0)
                
                # Create bar chart with exactly 2 bars
                bars = ax.bar([0, 1], [touching_counts[0], touching_counts[1]], 
                              color='steelblue', width=0.6, edgecolor='black', linewidth=1)
                
                ax.set_xlabel('Touching', fontsize=9)
                ax.set_ylabel('Count', fontsize=9)
                ax.set_xticks([0, 1])
                ax.set_xticklabels(['No', 'Yes'])
                ax.set_title(f"{sample}", fontsize=10)
                ax.set_xlim(-0.5, 1.5)
                
                # Force integer y-axis
                ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
                
                plot_idx += 1
            # Add a global title (dynamic based on contact_column)
            target_type = contact_column.replace('_contact', '').replace('_', ' ').title()
            fig.suptitle(f"Touching vs. Non-touching {target_type}", fontsize=14, fontweight="bold")

            # Adjust layout and save
            # plt.tight_layout(rect=[0, 0, 1, 0.96])
            fig.subplots_adjust(left=0.05, right=0.85, top=0.95, bottom=0.05)
            # plt.show(fig)
            pdf.savefig(fig, bbox_inches='tight', dpi=600)
            plt.close(fig)
            # plt.savefig("output.pdf", bbox_inches='tight', dpi=600)


def round_legend_ticks(value):
    import numpy as np
    if value <= 1.0:
        # return np.ceil(value * 10) / 10
        return 1.0
    elif value <= 100:
        return np.ceil(value / 10) * 10
    elif value <= 10000:
        return np.ceil(value / 500) * 500
    return value

