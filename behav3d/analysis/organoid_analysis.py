from pathlib import Path
import time

import json
import pandas as pd
import numpy as np
import scanpy as sc

from behav3d.core.utils import format_time
from behav3d.analysis import smooth_value_over_time


from behav3d.preprocessing import calc_z_projection
from behav3d.io.images import _ensure_zarr, load_image, load_zarr

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
import seaborn as sns

from matplotlib.colors import LinearSegmentedColormap
from matplotlib import colormaps
from matplotlib.lines import Line2D
from matplotlib.legend import Legend
from matplotlib.colors import LinearSegmentedColormap as _LSC
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle
from scipy.ndimage import binary_dilation, binary_erosion


from sklearn.preprocessing import MinMaxScaler


def _compute_general_death_dynamics(df_tracks):
    """Per-sample, per-timepoint organoid counts and cumulative fractions.

    Single source of truth shared by per-target Death Dynamics
    (``run_organoid_analysis``) and the cross-target comparison
    (``plot_multi_organoid_death_dynamics``) so the two can never diverge.

    ``df_tracks`` must carry ``sample_name``, ``TrackID``, ``position_t`` and a
    ``dead`` column (bool, or whatever feature extraction wrote it as). The
    ``dead`` column is normalised to boolean on a copy, so the caller's frame
    is left untouched and raw values are accepted.

    Returns ``(df_general, grid, dead_at_t)`` where ``df_general`` has one row
    per (``sample_name``, ``position_t``) with ``nr_alive`` / ``nr_dead`` /
    ``nr_disappeared`` / ``nr_organoids_t0`` and the matching ``percentage_*``
    columns. ``grid`` and ``dead_at_t`` are the per-track intermediates (grid
    carries ``TrackID``) reused by the per-target grouped plots; combined
    callers can ignore them.
    """
    df = df_tracks.copy()

    if pd.api.types.is_bool_dtype(df["dead"]):
        df["dead"] = df["dead"].fillna(False)
    else:
        df["dead"] = (
            df["dead"]
            .astype(str)
            .str.strip()
            .str.lower()
            .isin({"true", "1", "1.0", "yes"})
        )

    # Summarize tracks to their time of death
    summ = (
        df
        .groupby(["sample_name", "TrackID"])
        .agg(
            t_first=("position_t", "min"),
            t_last=("position_t", "max"),
            # first time it is ever dead (NaN if never dead)
            t_dead=("position_t", lambda s: s[df.loc[s.index, "dead"]].min()
                    if (df.loc[s.index, "dead"]).any() else np.nan),
        )
        .reset_index()
    )

    timeline = df[["sample_name", "position_t"]].drop_duplicates()

    grid = timeline.merge(summ, on="sample_name", how="left")
    t = grid["position_t"]

    dead_at_t = grid["t_dead"].notna() & (t >= grid["t_dead"])
    seen_at_t = (t >= grid["t_first"]) & (t <= grid["t_last"])
    never_dead = grid["t_dead"].isna()
    disappeared_by_t = (t > grid["t_last"]) & never_dead
    alive_at_t = seen_at_t & (~dead_at_t)

    counts = (
        grid.assign(alive=alive_at_t, dead=dead_at_t, disappeared=disappeared_by_t)
            .groupby(["sample_name", "position_t"])
            .agg(
                nr_alive=("alive", "sum"),
                nr_dead=("dead", "sum"),
                nr_disappeared=("disappeared", "sum"),
            )
            .reset_index()
    )

    df_t0 = (
        df[df["position_t"] == 0]
        .groupby("sample_name")
        .agg(nr_organoids_t0=("TrackID", "nunique"))
        .reset_index()
    )

    df_general = counts.merge(df_t0, on="sample_name", how="left")
    df_general["percentage_dead"]        = df_general["nr_dead"]        / df_general["nr_organoids_t0"]
    df_general["percentage_alive"]       = df_general["nr_alive"]       / df_general["nr_organoids_t0"]
    df_general["percentage_disappeared"] = df_general["nr_disappeared"] / df_general["nr_organoids_t0"]

    return df_general, grid, dead_at_t


def run_organoid_analysis(
    dead_perc_threshold=None,
    config=None,
    output_dir=None,
    metadata=None,
    df_tracks_path=None,
    org_type="organoid",
    group_cols=None,
    show_in_notebook=True,
    # df_tracks_summarized_path=None,
    ):
    print(f"--------------- Performing {org_type} Death Dynamics analysis ---------------")
    start_time = time.time()
    assert config is not None or all(
        [output_dir]
    ), "Either 'config' or 'output_dir, umap_minimal_distance, umap_n_neighbors, nr_of_clusters' parameters must be supplied"
        
    if output_dir is None:
        output_dir = config['output_dir']
        
    
    analysis_outdir = Path(output_dir, "analysis", org_type)
    feature_outdir = Path(analysis_outdir, "track_features")
    results_outdir = Path(analysis_outdir, "results")
    sample_outdir = Path(results_outdir, "per_sample")
    
    if not sample_outdir.exists():
        sample_outdir.mkdir(parents=True)
    if not analysis_outdir.exists():
        analysis_outdir.mkdir(parents=True)
    if not feature_outdir.exists():
        feature_outdir.mkdir(parents=True)    
    if not results_outdir.exists():
        results_outdir.mkdir(parents=True)    
        
    if df_tracks_path is None:
        df_tracks_path = Path(feature_outdir, f"BEHAV3D_{org_type}_combined_track_features_filtered.csv")
    # if df_tracks_summarized_path is None:
    #     df_tracks_summarized_path = Path(feature_outdir, f"BEHAV3D_organoid_combined_track_features_summarized.csv")
    
    df_tracks = pd.read_csv(df_tracks_path)
    # df_tracks_summarized = pd.read_csv(df_tracks_summarized_path)
    df_tracks=df_tracks.sort_values(by=["sample_name", "TrackID", "relative_time"])
    df_tracks["TrackID"] = df_tracks["TrackID"].astype(str)

    # Attach the user-selected condition-column grouping (e.g. organoid_line +
    # macrophage_line) once, so every downstream consumer (grouped mean/SEM,
    # dead-fraction pooling) sees it without re-merging metadata.
    df_tracks = add_condition_group_column(df_tracks, metadata, group_cols)
    has_condition_groups = group_cols and df_tracks["condition_group"].nunique() > 1
    group_sample_counts = (
        df_tracks.groupby("condition_group")["sample_name"].nunique().to_dict()
        if has_condition_groups else None
    )

    # Death dynamics now reads the final death classification from feature extraction.
    dead_columns = ["nr_dead_mask_pixels", "percentage_dead_mask", "dead"]
    has_dead_data = all(col in df_tracks.columns for col in dead_columns)
    
    if not has_dead_data:
        print(f"\n⚠️  WARNING: Dead channel data not found in features.")
        print(f"   Missing columns: {[col for col in dead_columns if col not in df_tracks.columns]}")
        print(f"   Re-run Feature Extraction and Track Filtering so the final 'dead' column is present in the track features CSV.")
        print(f"   Skipping death dynamics analysis.\n")
        
        # Cannot run death dynamics analysis without dead channel data
        end_time = time.time()
        h,m,s = format_time(start_time, end_time)
        print(f"### DONE (no death data) - elapsed time: {h}:{m:02}:{s:02}\n")
        return
    
    # death analysis (only runs if has_dead_data)
    df_tracks["smoothed_nr_dead_mask_pixels"] = smooth_value_over_time(
            df_tracks,
            column="nr_dead_mask_pixels",
            groupby=["TrackID", "sample_name"]
        )
    
    df_tracks["smoothed_percentage_dead_mask"] = smooth_value_over_time(
            df_tracks,
            column="percentage_dead_mask",
            groupby=["TrackID", "sample_name"]
        )
    
    if pd.api.types.is_bool_dtype(df_tracks["dead"]):
        df_tracks["dead"] = df_tracks["dead"].fillna(False)
    else:
        df_tracks["dead"] = (
            df_tracks["dead"]
            .astype(str)
            .str.strip()
            .str.lower()
            .isin({"true", "1", "1.0", "yes"})
        )

    if "mean_dead_dye" in df_tracks.columns:
        df_tracks["smoothed_mean_dead_dye"] = smooth_value_over_time(
                df_tracks,
                column="mean_dead_dye",
                groupby=["TrackID", "sample_name"]
            )
    else:
        df_tracks["smoothed_mean_dead_dye"] = np.nan

    # Per-sample cumulative dead/alive/disappeared dynamics. Computed by the
    # shared helper so this per-target output stays identical to what the
    # cross-target comparison recomputes from the same filtered tracks.
    # `grid` (carries TrackID) and `dead_at_t` are reused by the grouped
    # condition plot below.
    df_general, grid, dead_at_t = _compute_general_death_dynamics(df_tracks)

    df_general_outpath = Path(results_outdir, f"combined_general_{org_type}_dynamics_analysis.csv")
    df_general.to_csv(
        df_general_outpath,
        sep=",",
        index=False
    )

    # Track-level dead indicator pooled by condition group (page 1's grouped
    # % dead plot). `grid` already carries TrackID (broadcast from `summ`),
    # so no extra death-state computation is needed here.
    condition_color_map = None
    df_dead_at_t_grouped = None
    if has_condition_groups:
        df_dead_at_t = grid[["sample_name", "TrackID", "position_t"]].copy()
        df_dead_at_t["dead_at_t"] = dead_at_t.astype(float)
        group_map = df_tracks[["sample_name", "TrackID", "condition_group"]].drop_duplicates()
        df_dead_at_t = df_dead_at_t.merge(group_map, on=["sample_name", "TrackID"], how="left")
        df_dead_at_t_grouped = _compute_grouped_mean_sem_dynamics(
            df_dead_at_t, value_cols=["dead_at_t"], org_type=org_type
        )
        condition_color_map = _build_group_color_map(df_tracks["condition_group"].dropna().unique())

    # Per-sample AND per-condition-group mean +/- SEM dynamics across
    # organoids (raw + already-smoothed columns, absolute + baseline-from-0).
    # Saved as separate CSVs; plots are appended as extra pages to the PDF.
    raw_meansem_cols = [
        c for c in ("percentage_dead_mask", "nr_dead_mask_pixels", "mean_dead_dye")
        if c in df_tracks.columns
    ]
    meansem_value_cols = raw_meansem_cols + [f"smoothed_{c}" for c in raw_meansem_cols]

    df_meansem = _compute_per_sample_mean_sem_dynamics(df_tracks, value_cols=meansem_value_cols, org_type=org_type)
    df_meansem_grouped = (
        _compute_grouped_mean_sem_dynamics(df_tracks, value_cols=meansem_value_cols, org_type=org_type)
        if has_condition_groups else None
    )

    df_meansem_outpath = Path(
        results_outdir, f"combined_mean_sem_{org_type}_dynamics_analysis.csv"
    )
    df_meansem.to_csv(df_meansem_outpath, sep=",", index=False)
    if df_meansem_grouped is not None:
        df_meansem_grouped_outpath = Path(
            results_outdir, f"combined_mean_sem_grouped_{org_type}_dynamics_analysis.csv"
        )
        df_meansem_grouped.to_csv(df_meansem_grouped_outpath, sep=",", index=False)

    general_pdf_outpath = Path(results_outdir, f"combined_general_{org_type}_dynamics_analysis.pdf")

    plot_general_organoid_analysis(
            df_general=df_general,
            outpath=general_pdf_outpath,
            df_meansem=df_meansem,
            value_cols=raw_meansem_cols,
            org_type=org_type,
            figsize=(8.27, 11.69),
            df_meansem_grouped=df_meansem_grouped,
            df_dead_at_t_grouped=df_dead_at_t_grouped,
            group_cols=group_cols if has_condition_groups else None,
            color_map=condition_color_map,
            group_sample_counts=group_sample_counts,
            show_in_notebook=show_in_notebook,
            )
    
    ## TODO PLOT A STACKED BARPLOT OVER TIME WHERE THE TOTAL BAR IS 
    # THE DETECTED ORGANOIDS AT THAT TIMEPOINT. 
    # STACK ALIVE AND DEAD. 
    # THIS CAN SHOW IF ORGANOIDS DISSAPPEAR OR ACTUALLY DIE OVER TIME
    
    for sample_name in df_tracks["sample_name"].unique():
        """
        Create a analysis pdf for each sample separately
        """
        analysis_sample_outdir = Path(sample_outdir, sample_name, org_type)
        if not analysis_sample_outdir.exists():
            analysis_sample_outdir.mkdir(parents=True)
            
        df_tracks_sample = df_tracks[df_tracks["sample_name"] == sample_name]
        
        img_outdir = Path(output_dir, "images", sample_name)

        sample_metadata = pd.DataFrame()
        if metadata is not None:
            sample_metadata = metadata[metadata['sample_name'] == sample_name]

        df_tracks_sample_outpath = Path(analysis_sample_outdir, f"{sample_name}_{org_type}_track_analysis.csv")
        df_tracks_sample.to_csv(
            df_tracks_sample_outpath,
            sep=",",
            index=False
        )
        
        df_experiment_sample = df_tracks_sample.groupby(["sample_name", "position_t"]).agg(
            nr_organoids_t0=("TrackID", lambda x: x.nunique()),
            nr_dead=("TrackID", lambda x: x[df_tracks_sample.loc[x.index, "dead"]].nunique()),
        ).reset_index()
        df_experiment_sample["nr_alive"]=df_experiment_sample["nr_organoids_t0"] - df_experiment_sample["nr_dead"]
        df_experiment_sample["percentage_dead"]=df_experiment_sample["nr_dead"] / df_experiment_sample["nr_organoids_t0"]
        df_experiment_sample["percentage_alive"]= 1.0 - df_experiment_sample["percentage_dead"]
        
        df_experiment_outpath = Path(analysis_sample_outdir, f"{sample_name}_{org_type}_general_analysis.csv")
        df_experiment_sample.to_csv(
            df_experiment_outpath,
            sep=",",
            index=False
        )
        
        sample_pdf_outpath = Path(analysis_sample_outdir, f"{sample_name}_{org_type}_analysis.pdf")
    
        # raw_img_path: prefer metadata raw_image_path, fall back to constructed path
        raw_img_path = Path(img_outdir, f"{sample_name}.zarr")
        if not sample_metadata.empty and 'raw_image_path' in sample_metadata.columns:
            raw_path_val = sample_metadata['raw_image_path'].values[0]
            if pd.notna(raw_path_val) and str(raw_path_val).strip():
                raw_img_path = Path(raw_path_val)

        organoid_seg_path = Path(img_outdir, f"{sample_name}_{org_type}_tracked.zarr")

        # dead_mask_path: prefer metadata, fall back to constructed path
        dead_mask_path = None
        if not sample_metadata.empty and 'dead_mask_path' in sample_metadata.columns:
            mask_path_val = sample_metadata['dead_mask_path'].values[0]
            if pd.notna(mask_path_val) and str(mask_path_val).strip():
                dead_mask_path = Path(mask_path_val)
                if dead_mask_path.exists():
                    print(f"  Using dead_mask_path from metadata: {dead_mask_path}")
                else:
                    print(f"  ⚠️  dead_mask_path in metadata does not exist: {dead_mask_path}")
                    dead_mask_path = None

        # Fallback: construct path if not found in metadata
        if dead_mask_path is None:
            fallback_path = Path(img_outdir, f"{sample_name}_mask_dead.zarr")
            if fallback_path.exists():
                dead_mask_path = fallback_path
                print(f"  Using fallback dead_mask_path: {dead_mask_path}")
            else:
                print(f"  ⚠️  Fallback dead_mask_path does not exist: {fallback_path}")

        # Skip PDF plotting if dead_mask_path is still None
        if dead_mask_path is None:
            has_dead_channel = (
                not sample_metadata.empty
                and "dead_channel" in sample_metadata.columns
                and sample_metadata["dead_channel"].notna().any()
            )
            print(f"  ⚠️  Skipping PDF generation for {sample_name} - no dead_mask_path available")
            print(f"      (has_dead_channel={has_dead_channel}, metadata columns: {list(sample_metadata.columns) if not sample_metadata.empty else 'N/A'})")
            continue

        ### Plot analysis results per sample
        print(f"- Writing analysis pdf to {sample_pdf_outpath}")
        plot_sample_organoid_analysis(
            df_experiment_sample=df_experiment_sample,
            df_sample_tracks=df_tracks_sample,
            outpath=sample_pdf_outpath,
            raw_img_path=raw_img_path,
            organoid_seg_path=organoid_seg_path,
            dead_mask_path=dead_mask_path,
            dead_perc_threshold=dead_perc_threshold,
            figsize=(8.27, 11.69)
            )
    
    end_time = time.time()
    h,m,s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    return()




def plot_sample_organoid_analysis(
    df_experiment_sample,
    df_sample_tracks,
    outpath,
    raw_img_path,
    organoid_seg_path,
    dead_mask_path,
    dead_perc_threshold=None,
    figsize=(8.27, 11.69)
    ):
    _ensure_zarr(raw_img_path, label="Raw image")
    _ensure_zarr(organoid_seg_path, label="Organoid segmentation")
    _ensure_zarr(dead_mask_path, label="Dead mask")

    raw_img = load_image(raw_img_path)
    seg_img = load_image(organoid_seg_path)
    dead_mask = load_image(dead_mask_path)
    
    raw_img = calc_z_projection(raw_img, projection="max", z_axis=-3)
    seg_img = calc_z_projection(seg_img, projection="max", z_axis=-3)
    dead_mask = calc_z_projection(dead_mask, projection="max", z_axis=-3)
    
    # Use segmentation image time dimension as reference (since it was actually processed)
    seg_timepoints = seg_img.shape[0]
    
    # Truncate raw image to match segmentation timepoints
    raw_img = raw_img[:seg_timepoints]
    
    t_middle = raw_img.shape[0] // 2
    
    t_start_img = np.asarray(raw_img[0])
    t_mid_img = np.asarray(raw_img[t_middle])
    t_end_img = np.asarray(raw_img[-1])
    
    t_start_img_seg = np.asarray(seg_img[0])
    t_mid_img_seg = np.asarray(seg_img[t_middle])
    t_end_img_seg = np.asarray(seg_img[-1])

    t_start_img_dead = np.asarray(dead_mask[0])
    t_mid_img_dead = np.asarray(dead_mask[t_middle])
    t_end_img_dead = np.asarray(dead_mask[-1])
    
    colors = ['black'] + [colormaps.get_cmap('tab20')(i) for i in range(1, 20)]*10
    cmap = LinearSegmentedColormap.from_list('custom_tab20', colors, 200)
    
    with PdfPages(outpath) as pdf:
        fig = plt.figure(figsize=figsize)
        gs = GridSpec(4, 3, figure=fig, wspace=0.1, hspace=0.3)
        
        ax = fig.add_subplot(gs[0, 0])
        ax.imshow(t_start_img[1], cmap="grey", vmin=0, vmax=np.percentile(t_start_img[1], 99))
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_ylabel('Live Channel', fontsize=10)
        ax.set_title(f'T0')
        
        ax = fig.add_subplot(gs[0, 1])
        ax.imshow(t_mid_img[1], cmap="grey", vmin=0, vmax=np.percentile(t_mid_img[1], 99))
        ax.axis('off')
        ax.set_title(f'T{t_middle}')
        
        ax = fig.add_subplot(gs[0, 2])
        ax.imshow(t_end_img[1], cmap="gray", vmin=0, vmax=np.percentile(t_end_img[1], 99))
        ax.axis('off')
        ax.set_title(f'T{raw_img.shape[0]}')

        ax = fig.add_subplot(gs[1, 0])
        ax.imshow(t_start_img_dead, cmap="grey")
        ax.imshow(np.ma.masked_where(t_start_img_seg == 0, t_start_img_seg), cmap=cmap, alpha=0.5, interpolation='none')
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_ylabel('Dead Mask / Segmentation', fontsize=10)
        # ax.set_title(f'T0')
        
        ax = fig.add_subplot(gs[1, 1])
        ax.imshow(t_mid_img_dead, cmap="grey")
        ax.imshow(np.ma.masked_where(t_mid_img_seg == 0, t_mid_img_seg), cmap=cmap, alpha=0.5, interpolation='none')
        ax.axis('off')
        # ax.set_title(f'T{t_middle} - Dead channel/Segmentation')
        
        ax = fig.add_subplot(gs[1, 2])
        ax.imshow(t_end_img_dead, cmap="gray")
        ax.imshow(np.ma.masked_where(t_end_img_seg == 0, t_end_img_seg), cmap=cmap, alpha=0.5, interpolation='none')
        ax.axis('off')
        # ax.set_title(f'T{raw_img.shape[0]} - Dead channel/Segmentation')
        
        max_alive = df_experiment_sample["nr_alive"].max()
        ymax = max_alive + 1.0
        
        ax = fig.add_subplot(gs[2:, :])
        sns.lineplot(
        data = df_experiment_sample,
        x = 'position_t', 
        y = 'nr_alive', 
        ax = ax
        )
        ax.grid(True, linestyle=':', linewidth=1, alpha=0.5)
        ax.set_ylim(0, ymax)
        ax.set_xlim(0)
        ax.set_ylabel('')
        ax.set_title(f'Alive Organoids')
        
        # Secondary Y axis (percent scale)
        ax_percent = ax.twinx()
        ax_percent.set_ylim(0, 100 * ymax / max_alive)
        ax_percent.set_ylabel('% Alive')

        fig.subplots_adjust(left=0.05, right=0.85, top=0.95, bottom=0.05)
        pdf.savefig(fig)
        plt.close(fig)
        
        fig = plt.figure(figsize=figsize)
        gs = GridSpec(2, 2, figure=fig, wspace=0.5, hspace=0.3)
        
        ax = fig.add_subplot(gs[0, 0])
        ax = plot_feature_lineplot(
            df_sample_tracks,
            feature = 'smoothed_nr_dead_mask_pixels',
            title = f'# Dead Mask Pixels (smoothed)',
            ax = ax,
        )
        
        ax = fig.add_subplot(gs[0, 1])
        ax = plot_feature_lineplot(
            df_sample_tracks,
            feature = 'smoothed_percentage_dead_mask',
            title = f'# Percentage Dead Organoid (smoothed)',
            ax = ax,
        )
        if dead_perc_threshold is not None:
            ax.axhline(dead_perc_threshold, color='black', linestyle=':', linewidth=1)
            # Add small label on the right
            ax.text(
                x=ax.get_xlim()[1],  # far right of x-axis
                y=dead_perc_threshold,
                s=f'  Dead Threshold ({dead_perc_threshold})',
                va='center', ha='left',
                fontsize=5, color='black',
            )
            
        ax = fig.add_subplot(gs[1, 0])
        ax = plot_feature_lineplot(
            df_sample_tracks,
            feature = 'smoothed_mean_dead_dye',
            title = f'# Mean Dead Dye (smoothed)',
            ax = ax,
        )
        
        ax = fig.add_subplot(gs[1, 1])
        ax = plot_feature_lineplot(
            df_sample_tracks,
            feature = 'volume',
            title = f'Volume',
            ax = ax,
        )
        
        fig.subplots_adjust(left=0.05, right=0.85, top=0.95, bottom=0.05)
        
        # plt.show()
        pdf.savefig(fig)
        plt.close(fig)
    
def plot_feature_lineplot(
        df_tracks,
        feature,
        title,
        ax,
        t_colname = 'position_t',
        hue = "TrackID"
    ):
    sns.lineplot(
        data = df_tracks,
        x = t_colname, 
        y = feature, 
        hue=hue, 
        ax = ax
        )
    ax.set_ylabel('')
    ax.set_title(title)
    ax.legend(
        bbox_to_anchor=(1.05, 1),  # Position to the right
        loc='upper left',
        borderaxespad=0,
        prop={'size': 5},  # Smaller font size
        title=hue  # Optional: keep or remove the title
    ).set_title(hue, prop={'size': 6})
    return ax

def add_condition_group_column(df, metadata, group_cols, sample_col="sample_name", label_col="condition_group"):
    """
    Merge selected metadata condition column(s) onto ``df`` (by ``sample_col``)
    and combine them into a single grouping label column ``label_col``, joining
    values with " | ". Missing values become "(unknown)".

    Degrades to a constant "(all)" group when ``group_cols``/``metadata`` are
    not usable, so callers can always group by ``label_col``.
    """
    df = df.copy()
    group_cols = [c for c in (group_cols or []) if metadata is not None and c in getattr(metadata, "columns", [])]
    if not group_cols or metadata is None or sample_col not in metadata.columns:
        df[label_col] = "(all)"
        return df

    # Some track-feature CSVs (e.g. those produced by Filtering) already have
    # the grouping columns merged in from metadata. Only merge columns that
    # are missing to avoid pandas silently renaming collisions to _x/_y.
    missing_cols = [c for c in group_cols if c not in df.columns]
    if missing_cols:
        meta_subset = metadata[[sample_col] + missing_cols].drop_duplicates(subset=[sample_col])
        df = df.merge(meta_subset, on=sample_col, how="left")

    parts = [df[c].fillna("(unknown)").astype(str) for c in group_cols]
    label = parts[0]
    for p in parts[1:]:
        label = label + " | " + p
    df[label_col] = label
    return df


def _sem(x):
    """Standard error of the mean: std(ddof=1)/sqrt(n); 0 when n <= 1."""
    x = x.dropna()
    n = len(x)
    if n > 1:
        return x.std(ddof=1) / np.sqrt(n)
    return 0.0


def _compute_per_sample_mean_sem_dynamics(df_tracks, value_cols, org_type=None):
    """
    Aggregate per-sample mean +/- SEM dynamics across organoid tracks.

    For each column in ``value_cols`` present in ``df_tracks`` (pass both raw
    and ``smoothed_*`` column names to get both variants - smoothing already
    happened per-track upstream via ``smooth_value_over_time``, so no
    additional cohort-level smoothing is applied here):

    1. Aggregate across tracks per (sample_name, position_t):
        ``<col>_mean``  = mean of ``<col>`` across organoids at time t.
        ``<col>_sem``   = SEM of ``<col>`` across organoids at time t.

    2. Per-sample baseline-subtract (start each curve at 0). The SEM is
        unchanged by a constant shift, so we just copy it:
        ``<col>_from0_mean`` = ``<col>_mean`` - per-sample first-timepoint mean.
        ``<col>_from0_sem``  = ``<col>_sem``.

    Returns a tidy DataFrame keyed on (sample_name, position_t).
    Returns an empty DataFrame if no value columns are available.
    """
    time_col = "position_t"
    value_cols = [c for c in value_cols if c in df_tracks.columns]
    if not value_cols:
        return pd.DataFrame(columns=["sample_name", time_col])

    df = df_tracks.copy()
    df["TrackID"] = df["TrackID"].astype(str)

    # ---- 1. cohort aggregate per (sample, timepoint) -------------------
    g = (
        df.groupby(["sample_name", time_col])[value_cols]
          .agg(["mean", _sem])
          .reset_index()
    )
    g.columns = ["_".join([c for c in col if c]).strip("_") for col in g.columns]
    # function name "_sem" produces "<col>__sem"; collapse to "<col>_sem"
    g = g.rename(columns={f"{c}__sem": f"{c}_sem" for c in value_cols})
    g = g.sort_values(["sample_name", time_col]).reset_index(drop=True)

    # ---- 2. baseline-from-0 (per sample). SEM is constant-shift invariant
    for col in value_cols:
        base = g.groupby("sample_name")[f"{col}_mean"].transform("first")
        g[f"{col}_from0_mean"] = g[f"{col}_mean"] - base
        g[f"{col}_from0_sem"] = g[f"{col}_sem"]

    return g.sort_values(["sample_name", time_col]).reset_index(drop=True)


def _compute_grouped_mean_sem_dynamics(df, value_cols, group_col="condition_group", org_type=None):
    """
    Aggregate mean +/- SEM dynamics per (group_col, position_t), pooling
    across ALL rows (organoid tracks) that fall in each group - not
    per-sample averages. This matches the pooling convention used by the
    multi-organoid ``plot_grouped_dead_signal``: n at each timepoint is the
    number of contributing organoid tracks, and SEM = std / sqrt(n).

    Returns a tidy DataFrame keyed on (group_col, position_t). Returns an
    empty DataFrame if no value columns are available.
    """
    time_col = "position_t"
    value_cols = [c for c in value_cols if c in df.columns]
    if not value_cols or group_col not in df.columns:
        return pd.DataFrame(columns=[group_col, time_col])

    g = (
        df.groupby([group_col, time_col])[value_cols]
          .agg(["mean", _sem])
          .reset_index()
    )
    g.columns = ["_".join([c for c in col if c]).strip("_") for col in g.columns]
    g = g.rename(columns={f"{c}__sem": f"{c}_sem" for c in value_cols})
    return g.sort_values([group_col, time_col]).reset_index(drop=True)


def _build_group_color_map(levels):
    """Build a consistent tab10/husl color map for a set of group labels."""
    levels = sorted(set(str(l) for l in levels))
    n = len(levels)
    palette = sns.color_palette("tab10", n) if n <= 10 else sns.color_palette("husl", max(n, 1))
    return dict(zip(levels, palette))


def _plot_grouped_mean_sem_panel(
    ax,
    df_grouped,
    mean_col,
    sem_col,
    group_col="condition_group",
    title=None,
    ylabel=None,
    xlabel="Timepoint",
    legend_title=None,
    color_map=None,
    axhline_zero=False,
    group_sample_counts=None,
    suppress_single_sample_sem=True,
):
    """
    Draw one mean +/- SEM line per group into ``ax``, where each group's
    mean/SEM was computed by pooling across organoid tracks (see
    ``_compute_grouped_mean_sem_dynamics``) - the single-organoid analog of
    the multi-organoid grouped mean +/- SEM plots.

    ``group_sample_counts`` (optional dict of group -> number of distinct
    samples backing it) appends a "N sample(s)" suffix to each legend entry.
    When ``suppress_single_sample_sem`` is True (default), it also suppresses
    the SEM band for groups backed by only 1 sample - appropriate for metrics
    that collapse to a single value per sample (e.g. fraction dead), where a
    single biological replicate can't support a between-replicate error
    estimate. Metrics that vary per organoid (e.g. dead-dye intensity) still
    have a meaningful SEM across organoids within one sample, so callers for
    those should pass ``suppress_single_sample_sem=False``.
    """
    groups = sorted(df_grouped[group_col].dropna().astype(str).unique())
    if color_map is None:
        color_map = _build_group_color_map(groups)

    for grp in groups:
        s = df_grouped[df_grouped[group_col].astype(str) == grp].sort_values("position_t")
        if s.empty or s[mean_col].dropna().empty:
            continue
        color = color_map.get(grp, "C0")
        n_samples = (group_sample_counts or {}).get(grp)
        label = grp if n_samples is None else f"{grp} | {n_samples} sample{'s' if n_samples != 1 else ''}"
        ax.plot(s["position_t"], s[mean_col], linewidth=2.4, label=label, color=color, zorder=3)
        if not suppress_single_sample_sem or n_samples is None or n_samples >= 2:
            ax.fill_between(
                s["position_t"],
                s[mean_col] - s[sem_col],
                s[mean_col] + s[sem_col],
                alpha=0.2, color=color, linewidth=0, zorder=1,
            )

    if axhline_zero:
        ax.axhline(0, color="black", linewidth=0.6, alpha=0.35, zorder=0)

    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel or "", fontsize=11)
    ax.set_title(title or "", fontsize=12)
    ax.grid(True, linestyle=":", linewidth=0.7, alpha=0.35)
    ax.legend(
        bbox_to_anchor=(1.02, 1.0), loc="upper left", borderaxespad=0,
        prop={"size": 9}, title=legend_title or group_col, title_fontsize=10,
    )
    return ax


def _plot_per_sample_mean_sem_panel(
    ax,
    df_meansem,
    mean_col,
    sem_col,
    title=None,
    ylabel=None,
    xlabel="Timepoint",
    axhline_zero=False,
):
    """Draw one mean +/- SEM line per sample_name into ``ax``."""
    samples = sorted(df_meansem["sample_name"].dropna().unique().tolist())
    for sample in samples:
        s = df_meansem[df_meansem["sample_name"] == sample].sort_values("position_t")
        if s.empty or s[mean_col].dropna().empty:
            continue
        line, = ax.plot(s["position_t"], s[mean_col], linewidth=2, label=sample)
        ax.fill_between(
            s["position_t"], s[mean_col] - s[sem_col], s[mean_col] + s[sem_col],
            alpha=0.2, color=line.get_color(), linewidth=0,
        )

    if axhline_zero:
        ax.axhline(0, color="k", linewidth=0.6, alpha=0.4)

    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel or "", fontsize=11)
    ax.set_title(title or "", fontsize=12)
    ax.grid(alpha=0.3)
    ax.legend(
        bbox_to_anchor=(1.02, 1.0), loc="upper left", borderaxespad=0,
        prop={"size": 8}, title="Sample Name", title_fontsize=9,
    )
    return ax


_NICE_FEATURE_NAMES = {
    "percentage_dead_mask": "Dead-Mask Area Coverage",
    "nr_dead_mask_pixels": "Dead-Mask Pixel Count",
    "mean_dead_dye": "Dead-Dye Intensity",
}


def _nice_feature_name(col):
    """Human-readable label for a raw feature column name."""
    if col in _NICE_FEATURE_NAMES:
        return _NICE_FEATURE_NAMES[col]
    return col.replace("_", " ").strip().capitalize()


def _pretty_group_col_name(col):
    """Shorten a metadata condition-column name for legend titles, e.g.
    'or_organoid_line_condition' -> 'organoid', 'im_macrophage_line_condition'
    -> 'macrophage'."""
    name = col
    for prefix in ("or_", "im_", "ot_"):
        if name.startswith(prefix):
            name = name[len(prefix):]
            break
    if name.endswith("_line_condition"):
        name = name[: -len("_line_condition")]
    return name.replace("_", " ") or col


def _group_legend_title(group_cols):
    """Build a short legend title from the selected condition column(s)."""
    if not group_cols:
        return "Condition Group"
    return " / ".join(_pretty_group_col_name(c) for c in group_cols)


def _plot_feature_pages(
    pdf,
    df_meansem_persample,
    df_meansem_grouped,
    feature_cols,
    group_col="condition_group",
    group_label_title=None,
    color_map=None,
    org_type=None,
    figsize=(8.27, 11.69),
    group_sample_counts=None,
):
    """
    Append 1 A4 page per feature in ``feature_cols`` (RAW only - smoothed
    pages are hidden for now) to an open ``PdfPages`` object. Each page has
    3 stacked panels:
      (a) grouped mean +/- SEM, absolute (pooled across organoids per
          ``group_col``);
      (b) per-sample mean +/- SEM, absolute;
      (c) per-sample mean +/- SEM, baseline-shifted to start at 0.

    Unlike the page-1 fraction-dead panel, these features vary per organoid,
    so panel (a)'s SEM band is always shown (even for a single sample) as
    long as that sample has more than one organoid.

    Figures are constructed with ``matplotlib.figure.Figure`` +
    ``FigureCanvasAgg`` (NOT through ``pyplot``) so they never render inline
    in notebooks - only saved to the PDF.
    """
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    if not feature_cols:
        return

    org_label = f" [{org_type}]" if org_type else ""
    legend_title = group_label_title or group_col

    for col in feature_cols:
        nice_name = _nice_feature_name(col)
        for variant, source_col in (("raw", col),):
            mean_c, sem_c = f"{source_col}_mean", f"{source_col}_sem"
            from0_mean_c, from0_sem_c = f"{source_col}_from0_mean", f"{source_col}_from0_sem"

            has_persample = (
                df_meansem_persample is not None
                and mean_c in df_meansem_persample.columns
            )
            has_grouped = (
                df_meansem_grouped is not None
                and mean_c in df_meansem_grouped.columns
            )
            if not has_persample and not has_grouped:
                continue

            fig = Figure(figsize=figsize)
            FigureCanvasAgg(fig)
            gs = GridSpec(3, 1, figure=fig, hspace=0.65)
            variant_label = "Raw" if variant == "raw" else "Smoothed"
            fig.suptitle(
                f"{nice_name} Over Time (Per Organoid) — {variant_label}{org_label}",
                fontsize=13, y=0.98, fontweight="bold",
            )

            ax_a = fig.add_subplot(gs[0])
            if has_grouped:
                _plot_grouped_mean_sem_panel(
                    ax_a, df_meansem_grouped, mean_c, sem_c, group_col=group_col,
                    title="Grouped mean +/- SEM (pooled across organoids)", ylabel=col,
                    legend_title=legend_title, color_map=color_map, axhline_zero=False,
                    group_sample_counts=group_sample_counts,
                    suppress_single_sample_sem=False,
                )
            else:
                ax_a.axis("off")
                ax_a.text(0.5, 0.5, "No condition groups selected", ha="center", va="center", fontsize=10)

            ax_b = fig.add_subplot(gs[1])
            if has_persample:
                _plot_per_sample_mean_sem_panel(
                    ax_b, df_meansem_persample, mean_c, sem_c,
                    title="Per-sample mean +/- SEM (absolute)", ylabel=col, axhline_zero=False,
                )

            ax_c = fig.add_subplot(gs[2])
            if has_persample and from0_mean_c in df_meansem_persample.columns:
                _plot_per_sample_mean_sem_panel(
                    ax_c, df_meansem_persample, from0_mean_c, from0_sem_c,
                    title="Per-sample mean +/- SEM (from baseline)", ylabel=f"{col} - baseline", axhline_zero=True,
                )

            right_margin = 0.55 if has_grouped else 0.78
            fig.subplots_adjust(left=0.10, right=right_margin, top=0.91, bottom=0.06)
            pdf.savefig(fig)


def _draw_persample_percentage_dead_panel(ax, df_general):
    """Panel: % dead over time, one line per sample_name."""
    sns.lineplot(
        data=df_general,
        x='position_t',
        y='percentage_dead',
        hue='sample_name',
        ax=ax,
        )
    ax.grid(True, linestyle=':', linewidth=1, alpha=0.5)
    ax.set_ylim(0, 1.1)
    ax.set_xlim(0)
    ax.set_xlabel('Timepoint', fontsize=11)
    ax.set_ylabel('')
    ax.set_title('Percentage Dead Organoids (per sample)', fontsize=12)
    ax.legend(
        bbox_to_anchor=(1.02, 1.0),
        loc='upper left',
        borderaxespad=0,
        prop={'size': 9},
        title='Sample Name',
        title_fontsize=10,
        )


def _draw_end_state_barplot_panel(ax, df_general):
    """Panel: stacked bar of alive/dead/disappeared at each sample's endpoint."""
    idx = df_general.groupby("sample_name")["position_t"].idxmax()
    df_end = (
        df_general.loc[idx, ["sample_name", "percentage_alive", "percentage_dead", "percentage_disappeared"]]
        .copy()
        .sort_values("sample_name")
    )

    if "percentage_dead" not in df_end.columns:
        df_end["percentage_dead"] = 1.0 - df_end["percentage_alive"]

    if df_end.empty:
        ax.text(0.5, 0.5, "No organoids detected at experiment end",
                ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
        return

    plot_cols = ["percentage_alive", "percentage_dead", "percentage_disappeared"]

    df_end.set_index("sample_name")[plot_cols].plot(
        kind="bar",
        stacked=True,
        ax=ax,
        color=["#6699CC", "#CC6666", "#898989"],  # Blue=Alive, Red=Dead, Grey=Disappeared
        alpha=1.0,
        width=0.25,
        zorder=2,
    )

    ax.grid(True, linestyle=":", linewidth=1, alpha=0.5, zorder=1)
    ax.set_ylim(0, 1.1)
    ax.set_xlabel("")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", size=9)
    ax.set_ylabel("")
    ax.set_title("Percentage Alive Organoids at end of experiment", fontsize=12)

    handles, labels = ax.get_legend_handles_labels()
    order = [2, 0, 1]  # 0=Dead (red), 1=Alive (blue), 2=Disappeared (grey)
    nice_labels = {"percentage_alive": "Alive", "percentage_dead": "Dead", "percentage_disappeared": "Disappeared"}
    ax.legend(
        [handles[i] for i in order],
        [nice_labels.get(labels[i], labels[i]) for i in order],
        bbox_to_anchor=(1.02, 1.0),
        loc="upper left",
        borderaxespad=0,
        prop={"size": 9},
        title="",
    )


def plot_general_organoid_analysis(
    df_general,
    outpath,
    df_meansem=None,
    value_cols=None,
    org_type=None,
    figsize=(8.27, 11.69),
    df_meansem_grouped=None,
    df_dead_at_t_grouped=None,
    group_cols=None,
    color_map=None,
    group_sample_counts=None,
    show_in_notebook=True,
    ):

    has_grouped_dead = df_dead_at_t_grouped is not None and len(df_dead_at_t_grouped) > 0
    legend_title = _group_legend_title(group_cols)
    org_label = f" [{org_type}]" if org_type else ""

    n_groups = df_dead_at_t_grouped["condition_group"].nunique() if has_grouped_dead else 0
    n_samples = df_general["sample_name"].nunique()
    # Keep all 3 page-1 plots on one A4 page when legends fit; split onto a
    # second physical page once there are enough groups/samples that the
    # stacked legends would start crowding/overlapping.
    split_page1 = has_grouped_dead and (n_groups > 8 or n_samples > 15)

    with PdfPages(outpath) as pdf:
        if has_grouped_dead and split_page1:
            fig1 = plt.figure(figsize=figsize)
            fig1.suptitle(f"Death Dynamics Overview — Grouped by Condition{org_label}", fontsize=14, fontweight="bold", y=0.98)
            ax_top = fig1.add_subplot(GridSpec(1, 1, figure=fig1)[0])
            _plot_grouped_mean_sem_panel(
                ax_top, df_dead_at_t_grouped, "dead_at_t_mean", "dead_at_t_sem",
                title="Percentage Dead Organoids (grouped, pooled across organoids)",
                ylabel="Fraction Dead", legend_title=legend_title, color_map=color_map,
                group_sample_counts=group_sample_counts,
            )
            ax_top.set_ylim(0, 1.05)
            ax_top.set_xlim(0)
            fig1.subplots_adjust(left=0.10, right=0.55, top=0.90, bottom=0.08)
            if show_in_notebook:
                plt.show()
            pdf.savefig(fig1, bbox_inches="tight")
            plt.close(fig1)

            fig2 = plt.figure(figsize=figsize)
            fig2.suptitle(f"Death Dynamics Overview — Per Sample & End State{org_label}", fontsize=14, fontweight="bold", y=0.98)
            gs2 = GridSpec(2, 1, figure=fig2, height_ratios=[1.4, 1.0], hspace=0.5)
            _draw_persample_percentage_dead_panel(fig2.add_subplot(gs2[0]), df_general)
            _draw_end_state_barplot_panel(fig2.add_subplot(gs2[1]), df_general)
            fig2.subplots_adjust(left=0.10, right=0.78, top=0.90, bottom=0.10)
            if show_in_notebook:
                plt.show()
            pdf.savefig(fig2, bbox_inches="tight")
            plt.close(fig2)
        else:
            fig = plt.figure(figsize=figsize)
            fig.suptitle(f"Death Dynamics Overview{org_label}", fontsize=14, fontweight="bold", y=0.98)
            if has_grouped_dead:
                gs = GridSpec(3, 1, figure=fig, height_ratios=[1.3, 1.3, 1.0], hspace=0.6)
                _plot_grouped_mean_sem_panel(
                    fig.add_subplot(gs[0]), df_dead_at_t_grouped, "dead_at_t_mean", "dead_at_t_sem",
                    title="Percentage Dead Organoids (grouped, pooled across organoids)",
                    ylabel="Fraction Dead", legend_title=legend_title, color_map=color_map,
                    group_sample_counts=group_sample_counts,
                )
                fig.axes[-1].set_ylim(0, 1.05)
                fig.axes[-1].set_xlim(0)
                _draw_persample_percentage_dead_panel(fig.add_subplot(gs[1]), df_general)
                _draw_end_state_barplot_panel(fig.add_subplot(gs[2]), df_general)
                right_margin = 0.55
            else:
                gs = GridSpec(2, 1, figure=fig, height_ratios=[1.3, 1.0], hspace=0.5)
                _draw_persample_percentage_dead_panel(fig.add_subplot(gs[0]), df_general)
                _draw_end_state_barplot_panel(fig.add_subplot(gs[1]), df_general)
                right_margin = 0.78

            fig.subplots_adjust(left=0.10, right=right_margin, top=0.91, bottom=0.06)
            if show_in_notebook:
                plt.show()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        # Append the restructured feature pages (grouped + per-sample mean +/-
        # SEM) for percentage_dead_mask, nr_dead_mask_pixels, and mean_dead_dye
        # (whichever are available), in that order. Smoothed variants are
        # hidden for now (raw only).
        _plot_feature_pages(
            pdf=pdf,
            df_meansem_persample=df_meansem,
            df_meansem_grouped=df_meansem_grouped,
            feature_cols=value_cols,
            group_label_title=legend_title,
            color_map=color_map,
            org_type=org_type,
            group_sample_counts=group_sample_counts,
        )


def plot_dead_signal_per_organoid(
    df_combined,
    figsize=(12, 10),
    feature="smoothed_percentage_dead_mask",
    color_by="organoid_type",
    style_by=None,
    autoscale_y=True,
    y_quantile=0.99,
    min_ymax=None,
    title=None,
    ylabel=None,
    show_legends=True,
    linewidth=2.2,
    alpha=0.65,
    legend_max_items=40,
    color_map=None,
    linestyle_map=None,
    highlight_dead_phase=False,
    dead_phase_color="red",
    disappearance_markers=None,
    disappearance_marker_color="black",
    disappearance_marker_size=38,
    screen_show_scale=1.0,
    show_in_notebook=True,
    target_ax=None,
):
    """
    Plot individual organoid dead signal over time.

    - `color_by`: column used for colors (e.g. "organoid_type" or "condition_line")
    - `style_by`: optional column used for line styles (e.g. "organoid_type")
    - `autoscale_y`: when True, scales y-axis to the observed data range (not fixed 0–1),
      using the ``y_quantile`` upper bound so brief spikes do not flatten small signals.
    - `color_map`: optional dict mapping color_key levels to colors.
    - `linestyle_map`: optional dict mapping style_key levels to linestyles.
    - `show_in_notebook`: when False, skip the inline ``plt.show()`` (PDF saving is
      unaffected). Ignored when ``target_ax`` is provided (caller owns display).
    - `target_ax`: when provided, draw into this axis instead of creating a new
      figure. Figure-level adjustments (subplots_adjust, DPI scaling) and
      ``plt.show()`` are then the caller's responsibility.
    """
    if target_ax is None:
        fig, ax = plt.subplots(figsize=figsize)
        owns_figure = True
    else:
        ax = target_ax
        fig = ax.figure
        owns_figure = False

    dfp = df_combined.copy()

    # Allow passing multiple columns for color/style (builds a combined key)
    def _ensure_cols(cols, what):
        if cols is None:
            return []
        if isinstance(cols, (list, tuple)):
            cols_list = list(cols)
        else:
            cols_list = [cols]
        missing = [c for c in cols_list if c not in dfp.columns]
        if missing:
            raise ValueError(f"{what} column(s) missing in df_combined: {missing}")
        return cols_list

    color_cols = _ensure_cols(color_by, "color_by")
    style_cols = _ensure_cols(style_by, "style_by") if style_by is not None else []

    if feature not in dfp.columns:
        raise ValueError(f"feature='{feature}' not found in df_combined columns.")

    # Combined keys
    if len(color_cols) == 1:
        color_key_col = color_cols[0]
    else:
        color_key_col = "__color_key__"
        dfp[color_key_col] = dfp[color_cols].astype(str).agg(" | ".join, axis=1)

    if style_cols:
        if len(style_cols) == 1:
            style_key_col = style_cols[0]
        else:
            style_key_col = "__style_key__"
            dfp[style_key_col] = dfp[style_cols].astype(str).agg(" | ".join, axis=1)
    else:
        style_key_col = None
    
    # Create unique ID for each organoid if not present
    if "organoid_id" not in dfp.columns:
        dfp["organoid_id"] = (
            dfp["sample_name"].astype(str) + "_" +
            dfp.get("organoid_type", pd.Series(["unknown"] * len(dfp))).astype(str) + "_" +
            dfp["TrackID"].astype(str)
        )

    # Build palettes and linestyles if not provided
    if color_map is None:
        color_levels = sorted(dfp[color_key_col].dropna().astype(str).unique())
        n_colors = len(color_levels)
        palette = sns.color_palette("tab10", n_colors) if n_colors <= 10 else sns.color_palette("husl", n_colors)
        color_map = dict(zip(color_levels, palette))
    else:
        color_levels = sorted(color_map.keys())

    if linestyle_map is None:
        if style_key_col is not None:
            style_levels = sorted(dfp[style_key_col].dropna().astype(str).unique())
            # cycle through a small set of visually distinct linestyles
            linestyles = ["-", "--", ":", "-.", (0, (5, 2)), (0, (1, 2))]
            linestyle_map = {lvl: linestyles[i % len(linestyles)] for i, lvl in enumerate(style_levels)}
        else:
            style_levels = []
            linestyle_map = {}
    else:
        style_levels = sorted(linestyle_map.keys())

    # Plot each organoid track as its own line, inheriting group color/style
    group_cols = ["organoid_id", "position_t", feature, color_key_col] + ([style_key_col] if style_key_col else [])
    if highlight_dead_phase and "dead" in dfp.columns:
        group_cols.append("dead")
    df_plot = dfp[group_cols].copy()
    df_plot[color_key_col] = df_plot[color_key_col].astype(str)
    if style_key_col is not None:
        df_plot[style_key_col] = df_plot[style_key_col].astype(str)

    if disappearance_markers is not None and not disappearance_markers.empty:
        dm = disappearance_markers.copy()
        if "organoid_id" in dm.columns and "position_t" in dm.columns:
            dm["organoid_id"] = dm["organoid_id"].astype(str)
            dm["position_t"] = pd.to_numeric(dm["position_t"], errors="coerce")
            dm = dm.dropna(subset=["organoid_id", "position_t"])
            disappearance_t_map = (
                dm.groupby("organoid_id")["position_t"].first().to_dict()
            )
        else:
            disappearance_t_map = {}
    else:
        disappearance_t_map = {}

    for oid, df_oid in df_plot.groupby("organoid_id", sort=False):
        df_oid = df_oid.sort_values("position_t")
        cval = df_oid[color_key_col].iloc[0]
        if isinstance(cval, pd.Series):
            cval = " | ".join(cval.astype(str).tolist())
        else:
            cval = str(cval)

        sval = df_oid[style_key_col].iloc[0] if style_key_col is not None else None
        if isinstance(sval, pd.Series):
            sval = " | ".join(sval.astype(str).tolist())
        elif sval is not None:
            sval = str(sval)
        x = df_oid["position_t"].to_numpy()
        y = df_oid[feature].to_numpy()
        ls = linestyle_map.get(sval, "-") if style_key_col is not None else "-"
        base_color = color_map.get(cval, "C0")

        ax.plot(x, y, color=base_color, linestyle=ls, alpha=alpha, linewidth=linewidth, zorder=2)

        if highlight_dead_phase and "dead" in df_oid.columns:
            dead_bool = df_oid["dead"].fillna(False).astype(bool).to_numpy()
            if dead_bool.any():
                first_dead_idx = int(np.argmax(dead_bool))
                # Mark the first dead event with a "D" glyph instead of recoloring
                # the post-death segment.
                ax.scatter(
                    [x[first_dead_idx]],
                    [y[first_dead_idx]],
                    marker="$D$",
                    s=60,
                    color=dead_phase_color,
                    edgecolors=dead_phase_color,
                    linewidths=0.8,
                    alpha=0.95,
                    zorder=3.2,
                )

        t_dis = disappearance_t_map.get(str(oid))
        if t_dis is not None:
            t_dis = int(t_dis)
            hit = df_oid[df_oid["position_t"] == t_dis]
            if not hit.empty:
                y_dis = float(hit[feature].iloc[0])
                ax.scatter(
                    [t_dis],
                    [y_dis],
                    marker="x",
                    s=disappearance_marker_size,
                    linewidths=1.5,
                    color=disappearance_marker_color,
                    alpha=0.95,
                    zorder=3.5,
                )

    # Legends: colors (color_by) and optional styles (style_by)
    if show_legends:
        # Put legends outside the axes (right side) and reserve space for them.
        # If there are too many entries, show the first N (still better than none).
        shown_color_levels = color_levels[:legend_max_items]
        color_handles = [
            Line2D([0], [0], color=color_map[lvl], lw=3, linestyle="-", label=str(lvl))
            for lvl in shown_color_levels
        ]
        color_title = " / ".join(color_cols) if len(color_cols) > 1 else color_cols[0]
        leg1 = ax.legend(
            handles=color_handles,
            title=color_title,
            loc="upper left",
            bbox_to_anchor=(1.01, 1.00),
            borderaxespad=0,
            prop={"size": 8},
            frameon=True,
        )
        ax.add_artist(leg1)

        if style_key_col is not None and style_levels:
            shown_style_levels = style_levels[:legend_max_items]
            style_handles = [
                Line2D([0], [0], color="black", lw=3, linestyle=linestyle_map[lvl], label=str(lvl))
                for lvl in shown_style_levels
            ]
            style_title = " / ".join(style_cols) if len(style_cols) > 1 else style_cols[0]
            leg_style = ax.legend(
                handles=style_handles,
                title=style_title,
                loc="upper left",
                bbox_to_anchor=(1.01, 0.68),
                borderaxespad=0,
                prop={"size": 8},
                frameon=True,
                handlelength=3.5,
                handletextpad=0.8,
            )
            ax.add_artist(leg_style)

        if highlight_dead_phase and "dead" in dfp.columns:
            dead_handle = [
                Line2D(
                    [0], [0],
                    marker="$D$",
                    color=dead_phase_color,
                    markerfacecolor=dead_phase_color,
                    markeredgecolor=dead_phase_color,
                    markersize=9,
                    markeredgewidth=0.8,
                    lw=0,
                    label="First dead event (dead=True)",
                )
            ]
            leg_dead = ax.legend(
                handles=dead_handle,
                title="Overlay",
                loc="upper left",
                bbox_to_anchor=(1.01, 0.40),
                borderaxespad=0,
                prop={"size": 8},
                frameon=True,
            )
            ax.add_artist(leg_dead)

        if disappearance_t_map:
            dis_handle = [
                Line2D(
                    [0], [0],
                    marker="x",
                    color=disappearance_marker_color,
                    lw=0,
                    markersize=7,
                    markeredgewidth=1.6,
                    label="Track disappears (flat tail kept)",
                )
            ]
            leg_dis = ax.legend(
                handles=dis_handle,
                title="Markers",
                loc="upper left",
                bbox_to_anchor=(1.01, 0.20),
                borderaxespad=0,
                prop={"size": 8},
                frameon=True,
            )
            ax.add_artist(leg_dis)
    
    # Backwards-compatible defaults (keep old title/ylabel unless caller overrides)
    if title is None:
        title = "Dead Signal Increase Per Individual Organoid"
    if ylabel is None:
        if feature == "smoothed_percentage_dead_mask":
            ylabel = "% Dead Mask"
        elif feature == "smoothed_nr_dead_mask_pixels":
            ylabel = "# Dead Mask Pixels"
        else:
            ylabel = feature

    ax.set_title(title, fontsize=14)
    ax.set_xlabel("Timepoint", fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.grid(True, linestyle=":", alpha=0.5)
    
    # X-axis: give a small breathing room at the start and end of the data range.
    t_vals = pd.to_numeric(dfp["position_t"], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    if len(t_vals) > 0:
        t_min = float(t_vals.min())
        t_max = float(t_vals.max())
        x_pad = 0.03 * (t_max - t_min) if t_max > t_min else 1.0
        ax.set_xlim(t_min - x_pad, t_max + x_pad)
    else:
        ax.set_xlim(left=0)

    # Y scaling: small pad below 0 and above the top so traces don't touch the axes.
    if autoscale_y:
        y = pd.to_numeric(dfp[feature], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
        if len(y) > 0:
            ymax = float(y.max())
            if min_ymax is not None:
                ymax = max(ymax, float(min_ymax))
            if ymax <= 0:
                ymax = 1.0 if "percentage" in feature else 10.0
            y_pad = 0.05 * ymax
            ax.set_ylim(-y_pad, ymax + y_pad)
        else:
            ax.set_ylim(-0.05, 1.05 if "percentage" in feature else 1.0)
    else:
        if "percentage" in feature:
            ax.set_ylim(-0.05, 1.05)
        else:
            ax.set_ylim(bottom=-0.05)

    # Figure-level adjustments and inline display are skipped when drawing into
    # a caller-provided axis: the caller owns the figure and decides when (and
    # whether) to show it.
    if owns_figure:
        # Keep right margin for the stacked legends; tight_layout can clip outside legends.
        # Reserve a dedicated right gutter for legends (screen + PDF), no overlap with data.
        fig.subplots_adjust(right=0.64)
        # Screen-only downscaling for notebook display; keep original size for PDF saving.
        # Use DPI scaling so labels/legends shrink together with the plot area.
        orig_dpi = fig.get_dpi()
        if screen_show_scale != 1.0:
            fig.set_dpi(orig_dpi * screen_show_scale)
        if show_in_notebook:
            plt.show()
        if screen_show_scale != 1.0:
            fig.set_dpi(orig_dpi)
    return fig


def plot_grouped_dead_signal(
    df_processed,
    feature,
    figsize,
    color_key,
    style_key,
    condition_color_map,
    type_style_map,
    threshold_annotation,
    counts_annotation,
    screen_show_scale=1.0,
    band_alpha=0.18,
    show_in_notebook=True,
):
    """Plot grouped mean +/- SEM dead-signal dynamics (used for Plot 5a/5b).

    Tracks are pooled across samples within each (condition, organoid_type)
    group, so n at each timepoint is the number of contributing tracks and
    SEM = std / sqrt(n). This matches the uncertainty convention used by the
    interaction analysis, so bands are comparable across BEHAV3D plots.
    """
    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    nice_feat = "% Dead Mask" if "percentage" in feature else "# Dead Pixels"

    # Pool all tracks across samples within each (condition [x organoid_type], timepoint)
    # group. n is therefore the number of tracks contributing at that timepoint,
    # which is large and roughly constant over time (flat-tail extension keeps
    # disappearing tracks in the cohort), giving a stable, visually uniform band.
    # If color_key == style_key (e.g., fallback to organoid_type), avoid duplicate labels.
    group_keys = [color_key] if color_key == style_key else [color_key, style_key]
    group_cols = group_keys + ["position_t"]

    agg = (
        df_processed
        .groupby(group_cols)[feature]
        .agg(mean_val="mean", std_val="std", n_val="count")
        .reset_index()
    )

    agg["sem_val"] = agg["std_val"] / np.sqrt(agg["n_val"].clip(lower=1))

    for group_val, grp in agg.groupby(group_keys):
        grp = grp.sort_values("position_t")
        t = grp["position_t"].values
        mean = grp["mean_val"].values
        band_vals = grp["sem_val"].to_numpy()

        if color_key == style_key:
            cond_val = str(group_val)
            org_val = str(group_val)
        else:
            cond_val, org_val = group_val
            cond_val = str(cond_val)
            org_val = str(org_val)

        color = condition_color_map.get(cond_val, "C0")
        ls = type_style_map.get(str(org_val), "-")

        label = f"{cond_val} | {org_val}" if cond_val != org_val else cond_val
        ax.plot(t, mean, linewidth=2.8, label=label, color=color, linestyle=ls, zorder=3)
        ax.fill_between(
            t,
            mean - band_vals,
            mean + band_vals,
            alpha=band_alpha,
            color=color,
            linewidth=0,
            zorder=1,
        )

    ax.axhline(0, color="black", linewidth=1, alpha=0.35, zorder=0)

    # X-axis: give a small breathing room at the start and end of the data range.
    t_all = pd.to_numeric(agg["position_t"], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    if len(t_all) > 0:
        t_min = float(t_all.min())
        t_max = float(t_all.max())
        x_pad = 0.03 * (t_max - t_min) if t_max > t_min else 1.0
        ax.set_xlim(t_min - x_pad, t_max + x_pad)
    else:
        ax.set_xlim(left=0)
    ax.set_xlabel("Timepoint")
    ax.set_ylabel(f"Change from baseline ({nice_feat})")
    ax.set_title(
        f"Mean +/- SEM per Condition & Organoid Type (tracks pooled across samples)\n"
        f"({nice_feat}, smoothed, baseline-normalized, non-decreasing, disappearing tracks kept as flat tails)"
    )
    ax.grid(True, linestyle="--", linewidth=0.7, alpha=0.28)

    ax.text(
        0.02, 0.97,
        threshold_annotation,
        transform=ax.transAxes,
        fontsize=8,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lavender", alpha=0.8),
    )

    ax.text(
        0.02, 0.89,
        counts_annotation,
        transform=ax.transAxes,
        fontsize=8,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="honeydew", alpha=0.85),
    )

    # Y-axis: small pad below 0 and above the upper edge of the SEM band so the
    # traces don't touch the axes on either side.
    band_high = []
    for _, grp in agg.groupby(group_keys):
        mean = grp["mean_val"].to_numpy()
        band_vals = grp["sem_val"].to_numpy()
        high = mean + band_vals
        if np.isfinite(high).any():
            band_high.append(np.nanmax(high))
    if band_high:
        y_max = max(band_high)
        y_pad = 0.05 * (y_max if y_max > 0 else 1.0)
        ax.set_ylim(-y_pad, y_max + y_pad)
    else:
        ax.set_ylim(-0.05, 1.0)

    # Dual legends matching plots 3/4.
    color_handles = [
        Line2D([0], [0], color=condition_color_map[lvl], lw=3, linestyle="-", label=str(lvl))
        for lvl in sorted(condition_color_map.keys())
    ]
    leg1 = ax.legend(
        handles=color_handles,
        title=color_key.replace("_", " ").title(),
        loc="upper left",
        bbox_to_anchor=(1.01, 1.00),
        borderaxespad=0,
        prop={"size": 8},
        frameon=True,
    )
    ax.add_artist(leg1)

    style_handles = [
        Line2D([0], [0], color="black", lw=3, linestyle=type_style_map[lvl], label=str(lvl))
        for lvl in sorted(type_style_map.keys())
    ]
    leg_style = ax.legend(
        handles=style_handles,
        title=style_key.replace("_", " ").title(),
        loc="upper left",
        bbox_to_anchor=(1.01, 0.68),
        borderaxespad=0,
        prop={"size": 8},
        frameon=True,
        handlelength=3.5,
        handletextpad=0.8,
    )
    ax.add_artist(leg_style)

    fig.subplots_adjust(right=0.64)

    # Screen-only downscaling for notebook display; keep original size for PDF saving.
    # Use DPI scaling so labels/legends shrink together with the plot area.
    fig_orig_dpi = fig.get_dpi()
    if screen_show_scale != 1.0:
        fig.set_dpi(fig_orig_dpi * screen_show_scale)
    if show_in_notebook:
        plt.show()
    if screen_show_scale != 1.0:
        fig.set_dpi(fig_orig_dpi)

    return fig


def plot_multi_organoid_death_dynamics(
    output_dir,
    organoid_types,
    figsize=(12, 10),
    dead_perc_threshold_map=None,
    show_in_notebook=True,
    ):
    """
    Generate comparison plots for death dynamics across multiple organoid types.
    
    Parameters
    ----------
    output_dir : str or Path
        Base output directory containing analysis results
    organoid_types : list
        List of organoid type names (e.g., ['organoid1', 'organoid2'])
    figsize : tuple
        Figure size for plots
    
    Returns
    -------
    Path or None
        Path to generated PDF, or None if insufficient data
    """
    print(f"--------------- Multi-Organoid Death Dynamics Comparison ---------------")
    
    output_dir = Path(output_dir)

    # Threshold text for Plot 4a/4b annotations (from pipeline settings).
    # ``dead_perc_threshold_map`` values are fractions (0.0-1.0), matching
    # the ``dead_mask_percentage_threshold`` config scale; convert to
    # percent (0-100) for display here.
    threshold_annotation = "Dead threshold (% mask): not provided"
    if isinstance(dead_perc_threshold_map, dict) and dead_perc_threshold_map:
        cleaned_thr = {}
        for k, v in dead_perc_threshold_map.items():
            try:
                cleaned_thr[str(k)] = float(v) * 100.0
            except Exception:
                continue
        if cleaned_thr:
            uniq_thr = sorted(set(cleaned_thr.values()))
            if len(uniq_thr) == 1:
                threshold_annotation = f"Dead threshold (% mask): {uniq_thr[0]:.4g}"
            else:
                parts = [f"{k}={cleaned_thr[k]:.4g}" for k in sorted(cleaned_thr.keys())]
                threshold_annotation = "Dead threshold (% mask): " + "; ".join(parts)
    
    # Check which organoid types have death dynamics data available. Use each
    # target's *filtered* track features as the source, so this comparison
    # always reflects the current filtering — rather than a cached
    # combined_general_*.csv that only refreshes when per-target Death Dynamics
    # is re-run.
    dead_columns = ["nr_dead_mask_pixels", "percentage_dead_mask", "dead"]
    available_data = {}
    for org_type in organoid_types:
        csv_path = Path(output_dir, "analysis", org_type, "track_features", f"BEHAV3D_{org_type}_combined_track_features_filtered.csv")
        if csv_path.exists():
            available_data[org_type] = csv_path
            print(f"Found data for {org_type}")
        else:
            print(f"Missing data for {org_type}: {csv_path}")

    if len(available_data) < 1:
        print(f"\nNeed at least 1 organoid type with filtered track features.")
        print(f"   Available: {list(available_data.keys())}")
        print(f"   Run track filtering for each organoid type first.")
        return None

    # Recompute per-sample cumulative dynamics from the filtered tracks
    # (same helper the per-target step uses) and combine across targets.
    all_data = []
    for org_type, csv_path in available_data.items():
        df_src = pd.read_csv(csv_path)
        if not all(col in df_src.columns for col in dead_columns):
            print(f"Skipping {org_type}: filtered track features missing dead columns "
                  f"({[c for c in dead_columns if c not in df_src.columns]}).")
            continue
        df, _, _ = _compute_general_death_dynamics(df_src)
        df["organoid_type"] = org_type
        all_data.append(df)

    if not all_data:
        print(f"\nNo organoid types with usable dead data for the combined comparison.")
        return None

    df_combined = pd.concat(all_data, ignore_index=True)
    
    # Combined label for each sample + organoid_type combination
    df_combined["sample_organoid"] = df_combined["sample_name"] + "_" + df_combined["organoid_type"]
    
    # Output directory
    results_outdir = Path(output_dir, "analysis", "multi_organoid_comparison")
    results_outdir.mkdir(parents=True, exist_ok=True)
    
    # Save combined data
    combined_csv_path = results_outdir / "combined_multi_organoid_death_dynamics.csv"
    df_combined.to_csv(combined_csv_path, index=False)
    print(f"\nCombined data saved to: {combined_csv_path}")

    # Load and combine full track features for individual track plotting
    all_tracks_data = []
    skipped_types = []

    # Load metadata to fetch per-sample line condition info.
    # Prefer "metadata.csv", but also support "metadata1.csv", "metadata_*.csv", etc.
    metadata_df = None
    metadata_path = Path(output_dir, "metadata.csv")
    if not metadata_path.exists():
        candidates = sorted(Path(output_dir).glob("metadata*.csv"))
        metadata_path = candidates[0] if candidates else metadata_path

    if metadata_path.exists():
        try:
            metadata_df = pd.read_csv(metadata_path, low_memory=False)
            print(f"Using metadata file: {metadata_path}")
        except Exception as exc:
            print(f"\n⚠️  WARNING: Could not read metadata at {metadata_path}: {exc}")
            metadata_df = None
    else:
        print(
            f"\n⚠️  WARNING: No metadata*.csv found in {output_dir}. "
            f"Line-condition coloring will fall back to organoid_type."
        )
    
    for org_type in available_data.keys():
        track_csv_path = Path(output_dir, "analysis", org_type, "track_features", f"BEHAV3D_{org_type}_combined_track_features_filtered.csv")
        if track_csv_path.exists():
            df_t = pd.read_csv(track_csv_path)
            
            # Check if dead data columns exist
            dead_columns = ["nr_dead_mask_pixels", "percentage_dead_mask", "dead"]
            if not all(col in df_t.columns for col in dead_columns):
                skipped_types.append(org_type)
                continue
                
            df_t["organoid_type"] = org_type

            # Attach line condition from metadata (per sample) if available.
            # Metadata convention: columns like "or_{org_name}_line_condition".
            if metadata_df is not None and "sample_name" in metadata_df.columns:
                preferred = f"or_{org_type}_line_condition"
                if preferred in metadata_df.columns:
                    md_col = preferred
                else:
                    # fallback: any column that ends with "_line_condition" and contains org_type
                    candidates = [
                        c for c in metadata_df.columns
                        if c.endswith("_line_condition") and org_type.lower() in c.lower()
                    ]
                    md_col = candidates[0] if candidates else None

                if md_col is not None:
                    md_map = metadata_df.set_index("sample_name")[md_col].to_dict()
                    df_t["line_condition"] = df_t["sample_name"].map(md_map)
                else:
                    df_t["line_condition"] = np.nan
            else:
                df_t["line_condition"] = np.nan
            
            # Ensure smoothed_percentage_dead_mask is present and starts at 0 (min_periods=1)
            # Default smoothing helps reduce noise from flickering segmentation
            df_t["smoothed_percentage_dead_mask"] = smooth_value_over_time(
                df_t, column="percentage_dead_mask", groupby=["TrackID", "sample_name"], min_periods=1
            )

            # Smoothed dead pixel count (often more interpretable than percentages)
            df_t["smoothed_nr_dead_mask_pixels"] = smooth_value_over_time(
                df_t, column="nr_dead_mask_pixels", groupby=["TrackID", "sample_name"], min_periods=1
            )
            
            all_tracks_data.append(df_t)
            
    if skipped_types:
        print(f"\n⚠️  WARNING: Dead data not found for: {', '.join(skipped_types)}")
        print(f"   Run feature extraction and track filtering for these types first to enable individual dead signal plotting.\n")

    df_tracks_combined = pd.concat(all_tracks_data, ignore_index=True) if all_tracks_data else None

    # ------------------------------------------------------------------
    # Per-sample mean cumulative dead fraction per (line_condition x organoid_type).
    # `percentage_dead` in df_combined is already the per-sample cumulative dead
    # fraction (nr_dead / nr_organoids_t0, never decreasing). Meaning across
    # samples within each (line_condition, organoid_type) group keeps the curve
    # in [0, 1] regardless of group size, so groups are directly comparable.
    # Used by the small subpanel under Plots 4a/4b.
    # ------------------------------------------------------------------
    df_combined["line_condition"] = np.nan
    if metadata_df is not None and "sample_name" in metadata_df.columns:
        for org_type in df_combined["organoid_type"].unique():
            preferred = f"or_{org_type}_line_condition"
            if preferred in metadata_df.columns:
                md_col = preferred
            else:
                candidates_cols = [
                    c for c in metadata_df.columns
                    if c.endswith("_line_condition") and str(org_type).lower() in c.lower()
                ]
                md_col = candidates_cols[0] if candidates_cols else None
            if md_col is None:
                continue
            md_map = metadata_df.set_index("sample_name")[md_col].to_dict()
            mask = df_combined["organoid_type"] == org_type
            df_combined.loc[mask, "line_condition"] = (
                df_combined.loc[mask, "sample_name"].map(md_map)
            )

    df_cumdead = None
    if df_combined["line_condition"].notna().any():
        df_cumdead = (
            df_combined.dropna(subset=["line_condition"])
            .groupby(["line_condition", "organoid_type", "position_t"])["percentage_dead"]
            .mean()
            .reset_index(name="mean_cum_dead_fraction")
        )

    # Plots
    pdf_path = results_outdir / "multi_organoid_death_dynamics_comparison.pdf"
    
    # Unique combinations for color palette
    unique_combinations = df_combined["sample_organoid"].unique()
    n_combinations = len(unique_combinations)
    
    # Enough distinct colors
    if n_combinations <= 10:
        colors = sns.color_palette("tab10", n_combinations)
    elif n_combinations <= 20:
        colors = sns.color_palette("tab20", n_combinations)
    else:
        # Many combinations: continuous colormap
        colors = sns.color_palette("husl", n_combinations)
    
    color_map = dict(zip(unique_combinations, colors))
    
    screen_show_scale = 0.72  # Display only: smaller notebook rendering for plots 3/4/5.

    with PdfPages(pdf_path) as pdf:
        # Plot 1: Line plot - Percentage Dead per sample + organoid_type
        fig1, ax1 = plt.subplots(1, 1, figsize=figsize)
        
        sns.lineplot(
            data=df_combined,
            x='position_t',
            y='percentage_dead',
            hue='sample_organoid',
            ax=ax1,
            linewidth=2,
            palette=color_map
        )
        
        ax1.grid(True, linestyle=':', linewidth=1, alpha=0.5)
        ax1.set_ylim(0, 1.05)
        ax1.set_xlim(0)
        ax1.set_xlabel('Timepoint')
        ax1.set_ylabel('')
        ax1.set_title('Percentage of Dead Multi-Organoids')
        
        # Format y-axis as percentage
       #ax1.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x*100:.0f}%'))
        
        ax1.legend(
            bbox_to_anchor=(1.02, 1),
            loc='upper left',
            borderaxespad=0,
            prop={'size': 7},
            title='Sample Name'
        )
        
        plt.tight_layout()
        if show_in_notebook:
            plt.show()
        pdf.savefig(fig1, bbox_inches='tight')
        plt.close(fig1)
        
        # Plot 2: Bar plot - End state per sample + organoid_type
        fig2, ax2 = plt.subplots(1, 1, figsize=figsize)
        
        # Final timepoint for each sample and organoid type
        idx = df_combined.groupby(["sample_organoid"])["position_t"].idxmax()
        df_end = df_combined.loc[idx, ["organoid_type", "sample_name", "sample_organoid", 
                                        "percentage_dead", "percentage_alive", "percentage_disappeared"]].copy()
        
        # Sort by sample_name then organoid_type: consistent ordering
        df_end = df_end.sort_values(["sample_name", "organoid_type"]).reset_index(drop=True)
        
        # x-axis labels: "organoid_type_sample_name"
        df_end["bar_label"] = df_end["organoid_type"] + "_" + df_end["sample_name"]
        
        x = np.arange(len(df_end))
        width = 0.6
        
        # Stacked bar plot: Alive (bottom), Dead (middle), Disappeared (top)
        bars_alive = ax2.bar(x, df_end["percentage_alive"], width, 
                             label="percentage_alive", color="#6699CC", zorder=2)
        bars_dead = ax2.bar(x, df_end["percentage_dead"], width, 
                            bottom=df_end["percentage_alive"], 
                            label="percentage_dead", color="#CC6666", zorder=2)
        bars_disappeared = ax2.bar(x, df_end["percentage_disappeared"], width, 
                                   bottom=df_end["percentage_alive"] + df_end["percentage_dead"], 
                                   label="percentage_disappeared", color="#898989", zorder=2)
        
        ax2.set_xticks(x)
        ax2.set_xticklabels(df_end["bar_label"], rotation=45, ha='right', fontsize=8)
        ax2.grid(True, linestyle=':', linewidth=1, alpha=0.5, zorder=1)
        ax2.set_ylim(0, 1.05)
        ax2.set_xlabel('Organoid Name and Sample Name')
        ax2.set_ylabel('')
        ax2.set_title('Percentage of Multiple Organoids States at End of Experiment')
        
        # Format y-axis as percentage
        #x2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x*100:.0f}%'))
        
        # Legend with correct order: disappeared, alive, dead (matching single organoid plot)
        handles, labels = ax2.get_legend_handles_labels()
        order = [2, 0, 1]  # disappeared, alive, dead
        ax2.legend(
            [handles[i] for i in order],
            [labels[i] for i in order],
            bbox_to_anchor=(1.02, 1), 
            loc='upper left', 
            borderaxespad=0, 
            prop={'size': 8}
        )
        
        plt.tight_layout()
        if show_in_notebook:
            plt.show()
        pdf.savefig(fig2, bbox_inches='tight')
        plt.close(fig2)

        # Plot 3: Individual track signal over time
        if df_tracks_combined is not None:
            # Line style per organoid_type; color by line condition when metadata has values.
            df_indiv = df_tracks_combined.copy()
            has_line_cond = (
                "line_condition" in df_indiv.columns
                and df_indiv["line_condition"].notna().any()
            )
            if has_line_cond:
                df_indiv["line_condition"] = df_indiv["line_condition"].fillna("(no line condition)")
                color_key = "line_condition"
            else:
                color_key = "organoid_type"
            style_key = "organoid_type"

            # Pre-calculate consistent colors and linestyles using ALL tracks (before Plot 4/5 filtering)
            color_levels_global = sorted(df_indiv[color_key].dropna().astype(str).unique())
            n_colors_global = len(color_levels_global)
            palette_global = sns.color_palette("tab10", n_colors_global) if n_colors_global <= 10 else sns.color_palette("husl", n_colors_global)
            condition_color_map = dict(zip(color_levels_global, palette_global))

            style_levels_global = sorted(df_indiv[style_key].dropna().astype(str).unique())
            linestyles_cycle = ["-", "--", ":", "-.", (0, (5, 2)), (0, (1, 2))]
            type_style_map = {lvl: linestyles_cycle[i % len(linestyles_cycle)] for i, lvl in enumerate(style_levels_global)}

            # 3a) Smoothed percentage dead mask (auto-scaled so small values are visible)
            fig3a = plot_dead_signal_per_organoid(
                df_indiv,
                figsize=figsize,
                feature="smoothed_percentage_dead_mask",
                color_by=color_key,
                style_by=style_key,
                autoscale_y=True,
                title="Dead Signal Per Individual Organoid — smoothed (% Dead Mask)",
                ylabel=None,
                linewidth=2.6,
                color_map=condition_color_map,
                linestyle_map=type_style_map,
                screen_show_scale=screen_show_scale,
                show_in_notebook=False,
            )
            pdf.savefig(fig3a, bbox_inches="tight")
            plt.close(fig3a)

            # 3b) Smoothed number of dead pixels (auto-scaled)
            fig3b = plot_dead_signal_per_organoid(
                df_indiv,
                figsize=figsize,
                feature="smoothed_nr_dead_mask_pixels",
                color_by=color_key,
                style_by=style_key,
                autoscale_y=True,
                title="Dead Signal Per Individual Organoid — smoothed (# Dead Pixels)",
                ylabel=None,
                linewidth=2.6,
                color_map=condition_color_map,
                linestyle_map=type_style_map,
                screen_show_scale=screen_show_scale,
                show_in_notebook=False,
            )
            pdf.savefig(fig3b, bbox_inches="tight")
            plt.close(fig3b)

            # --- Plots 4a & 4b: baseline-normalized tracks with flat extension after disappearance ---

            # Identify tracks that disappear: a track "disappears" if its last recorded
            # timepoint is earlier than the last timepoint of its sample.
            track_groups = ["TrackID", "sample_name", "organoid_type"]

            # True sample endpoints must be computed on all tracks (before any filtering).
            sample_end_by_sample = (
                df_indiv.groupby("sample_name")["position_t"].max().to_dict()
            )

            n_total = df_indiv.groupby(track_groups).ngroups

            # Disappearance information for marker overlays (both dead and non-dead disappearing tracks).
            df_track_endpoints = (
                df_indiv.groupby(track_groups, as_index=False)["position_t"]
                .max()
                .rename(columns={"position_t": "disappearance_t"})
            )
            df_track_endpoints["sample_end_t"] = df_track_endpoints["sample_name"].map(sample_end_by_sample)
            df_track_endpoints["disappeared"] = df_track_endpoints["disappearance_t"] < df_track_endpoints["sample_end_t"]
            df_disappearance_markers = df_track_endpoints[df_track_endpoints["disappeared"]].copy()
            if not df_disappearance_markers.empty:
                df_disappearance_markers["organoid_id"] = (
                    df_disappearance_markers["sample_name"].astype(str)
                    + "_"
                    + df_disappearance_markers["organoid_type"].astype(str)
                    + "_"
                    + df_disappearance_markers["TrackID"].astype(str)
                )
                df_disappearance_markers = df_disappearance_markers[["organoid_id", "disappearance_t"]].rename(
                    columns={"disappearance_t": "position_t"}
                )

            def extend_disappearing_tracks_to_sample_end(df, group_cols, sample_end_map, time_column="position_t"):
                """
                Extend ALL disappearing tracks by repeating their last value until sample end.
                Keeps disappearing tracks in the cohort and avoids silently dropping them.
                """
                df = df.copy()
                extended_rows = []

                for _, track_group in df.groupby(group_cols):
                    track_max_t = track_group[time_column].max()
                    sample_name_val = str(track_group["sample_name"].iloc[0])
                    sample_max_t_val = sample_end_map.get(sample_name_val, track_max_t)

                    if track_max_t < sample_max_t_val:
                        last_row = track_group.sort_values(time_column).iloc[-1].copy()
                        missing_times = np.arange(track_max_t + 1, sample_max_t_val + 1)

                        for t in missing_times:
                            extended_row = last_row.copy()
                            extended_row[time_column] = t
                            extended_rows.append(extended_row)
                
                if extended_rows:
                    n_original_rows = len(df)
                    n_extended_rows = len(extended_rows)
                    df_extended = pd.concat([df, pd.DataFrame(extended_rows)], ignore_index=True)

                    print(f"  [Track extension] Extended {n_extended_rows} timepoint(s) for disappearing tracks")
                    print(f"    Original rows: {n_original_rows}, Extended rows: {n_extended_rows}, Total: {len(df_extended)}")
                    
                    return df_extended
                else:
                    print(f"  [Track extension] No disappearing tracks to extend")
                    return df

            print(f"\n[Track extension] Processing all disappearing tracks...")
            df_no_disapp = extend_disappearing_tracks_to_sample_end(
                df_indiv,
                track_groups,
                sample_end_by_sample,
                time_column="position_t",
            )
            n_disapp_total = int(df_disappearance_markers.shape[0])
            print(f"[Plot 4/5] Tracks with disappearance markers: {n_disapp_total}/{n_total}")

            if "dead" in df_no_disapp.columns:
                n_dead_tracks = int(
                    df_no_disapp.groupby(track_groups)["dead"].any().sum()
                )
            else:
                n_dead_tracks = 0

            counts_annotation = (
                f"Tracks: {n_total} | Dead: {n_dead_tracks} | Disappeared: {n_disapp_total}\n"
                f"(dead + non-dead; flat tail kept)"
            )

            def normalize_to_zero(df, feature, group_cols):
                """Subtract the first value of *feature* within each track group."""
                df = df.copy()
                first_val = df.sort_values("position_t").groupby(group_cols)[feature].transform("first")
                df[feature] = df[feature] - first_val
                return df

            def enforce_monotone_non_decreasing(df, features, group_cols, time_col="position_t"):
                """Force per-track non-decreasing signals (cummax over time)."""
                df = df.copy()
                df = df.sort_values(time_col)
                for feature in features:
                    if feature in df.columns:
                        df[feature] = df.groupby(group_cols)[feature].transform(lambda s: s.cummax())
                return df

            # Processing pipeline for plots 4 & 5:
            # 1. Normalize to zero (change from baseline)
            # 2. Enforce monotone non-decreasing (cummax on the change-from-baseline)
            # This ensures we see "cumulative increase from baseline" — lines never go down.
            processed_feats = [
                "smoothed_percentage_dead_mask",
                "smoothed_nr_dead_mask_pixels",
            ]

            df_processed = df_no_disapp.copy()
            for feat in processed_feats:
                df_processed = normalize_to_zero(df_processed, feat, track_groups)
            df_processed = enforce_monotone_non_decreasing(
                df_processed,
                features=processed_feats,
                group_cols=track_groups,
                time_col="position_t",
            )

            # Absolute-value counterpart of plots 4a/4b: same extension + cummax
            # pipeline but WITHOUT baseline normalization, so the raw dead-signal
            # scale of each track is visible. Drawn right before plots 4a/4b.
            df_absolute = enforce_monotone_non_decreasing(
                df_no_disapp.copy(),
                features=processed_feats,
                group_cols=track_groups,
                time_col="position_t",
            )

            for feat, label in [
                ("smoothed_percentage_dead_mask", "4a_abs"),
                ("smoothed_nr_dead_mask_pixels",  "4b_abs"),
            ]:
                nice_feat = "% Dead Mask" if "percentage" in feat else "# Dead Pixels"

                fig_abs = plot_dead_signal_per_organoid(
                    df_absolute,
                    figsize=figsize,
                    feature=feat,
                    color_by=color_key,
                    style_by=style_key,
                    autoscale_y=True,
                    title=(
                        f"Dead Signal Absolute ({nice_feat})\n"
                        f"(smoothed, non-decreasing, disappearing tracks kept as flat tails)"
                    ),
                    ylabel=f"Absolute ({nice_feat})",
                    linewidth=2.6,
                    color_map=condition_color_map,
                    linestyle_map=type_style_map,
                    highlight_dead_phase=True,
                    dead_phase_color="red",
                    disappearance_markers=df_disappearance_markers,
                    screen_show_scale=screen_show_scale,
                    show_in_notebook=False,
                )

                fig_abs.axes[0].text(
                    0.02, 0.98,
                    threshold_annotation,
                    transform=fig_abs.axes[0].transAxes,
                    fontsize=8,
                    verticalalignment="top",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="lavender", alpha=0.9),
                )

                fig_abs.axes[0].text(
                    0.02, 0.94,
                    counts_annotation,
                    transform=fig_abs.axes[0].transAxes,
                    fontsize=8,
                    verticalalignment="top",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="honeydew", alpha=0.9),
                )

                pdf.savefig(fig_abs, bbox_inches="tight")
                plt.close(fig_abs)

            # --- Plots 4a + 4b combined: single multipanel figure ---
            # Top row: 4a (% Dead Mask) on the left, 4b (# Dead Pixels) on the right.
            # Bottom row (full width): single cumulative-dead subpanel shared by both
            # metrics (same data, no point duplicating).
            # A single combined legend lives in the right gutter, positioned in
            # figure coordinates so it stays inside the canvas (in-axes
            # bbox_to_anchor sometimes gets clipped by the inline renderer).
            has_subpanel = df_cumdead is not None and not df_cumdead.empty

            # Sanity-check print: report per-(line_condition, organoid_type)
            # dead-event counts so the user can verify the dead-event computation.
            print(
                "\n[Plot 4 dead-event check] dead-event counts per "
                "(line_condition, organoid_type):"
            )
            if (
                "dead" in df_indiv.columns
                and "line_condition" in df_indiv.columns
                and df_indiv["line_condition"].notna().any()
            ):
                _check = df_indiv.copy()
                _check["dead"] = _check["dead"].fillna(False).astype(bool)
                _track_summary = (
                    _check.groupby(["line_condition", "organoid_type", "sample_name", "TrackID"])[
                        "dead"
                    ]
                    .any()
                    .reset_index(name="track_ever_dead")
                )
                _grouping_counts = (
                    _track_summary.groupby(["line_condition", "organoid_type"])
                    .agg(
                        n_tracks=("TrackID", "count"),
                        n_dead=("track_ever_dead", "sum"),
                        n_samples=("sample_name", "nunique"),
                    )
                    .reset_index()
                    .sort_values(["line_condition", "organoid_type"])
                )
                for _, _row in _grouping_counts.iterrows():
                    print(
                        f"  {_row['line_condition']} | {_row['organoid_type']} : "
                        f"{int(_row['n_dead'])} dead / {int(_row['n_tracks'])} tracks "
                        f"({int(_row['n_samples'])} sample(s))"
                    )
            else:
                print(
                    "  (no `dead` column or no metadata-mapped line_condition; "
                    "skipping check)"
                )

            # Wide figure: 2 main panels + a dedicated right-side legend gutter
            # that comfortably fits all stacked legends (including the per-group
            # subpanel legend so every condition is listed even when its curve
            # overlaps another at y=0).
            combo_figsize = (figsize[0] * 1.95, figsize[1] * (1.18 if has_subpanel else 1.0))
            fig_combo = plt.figure(figsize=combo_figsize)
            # The axes occupy ~70% of the width; the rest is reserved for the
            # legends in figure coordinates.
            legend_gutter_left = 0.73
            gs_kwargs = dict(
                hspace=0.14,
                wspace=0.18,
                left=0.05,
                right=legend_gutter_left - 0.02,
                top=0.92,
                bottom=0.09,
            )
            if has_subpanel:
                gs = fig_combo.add_gridspec(
                    2, 2,
                    height_ratios=[4, 1],
                    **gs_kwargs,
                )
                ax_4a = fig_combo.add_subplot(gs[0, 0])
                ax_4b = fig_combo.add_subplot(gs[0, 1])
                ax_sub = fig_combo.add_subplot(gs[1, :])
            else:
                gs = fig_combo.add_gridspec(1, 2, **gs_kwargs)
                ax_4a = fig_combo.add_subplot(gs[0, 0])
                ax_4b = fig_combo.add_subplot(gs[0, 1])
                ax_sub = None

            # Common context lives in the figure-level suptitle; each panel
            # only carries the short metric name as its title.
            fig_combo.suptitle(
                "Dead Signal"
                "(smoothed, baseline-normalized, non-decreasing, flat tails)",
                fontsize=14,
                y=0.97,
            )

            panel_specs = [
                (ax_4a, "smoothed_percentage_dead_mask", "% Dead Mask"),
                (ax_4b, "smoothed_nr_dead_mask_pixels",  "# Dead Pixels"),
            ]
            for ax_panel, feat, nice_feat in panel_specs:
                # show_legends=False on both: a single figure-level legend is
                # added once at the end so legend entries are not duplicated.
                # ylabel drops the "Change from baseline" prefix because the
                # figure suptitle already states the data is baseline-normalized
                # (avoids saying "from baseline" twice).
                plot_dead_signal_per_organoid(
                    df_processed,
                    target_ax=ax_panel,
                    feature=feat,
                    color_by=color_key,
                    style_by=style_key,
                    autoscale_y=True,
                    title=nice_feat,
                    ylabel=f"{nice_feat}",
                    linewidth=2.6,
                    color_map=condition_color_map,
                    linestyle_map=type_style_map,
                    highlight_dead_phase=True,
                    dead_phase_color="red",
                    disappearance_markers=df_disappearance_markers,
                    show_legends=False,
                )

                # Baseline reference line (0 = no change from baseline)
                ax_panel.axhline(0, color="black", linewidth=1, alpha=0.35, zorder=0)

            # ----- Subpanel: cumulative dead fraction per group -----
            if ax_sub is not None:
                # Order groups so the largest values are drawn LAST (on top).
                # This way smaller / flatter curves are not hidden underneath
                # larger ones when several conditions overlap at y=0.
                group_max = (
                    df_cumdead.groupby(["line_condition", "organoid_type"])[
                        "mean_cum_dead_fraction"
                    ]
                    .max()
                    .sort_values()
                )
                draw_order = list(group_max.index)

                print(
                    f"[Plot 4 subpanel] {len(draw_order)} (line_condition x organoid_type) "
                    f"group(s) in df_cumdead: "
                    f"{', '.join(f'{c}|{o}' for c, o in draw_order)}"
                )

                # Pre-compute y-axis scale so we can size the visual offset for
                # flat-at-zero curves relative to the overall y range.
                cur_max = float(
                    pd.to_numeric(df_cumdead["mean_cum_dead_fraction"], errors="coerce")
                    .replace([np.inf, -np.inf], np.nan)
                    .dropna()
                    .max()
                )
                if not np.isfinite(cur_max) or cur_max <= 0:
                    cur_max = 1.0

                # "Flat at zero" means the entire curve sits below 1% of the
                # overall y-axis max. Each such curve gets a small unique
                # upward offset so multiple stacked flat lines are visible.
                flat_zero_threshold = max(cur_max * 0.01, 1e-4)
                flat_offset_step = max(cur_max * 0.015, 1e-4)

                flat_groups = []
                for (cond_val, org_val) in draw_order:
                    grp_max = float(group_max.loc[(cond_val, org_val)])
                    if grp_max < flat_zero_threshold:
                        flat_groups.append((cond_val, org_val))

                flat_offsets = {
                    gk: (i + 1) * flat_offset_step
                    for i, gk in enumerate(flat_groups)
                }

                for (cond_val, org_val) in draw_order:
                    grp = (
                        df_cumdead[
                            (df_cumdead["line_condition"] == cond_val)
                            & (df_cumdead["organoid_type"] == org_val)
                        ]
                        .sort_values("position_t")
                    )
                    y_values = grp["mean_cum_dead_fraction"].to_numpy()
                    if (cond_val, org_val) in flat_offsets:
                        y_values = y_values + flat_offsets[(cond_val, org_val)]
                    ax_sub.plot(
                        grp["position_t"].to_numpy(),
                        y_values,
                        color=condition_color_map.get(str(cond_val), "C0"),
                        linestyle=type_style_map.get(str(org_val), "-"),
                        linewidth=2.2,
                        alpha=0.85,
                    )

                # Y range needs to include the largest flat-offset bump if any.
                max_offset = max(flat_offsets.values()) if flat_offsets else 0.0
                y_top = max(cur_max, max_offset) * 1.08
                ax_sub.set_ylim(-0.08 * y_top, y_top)

                ax_sub.set_ylabel(
                    "Mean cumulative dead fraction",
                    fontsize=10,
                )
                ax_sub.set_xlabel("Timepoint", fontsize=11)
                ax_sub.grid(True, linestyle=":", alpha=0.5)

                if flat_groups:
                    ax_sub.text(
                        0.01, 0.95,
                        "* flat-zero curves slightly offset for visibility",
                        transform=ax_sub.transAxes,
                        fontsize=7,
                        style="italic",
                        color="dimgray",
                        verticalalignment="top",
                    )

                # Align the subpanel's x-range with the top panels so timepoints
                # match visually across the full figure.
                top_xlims = [ax_4a.get_xlim(), ax_4b.get_xlim()]
                common_xlim = (
                    min(xl[0] for xl in top_xlims),
                    max(xl[1] for xl in top_xlims),
                )
                ax_4a.set_xlim(common_xlim)
                ax_4b.set_xlim(common_xlim)
                ax_sub.set_xlim(common_xlim)

                # When subpanel is present the top panels hand off the x label.
                ax_4a.set_xlabel("")
                ax_4b.set_xlabel("")
                ax_4a.tick_params(axis="x", labelbottom=False)
                ax_4b.tick_params(axis="x", labelbottom=False)

            # ----- Single combined figure-level legend in the right gutter -----
            # Anchored to fig.transFigure so legends are guaranteed to live
            # inside the figure canvas regardless of axes positioning.
            color_legend_title = (
                color_key.replace("_", " ").title() if isinstance(color_key, str) else "Color"
            )
            style_legend_title = (
                style_key.replace("_", " ").title() if isinstance(style_key, str) else "Style"
            )

            color_levels_combo = sorted(condition_color_map.keys())
            color_handles_combo = [
                Line2D([0], [0], color=condition_color_map[lvl], lw=3, linestyle="-", label=str(lvl))
                for lvl in color_levels_combo
            ]
            leg_color_combo = Legend(
                fig_combo,
                color_handles_combo,
                [h.get_label() for h in color_handles_combo],
                title=color_legend_title,
                loc="upper left",
                bbox_to_anchor=(legend_gutter_left, 0.92),
                bbox_transform=fig_combo.transFigure,
                prop={"size": 8},
                frameon=True,
            )
            fig_combo.add_artist(leg_color_combo)

            style_levels_combo = sorted(type_style_map.keys())
            style_handles_combo = [
                Line2D([0], [0], color="black", lw=3, linestyle=type_style_map[lvl], label=str(lvl))
                for lvl in style_levels_combo
            ]
            leg_style_combo = Legend(
                fig_combo,
                style_handles_combo,
                [h.get_label() for h in style_handles_combo],
                title=style_legend_title,
                loc="upper left",
                bbox_to_anchor=(legend_gutter_left, 0.70),
                bbox_transform=fig_combo.transFigure,
                prop={"size": 8},
                frameon=True,
                handlelength=3.5,
                handletextpad=0.8,
            )
            fig_combo.add_artist(leg_style_combo)

            dead_handle_combo = Line2D(
                [0], [0],
                marker="$D$", color="red", lw=0,
                markersize=9, markerfacecolor="red", markeredgecolor="red",
                markeredgewidth=0.8,
                label="First dead event (dead=True)",
            )
            leg_dead_combo = Legend(
                fig_combo,
                [dead_handle_combo],
                [dead_handle_combo.get_label()],
                title="Overlay",
                loc="upper left",
                bbox_to_anchor=(legend_gutter_left, 0.50),
                bbox_transform=fig_combo.transFigure,
                prop={"size": 8},
                frameon=True,
            )
            fig_combo.add_artist(leg_dead_combo)

            dis_handle_combo = Line2D(
                [0], [0],
                marker="x", color="black", lw=0,
                markersize=7, markeredgewidth=1.6,
                label="Track disappears (flat tail kept)",
            )
            leg_dis_combo = Legend(
                fig_combo,
                [dis_handle_combo],
                [dis_handle_combo.get_label()],
                title="Markers",
                loc="upper left",
                bbox_to_anchor=(legend_gutter_left, 0.38),
                bbox_transform=fig_combo.transFigure,
                prop={"size": 8},
                frameon=True,
            )
            fig_combo.add_artist(leg_dis_combo)

            # Single figure-level info box (threshold + track counts) in the
            # gutter, below all legends. Replaces the per-panel duplicates so
            # the same info isn't shown twice across the two top panels.
            # Anchored in figure coordinates and constrained to the gutter
            # width so it can't overlap legends above or get clipped at the
            # right edge when bbox_inches="tight" tightens the saved figure.
            info_text = f"{threshold_annotation}\n{counts_annotation}"
            fig_combo.text(
                legend_gutter_left,
                0.26,
                info_text,
                transform=fig_combo.transFigure,
                fontsize=8,
                verticalalignment="top",
                horizontalalignment="left",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="lavender", alpha=0.9),
            )

            # (No separate subpanel legend: each subpanel line's identity is
            # already fully described by the top color (line_condition) and
            # linestyle (organoid_type) legends above.)

            # Inline display: smaller notebook rendering, full size for PDF.
            orig_dpi = fig_combo.get_dpi()
            if screen_show_scale != 1.0:
                fig_combo.set_dpi(orig_dpi * screen_show_scale)
            if show_in_notebook:
                plt.show()
            if screen_show_scale != 1.0:
                fig_combo.set_dpi(orig_dpi)

            pdf.savefig(fig_combo, bbox_inches="tight")
            plt.close(fig_combo)


            # --- Plots 5a & 5b: mean +/- SEM per condition + organoid type ---
            # Uses the same fully-processed data as plot 4 (extend -> normalize -> cummax),
            # but aggregates by pooling all tracks across samples at each timepoint
            # (no per-sample averaging step), so n = number of contributing tracks.
            for feat, plot_label in [
                ("smoothed_percentage_dead_mask", "5a"),
                ("smoothed_nr_dead_mask_pixels", "5b"),
            ]:
                fig5 = plot_grouped_dead_signal(
                    df_processed=df_processed,
                    feature=feat,
                    figsize=figsize,
                    color_key=color_key,
                    style_key=style_key,
                    condition_color_map=condition_color_map,
                    type_style_map=type_style_map,
                    threshold_annotation=threshold_annotation,
                    counts_annotation=counts_annotation,
                    screen_show_scale=screen_show_scale,
                    show_in_notebook=False,
                )
                pdf.savefig(fig5, bbox_inches="tight")
                plt.close(fig5)

                
    print(f"PDF saved to: {pdf_path}")
    print(f"### DONE\n")
    
    return pdf_path

def run_organoid_morphology_dead_analysis(
    output_dir,
    organoid_types,
    initial_window=5,
    umap_n_neighbors=15,
    umap_min_dist=0.1,
    umap_spread=1.0,
    distance_metric = "euclidean",
    leiden_resolution=1.0,
    metadata=None,
    show_in_notebook=True,
):
    """
    Initial Morphology Related to Death: PCA + Leiden + UMAP

    Analyses whether the initial morphology of organoids (averaged over a
    selected time window) is predictive of whether they eventually die.

    Parameters
    ----------
    output_dir : str or Path
        Base output directory containing track features.
    organoid_types : list of str
        List of organoid types to analyze (e.g. ["organoid1"]).
    initial_window : int or tuple(int, int), default=5
        Window for baseline morphology. If int, uses 0 <= t <= N.
        If tuple, uses t_start <= t <= t_end. Features are averaged over this window.
    umap_n_neighbors : int, default=15
        Number of local neighbors used for graph construction.
    umap_min_dist : float, default=0.1
        Minimum distance between points in UMAP embedding.
    umap_spread : float, default=1.0
        The effective scale of embedded points.
    distance_metric : str, default="euclidean"
        Distance metric for UMAP/neighbors. Valid options include: euclidean,
        manhattan, chebyshev, cosine, correlation, canberra, braycurtis, etc.
        See umap-learn documentation for full list.
    leiden_resolution : float, default=1.0
        Controls cluster granularity. Higher values lead to more clusters.
    metadata : pd.DataFrame or None, default=None
        Metadata table (needed for organoid snapshot generation).
        If None, snapshots are skipped.

    Notes
    -----
    - Features: 4 key morphology descriptors (volume, sphericity,
      solidity, surface-to-volume ratio).
    - Normalization: Global Z-score across all features and samples.
    - Pipeline: PCA -> Neighborhood Graph (on top PCs) -> Leiden Clustering -> UMAP.
    - Data: All organoids from all types and samples are mixed for a shared morphological space.

    Outputs (in {output_dir}/analysis/morpho-dead/)
    ------------------------------------------------
    - morpho_dead_umap.pdf          : UMAP cluster + feature grid (T-cell style)
    - morpho_dead_heatmap.pdf       : Cluster × feature heatmap + violin plots
    - morpho_dead_dynamics.pdf      : Death dynamics per organoid, colored by cluster
    - morpho_dead_snapshots.pdf     : Organoid image crops per cluster at window start/end
    - morpho_dead_organoid_level.csv: Per-organoid results
    - morpho_dead_cluster_summary.csv: Cluster-level summary
    """
    # Hardcoded parameters for consistency
    random_state = 42
    #distance_metric = "euclidean"
    print("--------------- Initial Morphology Related to Death ---------------")

    output_dir = Path(output_dir)
    results_outdir = Path(output_dir, "analysis", "morpho-dead")
    results_outdir.mkdir(parents=True, exist_ok=True)

    # morphology features
    base_feature_cols = [
        'volume',
        'sphericity',
        'solidity',
        'surface_to_volume_ratio',
    ]

    required_cols = ["TrackID", "sample_name", "relative_time", "dead"] + base_feature_cols
    if not organoid_types:
        print("WARNING: No organoid types provided.")
        return None

    # Load data
    dfs = []
    missing = []
    for org_type in organoid_types:
        csv_path = Path(
            output_dir,
            "analysis",
            org_type,
            "track_features",
            f"BEHAV3D_{org_type}_combined_track_features_filtered.csv",
        )
        if not csv_path.exists():
            missing.append(org_type)
            continue
        df = pd.read_csv(csv_path, low_memory=False)
        df["organoid_type"] = org_type
        dfs.append(df)

    if missing:
        print(f"WARNING: Missing track features for: {', '.join(missing)}")
    if not dfs:
        print("WARNING: No track features found. Run track filtering/feature extraction first.")
        return None

    df_all = pd.concat(dfs, ignore_index=True)

    missing_cols = [c for c in required_cols if c not in df_all.columns]
    if missing_cols:
        print(f"WARNING: Missing required columns: {missing_cols}")
        return None

    if "SegmentID" not in df_all.columns:
        df_all["SegmentID"] = 0
    df_all["SegmentID"] = df_all["SegmentID"].fillna(0)

    # dead column boolean
    if df_all["dead"].dtype != bool:
        df_all["dead"] = (
            df_all["dead"]
            .astype(str)
            .str.lower()
            .isin(["true", "1", "t", "yes"])
        )

    df_all["relative_time"] = pd.to_numeric(df_all["relative_time"], errors="coerce")
    df_all = df_all[df_all["relative_time"].notna()].copy()

    df_all["TrackID"] = df_all["TrackID"].astype(str)
    df_all["SegmentID"] = df_all["SegmentID"].astype(str)
    df_all["organoid_id"] = (
        df_all["sample_name"].astype(str)
        + "_"
        + df_all["organoid_type"].astype(str)
        + "_"
        + df_all["TrackID"]
        + "_"
        + df_all["SegmentID"]
    )
    df_all = df_all.sort_values(by=["organoid_id", "relative_time"])
    df_all["dead"] = df_all.groupby("organoid_id")["dead"].transform(lambda x: x.cummax())

    # Select initial morphology window
    if isinstance(initial_window, (list, tuple)):
        w_start, w_end = initial_window
        print(f"\n  === Selecting Window Range ({w_start} <= t <= {w_end}) ===")
        df_init = df_all[(df_all["relative_time"] >= w_start) & (df_all["relative_time"] <= w_end)].copy()
        df_after = df_all[df_all["relative_time"] > w_end].copy()
    else:
        w_start = 0
        w_end = int(initial_window)
        print(f"\n  === Selecting Window Range ({w_start} <= t <= {w_end}) ===")
        df_init = df_all[df_all["relative_time"] <= w_end].copy()
        df_after = df_all[df_all["relative_time"] > w_end].copy()
    
    if df_init.empty:
        print("WARNING: No data within the initial morphology window.")
        return None

    # Organoids dead within initial window are excluded
    dead_in_init = df_init.groupby("organoid_id")["dead"].any()

    # Death after initial window
    if df_after.empty:
        print(f"WARNING: No data found after the initial window. "
              "All organoids will be labelled as 'alive'. Consider reducing the initial window.")
    dead_after = df_after.groupby("organoid_id")["dead"].any()

    # Feature extraction (mean morphology over window)
    print("\nCalculating baseline morphology (mean over window)")
    agg_dict = {
        "sample_name": ("sample_name", "first"),
        "organoid_type": ("organoid_type", "first"),
        "TrackID": ("TrackID", "first"),
        "SegmentID": ("SegmentID", "first"),
        "n_timepoints": ("relative_time", "nunique"),
    }
    
    # static features (mean)
    for col in base_feature_cols:
        agg_dict[f"{col}_mean"] = (col, "mean")
    
    df_feats = (
        df_init.groupby("organoid_id")
        .agg(**agg_dict)
        .reset_index()
    )

    feature_cols = [f"{col}_mean" for col in base_feature_cols]

    df_feats["dead_in_initial"] = df_feats["organoid_id"].map(dead_in_init).fillna(False)
    df_feats["dead_after"] = df_feats["organoid_id"].map(dead_after).fillna(False)
    
    total_orgs = len(df_feats)
    dead_init_count = int(df_feats["dead_in_initial"].sum())
    print(f"\n  Organoids in initial window: {total_orgs}")
    print(f"  Excluded (dead in initial window): {dead_init_count}")
    df_feats = df_feats[~df_feats["dead_in_initial"]].copy()

    # Drop NaN/inf feature rows
    df_feats[feature_cols] = df_feats[feature_cols].replace([np.inf, -np.inf], np.nan)
    before = len(df_feats)
    df_feats = df_feats.dropna(subset=feature_cols)
    if df_feats.empty:
        print("WARNING: No organoids left after filtering initial-dead and NaNs.")
        return None
    if before != len(df_feats):
        print(f"  Dropped {before - len(df_feats)} organoids due to NaN/inf features.")
    print(f"  Organoids used for clustering: {len(df_feats)}")

    n_obs = len(df_feats)
    if n_obs < 3:
        print("WARNING: Need at least 3 organoids for PCA/UMAP/Leiden clustering.")
        return None

    # Sample distribution
    print("\nSample distribution")
    sample_dist = df_feats.groupby(["organoid_type", "sample_name"]).size().reset_index(name="n_organoids")
    samples_per_type = df_feats.groupby("organoid_type")["sample_name"].nunique()
    
    for org_type in df_feats["organoid_type"].unique():
        df_type_dist = sample_dist[sample_dist["organoid_type"] == org_type]
        n_samples = samples_per_type[org_type]
        total_orgs = df_type_dist["n_organoids"].sum()
        print(f"  {org_type}: {n_samples} sample(s), {total_orgs} organoids")
        for _, row in df_type_dist.iterrows():
            print(f"    - {row['sample_name']}: {row['n_organoids']} organoids")

    #log transformation of any feature????
    '''print("\nLog transformation")
    df_feats[feature_cols] = np.log1p(df_feats[feature_cols])

    # Verify log transformation
    print("\nLog transformation verification")
    for col in feature_cols:
        mean_log = df_feats[col].mean()
        std_log = df_feats[col].std(ddof=0)
        status = "✓" if abs(mean_log) < 0.01 and abs(std_log - 1.0) < 0.1 else "⚠"
        print(f"    {status} {col}_log: μ={mean_log:10.7f}, σ={std_log:10.7f}")
    '''
    # Z-score standardization
    print("\nZ-score standardization")
    
    def _zscore(x):
        """Z-score normalization with safe handling of zero std."""
        std = x.std(ddof=0)
        if std == 0 or np.isnan(std):
            return x * 0.0
        return (x - x.mean()) / std

    df_scaled = df_feats.copy()
    
    for col in feature_cols:
        df_scaled[col + "_z"] = _zscore(df_scaled[col])

    # Verify standardization
    print("\nStandardization verification")
    for col in feature_cols:
        mean_z = df_scaled[col + "_z"].mean()
        std_z = df_scaled[col + "_z"].std(ddof=0)
        status = "✓" if abs(mean_z) < 0.01 and abs(std_z - 1.0) < 0.1 else "⚠"
        print(f"    {status} {col}_z: μ={mean_z:10.7f}, σ={std_z:10.7f}")

    # PCA > Neighbors > Leiden > UMAP
    z_cols = [c + "_z" for c in feature_cols]
    X = df_scaled[z_cols].to_numpy(dtype=float)

    obs_cols = ["organoid_id", "sample_name", "organoid_type", "TrackID", "SegmentID", "n_timepoints", "dead_after"]
    adata = sc.AnnData(X=X, obs=df_scaled[obs_cols].copy())
    adata.var_names = z_cols
    # Store raw (unscaled) features for reference
    adata.layers["raw"] = df_scaled[feature_cols].to_numpy(dtype=float)

    n_neighbors = int(umap_n_neighbors)
    if n_neighbors >= n_obs:
        n_neighbors = max(1, n_obs - 1)
        print(f"\n  Adjusted neighbors to {n_neighbors} due to small dataset.")

    # Step 1: PCA, compute all components for diagnostics
    print("\nRunning PCA -> Neighbors -> Leiden -> UMAP")
    n_comps = min(X.shape[0] - 1, X.shape[1] - 1, 50)
    sc.tl.pca(adata, n_comps=n_comps)

    # Determine n_pcs: enough PCs to explain >= 90% cumulative variance
    var_ratio = adata.uns["pca"]["variance_ratio"]
    cumvar = np.cumsum(var_ratio)
    n_pcs = int(np.searchsorted(cumvar, 0.90) + 1)
    n_pcs = max(2, min(n_pcs, n_comps))  # at least 2, at most all
    print(f"  PCA: {n_comps} components computed")
    print(f"  PCA Variance Ratio (top 5): {', '.join([f'{r:.1%}' for r in var_ratio[:5]])}")
    print(f"  Cumulative Variance: {cumvar[min(n_pcs-1, len(cumvar)-1)]:.1%} in {n_pcs} PCs (threshold: 90%)")

    # Step 2: Neighborhood graph, uses only top n_pcs components
    sc.pp.neighbors(
        adata,
        n_pcs=n_pcs,
        n_neighbors=n_neighbors,
        metric=str(distance_metric),
        random_state=int(random_state),
    )

    # Step 3: Leiden clustering, on the neighborhood graph
    sc.tl.leiden(
        adata,
        flavor="igraph",
        resolution=float(leiden_resolution),
        key_added="cluster",
        random_state=int(random_state),
    )

    # Step 4: UMAP, on the neighborhood graph (for visualization)
    sc.tl.umap(
        adata,
        min_dist=float(umap_min_dist),
        spread=float(umap_spread),
        random_state=int(random_state),
    )

    umap_coords = adata.obsm["X_umap"]
    obs = adata.obs.copy()
    obs["UMAP1"] = umap_coords[:, 0]
    obs["UMAP2"] = umap_coords[:, 1]
    # Leiden cluster labels
    obs["ClusterID"] = (obs["cluster"].astype(int) + 1).astype(str)

    df_out = df_scaled.merge(obs[["organoid_id", "ClusterID", "UMAP1", "UMAP2"]], on="organoid_id", how="left")
    df_out["ClusterID"] = df_out["ClusterID"].astype("category")
    df_out["status"] = df_out["dead_after"].map({True: "Dead", False: "Alive"})


    # Save results
    per_org_csv = results_outdir / "morpho_dead_organoid_level.csv"
    df_out.to_csv(per_org_csv, index=False)

    # Cluster summary with all features
    summary_agg = {
        "n": ("ClusterID", "size"),
        "pct_dead": ("dead_after", "mean"),
    }
    
    # Add mean and std for base morphology features
    for col in base_feature_cols:
        col_name = f"{col}_mean"  # actual column name in df_out
        summary_agg[f"{col}_mean"] = (col_name, "mean")
        summary_agg[f"{col}_std"] = (col_name, "std")
    
    summary = (
        df_out.groupby("ClusterID", observed=True)
        .agg(**summary_agg)
        .reset_index()
    )
    summary["pct_dead"] = summary["pct_dead"] * 100.0
    summary.to_csv(results_outdir / "morpho_dead_cluster_summary.csv", index=False)

    # Print summary table
    print("\n  Cluster Summary")
    for _, row in summary.iterrows():
        feat_str = ", ".join([f"{feat}={row[f'{feat}_mean']:.3f}±{row[f'{feat}_std']:.3f}" 
                              for feat in base_feature_cols])
        print(f"  Cluster {row['ClusterID']}: n={int(row['n'])}, "
              f"dead={row['pct_dead']:.1f}%, {feat_str}")
    print()

    # Composition tables
    comp_dead = df_out.groupby(["dead_after", "ClusterID"], observed=True).size().unstack(fill_value=0)
    comp_dead_pct = (comp_dead.div(comp_dead.sum(axis=1), axis=0) * 100.0).reset_index()
    comp_dead_pct.to_csv(results_outdir / "morpho_dead_cluster_composition_dead_alive.csv", index=False)

    comp_type = df_out.groupby(["ClusterID", "organoid_type"], observed=True).size().unstack(fill_value=0)
    comp_type_pct = (comp_type.div(comp_type.sum(axis=1), axis=0) * 100.0).reset_index()
    comp_type_pct.to_csv(results_outdir / "morpho_dead_cluster_composition_organoid_type.csv", index=False)

    comp_sample = df_out.groupby(["ClusterID", "sample_name"], observed=True).size().unstack(fill_value=0)
    comp_sample_pct = (comp_sample.div(comp_sample.sum(axis=1), axis=0) * 100.0).reset_index()
    comp_sample_pct.to_csv(results_outdir / "morpho_dead_cluster_composition_sample.csv", index=False)

    # Plots
    clusters = sorted(
        df_out["ClusterID"].dropna().unique().tolist(),
        key=lambda x: int(x) if str(x).isdigit() else x,
    )
    def _get_color_palette(n):
        if n <= 10: return sns.color_palette("Set1", n)
        elif n <= 20: return sns.color_palette("tab20", n)
        return sns.color_palette("husl", n)
    cluster_colors = _get_color_palette(len(clusters))
    cluster_color_map = dict(zip(clusters, cluster_colors))

    # Feature columns (continuous) vs metadata columns (categorical)
    umap_feature_cols = [f"{c}_mean" for c in base_feature_cols if f"{c}_mean" in df_out.columns]
    umap_meta_cols = [c for c in ["organoid_type", "sample_name", "status"] if c in df_out.columns]
    info_cols = umap_feature_cols + umap_meta_cols
    sample_cols_set = set(umap_meta_cols)

    print("\nGenerating Plots")

    # PDF 1: UMAP, cluster plot + feature/metadata subplots
    umap_pdf_path = results_outdir / "morpho_dead_umap.pdf"
    print("  Generating UMAP plots...")
    _rows_per_page = 4
    _nr_cols = 2
    _rows_first = 2
    n_sub_plots = len(info_cols)
    _n_rows_total = (n_sub_plots + _nr_cols - 1) // _nr_cols + _rows_first
    _nr_pages = max(1, (_n_rows_total + _rows_per_page - 1) // _rows_per_page)

    with PdfPages(umap_pdf_path) as pdf:
        plot_idx = 0
        for page in range(_nr_pages):
            fig = plt.figure(figsize=(8.27, 11.69))
            gs = GridSpec(_rows_per_page, _nr_cols, figure=fig, wspace=0.3)

            # Page 0: large cluster UMAP spanning top rows
            if page == 0:
                ax = fig.add_subplot(gs[:_rows_first, :])
                sns.scatterplot(
                    data=df_out, x="UMAP1", y="UMAP2", hue="ClusterID",
                    s=20, alpha=0.5, palette="Set1", ax=ax,
                )
                ax.legend(
                    loc="upper left", prop={"size": 8}, bbox_to_anchor=(1, 1),
                    borderpad=0.3, labelspacing=0.4, columnspacing=0.1, frameon=False,
                )
                ax.set_title("ClusterID", fontsize=10, loc="center")
                ax.set_xticks([]); ax.set_yticks([])
                ax.set_xlabel(""); ax.set_ylabel("")

            # Smaller subplots on remaining rows
            start_row = _rows_first if page == 0 else 0
            remaining_axes = [
                fig.add_subplot(gs[r, c])
                for r in range(start_row, _rows_per_page)
                for c in range(_nr_cols)
            ]

            for ax in remaining_axes:
                if plot_idx >= n_sub_plots:
                    ax.remove()
                    continue
                colorcol = info_cols[plot_idx]

                if colorcol in sample_cols_set or df_out.dtypes.get(colorcol) == bool:
                    # Categorical: colored points + legend
                    sns.scatterplot(
                        data=df_out, x="UMAP1", y="UMAP2", hue=colorcol,
                        s=20, alpha=0.5, palette="Set2", ax=ax,
                    )
                    ax.legend(
                        loc="upper left", prop={"size": 6}, bbox_to_anchor=(1, 1),
                        borderpad=0.3, labelspacing=0.3, columnspacing=0.1, frameon=False,
                    )
                else:
                    # Continuous: scatter + colorbar (clipped to 1st–99th percentile)
                    vals = pd.to_numeric(df_out[colorcol], errors="coerce")
                    vmin = np.nanpercentile(vals, 1)
                    vmax = np.nanpercentile(vals, 99)
                    sc_plot = ax.scatter(
                        df_out["UMAP1"], df_out["UMAP2"],
                        c=vals, cmap="viridis", s=20, alpha=0.5,
                        vmin=vmin, vmax=vmax,
                    )
                    cbar = plt.colorbar(sc_plot, ax=ax, fraction=0.046, pad=0.04)
                    cbar.ax.tick_params(labelsize=6)

                ax.set_title(colorcol, fontsize=10, loc="center")
                ax.set_xticks([]); ax.set_yticks([])
                ax.set_xlabel(""); ax.set_ylabel("")
                plot_idx += 1

            fig.subplots_adjust(left=0.05, right=0.85, top=0.95, bottom=0.05)
            pdf.savefig(fig, dpi=600)
            plt.close(fig)

    print(f"  UMAP PDF saved to: {umap_pdf_path}")
    print("  ClusterID")
    print(
        f"    PCs={n_pcs}, Leiden={leiden_resolution}, "
        f"UMAP(neighbors={umap_n_neighbors}, min_dist={umap_min_dist}, spread={umap_spread})"
    )

    # PDF 2: Heatmap + violin plots
    heatmap_pdf_path = results_outdir / "morpho_dead_heatmap.pdf"
    print("  Generating heatmap + violins...")

    # Cluster means (numeric features only)
    value_cols = [c for c in info_cols if c not in sample_cols_set]
    df_for_means = df_out[value_cols + ["ClusterID"]].copy()
    cluster_means = (
        df_for_means.groupby("ClusterID", observed=True)
        .mean(numeric_only=True).reset_index()
    )
    # Min-max scaling for heatmap
    cluster_means_scaled = cluster_means.copy()
    scale_columns = [c for c in cluster_means.columns if c != "ClusterID"]
    X_hm = cluster_means_scaled[scale_columns].apply(pd.to_numeric, errors="coerce")
    X_hm = X_hm.replace([np.inf, -np.inf], np.nan)
    all_nan_cols = X_hm.columns[X_hm.isna().all()].tolist()
    if all_nan_cols:
        X_hm = X_hm.drop(columns=all_nan_cols)
        scale_columns = [c for c in scale_columns if c not in all_nan_cols]
    if len(scale_columns) > 0:
        X_filled = X_hm.fillna(X_hm.median(numeric_only=True))
        cluster_means_scaled[scale_columns] = MinMaxScaler().fit_transform(X_filled[scale_columns])
        df_hm_long = cluster_means_scaled.melt(id_vars="ClusterID", var_name="var", value_name="AU")
        overall_heatmap_data = df_hm_long.pivot(index="var", columns="ClusterID", values="AU")
    else:
        overall_heatmap_data = pd.DataFrame()

    # Violin data
    df_violin = df_out[["ClusterID"] + value_cols].copy()
    for c in value_cols:
        df_violin[c] = pd.to_numeric(df_violin[c], errors="coerce")
    df_violin.replace([np.inf, -np.inf], np.nan, inplace=True)
    df_long = df_violin.melt(id_vars="ClusterID", var_name="var", value_name="value")
    cluster_order = sorted(df_violin["ClusterID"].dropna().unique().tolist())
    feat_names = [c for c in value_cols if c in df_long["var"].unique()]
    n_violin = len(feat_names)
    _vrows = 7
    _vcols = 2
    _v_rows_needed = (n_violin + _vcols - 1) // _vcols
    _v_pages = max(1, (_v_rows_needed + _vrows - 1) // _vrows)

    with PdfPages(heatmap_pdf_path) as pdf:
        # Page 1: full-page heatmap
        fig, ax = plt.subplots(figsize=(8.27, 11.69))
        if not overall_heatmap_data.empty:
            try:
                hm = sns.heatmap(
                    overall_heatmap_data, ax=ax, cmap="viridis",
                    cbar=True, yticklabels=True,
                )
                ax.set_title("Min–Max Scaled Cluster Means", fontsize=14, pad=14)
                ax.set_xlabel("ClusterID", fontsize=10)
                ax.set_ylabel("", fontsize=10)
                ax.tick_params(axis="y", labelsize=6)
                ax.tick_params(axis="x", labelsize=8, rotation=0)
                cbar_hm = hm.collections[0].colorbar
                cbar_hm.ax.tick_params(labelsize=8)
                fig.tight_layout(pad=2.0)
            except Exception:
                ax.text(0.5, 0.5, "Heatmap unavailable", ha="center", va="center")
                ax.axis("off")
        else:
            ax.text(0.5, 0.5, "No features for heatmap", ha="center", va="center")
            ax.axis("off")
        pdf.savefig(fig, dpi=600)
        if show_in_notebook:
            plt.show()
        plt.close(fig)

        # Other pages: violin plots
        v_idx = 0
        for vpage in range(_v_pages):
            fig = plt.figure(figsize=(8.27, 11.69))
            gs = GridSpec(_vrows, _vcols, figure=fig, hspace=1.5, wspace=0.3)
            vaxes = [
                fig.add_subplot(gs[r, c])
                for r in range(_vrows) for c in range(_vcols)
            ]
            for ax in vaxes:
                if v_idx >= n_violin:
                    ax.remove(); continue
                feat = feat_names[v_idx]
                sub = df_long.loc[df_long["var"] == feat, ["ClusterID", "value"]].dropna()
                if sub.empty:
                    ax.text(0.5, 0.5, f"{feat}\n(no data)", ha="center", va="center")
                    ax.axis("off"); v_idx += 1; continue
                try:
                    sns.violinplot(
                        data=sub, x="ClusterID", y="value",
                        order=cluster_order, inner=None, ax=ax, cut=0,
                        color="steelblue",
                    )
                except Exception:
                    ax.text(0.5, 0.5, f"{feat}\n(unavailable)", ha="center", va="center")
                    ax.axis("off"); v_idx += 1; continue
                # orange = mean
                means = sub.groupby("ClusterID", observed=False)["value"].mean().reindex(cluster_order)
                ax.scatter(
                    np.arange(len(cluster_order)), means.values,
                    s=60, color="orange", edgecolor="black", linewidths=0.8, zorder=3,
                )
                ax.set_title(feat, fontsize=9)
                ax.set_xlabel("ClusterID", fontsize=8)
                ax.set_ylabel("Value", fontsize=8)
                ax.tick_params(axis="x", rotation=0, labelsize=7)
                ax.tick_params(axis="y", labelsize=7)
                v_idx += 1

            fig.subplots_adjust(left=0.20, right=0.98, top=0.95, bottom=0.08)
            pdf.savefig(fig, dpi=600)
            plt.close(fig)

    print(f"  Heatmap PDF saved to: {heatmap_pdf_path}")

    # PDF 3: Death ---------------
    dynamics_pdf_path = results_outdir / "morpho_dead_dynamics.pdf"

    has_dead_feature = ("percentage_dead_mask" in df_all.columns
                        or "smoothed_percentage_dead_mask" in df_all.columns)
    has_dead_col = "dead" in df_all.columns

    if has_dead_feature or has_dead_col:
        print("  Generating death dynamics plots...")
        clustered_ids = set(df_out["organoid_id"].unique())
        df_dynamics = df_all[df_all["organoid_id"].isin(clustered_ids)].copy()
        id_to_cluster = df_out.set_index("organoid_id")["ClusterID"].to_dict()
        df_dynamics["ClusterID"] = df_dynamics["organoid_id"].map(id_to_cluster)

        # Individual curves: local smoothing with min_periods=1 so curves start at t=0
        if has_dead_feature and "percentage_dead_mask" in df_dynamics.columns:
            df_dynamics["viable_fraction"] = (
                df_dynamics
                .groupby("organoid_id")["percentage_dead_mask"]
                .transform(lambda s: s.rolling(window=20, min_periods=1, center=False).mean())
            )
            df_dynamics["viable_fraction"] = (1.0 - df_dynamics["viable_fraction"].clip(0, 1))

        # Cluster mean: binary alive from dead column
        if has_dead_col:
            df_dynamics["alive"] = (~df_dynamics["dead"]).astype(int)

        with PdfPages(dynamics_pdf_path) as pdf:
            # --- Individual organoid curves (viable fraction per cluster) ---
            if "viable_fraction" in df_dynamics.columns:
                for cl in clusters:
                    cl_data = df_dynamics[df_dynamics["ClusterID"] == cl]
                    n_orgs = cl_data["organoid_id"].nunique()
                    fig, ax = plt.subplots(figsize=(12, 6))
                    for oid in cl_data["organoid_id"].unique():
                        org_data = cl_data[cl_data["organoid_id"] == oid].dropna(subset=["viable_fraction"])
                        org_data = org_data.sort_values("relative_time")
                        ax.plot(
                            org_data["relative_time"], org_data["viable_fraction"],
                            color=cluster_color_map[cl], alpha=0.5, linewidth=0.8,
                        )
                    ax.axvspan(w_start, w_end, alpha=0.1, color="gray", label="Initial Window")
                    ax.set_xlabel("Relative Time", fontsize=11)
                    ax.set_ylabel("Viable Fraction (1 \u2212 dead mask coverage)", fontsize=11)
                    ax.set_title(f"Death Dynamics \u2014 Cluster {cl} ({n_orgs} organoids)", fontsize=13)
                    ax.set_ylim(-0.05, 1.05)
                    ax.set_xlim(left=0)
                    ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)
                    ax.legend(loc="upper right", prop={"size": 8})
                    plt.tight_layout()
                    pdf.savefig(fig, bbox_inches="tight")
                    plt.close(fig)

            # --- Cluster mean survival curve (binary alive, fixed cohort) ---
            if has_dead_col:
                fig, ax = plt.subplots(figsize=(12, 7))
                for cl in clusters:
                    cl_members = sorted(df_out.loc[df_out["ClusterID"] == cl, "organoid_id"].unique().tolist())
                    if not cl_members:
                        continue
                    cl_data = df_dynamics[df_dynamics["organoid_id"].isin(cl_members)].copy()
                    cl_data = cl_data[["organoid_id", "relative_time", "dead"]].dropna(subset=["relative_time"])
                    if cl_data.empty:
                        continue
                    # First timepoint where each organoid is dead (if never dead -> +inf)
                    first_dead_t = (
                        cl_data[cl_data["dead"] == True]
                        .groupby("organoid_id")["relative_time"]
                        .min()
                    )
                    first_dead_t = first_dead_t.reindex(cl_members).fillna(np.inf)
                    rel_times = np.sort(cl_data["relative_time"].unique())
                    frac_alive = [((first_dead_t > t).sum() / len(cl_members)) for t in rel_times]
                    ax.plot(
                        rel_times, frac_alive,
                        color=cluster_color_map[cl], linewidth=2.5, label=f"Cluster {cl}",
                    )

                ax.axvspan(w_start, w_end, alpha=0.1, color="gray", label="Initial Window")
                ax.set_xlabel("Relative Time", fontsize=11)
                ax.set_ylabel("% Alive", fontsize=11)
                ax.set_title("Death Dynamics per Morphology Cluster", fontsize=13)
                ax.set_ylim(-0.05, 1.05)
                ax.set_xlim(left=0)
                ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)
                ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0, prop={"size": 8})
                plt.tight_layout()
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)

        print(f"  Death dynamics PDF saved to: {dynamics_pdf_path}")
    else:
        print("  WARNING: No death feature found — skipping dynamics plot.")

    # ===== PDF 4: Organoid Snapshots — separate channels per organoid =====
    snapshot_pdf_path = results_outdir / "morpho_dead_snapshots.pdf"
    n_examples_per_cluster = 3
    t_col = "position_t" if "position_t" in df_init.columns else "relative_time"

    org_t_bounds = (
        df_init.groupby("organoid_id")[t_col]
        .agg(["min", "max"])
        .rename(columns={"min": "t_img_start", "max": "t_img_end"})
    )

    # Channel display colormaps (one per channel)
    _ch_cmaps = [
        _LSC.from_list("_g", ["black", "#00FF00"]),   # Ch 0: green
        _LSC.from_list("_m", ["black", "#FF00FF"]),   # Ch 1: magenta
        _LSC.from_list("_c", ["black", "#00FFFF"]),   # Ch 2: cyan
        _LSC.from_list("_y", ["black", "#FFFF00"]),   # Ch 3: yellow
        _LSC.from_list("_o", ["black", "#FF8000"]),   # Ch 4: orange
        _LSC.from_list("_b", ["black", "#8080FF"]),   # Ch 5: light blue
    ]
    _ch_color_names = ["green", "magenta", "cyan", "yellow", "orange", "light blue"]

    _zarr_cache = {}

    def _load_zarr_cached(zpath, cache_key):
        if cache_key in _zarr_cache:
            return _zarr_cache[cache_key]
        if not Path(zpath).exists():
            _zarr_cache[cache_key] = None
            return None
        try:
            z = load_zarr(zpath)
            _zarr_cache[cache_key] = z
            return z
        except Exception as e:
            print(f"    WARNING: Could not load {zpath}: {e}")
            _zarr_cache[cache_key] = None
            return None

    def _crop_per_channel(raw_zarr_data, seg_zarr_data, t_idx, track_id, margin=30):
        """Return (crops_list, contour_mask) for one organoid at one timepoint."""
        t_idx = int(t_idx)
        if seg_zarr_data is None or raw_zarr_data is None:
            return None
        if t_idx < 0 or t_idx >= seg_zarr_data.shape[0]:
            return None

        seg_vol = np.asarray(seg_zarr_data[t_idx])
        seg_2d = seg_vol.max(axis=0) if seg_vol.ndim == 3 else seg_vol
        mask = (seg_2d == track_id)
        if not mask.any():
            return None
        ys, xs = np.where(mask)
        y0 = max(0, ys.min() - margin)
        y1 = min(seg_2d.shape[0], ys.max() + margin + 1)
        x0 = max(0, xs.min() - margin)
        x1 = min(seg_2d.shape[1], xs.max() + margin + 1)

        # Contour: dilated − eroded boundary of the mask
        mask_crop = mask[y0:y1, x0:x1]
        contour = binary_dilation(mask_crop, iterations=1) ^ binary_erosion(mask_crop, iterations=1)

        has_ch = (raw_zarr_data.ndim == 5)
        n_ch = raw_zarr_data.shape[1] if has_ch else 1
        crops = []
        for ci in range(n_ch):
            if has_ch:
                ch_vol = np.asarray(raw_zarr_data[t_idx, ci])
            else:
                ch_vol = np.asarray(raw_zarr_data[t_idx])
            ch_2d = ch_vol.max(axis=0) if ch_vol.ndim == 3 else ch_vol
            crops.append(ch_2d[y0:y1, x0:x1])
        return crops, contour

    # Build organoid_type → channel index map from metadata
    _type_to_channel = {}
    _channel_cols = []
    if metadata is not None:
        _channel_cols = sorted([c for c in metadata.columns if c.startswith("channel") and c[7:].isdigit()])
        if _channel_cols:
            for _, mrow in metadata.iterrows():
                for col in _channel_cols:
                    ch_idx = int(col.replace("channel", ""))
                    label = str(mrow[col]).strip() if pd.notna(mrow[col]) else ""
                    if label and label not in _type_to_channel:
                        _type_to_channel[label] = ch_idx

    def _get_organoid_channel(sample_name, organoid_type):
        """Resolve organoid channel per sample, fallback to global mapping."""
        if metadata is None or not _channel_cols:
            return _type_to_channel.get(organoid_type, -1)
        row = metadata[metadata["sample_name"] == sample_name]
        if row.empty:
            return _type_to_channel.get(organoid_type, -1)
        row0 = row.iloc[0]
        for col in _channel_cols:
            label = str(row0[col]).strip() if pd.notna(row0[col]) else ""
            if label == organoid_type:
                return int(col.replace("channel", ""))
        return _type_to_channel.get(organoid_type, -1)

    # Collect snapshot data per cluster — diverse sampling across types & samples
    snapshot_data = {}
    n_channels_found = 1
    try:
        print("  Generating organoid snapshots (separate channels)...")
        for cl in clusters:
            cl_orgs = df_out[df_out["ClusterID"] == cl]
            # Diverse sampling: one per (organoid_type, sample_name), then fill randomly
            selected = (
                cl_orgs
                .groupby(["organoid_type", "sample_name"], group_keys=False)
                .sample(n=1, random_state=random_state)
                .reset_index(drop=True)
            )
            selected = selected.drop_duplicates(subset=["organoid_id"]).copy()
            remaining = cl_orgs[~cl_orgs["organoid_id"].isin(selected["organoid_id"])]
            n_extra = max(0, n_examples_per_cluster - len(selected))
            if n_extra > 0 and len(remaining) > 0:
                extra = remaining.sample(min(n_extra, len(remaining)), random_state=random_state)
                selected = pd.concat([selected, extra], ignore_index=True)
            selected = selected.drop_duplicates(subset=["organoid_id"]).copy()
            selected = selected.head(n_examples_per_cluster)

            snapshot_data[cl] = []
            for _, org in selected.iterrows():
                oid = org["organoid_id"]
                if oid not in org_t_bounds.index:
                    continue
                ts = int(org_t_bounds.loc[oid, "t_img_start"])
                te = int(org_t_bounds.loc[oid, "t_img_end"])
                try:
                    track_id = int(float(org["TrackID"]))
                except (ValueError, TypeError):
                    continue
                sname = org["sample_name"]
                otype = org["organoid_type"]
                seg_path = output_dir / "images" / sname / f"{sname}_{otype}_tracked.zarr"
                raw_path = output_dir / "images" / sname / f"{sname}.zarr"
                seg_z = _load_zarr_cached(seg_path, f"seg_{sname}_{otype}")
                raw_z = _load_zarr_cached(raw_path, f"raw_{sname}")
                if seg_z is None or raw_z is None:
                    continue
                result_s = _crop_per_channel(raw_z, seg_z, ts, track_id)
                result_e = _crop_per_channel(raw_z, seg_z, te, track_id)
                crops_s, contour_s = result_s if result_s else (None, None)
                crops_e, contour_e = result_e if result_e else (None, None)
                if crops_s is not None:
                    n_channels_found = len(crops_s)
                elif crops_e is not None:
                    n_channels_found = len(crops_e)
                org_ch_idx = _get_organoid_channel(sname, otype)
                if crops_s is not None or crops_e is not None:
                    snapshot_data[cl].append({
                        "crops_s": crops_s, "crops_e": crops_e,
                        "contour_s": contour_s, "contour_e": contour_e,
                        "oid": oid, "otype": otype, "sname": sname,
                        "ts": ts, "te": te, "org_ch": org_ch_idx,
                    })
    except Exception as e:
        print(f"  WARNING: Snapshot generation failed: {e}")
        import traceback; traceback.print_exc()
        snapshot_data = {}
    finally:
        _zarr_cache.clear()

    clusters_with_data = [cl for cl in clusters if snapshot_data.get(cl)]
    if clusters_with_data:
        n_ch = n_channels_found

        def _show_crop(ax, img, cmap, contour=None, highlight=False):
            """Display a single channel crop. Red border if highlight=True."""
            if img is None:
                ax.text(0.5, 0.5, "N/A", ha="center", va="center",
                        transform=ax.transAxes, fontsize=9, color="white")
                ax.set_facecolor("black")
                return
            img_f = img.astype(float)
            pos = img_f[img_f > 0]
            vmin = float(np.percentile(pos, 1)) if pos.size > 0 else 0
            vmax = float(np.percentile(pos, 99)) if pos.size > 0 else max(1, img_f.max())
            ax.imshow(img_f, cmap=cmap, vmin=vmin, vmax=vmax)
            if highlight and contour is not None:
                red_cmap = ListedColormap(["red"])
                contour_overlay = np.ma.masked_where(~contour, contour.astype(float))
                ax.imshow(contour_overlay, cmap=red_cmap, vmin=0, vmax=1, alpha=0.9)

        # Build channel labels from metadata for column headers
        _ch_labels = {}
        if metadata is not None:
            ch_cols = sorted([c for c in metadata.columns if c.startswith("channel") and c[7:].isdigit()])
            if ch_cols:
                first_row = metadata.iloc[0]
                for col in ch_cols:
                    ci = int(col.replace("channel", ""))
                    label = str(first_row[col]).strip() if pd.notna(first_row[col]) else ""
                    if label:
                        _ch_labels[ci] = label

        with PdfPages(snapshot_pdf_path) as pdf:
            for cl in clusters_with_data:
                examples = snapshot_data[cl]
                n_ex = len(examples)
                n_cols = n_ch * 2
                fig, axes = plt.subplots(
                    n_ex, n_cols,
                    figsize=(2.2 * n_cols + 1.5, 2.5 * n_ex + 1.5),
                    squeeze=False,
                )
                fig.suptitle(
                    f"Cluster {cl} \u2014 Organoid Snapshots",
                    fontsize=13, y=0.99, fontweight="bold",
                    color=cluster_color_map.get(cl, "black"),
                )
                fig.subplots_adjust(left=0.14, top=0.88, bottom=0.04, hspace=0.3, wspace=0.08)

                # Column headers (black text) — channel index + label if available
                for ci in range(n_ch):
                    ch_label = _ch_labels.get(ci, "")
                    header = f"Ch {ci}" + (f" ({ch_label})" if ch_label else "")
                    axes[0, ci].set_title(header, fontsize=8, color="black")
                    axes[0, n_ch + ci].set_title(header, fontsize=8, color="black")

                # Timepoint group headers
                fig.text(0.14 + (0.86 - 0.14) * 0.25, 0.92,
                         f"Window Start (t={w_start})", ha="center", fontsize=10, style="italic")
                fig.text(0.14 + (0.86 - 0.14) * 0.75, 0.92,
                         f"Window End (t={w_end})", ha="center", fontsize=10, style="italic")

                for ex_i, snap in enumerate(examples):
                    crops_s = snap["crops_s"]
                    crops_e = snap["crops_e"]
                    contour_s = snap["contour_s"]
                    contour_e = snap["contour_e"]
                    otype = snap["otype"]
                    sname = snap["sname"]
                    org_ch = snap["org_ch"]

                    for ci in range(n_ch):
                        cmap = _ch_cmaps[ci % len(_ch_cmaps)]
                        is_org_ch = (ci == org_ch)
                        img_s = crops_s[ci] if (crops_s and ci < len(crops_s)) else None
                        img_e = crops_e[ci] if (crops_e and ci < len(crops_e)) else None
                        _show_crop(axes[ex_i, ci], img_s, cmap,
                                   contour=contour_s, highlight=is_org_ch)
                        _show_crop(axes[ex_i, n_ch + ci], img_e, cmap,
                                   contour=contour_e, highlight=is_org_ch)
                        axes[ex_i, ci].axis("off")
                        axes[ex_i, n_ch + ci].axis("off")

                        # Red border around the organoid's channel
                        if is_org_ch:
                            for ax, img in [
                                (axes[ex_i, ci], img_s),
                                (axes[ex_i, n_ch + ci], img_e),
                            ]:
                                if img is not None:
                                    h, w = img.shape[:2]
                                    ax.add_patch(
                                        Rectangle(
                                            (-1.0, -1.0), w + 2.0, h + 2.0,
                                            transform=ax.transData,
                                            fill=False, edgecolor="red",
                                            linewidth=2.5, zorder=30, clip_on=False
                                        )
                                    )

                    # Row label
                    row_y = axes[ex_i, 0].get_position().y0 + axes[ex_i, 0].get_position().height / 2
                    fig.text(
                        0.01, row_y, f"{otype}\n{sname}",
                        fontsize=7, va="center", ha="left", style="italic",
                    )

                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)

        print(f"  Snapshots PDF saved to: {snapshot_pdf_path}")
    else:
        print("  WARNING: Could not generate snapshots — no images found.")

    # ===== PDF 5: Composition bar plots =====
    comp_pdf_path = results_outdir / "morpho_dead_composition.pdf"
    print("  Generating composition bar plots...")

    with PdfPages(comp_pdf_path) as pdf:
        # Bar 1: Dead vs Alive composition by cluster
        comp_dead_plot = comp_dead.div(comp_dead.sum(axis=1), axis=0) * 100.0
        comp_dead_plot = comp_dead_plot.reindex([True, False]).fillna(0.0)
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        x_labels = ["Dead", "Alive"]
        x = np.arange(len(x_labels))
        bottom = np.zeros(len(x_labels))
        for cl in clusters:
            vals = comp_dead_plot.get(cl, pd.Series([0.0, 0.0])).values
            ax.bar(x, vals, bottom=bottom, label=f"Cluster {cl}", color=cluster_color_map.get(cl))
            bottom += vals
        ax.set_xticks(x); ax.set_xticklabels(x_labels)
        ax.set_ylabel("Percent of group"); ax.set_ylim(0, 100)
        ax.set_title("Cluster Composition: Dead vs Alive")
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0, prop={"size": 7})
        ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # Bar 2: Organoid type composition per cluster
        comp_type_plot = comp_type.div(comp_type.sum(axis=1), axis=0) * 100.0
        comp_type_plot = comp_type_plot.reindex(clusters).fillna(0.0)
        types = comp_type_plot.columns.tolist()
        type_colors = _get_color_palette(len(types))
        type_color_map = dict(zip(types, type_colors))
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        x = np.arange(len(comp_type_plot.index))
        bottom = np.zeros(len(comp_type_plot.index))
        for t in types:
            vals = comp_type_plot[t].values
            ax.bar(x, vals, bottom=bottom, label=t, color=type_color_map[t])
            bottom += vals
        ax.set_xticks(x); ax.set_xticklabels([str(c) for c in comp_type_plot.index])
        ax.set_ylabel("Percent in cluster"); ax.set_ylim(0, 100)
        ax.set_title("Organoid Type Composition per Cluster")
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0, prop={"size": 7})
        ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # Bar 3: Sample composition per cluster
        comp_sample_plot = comp_sample.div(comp_sample.sum(axis=1), axis=0) * 100.0
        comp_sample_plot = comp_sample_plot.reindex(clusters).fillna(0.0)
        samples = comp_sample_plot.columns.tolist()
        sample_colors = _get_color_palette(len(samples))
        sample_color_map = dict(zip(samples, sample_colors))
        fig, ax = plt.subplots(1, 1, figsize=(12, 6))
        x = np.arange(len(comp_sample_plot.index))
        bottom = np.zeros(len(comp_sample_plot.index))
        for s in samples:
            vals = comp_sample_plot[s].values
            ax.bar(x, vals, bottom=bottom, label=s, color=sample_color_map[s])
            bottom += vals
        ax.set_xticks(x); ax.set_xticklabels([str(c) for c in comp_sample_plot.index])
        ax.set_ylabel("Percent in cluster"); ax.set_ylim(0, 100)
        ax.set_title("Sample Composition per Cluster")
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0, prop={"size": 6})
        ax.grid(True, linestyle=":", linewidth=0.5, alpha=0.5)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # Box 4: Time of death per cluster (boxplot + individual organoid points)
        dead_orgs = df_out[df_out["dead_after"] == True].copy()
        if not dead_orgs.empty:
            dead_ids = set(dead_orgs["organoid_id"].unique())
            df_death_times = (
                df_all[df_all["organoid_id"].isin(dead_ids) & (df_all["dead"] == True)]
                .groupby("organoid_id")["relative_time"]
                .min()
                .reset_index()
                .rename(columns={"relative_time": "death_time"})
            )
            df_death_times = df_death_times.merge(
                df_out[["organoid_id", "ClusterID"]].drop_duplicates(),
                on="organoid_id", how="left",
            )
            df_death_times = df_death_times.dropna(subset=["ClusterID", "death_time"])

            if not df_death_times.empty:
                fig, ax = plt.subplots(1, 1, figsize=(10, 7))
                sns.boxplot(
                    data=df_death_times, x="ClusterID", y="death_time",
                    order=clusters, color="white", fliersize=0,
                    boxprops=dict(edgecolor="black"),
                    medianprops=dict(color="black", linewidth=2),
                    whiskerprops=dict(color="black"),
                    capprops=dict(color="black"),
                    ax=ax,
                )
                sns.stripplot(
                    data=df_death_times, x="ClusterID", y="death_time",
                    order=clusters, ax=ax, size=5, alpha=0.6, jitter=0.2,
                    palette=cluster_color_map, hue="ClusterID", hue_order=clusters,
                    legend=False,
                )
                ax.set_xlabel("Morphology Cluster", fontsize=11)
                ax.set_ylabel("Death Timepoint (relative time)", fontsize=11)
                ax.set_title("Organoid Death Timepoint per Morphology Cluster", fontsize=13)
                ax.grid(True, axis="y", linestyle=":", linewidth=0.5, alpha=0.5)
                plt.tight_layout()
                pdf.savefig(fig, bbox_inches="tight")
                if show_in_notebook:
                    plt.show()
                plt.close(fig)

    print(f"  Composition PDF saved to: {comp_pdf_path}")

    print(f"\n  Results saved to: {results_outdir}")
    print("\n### ANALYSIS COMPLETE ###\n")
    return results_outdir
