import time
import random
import pickle
from contextlib import nullcontext
from types import SimpleNamespace
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
from matplotlib.colors import to_hex
import seaborn as sns
import umap
import scanpy as sc

from sklearn.cluster import KMeans, HDBSCAN, AgglomerativeClustering
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    confusion_matrix,
    classification_report,
)
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from sklearn.decomposition import PCA
from sklearn.feature_selection import VarianceThreshold


import re
import anndata as ad
import yaml
try:
    from dtaidistance import dtw_ndim
except Exception:
    dtw_ndim = None
from sklearn.metrics import silhouette_score
from behav3d.analysis.behavior.state.classification import FULL_STATE_COL
from behav3d.analysis.behavior.track.feature_dtw import run_tcell_analysis
from behav3d.analysis.behavior.track.visualization.plots.feature_dtw import (
    plot_cluster_percentage_bars,
    plot_clustering_feature_heatmap,
    plot_feature_umap,
)
from behav3d.core.metadata import load_behav3d_metadata, check_behav3d_metadata
from behav3d.analysis.behavior.general import relabel_cluster_ids
from behav3d.features.state_descriptive_features import (
    extract_descibing_track_state_features, 
    scale_feature_blocks, 
    l2_normalize_features_blocks,
    drop_highly_correlated_features,
    drop_low_variance_features
)
from behav3d.analysis.filtering import filter_and_truncate_tracks_anndata
from behav3d.analysis.behavior.general.leiden import (
    run_pca, 
    run_leiden_clustering
)
from behav3d.analysis.behavior.state.visualization.backprojection import (
    export_behavioral_state_backprojection_zarrs,
    show_behavioral_state_backprojection,
)
from behav3d.analysis.behavior.track.visualization.plots.exemplar_track_per_cluster import (
    plot_exemplar_tracks_by_cluster,
    save_exemplar_statebar_track_pdf_per_cluster,
    save_exemplar_statebar_backprojection_pdf,
    save_exemplar_statebar_backprojection_video_per_cluster,
    select_exemplar_tracks_by_cluster,
    _build_state_color_map,
    _compute_state_bar_segments,
    _plot_statebar_segments_on_ax,
    _prepare_exemplar_backprojection_data,
)
from behav3d.analysis.behavior.track.visualization.plots.exemplar_coordinate_utils import (
    ensure_exemplar_coordinate_columns as _shared_ensure_exemplar_coordinate_columns,
    merge_coordinate_columns_into_obs as _shared_merge_coordinate_columns_into_obs,
    resolve_exemplar_positions_csv_path as _shared_resolve_exemplar_positions_csv_path,
)
from behav3d.analysis.behavior.general.visualization.plots import plot_top_ranking_features
from behav3d.analysis.behavior.general.visualization.plots.proportion_bars import (
    draw_thin_stacked_proportion_barh,
    draw_stacked_proportion_barv,
    compute_class_by_stack_proportions,
    compute_condition_diff_stats_pairwise,
    legend_layout,
    plot_condition_diff_grid,
    plot_condition_diff_grid_2d,
    plot_page_stacked_proportion_barh_grid,
    stacked_proportion_barh_rows_per_page,
    _resolve_effective_group_cols,
    _make_group_label,
    _chunk_list,
    _wrap_row_label,
)
from behav3d.analysis.behavior.general.visualization.plots.feature_violin import (
    plot_feature_violin_by_group,
)
from behav3d.analysis.behavior.track.contact_grouping import (
    compute_track_contact_features,
    merge_track_contact_features_into_obs,
    _contact_group_col_name,
    _contact_mean_col_name,
    _contact_max_bout_col_name,
)
from behav3d.analysis.behavior.state.utils import (
    _apply_state_order,
    _get_classification_state_colors,
    _get_classification_state_order,
    _normalize_label_color_map,
)
from behav3d.analysis.behavior.utils import (
    _mixed_label_sort_key,
    _resolve_output_dir,
    _sanitize_filename_token,
    _to_numpy_2d,
    _vdone,
    _vinfo,
    _vsave,
    _vstart,
)


_GRID_PANEL_ASPECT = 1.0 / 6.0   # target panel height:width ratio -> thin bars
_GRID_PANEL_W_IN = 1.6           # comfortable width per grouped-grid column
_GRID_HEADER_H_IN = 0.3          # space for column-header labels (2D grid)
_GRID_HEADER_W_IN = 1.2          # space for row-header labels (2D grid) — horizontal, up to 2 wrapped lines
_GRID_ROW_LABEL_MAX_CHARS = 13   # chars/line at _GRID_HEADER_W_IN, fontsize=9 bold
_GRID_TOP_MARGIN_IN = 0.5        # space reserved above the grid for the suptitle
_GRID_BOTTOM_MARGIN_IN = 0.55    # space reserved below the grid for the legend


def archive_track_clustering_pdfs(outfolder, archive_dir_name="clustering_originals"):
    outfolder = Path(outfolder)
    outfolder.mkdir(parents=True, exist_ok=True)

    pdf_paths = sorted([p for p in outfolder.glob("*.pdf") if p.is_file()])
    if len(pdf_paths) == 0:
        return {"archive_dir": None, "archived_paths": []}

    ts = time.strftime("%Y%m%d_%H%M%S")
    archive_dir = outfolder / str(archive_dir_name) / ts
    archive_dir.mkdir(parents=True, exist_ok=True)

    archived_paths = []
    for src in pdf_paths:
        dst = archive_dir / src.name
        suffix_idx = 1
        while dst.exists():
            dst = archive_dir / f"{src.stem}_{suffix_idx}{src.suffix}"
            suffix_idx += 1
        src.replace(dst)
        archived_paths.append(dst)

    return {"archive_dir": archive_dir, "archived_paths": archived_paths}


def _apply_best_pdf_orientation(fig, default_orientation="landscape"):
    """
    Keep figure orientation matched to content shape before PDF export.
    """
    if fig is None:
        return None
    default_orientation = str(default_orientation).strip().lower()
    width, height = fig.get_size_inches()
    if abs(width - height) < 0.05:
        if default_orientation == "portrait":
            fig.set_size_inches(height, max(width, height * 1.25), forward=True)
        else:
            fig.set_size_inches(max(width, height * 1.25), height, forward=True)
    width, height = fig.get_size_inches()
    return "landscape" if float(width) >= float(height) else "portrait"


def _rank_cluster_counts(adata_tracks, cluster_key, cluster_order=None):
    if cluster_key not in adata_tracks.obs.columns:
        raise ValueError(f"Missing '{cluster_key}' in adata_tracks.obs.")

    labels = (
        pd.Series(adata_tracks.obs[cluster_key])
        .astype("string")
        .dropna()
        .astype(str)
    )
    counts = labels.value_counts(dropna=False).to_dict()
    if cluster_order:
        ranked = [(str(label), int(counts.get(str(label), 0))) for label in cluster_order]
        leftovers = sorted(
            [(label, n) for label, n in counts.items() if str(label) not in {str(o) for o in cluster_order}],
            key=lambda x: (-int(x[1]), _mixed_label_sort_key(x[0])),
        )
        ranked = ranked + leftovers
    else:
        ranked = sorted(
            [(str(label), int(n)) for label, n in counts.items()],
            key=lambda x: (-int(x[1]), _mixed_label_sort_key(x[0])),
        )
    return ranked


def _plot_cluster_rankings(adata_tracks, cluster_key, cluster_order=None, cluster_colors=None):
    ranked = _rank_cluster_counts(adata_tracks=adata_tracks, cluster_key=cluster_key, cluster_order=cluster_order)
    cluster_labels = [x[0] for x in ranked]
    cluster_counts = [x[1] for x in ranked]
    bar_colors = (
        [cluster_colors.get(label, "#2E6FBA") for label in cluster_labels]
        if cluster_colors
        else "#2E6FBA"
    )

    fig_h = max(4.5, 0.45 * len(cluster_labels) + 2.0)
    fig, ax = plt.subplots(figsize=(10.5, fig_h))
    y = np.arange(len(cluster_labels))
    bars = ax.barh(y, cluster_counts, color=bar_colors)
    ax.set_yticks(y)
    ax.set_yticklabels(cluster_labels)
    ax.invert_yaxis()
    ax.set_xlabel("N tracks")
    ax.set_ylabel(str(cluster_key))
    ax.set_title(f"Cluster occurrence ranking ({cluster_key})")
    ax.grid(axis="x", alpha=0.2)
    for bar, n in zip(bars, cluster_counts):
        ax.text(
            bar.get_width(),
            bar.get_y() + bar.get_height() / 2.0,
            f" {int(n)}",
            va="center",
            ha="left",
            fontsize=9,
        )
    fig.tight_layout()
    return fig




def _compute_grouped_class_proportions(df, *, group_cols, class_col, class_order):
    """
    Bucket rows by group label; return per-group class proportions, per-group
    raw class counts, and per-group total track counts.
    """
    group_labels = _make_group_label(df, group_cols)
    tmp = df.copy()
    tmp["_group_label"] = group_labels.values
    unique_groups = tmp["_group_label"].dropna().unique().tolist()

    proportions_by_group = {}
    counts_by_group = {}
    n_by_group = {}
    for grp in unique_groups:
        panel = tmp[tmp["_group_label"] == grp]
        counts = panel[class_col].value_counts().reindex(class_order, fill_value=0).astype(int)
        total = int(counts.sum())
        proportions = (counts / total) if total > 0 else counts.astype(float)
        proportions_by_group[str(grp)] = proportions
        counts_by_group[str(grp)] = counts
        n_by_group[str(grp)] = total

    return proportions_by_group, counts_by_group, n_by_group


def _class_legend_handles(class_order, colors):
    return (
        [
            plt.Line2D([0], [0], marker="s", linestyle="none", color=colors[c], label=str(c), markersize=7)
            for c in class_order
        ],
        [str(c) for c in class_order],
    )


def _plot_page_grouped_class_proportions_2d_grid(
    data_by_group,
    *,
    group_cols,
    unique_vals_per_col,
    class_order,
    class_colors,
    group_label_title,
    axis_cols=None,
    panel_w_in=_GRID_PANEL_W_IN,
    group_title_fontsize=9,
):
    """
    One page: grouped track-class proportions arranged in a 2D grid.

    For 1 group_col: single column of panels, one row per unique value.
    For 2 group_cols: by default the column with more unique values goes to
    rows (Y) and the column with fewer unique values goes to columns (X);
    pass ``axis_cols=(col_y, col_x)`` to pick the axes explicitly instead.
    Header labels are drawn on each axis.

    Figure size scales with the number of X/Y conditions (rather than being
    clipped to a fixed page), and each panel keeps a thin, consistent
    ``_GRID_PANEL_ASPECT`` (~1:6) height:width ratio regardless of how many
    conditions are present.
    """
    panel_h_in = panel_w_in * _GRID_PANEL_ASPECT
    legend_ncol, _, legend_margin_in = legend_layout(len(class_order), base_margin_in=_GRID_BOTTOM_MARGIN_IN)

    if len(group_cols) == 2:
        if axis_cols is not None:
            col_y, col_x = axis_cols
        else:
            col_y, col_x = sorted(group_cols, key=lambda c: -len(unique_vals_per_col[c]))
        col_x_vals = unique_vals_per_col[col_x]
        col_y_vals = unique_vals_per_col[col_y]

        nrows_data = max(1, len(col_y_vals))
        ncols_data = max(1, len(col_x_vals))

        fig_w = _GRID_HEADER_W_IN + ncols_data * panel_w_in
        fig_h = _GRID_TOP_MARGIN_IN + _GRID_HEADER_H_IN + nrows_data * panel_h_in + legend_margin_in
        fig = plt.figure(figsize=(fig_w, fig_h))
        outer = GridSpec(
            nrows=nrows_data + 1,
            ncols=ncols_data + 1,
            figure=fig,
            height_ratios=[_GRID_HEADER_H_IN] + [panel_h_in] * nrows_data,
            width_ratios=[_GRID_HEADER_W_IN] + [panel_w_in] * ncols_data,
            hspace=0.6,
            wspace=0.3,
            top=1.0 - _GRID_TOP_MARGIN_IN / fig_h,
            bottom=legend_margin_in / fig_h,
        )

        ax_corner = fig.add_subplot(outer[0, 0])
        ax_corner.axis("off")

        for c, col_x_val in enumerate(col_x_vals):
            ax_ch = fig.add_subplot(outer[0, c + 1])
            ax_ch.axis("off")
            ax_ch.text(
                0.5, 0.5, str(col_x_val),
                ha="center", va="center", fontsize=group_title_fontsize, fontweight="bold",
                transform=ax_ch.transAxes,
            )

        for r, col_y_val in enumerate(col_y_vals):
            ax_rh = fig.add_subplot(outer[r + 1, 0])
            ax_rh.axis("off")
            ax_rh.text(
                0.93, 0.5, _wrap_row_label(str(col_y_val), max_chars=_GRID_ROW_LABEL_MAX_CHARS),
                ha="right", va="center", ma="right", fontsize=group_title_fontsize, fontweight="bold",
                transform=ax_rh.transAxes,
            )

            for c, col_x_val in enumerate(col_x_vals):
                vals_map = {col_y: col_y_val, col_x: col_x_val}
                key = " | ".join(str(vals_map[gc]) for gc in group_cols)

                ax = fig.add_subplot(outer[r + 1, c + 1])

                if key not in data_by_group:
                    ax.text(0.5, 0.5, "—", ha="center", va="center", fontsize=10, transform=ax.transAxes)
                    ax.axis("off")
                    continue

                draw_thin_stacked_proportion_barh(ax, data_by_group[key], class_order, class_colors, xmax=1.0)

        title_str = f"Grouped Track-Class Proportions — {col_x} (columns) × {col_y} (rows)"

    else:
        col_y = group_cols[0]
        col_y_vals = unique_vals_per_col[col_y]
        nrows_data = max(1, len(col_y_vals))

        fig_w = max(4.0 * panel_w_in, 6.0)
        fig_h = _GRID_TOP_MARGIN_IN + nrows_data * panel_h_in + legend_margin_in
        fig = plt.figure(figsize=(fig_w, fig_h))
        outer = GridSpec(
            nrows=nrows_data,
            ncols=1,
            figure=fig,
            height_ratios=[panel_h_in] * nrows_data,
            hspace=0.5,
            top=1.0 - _GRID_TOP_MARGIN_IN / fig_h,
            bottom=legend_margin_in / fig_h,
        )

        for r, col_y_val in enumerate(col_y_vals):
            key = str(col_y_val)
            ax = fig.add_subplot(outer[r])

            if key not in data_by_group:
                ax.text(0.5, 0.5, "—", ha="center", va="center", fontsize=10, transform=ax.transAxes)
                ax.axis("off")
                continue

            draw_thin_stacked_proportion_barh(ax, data_by_group[key], class_order, class_colors, xmax=1.0)
            ax.set_title(f"{col_y}: {col_y_val}", fontsize=group_title_fontsize, fontweight="bold", pad=1)

        title_str = f"Grouped Track-Class Proportions — {group_label_title}"

    legend_handles, legend_labels = _class_legend_handles(list(reversed(class_order)), class_colors)
    fig.legend(
        handles=legend_handles, labels=legend_labels,
        loc="lower center", ncol=legend_ncol, frameon=False, fontsize=7,
    )
    fig.suptitle(title_str, y=0.99, fontsize=10, fontweight="bold", wrap=True)
    return fig


def _plot_page_grouped_class_proportions_flat_grid(
    data_by_group,
    *,
    class_order,
    class_colors,
    group_label_title,
    grid_ncols,
    panel_w_in=_GRID_PANEL_W_IN,
    group_title_fontsize=8,
):
    """One flat wrapped grid page (for 3+ group_cols): one panel per group-label combo.

    Figure size scales with the resulting grid shape (rather than being clipped to a
    fixed page), and each panel keeps a thin, consistent ``_GRID_PANEL_ASPECT`` (~1:6)
    height:width ratio regardless of how many group combinations exist.
    """
    groups = list(data_by_group.keys())
    ncols = max(1, int(grid_ncols))
    nrows = max(1, int(np.ceil(float(len(groups)) / float(ncols))))
    panel_h_in = panel_w_in * _GRID_PANEL_ASPECT
    legend_ncol, _, legend_margin_in = legend_layout(len(class_order), base_margin_in=_GRID_BOTTOM_MARGIN_IN)

    fig_w = ncols * panel_w_in
    fig_h = _GRID_TOP_MARGIN_IN + nrows * panel_h_in + legend_margin_in
    fig = plt.figure(figsize=(fig_w, fig_h))
    outer = GridSpec(
        nrows=nrows,
        ncols=ncols,
        figure=fig,
        hspace=0.6,
        wspace=0.3,
        top=1.0 - _GRID_TOP_MARGIN_IN / fig_h,
        bottom=legend_margin_in / fig_h,
    )

    for i, grp in enumerate(groups):
        ax = fig.add_subplot(outer[i])
        draw_thin_stacked_proportion_barh(ax, data_by_group[grp], class_order, class_colors, xmax=1.0)
        ax.set_title(grp, fontsize=group_title_fontsize, fontweight="bold", pad=1)

    for i in range(len(groups), nrows * ncols):
        fig.add_subplot(outer[i]).axis("off")

    legend_handles, legend_labels = _class_legend_handles(list(reversed(class_order)), class_colors)
    fig.legend(
        handles=legend_handles, labels=legend_labels,
        loc="lower center", ncol=legend_ncol, frameon=False, fontsize=7,
    )
    fig.suptitle(
        f"Grouped Track-Class Proportions — {group_label_title}",
        y=0.99, fontsize=10, fontweight="bold", wrap=True,
    )
    return fig


_STACK_GRID_PANEL_H_IN = 1.8         # taller than the thin-bar grid: room for a full multi-class bar chart
_STACK_GRID_BOTTOM_MARGIN_IN = 1.0   # extra room below panels for rotated x-tick class labels + legend
_STACK_GRID_MIN_PANEL_W_IN = 2.0


def _plot_page_class_stack_grid(
    data_by_facet,
    *,
    facet_cols,
    unique_vals_per_col,
    class_order,
    stack_order,
    colors,
    title,
    axis_cols=None,
    panel_w_in=None,
    group_title_fontsize=9,
):
    """One page: per-facet vertical stacked bar charts (all `class_order` classes on the x-axis,
    `stack_order` stacked on the y-axis) arranged in a 2D grid (or a single column for 1 facet
    column). Mirrors `_plot_page_grouped_class_proportions_2d_grid`'s row/column-header layout,
    but each grid cell holds a full multi-class bar chart instead of a single thin bar, so panel
    size scales with the number of classes rather than a fixed thin-bar aspect ratio.
    """
    panel_w_in = panel_w_in if panel_w_in is not None else max(0.45 * len(class_order), _STACK_GRID_MIN_PANEL_W_IN)
    panel_h_in = _STACK_GRID_PANEL_H_IN
    legend_ncol, _, legend_margin_in = legend_layout(len(stack_order), base_margin_in=_STACK_GRID_BOTTOM_MARGIN_IN)

    if len(facet_cols) == 2:
        if axis_cols is not None:
            col_y, col_x = axis_cols
        else:
            col_y, col_x = sorted(facet_cols, key=lambda c: -len(unique_vals_per_col[c]))
        col_x_vals = unique_vals_per_col[col_x]
        col_y_vals = unique_vals_per_col[col_y]

        nrows_data = max(1, len(col_y_vals))
        ncols_data = max(1, len(col_x_vals))

        fig_w = _GRID_HEADER_W_IN + ncols_data * panel_w_in
        fig_h = _GRID_TOP_MARGIN_IN + _GRID_HEADER_H_IN + nrows_data * panel_h_in + legend_margin_in
        fig = plt.figure(figsize=(fig_w, fig_h))
        outer = GridSpec(
            nrows=nrows_data + 1,
            ncols=ncols_data + 1,
            figure=fig,
            height_ratios=[_GRID_HEADER_H_IN] + [panel_h_in] * nrows_data,
            width_ratios=[_GRID_HEADER_W_IN] + [panel_w_in] * ncols_data,
            hspace=1.1,
            wspace=0.35,
            top=1.0 - _GRID_TOP_MARGIN_IN / fig_h,
            bottom=legend_margin_in / fig_h,
        )

        ax_corner = fig.add_subplot(outer[0, 0])
        ax_corner.axis("off")

        for c, col_x_val in enumerate(col_x_vals):
            ax_ch = fig.add_subplot(outer[0, c + 1])
            ax_ch.axis("off")
            ax_ch.text(
                0.5, 0.5, str(col_x_val),
                ha="center", va="center", fontsize=group_title_fontsize, fontweight="bold",
                transform=ax_ch.transAxes,
            )

        for r, col_y_val in enumerate(col_y_vals):
            ax_rh = fig.add_subplot(outer[r + 1, 0])
            ax_rh.axis("off")
            ax_rh.text(
                0.93, 0.5, _wrap_row_label(str(col_y_val), max_chars=_GRID_ROW_LABEL_MAX_CHARS),
                ha="right", va="center", ma="right", fontsize=group_title_fontsize, fontweight="bold",
                transform=ax_rh.transAxes,
            )

            for c, col_x_val in enumerate(col_x_vals):
                vals_map = {col_y: col_y_val, col_x: col_x_val}
                key = " | ".join(str(vals_map[fc]) for fc in facet_cols)

                ax = fig.add_subplot(outer[r + 1, c + 1])

                if key not in data_by_facet:
                    ax.text(0.5, 0.5, "—", ha="center", va="center", fontsize=10, transform=ax.transAxes)
                    ax.axis("off")
                    continue

                draw_stacked_proportion_barv(ax, data_by_facet[key], class_order, stack_order, colors, ymax=1.0)
                if c > 0:
                    ax.set_yticklabels([])

        title_str = f"{title} — {col_x} (columns) × {col_y} (rows)"

    else:
        col_y = facet_cols[0]
        col_y_vals = unique_vals_per_col[col_y]
        nrows_data = max(1, len(col_y_vals))

        fig_w = max(panel_w_in, _STACK_GRID_MIN_PANEL_W_IN)
        fig_h = _GRID_TOP_MARGIN_IN + nrows_data * panel_h_in + legend_margin_in
        fig = plt.figure(figsize=(fig_w, fig_h))
        outer = GridSpec(
            nrows=nrows_data,
            ncols=1,
            figure=fig,
            height_ratios=[panel_h_in] * nrows_data,
            hspace=1.1,
            top=1.0 - _GRID_TOP_MARGIN_IN / fig_h,
            bottom=legend_margin_in / fig_h,
        )

        for r, col_y_val in enumerate(col_y_vals):
            key = str(col_y_val)
            ax = fig.add_subplot(outer[r])

            if key not in data_by_facet:
                ax.text(0.5, 0.5, "—", ha="center", va="center", fontsize=10, transform=ax.transAxes)
                ax.axis("off")
                continue

            draw_stacked_proportion_barv(ax, data_by_facet[key], class_order, stack_order, colors, ymax=1.0)
            ax.set_title(f"{col_y}: {col_y_val}", fontsize=group_title_fontsize, fontweight="bold", pad=4)

        title_str = f"{title} — {col_y}"

    legend_handles, legend_labels = _class_legend_handles(list(reversed(stack_order)), colors)
    fig.legend(
        handles=legend_handles, labels=legend_labels,
        loc="lower center", ncol=legend_ncol, frameon=False, fontsize=7,
    )
    fig.suptitle(title_str, y=0.99, fontsize=10, fontweight="bold", wrap=True)
    return fig


def _plot_single_class_stack_panel(props_df, *, class_order, stack_order, colors, title):
    """A single vertical stacked-bar panel (no facet grid) — used when only page-pooling
    columns are selected (no group_x/group_y), so each page is one plain bar chart."""
    panel_w_in = max(0.45 * len(class_order), 4.0)
    legend_ncol, _, legend_margin_in = legend_layout(len(stack_order), base_margin_in=0.42)
    fig_h = 3.78 + legend_margin_in
    fig, ax = plt.subplots(figsize=(panel_w_in, fig_h))
    draw_stacked_proportion_barv(ax, props_df, class_order, stack_order, colors, ymax=1.0, xtick_fontsize=8)
    ax.set_ylabel("Proportion")
    legend_handles, legend_labels = _class_legend_handles(list(reversed(stack_order)), colors)
    fig.legend(
        handles=legend_handles, labels=legend_labels,
        loc="lower center", ncol=legend_ncol, frameon=False, fontsize=8,
    )
    fig.suptitle(title, fontsize=11, fontweight="bold", wrap=True)
    fig.tight_layout(rect=(0.0, legend_margin_in / fig_h, 1.0, 0.92))
    return fig


_CONTACT_STACK_ORDER = ["no_contact", "contact"]
_CONTACT_STACK_COLORS = {"no_contact": "#B0B0B0", "contact": "#D1495B"}


def _save_contact_cluster_stack_grid_page(
    pdf,
    adata_tracks,
    *,
    class_col,
    group_col,
    extra_group_cols,
    group_x,
    group_y,
    class_order,
    csv_dir,
    dpi=300,
    verbose=False,
):
    """Add the `class_col`-on-x-axis contact/no_contact stacked-bar grid page(s) to `pdf`,
    faceted by `group_x`/`group_y` as a true 2D grid; `extra_group_cols` ("group per page")
    paginate into separate grid pages instead of adding a third grid dimension.

    Skipped entirely (returns ``n_pages=0``) when none of `group_x`/`group_y`/`extra_group_cols`
    are set, matching the plain per-sample contact-rate page which is always produced separately.
    """
    obs_cols = list(adata_tracks.obs.columns)
    extra_group_cols = [c for c in (extra_group_cols or []) if c in obs_cols]
    grid_axis_cols = [c for c in (group_x, group_y) if c and c in obs_cols]

    if not grid_axis_cols and not extra_group_cols:
        if bool(verbose):
            print("  Note: no group_x/group_y/extra_group_cols set — skipping cluster-stack grid page.")
        return {"n_pages": 0, "csv_path": None}

    stack_order = list(_CONTACT_STACK_ORDER)
    stack_colors = _normalize_label_color_map(stack_order, colors=_CONTACT_STACK_COLORS, cmap_name="tab10")

    needed_cols = [class_col, group_col] + list(dict.fromkeys(grid_axis_cols + extra_group_cols))
    plot_df = adata_tracks.obs[needed_cols].copy()
    plot_df[class_col] = plot_df[class_col].astype("string").fillna("unassigned").astype(str)
    plot_df[group_col] = plot_df[group_col].astype("string").fillna("(unknown)").astype(str)
    for gc in dict.fromkeys(grid_axis_cols + extra_group_cols):
        plot_df[gc] = plot_df[gc].astype("string").fillna("(unknown)").astype(str)

    if extra_group_cols:
        page_labels = _make_group_label(plot_df, extra_group_cols)
        page_keys = sorted(page_labels.dropna().unique().tolist(), key=_mixed_label_sort_key)
    else:
        page_labels = pd.Series(["(all)"] * len(plot_df), index=plot_df.index)
        page_keys = ["(all)"]

    long_rows = []
    n_pages = 0
    for page_key in page_keys:
        page_df = plot_df[page_labels == page_key]
        if len(page_df) == 0:
            continue
        extra_vals = str(page_key).split(" | ") if extra_group_cols else []

        if grid_axis_cols:
            unique_vals_per_col = {
                c: sorted(page_df[c].dropna().unique().tolist(), key=_mixed_label_sort_key)
                for c in grid_axis_cols
            }
            facet_labels = _make_group_label(page_df, grid_axis_cols)
            data_by_facet = {}
            for facet_key in facet_labels.dropna().unique().tolist():
                facet_df = page_df[facet_labels == facet_key]
                props = compute_class_by_stack_proportions(
                    facet_df, class_col=class_col, stack_col=group_col,
                    class_order=class_order, stack_order=stack_order,
                )
                data_by_facet[facet_key] = props
                facet_vals = str(facet_key).split(" | ")
                for cls in class_order:
                    row = dict(zip(grid_axis_cols, facet_vals))
                    row.update(dict(zip(extra_group_cols, extra_vals)))
                    row[class_col] = cls
                    for stack_val in stack_order:
                        row[f"{group_col}={stack_val}"] = float(props.loc[cls, stack_val]) if cls in props.index else 0.0
                    long_rows.append(row)

            axis_cols = (
                (group_y, group_x)
                if (group_x and group_y and group_x in grid_axis_cols and group_y in grid_axis_cols)
                else None
            )
            title = f"Contact-Group Composition by {class_col}"
            if extra_group_cols:
                title = f"{title} — {', '.join(extra_group_cols)}: {page_key}"
            fig = _plot_page_class_stack_grid(
                data_by_facet,
                facet_cols=grid_axis_cols,
                unique_vals_per_col=unique_vals_per_col,
                class_order=class_order,
                stack_order=stack_order,
                colors=stack_colors,
                title=title,
                axis_cols=axis_cols,
            )
        else:
            props = compute_class_by_stack_proportions(
                page_df, class_col=class_col, stack_col=group_col,
                class_order=class_order, stack_order=stack_order,
            )
            for cls in class_order:
                row = dict(zip(extra_group_cols, extra_vals))
                row[class_col] = cls
                for stack_val in stack_order:
                    row[f"{group_col}={stack_val}"] = float(props.loc[cls, stack_val]) if cls in props.index else 0.0
                long_rows.append(row)

            title = f"Contact-Group Composition by {class_col}"
            if extra_group_cols:
                title = f"{title} — {', '.join(extra_group_cols)}: {page_key}"
            fig = _plot_single_class_stack_panel(
                props, class_order=class_order, stack_order=stack_order, colors=stack_colors, title=title,
            )

        pdf.savefig(fig, dpi=int(dpi))
        plt.close(fig)
        n_pages += 1

    csv_path = None
    if long_rows:
        class_token = _sanitize_filename_token(class_col, fallback="ClusterID")
        csv_path = Path(csv_dir) / f"contact_cluster_stack_grid_{class_token}.csv"
        pd.DataFrame(long_rows).to_csv(csv_path, index=False)

    return {"n_pages": n_pages, "csv_path": str(csv_path) if csv_path is not None else None}


def save_track_class_proportions_by_sample_plot(
    adata_tracks,
    out_dir,
    *,
    sample_col="sample_name",
    class_col="ClusterID",
    group_cols=None,
    group_x=None,
    group_y=None,
    grid_ncols=3,
    dpi=300,
    cmap_name="tab20",
    verbose=False,
    class_colors=None,
    pdf_pages=None,
    csv_dir=None,
):
    """
    Save one horizontal stacked bar per sample showing track-class proportions,
    plus optional grouped grid pages (1-2 effective group columns: true 2D grid;
    3+: flat grid). ``group_x``/``group_y`` explicitly pick the 2D grid's axes;
    ``group_cols`` is the "group per page" column list (unaffected in meaning).

    If ``pdf_pages`` (an open ``PdfPages``) is given, pages are appended to it instead of a
    standalone PDF being created. ``csv_dir`` overrides where CSVs are written (defaults to
    ``out_dir``).
    """
    if sample_col not in adata_tracks.obs.columns:
        raise ValueError(f"Missing '{sample_col}' in adata_tracks.obs.")
    if class_col not in adata_tracks.obs.columns:
        raise ValueError(f"Missing '{class_col}' in adata_tracks.obs.")

    effective_group_cols, requested_axis_cols = _resolve_effective_group_cols(
        group_cols, group_x, group_y,
    )

    valid_group_cols = []
    if effective_group_cols:
        obs_cols_available = list(adata_tracks.obs.columns)
        for gc in effective_group_cols:
            if gc in obs_cols_available:
                valid_group_cols.append(gc)
            elif verbose:
                print(f"  Warning: group_col '{gc}' not found in adata_tracks.obs — skipping.")

    axis_cols = (
        requested_axis_cols
        if requested_axis_cols is not None
        and requested_axis_cols[0] in valid_group_cols
        and requested_axis_cols[1] in valid_group_cols
        else None
    )

    plot_df = adata_tracks.obs[[sample_col, class_col] + valid_group_cols].copy()
    plot_df[sample_col] = (
        plot_df[sample_col]
        .astype("string")
        .fillna("unassigned")
        .astype(str)
        .str.strip()
        .replace("", "unassigned")
    )
    plot_df[class_col] = (
        plot_df[class_col]
        .astype("string")
        .fillna("unassigned")
        .astype(str)
        .str.strip()
        .replace("", "unassigned")
    )
    for gc in valid_group_cols:
        plot_df[gc] = plot_df[gc].astype("string").fillna("(unknown)").astype(str)

    if len(plot_df) == 0:
        raise ValueError("No rows available to plot track class proportions.")

    sample_order = (
        pd.Series(plot_df[sample_col], dtype="string")
        .dropna()
        .astype(str)
        .drop_duplicates()
        .tolist()
    )
    class_order = sorted(
        pd.Series(plot_df[class_col], dtype="string")
        .dropna()
        .astype(str)
        .unique()
        .tolist(),
        key=_mixed_label_sort_key,
    )
    class_order = _apply_state_order(class_order, _get_classification_state_order(adata_tracks, class_col))
    if len(sample_order) == 0:
        raise ValueError("No samples available to plot track class proportions.")
    if len(class_order) == 0:
        raise ValueError("No classes available to plot track class proportions.")

    proportion_table = pd.crosstab(
        plot_df[sample_col],
        plot_df[class_col],
        normalize="index",
        dropna=False,
    ).reindex(index=sample_order, columns=class_order, fill_value=0.0)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    class_token = _sanitize_filename_token(class_col, fallback="ClusterID")
    pdf_path = out_dir / f"track_class_proportions_by_sample_{class_token}.pdf"
    csv_dir = Path(csv_dir) if csv_dir is not None else out_dir
    csv_dir.mkdir(parents=True, exist_ok=True)
    csv_path = csv_dir / pdf_path.with_suffix(".csv").name
    proportion_table.to_csv(csv_path)

    resolved_colors = dict(class_colors) if class_colors else _get_classification_state_colors(adata_tracks, class_col)
    colors = _normalize_label_color_map(class_order, colors=resolved_colors, cmap_name=cmap_name)
    grouped_csv_path = None
    n_grouped_pages = 0

    with (nullcontext(pdf_pages) if pdf_pages is not None else PdfPages(pdf_path)) as pdf:
        rows_per_page = stacked_proportion_barh_rows_per_page()
        data_by_sample = {s: proportion_table.loc[s] for s in sample_order}
        for page_samples in _chunk_list(sample_order, rows_per_page):
            fig = plot_page_stacked_proportion_barh_grid(
                data_by_sample,
                row_order=page_samples,
                class_order=class_order,
                colors=colors,
                title=f"Track-Class Proportions by Sample — {class_col}",
            )
            pdf.savefig(fig, dpi=int(dpi))
            plt.close(fig)

        if valid_group_cols:
            proportions_by_group, counts_by_group, _n_by_group = _compute_grouped_class_proportions(
                plot_df, group_cols=valid_group_cols, class_col=class_col, class_order=class_order,
            )
            group_label_title = ", ".join(valid_group_cols)
            unique_vals_per_col = {}
            if len(valid_group_cols) in (1, 2):
                for col in valid_group_cols:
                    unique_vals_per_col[col] = sorted(plot_df[col].astype(str).dropna().unique().tolist())

            if len(valid_group_cols) in (1, 2):
                fig_rel = _plot_page_grouped_class_proportions_2d_grid(
                    proportions_by_group,
                    group_cols=valid_group_cols,
                    unique_vals_per_col=unique_vals_per_col,
                    class_order=class_order,
                    class_colors=colors,
                    group_label_title=group_label_title,
                    axis_cols=axis_cols,
                )
                pdf.savefig(fig_rel, dpi=dpi, bbox_inches="tight")
                plt.close(fig_rel)
                n_grouped_pages += 1
            else:
                if verbose:
                    print(
                        f"  Note: {len(valid_group_cols)} group_cols given — "
                        "using a flat grid (not a true 2D layout)."
                    )
                fig_rel = _plot_page_grouped_class_proportions_flat_grid(
                    proportions_by_group,
                    class_order=class_order,
                    class_colors=colors,
                    group_label_title=group_label_title,
                    grid_ncols=grid_ncols,
                )
                pdf.savefig(fig_rel, dpi=dpi)
                plt.close(fig_rel)
                n_grouped_pages += 1

            grouped_long_rows = []
            for grp_label, prop_series in proportions_by_group.items():
                cnt_series = counts_by_group[grp_label]
                group_vals = grp_label.split(" | ")
                for class_name in class_order:
                    row = dict(zip(valid_group_cols, group_vals))
                    row["class"] = class_name
                    row["proportion"] = float(prop_series.get(class_name, 0.0))
                    row["n_tracks"] = int(cnt_series.get(class_name, 0))
                    grouped_long_rows.append(row)
            grouped_csv_path = csv_dir / f"track_class_proportions_by_group_{class_token}.csv"
            pd.DataFrame(grouped_long_rows).to_csv(grouped_csv_path, index=False)

    color_hex = {str(k): str(to_hex(v)) for k, v in colors.items()}
    return {
        "pdf_path": str(pdf_path),
        "csv_path": str(csv_path),
        "sample_order": [str(x) for x in sample_order],
        "class_order": [str(x) for x in class_order],
        "colors": color_hex,
        "group_cols": valid_group_cols,
        "grouped_csv_path": str(grouped_csv_path) if grouped_csv_path is not None else None,
        "n_grouped_pages": n_grouped_pages,
    }


def save_track_condition_comparison_report(
    adata_tracks,
    out_dir,
    *,
    sample_col="sample_name",
    class_col="ClusterID",
    condition_col,
    group_cols=None,
    group_x=None,
    group_y=None,
    class_order=None,
    class_colors=None,
    verbose=False,
    pdf_pages=None,
    csv_dir=None,
):
    """Per-cluster overall (pooled) track-class proportion difference between every pairwise
    combination of `condition_col`'s levels, with Welch's two-sided unpaired t-test significance
    stars — one row per pairwise comparison, one column per group, optionally split into
    multiple groups (columns) via `group_x`/`group_cols`.

    diff = mean of the second level minus mean of the first level in each pair.

    When `condition_col` is binary (exactly 2 levels — e.g. "contact"/"no_contact") and both
    `group_x` and `group_y` are given, the report switches to a true 2D grid (columns =
    `group_x`, rows = `group_y`) instead of pooling everything into columns, since the pairwise-
    comparison row axis collapses to a single row anyway in that case; `group_cols` still
    paginates into separate grid pages. Otherwise (non-binary `condition_col`, or only one of
    `group_x`/`group_y` set), behavior is unchanged: `group_x`/`group_cols` pool into columns and
    the pairwise comparisons are the rows.

    If ``pdf_pages`` (an open ``PdfPages``) is given, the figure is appended to it instead of
    being written to its own PDF file. ``csv_dir`` overrides where the CSV is written (defaults
    to ``out_dir``).
    """
    obs = adata_tracks.obs
    extra_group_cols = [c for c in (group_cols or []) if c in obs.columns]
    grid_axis_cols = [c for c in (group_x, group_y) if c and c in obs.columns]
    valid_group_cols = list(dict.fromkeys(extra_group_cols + grid_axis_cols))
    required = [sample_col, class_col, condition_col] + valid_group_cols
    missing = [c for c in required if c not in obs.columns]
    if len(missing) > 0:
        raise KeyError(f"Missing required columns in adata_tracks.obs: {missing}")

    df = obs[required].dropna(subset=[sample_col, class_col, condition_col]).copy()
    df[sample_col] = df[sample_col].astype(str)
    df[class_col] = df[class_col].astype(str)
    df[condition_col] = df[condition_col].astype(str)
    for gc in valid_group_cols:
        df[gc] = df[gc].astype(str)
    if len(df) == 0:
        raise ValueError("No valid rows remain after filtering NaNs in required columns.")

    if class_order is not None:
        resolved_class_order = [str(c) for c in class_order]
    else:
        resolved_class_order = sorted(df[class_col].dropna().unique().tolist(), key=_mixed_label_sort_key)
    resolved_class_order = _apply_state_order(
        resolved_class_order, _get_classification_state_order(adata_tracks, class_col)
    )

    # The unit compared below is (sample, condition_col level), not the sample alone.
    # condition_col is usually constant within a sample (e.g. a per-well treatment), in which
    # case this is a no-op relabeling. But for a contact/no_contact style grouping, condition_col
    # is a per-track tag that varies *within* a sample. Collapsing straight to one row per sample
    # would arbitrarily keep whichever level happened to appear first and silently leave a single
    # degenerate condition level with zero pairs left to compare (empty report). Using this
    # composite unit lets a sample that contains multiple condition levels contribute one row per
    # level. Mirrors the same fix in state_composition.save_state_condition_comparison_report.
    unit_col = "__comparison_unit__"
    df[unit_col] = df[sample_col] + " | " + df[condition_col]

    proportions_by_unit, _, _ = _compute_grouped_class_proportions(
        df, group_cols=[unit_col], class_col=class_col, class_order=resolved_class_order,
    )
    per_sample_df = pd.DataFrame(proportions_by_unit).T.reindex(columns=resolved_class_order, fill_value=0.0)

    metadata_cols = [unit_col, condition_col] + valid_group_cols
    sample_metadata = df[metadata_cols].drop_duplicates(subset=[unit_col]).set_index(unit_col)

    resolved_colors = dict(class_colors) if class_colors else _get_classification_state_colors(adata_tracks, class_col)
    resolved_colors = _normalize_label_color_map(resolved_class_order, colors=resolved_colors, cmap_name="tab20")

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    condition_token = _sanitize_filename_token(condition_col, fallback="condition")
    out_pdf = out_dir / f"condition_comparison_{condition_token}.pdf"
    csv_dir = Path(csv_dir) if csv_dir is not None else out_dir
    csv_dir.mkdir(parents=True, exist_ok=True)
    out_csv = csv_dir / out_pdf.with_suffix(".csv").name

    title = f"{class_col} — {condition_col} pairwise comparison"
    n_condition_levels = sample_metadata[condition_col].astype(str).nunique()
    use_2d_grid = bool(group_x) and bool(group_y) and n_condition_levels == 2

    if use_2d_grid:
        result = plot_condition_diff_grid_2d(
            per_sample_df,
            sample_metadata,
            class_order=resolved_class_order,
            colors=resolved_colors,
            condition_col=condition_col,
            grid_axis_cols=grid_axis_cols,
            extra_group_cols=extra_group_cols,
            title=title,
            out_pdf=out_pdf,
            out_csv=out_csv,
            pdf_pages=pdf_pages,
        )
    else:
        diff_stats_by_group = compute_condition_diff_stats_pairwise(
            per_sample_df,
            sample_metadata,
            class_order=resolved_class_order,
            condition_col=condition_col,
            group_cols=valid_group_cols,
        )
        result = plot_condition_diff_grid(
            diff_stats_by_group,
            class_order=resolved_class_order,
            colors=resolved_colors,
            title=title,
            out_pdf=out_pdf,
            out_csv=out_csv,
            pdf_pages=pdf_pages,
        )
    if bool(verbose):
        print(f"Saved track condition comparison report: {result['pdf_path']}")
    return result


def save_track_contact_group_analysis(
    adata_tracks,
    df_timepoints,
    out_dir,
    *,
    contact_col,
    min_bout_length,
    sample_col="sample_name",
    class_col="ClusterID",
    class_order=None,
    class_colors=None,
    extra_group_cols=None,
    group_x=None,
    group_y=None,
    verbose=False,
):
    """Group classified tracks by whether they had a sufficiently long contiguous bout of
    ``contact_col`` (>= ``min_bout_length`` timepoints, within each track's classified time
    window) versus not, and produce the full comparison bundle:

    1. Contact-rate proportions per sample (one thin bar per sample).
    2. A `class_col`-on-x-axis grid of contact/no_contact composition, faceted by `group_x`/
       `group_y` (true 2D grid) with `extra_group_cols` paginating into separate grid pages;
       skipped if none of `group_x`/`group_y`/`extra_group_cols` are set.
    3. A condition-comparison report (Welch's t-test) between the two contact groups — a true
       `group_x` × `group_y` 2D grid (since the comparison is always binary, both axes are free
       to be real grid axes); `extra_group_cols` paginate into separate grid pages.
    4. Violin plots of mean contact fraction and max contact-bout length per cluster.

    Returns a dict of all output artifact paths.
    """
    extra_group_cols = list(extra_group_cols) if extra_group_cols else []
    group_col = _contact_group_col_name(contact_col)
    mean_col = _contact_mean_col_name(contact_col)
    max_bout_col = _contact_max_bout_col_name(contact_col)

    contact_features = compute_track_contact_features(
        df_timepoints,
        adata_tracks,
        contact_col=contact_col,
        min_bout_length=min_bout_length,
        verbose=verbose,
    )
    merge_track_contact_features_into_obs(
        adata_tracks,
        contact_features,
        contact_col=contact_col,
        min_bout_length=min_bout_length,
    )
    if bool(verbose):
        n_contact = int((adata_tracks.obs[group_col] == "contact").sum())
        n_total = len(adata_tracks.obs)
        print(f"Contact grouping on '{contact_col}' (min_bout_length={min_bout_length}): {n_contact}/{n_total} tracks labeled 'contact'.")

    # Own subfolder, named after the raw contact column (e.g. "macrophage_contact"), holding a
    # single combined report PDF and a "csv" subfolder for all underlying CSVs.
    out_dir = Path(out_dir) / "contact_analysis" / str(contact_col)
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_dir = out_dir / "csv"
    csv_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / "contact_analysis.pdf"

    obs = adata_tracks.obs
    resolved_class_order = (
        [str(c) for c in class_order] if class_order is not None
        else sorted(obs[class_col].dropna().astype(str).unique().tolist(), key=_mixed_label_sort_key)
    )
    resolved_class_order = _apply_state_order(resolved_class_order, _get_classification_state_order(adata_tracks, class_col))
    resolved_colors = dict(class_colors) if class_colors else _get_classification_state_colors(adata_tracks, class_col)
    resolved_colors = _normalize_label_color_map(resolved_class_order, colors=resolved_colors, cmap_name="tab20")

    with PdfPages(pdf_path) as pdf:
        contact_rate_by_cluster = save_track_class_proportions_by_sample_plot(
            adata_tracks,
            out_dir,
            sample_col=sample_col,
            class_col=group_col,
            group_cols=None,
            verbose=verbose,
            pdf_pages=pdf,
            csv_dir=csv_dir,
        )

        cluster_stack_grid = _save_contact_cluster_stack_grid_page(
            pdf,
            adata_tracks,
            class_col=class_col,
            group_col=group_col,
            extra_group_cols=extra_group_cols,
            group_x=group_x,
            group_y=group_y,
            class_order=resolved_class_order,
            csv_dir=csv_dir,
            verbose=verbose,
        )

        condition_comparison = save_track_condition_comparison_report(
            adata_tracks,
            out_dir,
            sample_col=sample_col,
            class_col=class_col,
            condition_col=group_col,
            group_x=group_x,
            group_y=group_y,
            group_cols=list(extra_group_cols),
            class_order=resolved_class_order,
            class_colors=class_colors,
            verbose=verbose,
            pdf_pages=pdf,
            csv_dir=csv_dir,
        )

        fig = plot_feature_violin_by_group(
            obs, mean_col, class_col,
            group_order=resolved_class_order, colors=resolved_colors,
            ylabel=f"Mean fraction of timepoints in contact ({contact_col})",
            title=f"Mean contact fraction by {class_col}",
        )
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        fig = plot_feature_violin_by_group(
            obs, max_bout_col, class_col,
            group_order=resolved_class_order, colors=resolved_colors,
            ylabel=f"Longest contiguous contact bout, timepoints ({contact_col})",
            title=f"Max contact-bout length by {class_col}",
        )
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    # The per-plot helpers above report their own nominal (standalone) pdf paths; since they all
    # actually wrote into the single combined `pdf_path` here, point their "pdf_path" entries at
    # that combined file instead so downstream consumers reference something that exists on disk.
    contact_rate_by_cluster["pdf_path"] = str(pdf_path)
    condition_comparison["pdf_path"] = str(pdf_path)

    if bool(verbose):
        print(f"Saved contact-group analysis bundle for '{contact_col}' to: {pdf_path}")

    return {
        "contact_col": str(contact_col),
        "min_bout_length": int(min_bout_length),
        "group_col": group_col,
        "extra_group_cols": extra_group_cols,
        "pdf_path": str(pdf_path),
        "csv_dir": str(csv_dir),
        "contact_rate_by_cluster": contact_rate_by_cluster,
        "cluster_stack_grid": cluster_stack_grid,
        "condition_comparison": condition_comparison,
        "mean_fraction_violin_pdf": str(pdf_path),
        "max_bout_length_violin_pdf": str(pdf_path),
    }


def generate_track_clustering_report_pdfs(
    adata_tracks,
    outfolder,
    cluster_key="ClusterID",
    heatmap_figsize=(25, 20),
    matrixplot_figsize=(20, 40),
    umap_size=1,
    umap_alpha=0.5,
    plot_dpi=300,
    filename_suffix="",
    verbose=False,
):
    """
    Write track clustering diagnostics as one multi-page PDF:
    - page 1: UMAP + correlation matrix (side-by-side)
    - page 2: heatmap
    - page 3: matrixplot
    - page 4: cluster occurrence ranking
    """
    outfolder = Path(outfolder)
    outfolder.mkdir(parents=True, exist_ok=True)

    if cluster_key not in adata_tracks.obs.columns:
        raise ValueError(f"Missing '{cluster_key}' in adata_tracks.obs.")

    resolved_order = _apply_state_order(
        sorted(adata_tracks.obs[cluster_key].astype(str).unique().tolist(), key=_mixed_label_sort_key),
        _get_classification_state_order(adata_tracks, cluster_key),
    )
    adata_tracks.obs[cluster_key] = pd.Categorical(
        adata_tracks.obs[cluster_key].astype(str), categories=resolved_order,
    )
    resolved_colors = _normalize_label_color_map(
        resolved_order, colors=_get_classification_state_colors(adata_tracks, cluster_key),
    )
    adata_tracks.uns[f"{cluster_key}_colors"] = [resolved_colors[c] for c in resolved_order]

    if "X_umap" not in adata_tracks.obsm:
        sc.tl.umap(adata_tracks, random_state=0)

    sc.tl.dendrogram(adata_tracks, groupby=cluster_key)

    suffix = str(filename_suffix).strip()
    diagnostics_pdf = Path(
        outfolder,
        f"clustering_diagnostics_{_sanitize_filename_token(cluster_key)}{suffix}.pdf",
    )
    diagnostics_started = _vstart(verbose, "trajectory-clustering", "write clustering diagnostics PDF")

    with PdfPages(diagnostics_pdf) as pdf:
        fig, (ax_umap, ax_corr) = plt.subplots(
            1,
            2,
            figsize=(14.5, 8.5),
            gridspec_kw={"width_ratios": [1.25, 0.75]},
        )
        sc.pl.umap(
            adata_tracks,
            color=cluster_key,
            size=umap_size,
            alpha=umap_alpha,
            ax=ax_umap,
            show=False,
            title="UMAP",
        )
        ax_umap.set_aspect("equal", adjustable="box")

        sc.pl.correlation_matrix(
            adata_tracks,
            groupby=cluster_key,
            show_correlation_numbers=True,
            ax=ax_corr,
            show=False,
        )
        ax_corr.set_title("Cluster correlations", fontsize=10)
        fig.suptitle("Track clustering diagnostics | UMAP + Correlation", fontsize=12, fontweight="bold", y=0.995)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        umap_bbox = ax_umap.get_position().frozen()
        bbox = ax_corr.get_position()
        extra_gap = 0.03
        new_x = min(0.98 - bbox.width, bbox.x0 + extra_gap)
        height_shrink = 0.5
        new_h = bbox.height * height_shrink
        ax_corr.set_position([new_x, bbox.y0 + (bbox.height - new_h) / 2, bbox.width, new_h])
        ax_umap.set_position(umap_bbox)
        pdf.savefig(fig, dpi=int(plot_dpi), bbox_inches="tight")
        plt.close(fig)

        sc.pl.heatmap(
            adata_tracks,
            var_names=adata_tracks.var_names,
            groupby=cluster_key,
            standard_scale="var",
            figsize=heatmap_figsize,
            swap_axes=True,
            dendrogram=True,
            show_gene_labels=True,
            show=False,
        )
        fig = plt.gcf()
        default_orientation = "landscape" if float(heatmap_figsize[0]) >= float(heatmap_figsize[1]) else "portrait"
        _apply_best_pdf_orientation(fig, default_orientation=default_orientation)
        pdf.savefig(fig, dpi=int(plot_dpi), bbox_inches="tight")
        plt.close(fig)

        sc.pl.matrixplot(
            adata_tracks,
            var_names=adata_tracks.var_names,
            groupby=cluster_key,
            standard_scale="var",
            figsize=matrixplot_figsize,
            swap_axes=True,
            dendrogram=True,
            show=False,
        )
        fig = plt.gcf()
        default_orientation = "landscape" if float(matrixplot_figsize[0]) >= float(matrixplot_figsize[1]) else "portrait"
        _apply_best_pdf_orientation(fig, default_orientation=default_orientation)
        pdf.savefig(fig, dpi=int(plot_dpi), bbox_inches="tight")
        plt.close(fig)

        fig = _plot_cluster_rankings(
            adata_tracks=adata_tracks,
            cluster_key=cluster_key,
            cluster_order=resolved_order,
            cluster_colors=resolved_colors,
        )
        _apply_best_pdf_orientation(fig, default_orientation="landscape")
        pdf.savefig(fig, dpi=int(plot_dpi), bbox_inches="tight")
        plt.close(fig)

    _vsave(verbose, "trajectory-clustering", "clustering diagnostics PDF", diagnostics_pdf)
    _vdone(verbose, "trajectory-clustering", "write clustering diagnostics PDF", diagnostics_started)
    return {
        "diagnostics_pdf": diagnostics_pdf,
    }
