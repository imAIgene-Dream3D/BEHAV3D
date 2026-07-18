import time
import random
import pickle
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
    compute_condition_diff_stats,
    plot_condition_diff_composite,
    plot_page_stacked_proportion_barh_grid,
    stacked_proportion_barh_rows_per_page,
    _chunk_list,
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


A4_PORTRAIT = (8.27, 11.69)
_A4_CONTENT_H = 9.0    # usable height (in) after title + legend + margins
_A4_CONTENT_W = 7.2    # usable width (in) after side margins
_MIN_PANEL_H = 2.2     # minimum comfortable panel height -> max 4 rows on A4
_MIN_PANEL_W = 2.0     # minimum comfortable panel width  -> max 3 cols on A4


def _panels_per_a4_page(grid_ncols):
    """Return (panels_per_page, ncols, max_rows) fitting comfortably on A4."""
    ncols = max(1, min(int(grid_ncols), int(_A4_CONTENT_W / _MIN_PANEL_W)))
    max_rows = max(1, int(_A4_CONTENT_H / _MIN_PANEL_H))
    return ncols * max_rows, ncols, max_rows


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


def _rank_cluster_counts(adata_tracks, cluster_key):
    if cluster_key not in adata_tracks.obs.columns:
        raise ValueError(f"Missing '{cluster_key}' in adata_tracks.obs.")

    labels = (
        pd.Series(adata_tracks.obs[cluster_key])
        .astype("string")
        .dropna()
        .astype(str)
    )
    counts = labels.value_counts(dropna=False).to_dict()
    ranked = sorted(
        [(str(label), int(n)) for label, n in counts.items()],
        key=lambda x: (-int(x[1]), _mixed_label_sort_key(x[0])),
    )
    return ranked


def _plot_cluster_rankings(adata_tracks, cluster_key):
    ranked = _rank_cluster_counts(adata_tracks=adata_tracks, cluster_key=cluster_key)
    cluster_labels = [x[0] for x in ranked]
    cluster_counts = [x[1] for x in ranked]

    fig_h = max(4.5, 0.45 * len(cluster_labels) + 2.0)
    fig, ax = plt.subplots(figsize=(10.5, fig_h))
    y = np.arange(len(cluster_labels))
    bars = ax.barh(y, cluster_counts, color="#2E6FBA")
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




def _make_group_label(df, group_cols):
    """Return a Series of concatenated group-column values (index aligned with df)."""
    parts = []
    for col in group_cols:
        if col in df.columns:
            parts.append(df[col].fillna("(unknown)").astype(str))
        else:
            parts.append(pd.Series(["(unknown)"] * len(df), index=df.index))
    if len(parts) == 0:
        return pd.Series(["(all)"] * len(df), index=df.index)
    result = parts[0]
    for p in parts[1:]:
        result = result + " | " + p
    return result


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
    row_slice=None,
    page_size=A4_PORTRAIT,
    group_title_fontsize=9,
):
    """
    One page: grouped track-class proportions arranged in a 2D grid.

    For 1 group_col: single column of panels, one row per unique value.
    For 2 group_cols: the column with more unique values goes to rows (Y),
    the column with fewer unique values goes to columns (X), with header
    labels on each axis. row_slice=(start, end) selects a subset of Y rows
    for pagination.
    """
    if len(group_cols) == 2:
        col_y, col_x = sorted(group_cols, key=lambda c: -len(unique_vals_per_col[c]))
        col_x_vals = unique_vals_per_col[col_x]
        col_y_vals_all = unique_vals_per_col[col_y]
        page_col_y_vals = col_y_vals_all[row_slice[0] : row_slice[1]] if row_slice is not None else col_y_vals_all

        nrows_data = max(1, len(page_col_y_vals))
        ncols_data = max(1, len(col_x_vals))

        header_h = 0.06
        header_w = 0.06
        fig = plt.figure(figsize=page_size)
        outer = GridSpec(
            nrows=nrows_data + 1,
            ncols=ncols_data + 1,
            figure=fig,
            height_ratios=[header_h] + [1.0] * nrows_data,
            width_ratios=[header_w] + [1.0] * ncols_data,
            hspace=0.55,
            wspace=0.3,
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

        for r, col_y_val in enumerate(page_col_y_vals):
            ax_rh = fig.add_subplot(outer[r + 1, 0])
            ax_rh.axis("off")
            ax_rh.text(
                0.5, 0.5, str(col_y_val),
                ha="center", va="center", fontsize=group_title_fontsize, fontweight="bold",
                rotation=90, transform=ax_rh.transAxes,
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
        if row_slice is not None and row_slice[0] > 0:
            title_str += f"\n(rows {row_slice[0] + 1}–{row_slice[1]} of {len(col_y_vals_all)})"

    else:
        col_y = group_cols[0]
        col_y_vals_all = unique_vals_per_col[col_y]
        page_col_y_vals = col_y_vals_all[row_slice[0] : row_slice[1]] if row_slice is not None else col_y_vals_all

        nrows_data = max(1, len(page_col_y_vals))
        fig = plt.figure(figsize=page_size)
        outer = GridSpec(nrows=nrows_data, ncols=1, figure=fig, hspace=0.35)

        for r, col_y_val in enumerate(page_col_y_vals):
            key = str(col_y_val)
            ax = fig.add_subplot(outer[r])

            if key not in data_by_group:
                ax.text(0.5, 0.5, "—", ha="center", va="center", fontsize=10, transform=ax.transAxes)
                ax.axis("off")
                continue

            draw_thin_stacked_proportion_barh(ax, data_by_group[key], class_order, class_colors, xmax=1.0)
            ax.set_title(f"{col_y}: {col_y_val}", fontsize=group_title_fontsize, fontweight="bold", pad=1)

        title_str = f"Grouped Track-Class Proportions — {group_label_title}"
        if row_slice is not None and row_slice[0] > 0:
            title_str += f"\n(rows {row_slice[0] + 1}–{row_slice[1]} of {len(col_y_vals_all)})"

    legend_handles, legend_labels = _class_legend_handles(class_order, class_colors)
    fig.legend(
        handles=legend_handles, labels=legend_labels,
        loc="lower center", ncol=min(len(legend_labels), 8), frameon=False, fontsize=7,
    )
    fig.suptitle(title_str, y=0.97, fontsize=10, fontweight="bold", wrap=True)
    fig.tight_layout(rect=(0.02, 0.06, 0.98, 0.93))
    return fig


def _plot_page_grouped_class_proportions_flat_grid(
    data_by_group,
    *,
    class_order,
    class_colors,
    group_label_title,
    grid_ncols,
    group_title_fontsize=8,
):
    """One flat wrapped grid page (for 3+ group_cols): one panel per group-label combo."""
    groups = list(data_by_group.keys())
    ncols = max(1, int(grid_ncols))
    nrows = max(1, int(np.ceil(float(len(groups)) / float(ncols))))

    fig = plt.figure(figsize=A4_PORTRAIT)
    outer = GridSpec(nrows=nrows, ncols=ncols, figure=fig, hspace=0.35, wspace=0.3)

    for i, grp in enumerate(groups):
        ax = fig.add_subplot(outer[i])
        draw_thin_stacked_proportion_barh(ax, data_by_group[grp], class_order, class_colors, xmax=1.0)
        ax.set_title(grp, fontsize=group_title_fontsize, fontweight="bold", pad=1)

    for i in range(len(groups), nrows * ncols):
        fig.add_subplot(outer[i]).axis("off")

    legend_handles, legend_labels = _class_legend_handles(class_order, class_colors)
    fig.legend(
        handles=legend_handles, labels=legend_labels,
        loc="lower center", ncol=min(len(legend_labels), 8), frameon=False, fontsize=7,
    )
    fig.suptitle(
        f"Grouped Track-Class Proportions — {group_label_title}",
        y=0.97, fontsize=10, fontweight="bold", wrap=True,
    )
    fig.tight_layout(rect=(0.02, 0.07, 0.98, 0.93))
    return fig


def save_track_class_proportions_by_sample_plot(
    adata_tracks,
    out_dir,
    *,
    sample_col="sample_name",
    class_col="ClusterID",
    group_cols=None,
    grid_ncols=3,
    dpi=300,
    cmap_name="tab20",
    verbose=False,
    class_colors=None,
):
    """
    Save one horizontal stacked bar per sample showing track-class proportions,
    plus optional grouped grid pages (1-2 group_cols: true 2D grid; 3+: flat grid).
    """
    if sample_col not in adata_tracks.obs.columns:
        raise ValueError(f"Missing '{sample_col}' in adata_tracks.obs.")
    if class_col not in adata_tracks.obs.columns:
        raise ValueError(f"Missing '{class_col}' in adata_tracks.obs.")

    valid_group_cols = []
    if group_cols:
        obs_cols_available = list(adata_tracks.obs.columns)
        for gc in group_cols:
            if gc in obs_cols_available:
                valid_group_cols.append(gc)
            elif verbose:
                print(f"  Warning: group_col '{gc}' not found in adata_tracks.obs — skipping.")

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
    csv_path = pdf_path.with_suffix(".csv")
    proportion_table.to_csv(csv_path)

    resolved_colors = dict(class_colors) if class_colors else _get_classification_state_colors(adata_tracks, class_col)
    colors = _normalize_label_color_map(class_order, colors=resolved_colors, cmap_name=cmap_name)
    grouped_csv_path = None
    n_grouped_pages = 0

    with PdfPages(pdf_path) as pdf:
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

            _, ncols_eff, max_rows = _panels_per_a4_page(grid_ncols)

            if len(valid_group_cols) in (1, 2):
                col_y_vals = (
                    unique_vals_per_col[max(valid_group_cols, key=lambda c: len(unique_vals_per_col[c]))]
                    if len(valid_group_cols) == 2
                    else unique_vals_per_col[valid_group_cols[0]]
                )
                for row_start in range(0, max(1, len(col_y_vals)), max_rows):
                    row_end = min(row_start + max_rows, len(col_y_vals))
                    _rs = (row_start, row_end)

                    fig_rel = _plot_page_grouped_class_proportions_2d_grid(
                        proportions_by_group,
                        group_cols=valid_group_cols,
                        unique_vals_per_col=unique_vals_per_col,
                        class_order=class_order,
                        class_colors=colors,
                        group_label_title=group_label_title,
                        row_slice=_rs,
                    )
                    pdf.savefig(fig_rel, dpi=dpi)
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
                    grid_ncols=ncols_eff,
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
            grouped_csv_path = pdf_path.parent / f"track_class_proportions_by_group_{class_token}.csv"
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
    level_a,
    level_b,
    facet_col=None,
    class_order=None,
    class_colors=None,
    verbose=False,
):
    """Per-cluster overall (pooled) track-class proportion difference between two condition
    levels, with Welch's two-sided unpaired t-test significance stars — optionally faceted
    into side-by-side panels by the unique values of `facet_col`.

    diff = mean(level_b) - mean(level_a), computed across per-sample proportions.
    """
    obs = adata_tracks.obs
    required = [sample_col, class_col, condition_col] + ([facet_col] if facet_col else [])
    missing = [c for c in required if c not in obs.columns]
    if len(missing) > 0:
        raise KeyError(f"Missing required columns in adata_tracks.obs: {missing}")

    df = obs[required].dropna(subset=[sample_col, class_col, condition_col]).copy()
    df[sample_col] = df[sample_col].astype(str)
    df[class_col] = df[class_col].astype(str)
    df[condition_col] = df[condition_col].astype(str)
    if facet_col is not None:
        df[facet_col] = df[facet_col].astype(str)
    if len(df) == 0:
        raise ValueError("No valid rows remain after filtering NaNs in required columns.")

    if class_order is not None:
        resolved_class_order = [str(c) for c in class_order]
    else:
        resolved_class_order = sorted(df[class_col].dropna().unique().tolist(), key=_mixed_label_sort_key)
        resolved_class_order = _apply_state_order(
            resolved_class_order, _get_classification_state_order(adata_tracks, class_col)
        )

    proportions_by_sample, _, _ = _compute_grouped_class_proportions(
        df, group_cols=[sample_col], class_col=class_col, class_order=resolved_class_order,
    )
    per_sample_df = pd.DataFrame(proportions_by_sample).T.reindex(columns=resolved_class_order, fill_value=0.0)

    metadata_cols = [sample_col, condition_col] + ([facet_col] if facet_col is not None else [])
    sample_metadata = df[metadata_cols].drop_duplicates(subset=[sample_col]).set_index(sample_col)

    resolved_colors = dict(class_colors) if class_colors else _get_classification_state_colors(adata_tracks, class_col)
    resolved_colors = _normalize_label_color_map(resolved_class_order, colors=resolved_colors, cmap_name="tab20")

    diff_stats_by_facet = compute_condition_diff_stats(
        per_sample_df,
        sample_metadata,
        class_order=resolved_class_order,
        condition_col=condition_col,
        level_a=level_a,
        level_b=level_b,
        facet_col=facet_col,
    )

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    condition_token = _sanitize_filename_token(condition_col, fallback="condition")
    level_a_token = _sanitize_filename_token(str(level_a), fallback="a")
    level_b_token = _sanitize_filename_token(str(level_b), fallback="b")
    out_pdf = out_dir / f"condition_comparison_{condition_token}_{level_a_token}_vs_{level_b_token}.pdf"
    out_csv = out_pdf.with_suffix(".csv")

    title = f"{class_col}: {level_b} vs {level_a} ({condition_col})"
    result = plot_condition_diff_composite(
        diff_stats_by_facet,
        class_order=resolved_class_order,
        colors=resolved_colors,
        level_a_label=str(level_a),
        level_b_label=str(level_b),
        title=title,
        out_pdf=out_pdf,
        out_csv=out_csv,
    )
    if bool(verbose):
        print(f"Saved track condition comparison report: {result['pdf_path']}")
    return result


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

        fig = _plot_cluster_rankings(adata_tracks=adata_tracks, cluster_key=cluster_key)
        _apply_best_pdf_orientation(fig, default_orientation="landscape")
        pdf.savefig(fig, dpi=int(plot_dpi), bbox_inches="tight")
        plt.close(fig)

    _vsave(verbose, "trajectory-clustering", "clustering diagnostics PDF", diagnostics_pdf)
    _vdone(verbose, "trajectory-clustering", "write clustering diagnostics PDF", diagnostics_started)
    return {
        "diagnostics_pdf": diagnostics_pdf,
    }
