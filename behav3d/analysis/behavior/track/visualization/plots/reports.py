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


def _build_label_color_map(labels, cmap_name="tab20"):
    """Create deterministic label->color mapping for stacked plots."""
    labels = [str(x) for x in list(labels)]
    if len(labels) == 0:
        return {}
    cmap = plt.get_cmap(cmap_name)
    if len(labels) <= getattr(cmap, "N", 256):
        denom = max(len(labels) - 1, 1)
        color_values = [cmap(i / denom) for i in range(len(labels))]
    else:
        hsv = plt.get_cmap("hsv")
        color_values = [hsv(i / max(len(labels), 1)) for i in range(len(labels))]
    return {lab: color_values[i] for i, lab in enumerate(labels)}


def save_track_class_proportions_by_sample_plot(
    adata_tracks,
    out_dir,
    *,
    sample_col="sample_name",
    class_col="ClusterID",
    dpi=300,
    cmap_name="tab20",
):
    """
    Save one horizontal stacked bar per sample showing track-class proportions.
    """
    if sample_col not in adata_tracks.obs.columns:
        raise ValueError(f"Missing '{sample_col}' in adata_tracks.obs.")
    if class_col not in adata_tracks.obs.columns:
        raise ValueError(f"Missing '{class_col}' in adata_tracks.obs.")

    plot_df = adata_tracks.obs[[sample_col, class_col]].copy()
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

    colors = _build_label_color_map(class_order, cmap_name=cmap_name)
    fig_h = max(2.0, 1.2 * len(sample_order))
    fig, axes = plt.subplots(
        nrows=len(sample_order),
        ncols=1,
        figsize=(10.0, fig_h),
        squeeze=False,
    )
    axes = axes.flatten()
    for i, sample_name in enumerate(sample_order):
        ax = axes[i]
        left = 0.0
        for class_name in class_order:
            val = float(proportion_table.loc[sample_name, class_name])
            if val <= 0.0:
                continue
            ax.barh(
                [0.0],
                [val],
                left=left,
                color=colors[class_name],
                height=0.92,
                edgecolor="none",
                linewidth=0.0,
            )
            left += val
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(-0.7, 0.7)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_title(str(sample_name), loc="left", fontsize=10, pad=2.0)

    legend_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="s",
            linestyle="none",
            color=colors[class_name],
            label=str(class_name),
            markersize=7,
        )
        for class_name in class_order
    ]
    fig.legend(
        handles=legend_handles,
        labels=[str(class_name) for class_name in class_order],
        title=str(class_col),
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        borderaxespad=0.0,
        frameon=False,
    )
    fig.tight_layout(rect=[0, 0, 0.82, 1.0])
    fig.savefig(pdf_path, dpi=int(dpi), bbox_inches="tight")
    plt.close(fig)

    color_hex = {str(k): str(to_hex(v)) for k, v in colors.items()}
    return {
        "pdf_path": str(pdf_path),
        "csv_path": str(csv_path),
        "sample_order": [str(x) for x in sample_order],
        "class_order": [str(x) for x in class_order],
        "colors": color_hex,
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
