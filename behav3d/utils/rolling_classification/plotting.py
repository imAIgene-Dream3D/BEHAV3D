import numpy as np
import pandas as pd
from typing import Optional, Tuple, Dict, Literal, List, Iterable
from pandas.api.types import is_numeric_dtype
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg

import math

from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from dtaidistance import dtw, dtw_ndim
import seaborn as sns
from sklearn.cluster import KMeans, HDBSCAN, AgglomerativeClustering
import umap
from tqdm import tqdm

from sklearn.metrics import adjusted_rand_score
from sklearn.cluster import AgglomerativeClustering
from itertools import combinations

import numpy as np
from scipy import sparse
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
import igraph as ig
import leidenalg as la

from typing import Dict, Iterable, Optional
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy import sparse

import scanpy

from collections import defaultdict
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import TwoSlopeNorm
from behav3d.utils.fileio import load_zarr
from behav3d.utils.preprocessing import calc_z_projection

import imageio_ffmpeg as iioff

from pathlib import Path
from typing import Tuple, Optional, List, Dict

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg

import imageio
import imageio_ffmpeg
import scanpy as sc

def plot_cluster_feature_heatmap(
    df: pd.DataFrame,
    cluster_col: str = "ClusterID",
    feature_cols=None,
    agg: str = "mean",              # "mean", "median", or callable
    sort_clusters: str = "size",    # "size", "index", "hierarchical", or None
    linkage_method: str = "ward",   # used if sort_clusters="hierarchical"
    metric: str = "euclidean",      # used if sort_clusters="hierarchical"
    figsize_scale=(0.35, 0.4),      # (per_cluster_width, per_feature_height)
    title: str | None = None
):
    """
    Heatmap of aggregated feature values per cluster.

    Output axes:
      X = clusters
      Y = features

    Colors:
      blue  = values < 0
      white = 0
      red   = values > 0
    """
    if feature_cols is None:
        raise ValueError("Pass feature_cols (e.g., descriptive_feature_cols).")

    sub = df[[cluster_col] + feature_cols].dropna(subset=[cluster_col]).copy()
    grouped = sub.groupby(cluster_col)[feature_cols].agg(agg)

    # ---- cluster sorting options ----
    if sort_clusters == "size":
        order = sub[cluster_col].value_counts().index
        grouped = grouped.loc[order]

    elif sort_clusters == "index":
        grouped = grouped.sort_index()

    elif sort_clusters == "hierarchical":
        X = grouped.values
        if linkage_method.lower() == "ward":
            Z = linkage(X, method="ward")
        else:
            d = pdist(X, metric=metric)
            Z = linkage(d, method=linkage_method)
        grouped = grouped.iloc[leaves_list(Z)]

    # enforce feature order
    grouped = grouped[feature_cols]

    # TRANSPOSE so rows=features (Y), cols=clusters (X)
    plot_data = grouped.T  # shape: (features, clusters)

    per_clust_w, per_feat_h = figsize_scale

    # --- tweak row height here ---
    per_feat_h = min(per_feat_h, 0.14)   # cap default row height
    # you can try 0.15–0.22 depending on taste

    fig_w = per_clust_w * plot_data.shape[1] + 4
    fig_h = max(2.0, per_feat_h * plot_data.shape[0] + 2)

    # ----- diverging cmap centered at 0 -----
    # symmetric scaling around 0 (best for z-scored inputs)
    max_abs = np.nanmax(np.abs(plot_data.values))
    norm = TwoSlopeNorm(vmin=-max_abs, vcenter=0.0, vmax=max_abs)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(
        plot_data.values,
        aspect="auto",
        interpolation="nearest",
        cmap="RdBu_r",   # blue for <0, red for >0
        norm=norm
    )

    cbar = plt.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label(f"{agg} (input already z-normalized)")

    # X axis = clusters
    ax.set_xticks(np.arange(plot_data.shape[1]))
    ax.set_xticklabels(plot_data.columns, rotation=45, ha="right", fontsize=9)

    # Y axis = features
    ax.set_yticks(np.arange(plot_data.shape[0]))
    ax.set_yticklabels(plot_data.index, fontsize=8)

    ax.set_xlabel(cluster_col)
    ax.set_ylabel("Features")
    ax.set_title(title or f"{agg} feature values per {cluster_col} ({sort_clusters} sorted)")

    plt.tight_layout()
    plt.show()

    return grouped

def plot_feature_cluster_heatmap(
    df_analysis: pd.DataFrame,
    feature_cols,
    cluster_col: str = "ClusterID",
    drop_noise_label: int | None = -1,            # set None to keep all clusters
    feature_category_map: dict[str, str] | None = None,
    zscore_rows: bool = True,                     # z-score features across clusters
    row_distance: str = "abs_correlation",        # {"abs_correlation","correlation","euclidean"}
    col_distance: str = "correlation",            # {"correlation","euclidean"}
    row_linkage: str = "average",
    col_linkage: str = "average",
    cmap: str = "vlag",
    figsize=(12, 14),
    savepath: str | None = None,
):
    """
    RNA-style clustermap:
      - rows = features clustered by (absolute) correlation
      - cols = clusters ordered by similarity
      - cells = mean(feature) within each cluster
    """
    # --- safety checks
    assert cluster_col in df_analysis.columns, f"{cluster_col} not found in df_analysis"
    for c in feature_cols:
        if c not in df_analysis.columns:
            raise ValueError(f"Feature '{c}' not found in df_analysis")

    # --- optionally remove HDBSCAN noise
    df = df_analysis.copy()
    if drop_noise_label is not None:
        df = df[df[cluster_col] != drop_noise_label]

    # --- compute feature x cluster matrix of means
    cluster_means = df.groupby(cluster_col)[list(feature_cols)].mean().T
    cluster_means = cluster_means.replace([np.inf, -np.inf], np.nan).dropna(how="all")

    # remove zero-variance rows (cannot compute correlation)
    nonconst = cluster_means.loc[cluster_means.std(axis=1) > 0]
    if nonconst.empty:
        raise ValueError("All features have zero variance across clusters.")
    M = nonconst.copy()

    # optional row z-scoring (emphasize patterns)
    if zscore_rows:
        M = (M - M.mean(axis=1).to_numpy()[:, None]) / (M.std(axis=1, ddof=0).to_numpy()[:, None] + 1e-9)

    # --- row (feature) distances
    if row_distance == "euclidean":
        d_rows = pdist(M.values, metric="euclidean")
    elif row_distance == "correlation":
        d_rows = pdist(M.values, metric="correlation")  # 1 - corr
    elif row_distance == "abs_correlation":
        C = np.corrcoef(M.values)                       # (n_features x n_features)
        D = 1.0 - np.abs(C)
        np.fill_diagonal(D, 0.0)
        d_rows = squareform(D, checks=False)
    else:
        raise ValueError(f"Unsupported row_distance='{row_distance}'")
    row_link = linkage(d_rows, method=row_linkage)
    row_order = leaves_list(row_link)

    # --- column (cluster) distances
    if col_distance == "euclidean":
        d_cols = pdist(M.T.values, metric="euclidean")
    elif col_distance == "correlation":
        d_cols = pdist(M.T.values, metric="correlation")
    else:
        raise ValueError(f"Unsupported col_distance='{col_distance}'")
    col_link = linkage(d_cols, method=col_linkage)
    col_order = leaves_list(col_link)

    # --- annotation bars
    row_colors = None
    if feature_category_map is not None:
        cats = pd.Series([feature_category_map.get(f, "other") for f in M.index], index=M.index)
        palette = sns.color_palette("tab20", n_colors=cats.nunique())
        lut = dict(zip(cats.unique(), palette))
        row_colors = cats.map(lut)

    cluster_counts = df[cluster_col].value_counts().reindex(M.columns).fillna(0).astype(int)
    norm = Normalize(vmin=cluster_counts.min(), vmax=max(cluster_counts.max(), 1))
    col_colors = [plt.cm.Blues(norm(v)) for v in cluster_counts.values]

    n_rows = M.shape[0]
    row_height = 0.25  # inches per feature row (adjust to taste)
    fig_height = max(6, n_rows * row_height)
    fig_width = figsize[0]
    # --- plot
    g = sns.clustermap(
        M,
        row_linkage=row_link,
        col_linkage=col_link,
        row_colors=row_colors,
        col_colors=col_colors,
        cmap=cmap,
        figsize=(fig_width, fig_height),
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": "mean (row z-score)" if zscore_rows else "mean"},
        dendrogram_ratio=(0.12, 0.12),
        colors_ratio=(0.02, 0.04),
    )
    g.ax_heatmap.set_xlabel("Cluster")
    g.ax_heatmap.set_ylabel("Feature")

    # row legend (feature categories), if provided
    if feature_category_map is not None:
        for cat, color in lut.items():
            g.ax_col_dendrogram.bar(0, 0, color=color, label=cat, linewidth=0)
        g.ax_col_dendrogram.legend(title="Feature category", ncols=min(3, len(lut)), loc="center",
                                   bbox_to_anchor=(0.5, 1.2), frameon=False)

    # column legend (cluster sizes)
    import matplotlib.patches as mpatches
    ticks = np.unique(np.linspace(cluster_counts.min(),
                                  max(cluster_counts.max(), 1), 3).astype(int))
    if len(ticks):
        handles = [mpatches.Patch(color=plt.cm.Blues(norm(v)), label=f"N={v}") for v in ticks]
        g.ax_row_dendrogram.legend(handles=handles, title="Cluster size",
                                   loc="center", bbox_to_anchor=(0.5, 1.15), frameon=False)

    plt.tight_layout()
    if savepath:
        g.savefig(savepath, dpi=300, bbox_inches="tight")

    return g, {
        "matrix": M,
        "row_order": row_order,
        "col_order": col_order,
        "row_linkage": row_link,
        "col_linkage": col_link,
        "cluster_counts": cluster_counts,
    }
    
def plot_hmm_transition_matrix(hmm_model):
    K = hmm_model.n_components
    state_names = [f"S{i}" for i in range(K)]
    T_df = pd.DataFrame(hmm_model.transmat_, index=state_names, columns=state_names)
    T_df

    plt.figure(figsize=(5,4))
    plt.imshow(hmm_model.transmat_, aspect="auto")
    plt.colorbar(label="P(next state)")
    plt.xticks(range(K), state_names)
    plt.yticks(range(K), state_names)
    plt.xlabel("to state")
    plt.ylabel("from state")
    plt.title("HMM transition matrix")
    plt.show()
    
def _pca_rotation(points):
    """
    Get rotation matrix R from PCA (via SVD) WITHOUT changing the origin.
    We compute R on mean-centered data but apply it to the original points.
    Returns R such that rotated = original @ R.T
    """
    P0 = points - points.mean(axis=0, keepdims=True)  # only for orientation
    _, _, Vt = np.linalg.svd(P0, full_matrices=False)
    return Vt  # rows are principal axes

def _align_keep_origin(points):
    """
    Apply PCA rotation while preserving the current origin of `points`.
    """
    R = _pca_rotation(points)
    return points @ R.T

def plot_cluster_window_max_projection(
    df_positions,
    df_windows,
    cluster_id,
    cluster_col="cluster_label_hdbscan",
    # identifiers
    sample_col="sample_name",
    track_col="TrackID",
    time_col="position_t",
    window_start_col="window_start_position_t",
    window_end_col="window_end_position_t",
    # coordinates
    x_col="position_x", 
    y_col="position_y", 
    z_col="position_z",
    # subsampling
    sample_frac=None,          # e.g. 0.2 keeps 20% of windows
    max_windows=None,          # hard cap after (optional) frac sampling
    random_state=123,
    # aesthetics
    line_alpha=0.15,
    line_width=0.8,
    figsize=(12, 4),
    # axis control
    axis_limits=None,          # dict {"pc1": L1, "pc2": L2, "pc3": L3} for fixed limits
    equal_axes=True,
    title_prefix="Max projection (PCA-aligned) for windows in cluster",
):
    """
    For all windows with df_windows[cluster_col] == cluster_id:
      - slice df_positions by (sample, track) and time ∈ [start, end],
      - translate each window so its first point (t0) is at the origin,
      - rotate each window so PC1 is horizontal (X),
      - draw three 2D projections: (PC1, PC2), (PC1, PC3), (PC2, PC3).

    Distances remain in the same units as the coordinate columns.
    """
    # checks
    req_win_cols = {sample_col, track_col, window_start_col, window_end_col, cluster_col}
    missing = req_win_cols - set(df_windows.columns)
    if missing:
        raise ValueError(f"df_windows missing columns: {missing}")
    
    # select cluster's windows
    wins = df_windows.loc[
        df_windows[cluster_col] == cluster_id,
        [sample_col, track_col, window_start_col, window_end_col]
    ].copy()
    if wins.empty:
        raise ValueError(f"No windows found for {cluster_col} == {cluster_id}")

    # subsample windows
    if sample_frac is not None and 0 < sample_frac < 1:
        wins = wins.sample(frac=sample_frac, random_state=random_state)
    if max_windows is not None and len(wins) > max_windows:
        wins = wins.sample(n=max_windows, random_state=random_state)

    # prep positions grouped by (sample, track)
    needed_cols = [sample_col, track_col, time_col, x_col, y_col, z_col]
    for c in needed_cols:
        if c not in df_positions.columns:
            raise ValueError(f"df_positions missing column: {c}")
    pos_sorted = df_positions[needed_cols].sort_values([sample_col, track_col, time_col])

    groups = defaultdict(lambda: None)
    for (s, t), sub in pos_sorted.groupby([sample_col, track_col], sort=False):
        groups[(s, t)] = {
            "t": sub[time_col].to_numpy(),
            "xyz": sub[[x_col, y_col, z_col]].to_numpy(dtype=float)
        }

    # figure
    fig, axes = plt.subplots(1, 3, figsize=figsize, constrained_layout=True)
    ax_xy, ax_xz, ax_yz = axes
    all_xy, all_xz, all_yz = [], [], []
    n_plotted = 0

    # iterate windows
    for s, t, t0, t1 in wins.itertuples(index=False, name=None):
        buf = groups.get((s, t))
        if buf is None:
            continue
        tt = buf["t"]
        sel = (tt >= t0) & (tt <= t1)
        if not np.any(sel):
            continue

        coords = buf["xyz"][sel]
        if coords.shape[0] < 2:
            continue

        coords_centered = coords - coords[0]
        coords_aligned  = _align_keep_origin(coords_centered) 

        all_xy.append(coords_aligned[:, [0, 1]])
        all_xz.append(coords_aligned[:, [0, 2]])
        all_yz.append(coords_aligned[:, [1, 2]])

        ax_xy.plot(coords_aligned[:, 0], coords_aligned[:, 1], linewidth=line_width, alpha=line_alpha)
        ax_xz.plot(coords_aligned[:, 0], coords_aligned[:, 2], linewidth=line_width, alpha=line_alpha)
        ax_yz.plot(coords_aligned[:, 1], coords_aligned[:, 2], linewidth=line_width, alpha=line_alpha)

        n_plotted += 1

    # labels
    ax_xy.set_xlabel("PC1 (max displacement axis)"); ax_xy.set_ylabel("PC2"); ax_xy.set_title("PC1 vs PC2")
    ax_xz.set_xlabel("PC1 (max displacement axis)"); ax_xz.set_ylabel("PC3"); ax_xz.set_title("PC1 vs PC3")
    ax_yz.set_xlabel("PC2"); ax_yz.set_ylabel("PC3"); ax_yz.set_title("PC2 vs PC3")

    # axis limits
    if axis_limits is not None:
        L1, L2, L3 = axis_limits["pc1"], axis_limits["pc2"], axis_limits["pc3"]
        ax_xy.set_xlim(-L1, L1); ax_xy.set_ylim(-L2, L2); ax_xy.set_aspect("equal", adjustable="box")
        ax_xz.set_xlim(-L1, L1); ax_xz.set_ylim(-L3, L3); ax_xz.set_aspect("equal", adjustable="box")
        ax_yz.set_xlim(-L2, L2); ax_yz.set_ylim(-L3, L3); ax_yz.set_aspect("equal", adjustable="box")
    elif equal_axes:
        def _set_equal(ax, pairs):
            if len(pairs) == 0: return
            data = np.vstack(pairs)
            xmin, ymin = data.min(axis=0); xmax, ymax = data.max(axis=0)
            lim = max(abs(xmin), abs(xmax), abs(ymin), abs(ymax))
            lim = 1.0 if lim == 0 else lim
            pad = 0.05 * lim
            ax.set_xlim(-lim - pad, lim + pad); ax.set_ylim(-lim - pad, lim + pad)
            ax.set_aspect("equal", adjustable="box")
        _set_equal(ax_xy, all_xy); _set_equal(ax_xz, all_xz); _set_equal(ax_yz, all_yz)

    fig.suptitle(f"{title_prefix} {cluster_id}  |  windows plotted: {n_plotted}", y=1.02, fontsize=12)
    return fig, axes

def compute_global_pc_axis_limits_for_windows(
    df_positions,
    df_windows,
    cluster_col="cluster_label_hdbscan",
    cluster_ids=None,          # if None, use all unique in df_windows
    include_noise=True,        # include -1 for HDBSCAN noise if present
    # identifiers
    sample_col="sample_name",
    track_col="TrackID",
    time_col="position_t",
    window_start_col="window_start_position_t",
    window_end_col="window_end_position_t",
    # coordinates
    x_col="position_x", 
    y_col="position_y", 
    z_col="position_z",
    # optional subsampling for speed (applied per cluster)
    sample_frac=None,
    max_windows=None,
    random_state=123,
    # alignment function (leave None to use the same as the plotter)
    _pca_align_fn=None,
):
    """
    Returns {'pc1': L1, 'pc2': L2, 'pc3': L3}, each a half-range so you can use [-L, +L].
    """
    if _pca_align_fn is None:
        def _pca_align_fn(P):
            P = P - P.mean(axis=0, keepdims=True)
            U, S, Vt = np.linalg.svd(P, full_matrices=False)
            return P @ Vt.T

    # choose clusters
    clust_series = df_windows[cluster_col].dropna()
    chosen = sorted(clust_series.unique().tolist())
    if not include_noise and -1 in chosen:
        chosen.remove(-1)
    if cluster_ids is not None:
        chosen = [c for c in chosen if c in set(cluster_ids)]
    if not chosen:
        raise ValueError("No clusters to compute limits on.")

    # prep positions grouped by (sample, track)
    needed_cols = [sample_col, track_col, time_col, x_col, y_col, z_col]
    for c in needed_cols:
        if c not in df_positions.columns:
            raise ValueError(f"df_positions missing column: {c}")
    pos_sorted = df_positions[needed_cols].sort_values([sample_col, track_col, time_col])

    groups = defaultdict(lambda: None)
    for (s, t), sub in pos_sorted.groupby([sample_col, track_col], sort=False):
        groups[(s, t)] = {
            "t": sub[time_col].to_numpy(),
            "xyz": sub[[x_col, y_col, z_col]].to_numpy(dtype=float)
        }

    rng = np.random.RandomState(random_state)
    max_abs_pc = np.zeros(3, dtype=float)

    for cid in chosen:
        wins = df_windows.loc[
            df_windows[cluster_col] == cid,
            [sample_col, track_col, window_start_col, window_end_col]
        ]
        if wins.empty:
            continue

        # subsampling for bounds (optional)
        if sample_frac is not None and 0 < sample_frac < 1:
            wins = wins.sample(frac=sample_frac, random_state=random_state)
        if max_windows is not None and len(wins) > max_windows:
            wins = wins.sample(n=max_windows, random_state=random_state)

        for s, t, t0, t1 in wins.itertuples(index=False, name=None):
            buf = groups.get((s, t))
            if buf is None:
                continue
            tt = buf["t"]
            sel = (tt >= t0) & (tt <= t1)
            if not np.any(sel):
                continue

            coords = buf["xyz"][sel]
            if coords.shape[0] < 2:
                continue

            # center at window start; PCA-align (no scaling)
            coords_centered = coords - coords[0]
            R = _pca_rotation(coords_centered)
            coords_aligned = coords_centered @ R.T

            # update global maxima per principal axis
            max_abs_pc = np.maximum(max_abs_pc, np.max(np.abs(coords_aligned), axis=0))

    # padding to keep lines off the border
    pad = np.maximum(1e-9, 0.05 * max_abs_pc)
    max_abs_pc = max_abs_pc + pad

    return {"pc1": float(max_abs_pc[0]), "pc2": float(max_abs_pc[1]), "pc3": float(max_abs_pc[2])}

def plot_all_clusters_window_max_projection(
    df_positions,
    df_windows,
    cluster_col="cluster_label_hdbscan",
    cluster_ids=None,          # e.g. [0,1,2]; if None, uses all unique in df_windows
    include_noise=True,        # set False to skip -1 (HDBSCAN noise)
    # subsampling (applied independently per cluster)
    sample_frac=None,
    max_windows=None,
    random_state=123,
    # fixed axes (pass output of compute_global_pc_axis_limits_for_windows)
    axis_limits=None,
    # saving
    save_dir=None,             # e.g. r"C:\plots\clusters"
    save_pdf_path=None,        # e.g. r"C:\plots\clusters\all_clusters.pdf"
    file_prefix="cluster_windows_",  # filename prefix for per-cluster images
    image_dpi=300,
    # aesthetics forwarded to single-plotter
    line_alpha=0.15,
    line_width=0.8,
    figsize=(12, 4),
    equal_axes=True,
    title_prefix="Max projection (PCA-aligned) for windows in cluster",
    # id/coord column names (forwarded)
    sample_col="sample_name",
    track_col="TrackID",
    time_col="position_t",
    window_start_col="window_start_position_t",
    window_end_col="window_end_position_t",
    x_col="position_x", 
    y_col="position_y", 
    z_col="position_z",
):
    """
    Creates one window-level max-projection plot per cluster_id.
    Returns: {cluster_id: (fig, axes)}
    """
    # cluster list
    clust_series = df_windows[cluster_col].dropna()
    unique_ids = sorted(clust_series.unique().tolist())
    if not include_noise and -1 in unique_ids:
        unique_ids.remove(-1)
    if cluster_ids is not None:
        keep = set(cluster_ids)
        unique_ids = [c for c in unique_ids if c in keep]
    if len(unique_ids) == 0:
        raise ValueError(f"No cluster IDs found in '{cluster_col}' with the current filters.")

    # output dirs / pdf
    if save_dir is not None:
        Path(save_dir).mkdir(parents=True, exist_ok=True)
    pdf = PdfPages(save_pdf_path) if save_pdf_path else None

    axis_limits = compute_global_pc_axis_limits_for_windows(
        df_positions=df_positions,
        df_windows=df_windows,
        cluster_col=cluster_col,
        include_noise=False,
        x_col=x_col, 
        y_col=y_col, 
        z_col=z_col,
    )
    
    out = {}
    try:
        for cid in unique_ids:
            fig, axes = plot_cluster_window_max_projection(
                df_positions=df_positions,
                df_windows=df_windows,
                cluster_id=cid,
                cluster_col=cluster_col,
                sample_col=sample_col,
                track_col=track_col,
                time_col=time_col,
                window_start_col=window_start_col,
                window_end_col=window_end_col,
                x_col=x_col, y_col=y_col, z_col=z_col,
                sample_frac=sample_frac,
                max_windows=max_windows,
                random_state=random_state,
                line_alpha=line_alpha,
                line_width=line_width,
                figsize=figsize,
                axis_limits=axis_limits,   # ensure identical axes across figures (if provided)
                equal_axes=equal_axes,
                title_prefix=title_prefix,
            )
            out[cid] = (fig, axes)

            if save_dir is not None:
                fname = f"{file_prefix}{cluster_col}_{cid}.png"
                fig.savefig(Path(save_dir) / fname, dpi=image_dpi, bbox_inches="tight")
            if pdf is not None:
                pdf.savefig(fig, bbox_inches="tight")
    finally:
        if pdf is not None:
            pdf.close()

    return out

def plot_per_cluster_proportions(
    df,
    groupby = ["ClusterID", "sample_name"],
    show=True
    ):
    prop_df = (
        df
        .groupby(groupby)
        .size()
        .groupby(level=0)
        .apply(lambda x: x / x.sum())  # normalize within each cluster
        .unstack(fill_value=0)         # make sample_name columns
    )
    # Plot stacked bar chart
    prop_df.index = prop_df.index.get_level_values(0)
    # Create stacked bar plot
    ax = prop_df.plot(kind='bar', stacked=True, figsize=(10, 6))

    plt.title("Proportion of sample_name per ClusterID")
    plt.xlabel("ClusterID")
    plt.ylabel("Proportion")
    plt.legend(title="Sample Name", bbox_to_anchor=(1.05, 1), loc='upper left')

    # Ensure upright tick labels
    plt.xticks(ticks=range(len(prop_df.index)), labels=prop_df.index.astype(str), rotation=0, ha='center')

    plt.tight_layout()

    if show:
        plt.show()
    else:
        return ax
 
def plot_number_per_clusters(
    df,
    cluster_col="ClusterID",
    show=True
    ):
    plt.figure(figsize=(6, 3))
    h_counts = df[cluster_col].value_counts().sort_index()
    ax = sns.barplot(x=h_counts.index.astype(str), y=h_counts.values, color="tab:green")
    ax.bar_label(ax.containers[0], padding=3)
    plt.title("Cluster sizes")
    plt.xlabel("Cluster")
    plt.ylabel("Count")
    ax.margins(y=0.15)
    plt.tight_layout()
    
    if show:
        plt.show()
    else:
        return ax
     
def create_max_projection_cutout(
    df_window,
    df_positions,
    output_folder,
    margin = 10,
    pmin = 0,
    pmax = 99.99,
    mask_margin = False
    ):
    
    sample_name = df_window["sample_name"]
    start_t = int(df_window["window_start_position_t"])
    end_t = int(df_window["window_end_position_t"])
    
    df_track = df_positions[
        (df_positions["sample_name"] == sample_name) &
        (df_positions["TrackID"] == df_window["TrackID"])
    ]
    
    zarr_path = Path(output_folder, "images", sample_name, f"{sample_name}.zarr")
    zarr = load_zarr(zarr_path)
    T, C, Z, Y, X = zarr.shape
    
    p_img = np.asarray(zarr[-1])
    percentiles = {}
    for c in range(C):
        ch = p_img[c]
        lo = np.percentile(ch, pmin)
        hi = np.percentile(ch, pmax)

        if hi <= lo:
            hi = lo + 1e-6  # avoid div0
        percentiles[c] = (float(lo), float(hi))
    
    if mask_margin:
        masked = np.zeros_like(zarr)
        for t in range(start_t, end_t+1):
            df_track_t = df_track[df_track["position_t"] == t]
            pos_x = int(df_track_t["pixel_position_x"].values[0])
            pos_y = int(df_track_t["pixel_position_y"].values[0])
            pos_z = int(df_track_t["pixel_position_z"].values[0])
            
            x0 = max(0, pos_x - margin)
            x1 = min(X, pos_x + margin + 1)
            y0 = max(0, pos_y - margin)
            y1 = min(Y, pos_y + margin + 1)
            z0 = max(0, pos_z - margin)
            z1 = min(Z, pos_z + margin + 1)
            
            masked[t, :, z0:z1, y0:y1, x0:x1] = zarr[t, :, z0:z1, y0:y1, x0:x1]
    else:  
        masked = zarr 
    
    z_min = int(df_track["pixel_position_z"].min())
    z_max = int(df_track["pixel_position_z"].max())
    y_min = int(df_track["pixel_position_y"].min())
    y_max = int(df_track["pixel_position_y"].max())
    x_min = int(df_track["pixel_position_x"].min())
    x_max = int(df_track["pixel_position_x"].max())
    
    z_min = max(0, int(z_min - margin))
    y_min = max(0, int(y_min - margin))
    x_min = max(0, int(x_min - margin))

    z_max = min(Z - 1, int(z_max + margin))
    y_max = min(Y - 1, int(y_max + margin))
    x_max = min(X - 1, int(x_max + margin))
    
    cropped_img = masked[
        start_t:end_t,
        :,
        z_min:z_max,
        y_min:y_max,
        x_min:x_max
        ]
    
    cropped_img = np.asarray(cropped_img)
    xy_proj = calc_z_projection(cropped_img, z_axis=-3, projection='max')
    xz_proj = calc_z_projection(cropped_img, z_axis=-2, projection='max')
    yz_proj = calc_z_projection(cropped_img, z_axis=-1, projection='max')
    
    xy_proj_rgb = colorize_channels_to_rgb(xy_proj, percentiles=percentiles)
    xz_proj_rgb = colorize_channels_to_rgb(xz_proj, percentiles=percentiles)
    yz_proj_rgb = colorize_channels_to_rgb(yz_proj, percentiles=percentiles)
    
    return(xy_proj_rgb, xz_proj_rgb, yz_proj_rgb)

def create_fulltrack_max_projection_stacks_with_track(
    df_window,
    df_positions,
    output_folder,
    margin=10,
    pmin=0,
    pmax=99,
    mask_margin=False,
    normalize_per_channel=True,
):
    """
      * crops a fixed box that covers the *full track* (min/max over time + margin)
      * computes per-time max-projection stacks (XY, XZ, YZ)
      * returns 2D track coordinates in crop space for each projection so that
        we can overlay the track on top of the max-projection videos.
    """
    sample_name = df_window["sample_name"]
    start_t = int(df_window["window_start_position_t"])
    end_t = int(df_window["window_end_position_t"])

    # track over all times (we'll index into it per t)
    df_track = df_positions[
        (df_positions["sample_name"] == sample_name)
        & (df_positions["TrackID"].astype(int) == int(df_window["TrackID"]))
    ]

    zarr_path = Path(output_folder, "images", sample_name, f"{sample_name}.zarr")
    zarr = load_zarr(zarr_path)  # (T, C, Z, Y, X)
    T, C, Z, Y, X = zarr.shape

    # Normalize per channel from last timepoint like in create_max_projection_cutout
    if normalize_per_channel:
        p_img = np.asarray(zarr[-1])
        percentiles = {}
        for c in range(C):
            ch = p_img[c]
            lo = np.percentile(ch, pmin)
            hi = np.percentile(ch, pmax)
            hi_floor = 30000
            hi = max(hi, hi_floor)
            # hi = np.percentile(ch, pmax)
            if hi <= lo:
                hi = lo + 1e-6
            percentiles[c] = (float(lo), float(hi))
    else:
        percentiles = None

    # Optionally mask outside a small margin around the track before cropping
    if mask_margin:
        masked = np.zeros_like(zarr)
        for t in range(start_t, end_t + 1):
            df_track_t = df_track[df_track["position_t"] == t]
            if len(df_track_t) == 0:
                continue
            pos_x = int(df_track_t["pixel_position_x"].values[0])
            pos_y = int(df_track_t["pixel_position_y"].values[0])
            pos_z = int(df_track_t["pixel_position_z"].values[0])

            x0 = max(0, pos_x - margin)
            x1 = min(X, pos_x + margin + 1)
            y0 = max(0, pos_y - margin)
            y1 = min(Y, pos_y + margin + 1)
            z0 = max(0, pos_z - margin)
            z1 = min(Z, pos_z + margin + 1)

            masked[t, :, z0:z1, y0:y1, x0:x1] = zarr[t, :, z0:z1, y0:y1, x0:x1]
    else:
        masked = zarr

        # Compute bounding box of the track *within the current time window*
    df_track_window = df_track[
        (df_track["position_t"] >= start_t) & (df_track["position_t"] <= end_t)
    ]
    if len(df_track_window) == 0:
        # Fallback: use the full track if, for some reason, there are no
        # detections inside the window.
        df_bbox = df_track
    else:
        df_bbox = df_track_window
        
    z_min = int(df_bbox["pixel_position_z"].min())
    z_max = int(df_bbox["pixel_position_z"].max())
    y_min = int(df_bbox["pixel_position_y"].min())
    y_max = int(df_bbox["pixel_position_y"].max())
    x_min = int(df_bbox["pixel_position_x"].min())
    x_max = int(df_bbox["pixel_position_x"].max())

    # lower bounds: inclusive
    z_min = max(0, int(z_min - margin))
    y_min = max(0, int(y_min - margin))
    x_min = max(0, int(x_min - margin))

    # upper bounds: make them *exclusive* and clip to array size
    z_max = min(Z, int(z_max + margin + 1))
    y_max = min(Y, int(y_max + margin + 1))
    x_max = min(X, int(x_max + margin + 1))

    # Crop volumes over the time window, fixed spatial box
# Crop volumes over the time window, fixed spatial box
    cropped_img = masked[
        start_t : end_t + 1,
        :,
        z_min:z_max,
        y_min:y_max,
        x_min:x_max,
    ]
    cropped_img = np.asarray(cropped_img)
    # Shapes now: (T_window, C, Zc, Yc, Xc)
    T_window = cropped_img.shape[0]

    # Max-projections per time point
    xy_proj = calc_z_projection(cropped_img, z_axis=-3, projection="max")  # (T_w, C, Yc, Xc)
    xz_proj = calc_z_projection(cropped_img, z_axis=-2, projection="max")  # (T_w, C, Zc, Xc)
    yz_proj = calc_z_projection(cropped_img, z_axis=-1, projection="max")  # (T_w, C, Zc, Yc)

    xy_proj_rgb = colorize_channels_to_rgb(xy_proj, percentiles=percentiles)
    xz_proj_rgb = colorize_channels_to_rgb(xz_proj, percentiles=percentiles)
    yz_proj_rgb = colorize_channels_to_rgb(yz_proj, percentiles=percentiles)

    # Build track coords in crop-space for each timepoint
    track_xy = np.full((T_window, 2), np.nan, dtype=float)
    track_xz = np.full((T_window, 2), np.nan, dtype=float)
    track_yz = np.full((T_window, 2), np.nan, dtype=float)

    for idx, t in enumerate(range(start_t, end_t + 1)):
        df_t = df_track[df_track["position_t"] == t]
        if len(df_t) == 0:
            continue
        pos_x = float(df_t["pixel_position_x"].values[0])
        pos_y = float(df_t["pixel_position_y"].values[0])
        pos_z = float(df_t["pixel_position_z"].values[0])

        # local coords within the cropped box
        x_local = pos_x - x_min
        y_local = pos_y - y_min
        z_local = pos_z - z_min

        # XY projection: axes (Y, X)
        track_xy[idx, 0] = x_local
        track_xy[idx, 1] = y_local

        # XZ projection: axes (Z, X)
        track_xz[idx, 0] = x_local
        track_xz[idx, 1] = z_local

        # YZ projection: axes (Z, Y)
        track_yz[idx, 0] = y_local
        track_yz[idx, 1] = z_local

    return xy_proj_rgb, xz_proj_rgb, yz_proj_rgb, track_xy, track_xz, track_yz

def plot_umap_feature_grid(
    df: pd.DataFrame,
    feature_cols: list[str],
    x_col: str = "UMAP1",
    y_col: str = "UMAP2",
    ncols: int = 4,
    max_plots: int | None = None,
    point_size: int = 5,
    alpha: float = 0.5,
    numeric_cmap: str = "viridis",
    categorical_palette: str = "tab20",
    add_colorbar: bool = True,
    page: int = 0,   # for pagination: 0-based page index
    ):
    """
    Creates a multi-row, multi-column grid of UMAP scatterplots colored by each feature in feature_cols.
    Filters out non-scalar or missing features automatically. Supports pagination via `page`.
    """
    # Filter valid features (exist, scalar, not all NaN)
    valid = []
    for c in feature_cols:
        if c in df.columns and df[c].notna().any():
            valid.append(c)
    if max_plots is not None:
        valid = valid[:max_plots]

    if len(valid) == 0:
        raise ValueError("No valid features to plot.")

    n = len(valid)
    nrows = math.ceil(n / ncols)

    # Pagination support: choose a slice of features per page
    per_page = nrows * ncols
    start = page * per_page
    end = min(start + per_page, len(valid))
    feats = valid[start:end]
    if len(feats) == 0:
        raise ValueError(f"No features to plot on page {page} (only {math.ceil(len(valid)/per_page)} page(s) available).")

    # Axes limits shared across panels
    x_min, x_max = df[x_col].min(), df[x_col].max()
    y_min, y_max = df[y_col].min(), df[y_col].max()

    # Build grid
    fig, axes = plt.subplots(
        nrows=math.ceil(len(feats)/ncols),
        ncols=ncols,
        figsize=(4*ncols, 3.5*math.ceil(len(feats)/ncols)),
        squeeze=False,
        constrained_layout=True
    )

    for i, feat in enumerate(feats):
        r, c = divmod(i, ncols)
        ax = axes[r, c]

        s = df[feat]
        # Numeric vs categorical handling
        if is_numeric_dtype(s):
            # Numeric: use matplotlib scatter for easy colorbar handling
            scanpy = ax.scatter(
                df[x_col], df[y_col],
                s=point_size, alpha=alpha,
                c=s, cmap=numeric_cmap, edgecolors="none"
            )
            if add_colorbar:
                cb = plt.colorbar(scanpy, ax=ax, fraction=0.046, pad=0.04)
                cb.ax.tick_params(labelsize=8)
        else:
            # Categorical: enforce category dtype and use seaborn palette
            s_cat = s.astype("category")
            tmp = df.copy()
            tmp[feat] = s_cat
            sns.scatterplot(
                data=tmp, x=x_col, y=y_col, hue=feat,
                palette=categorical_palette, s=point_size, alpha=alpha,
                legend=False, ax=ax
            )

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_title(feat, fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("")

    # Hide any unused axes (if grid not full)
    total_cells = axes.size
    for j in range(len(feats), total_cells):
        r, c = divmod(j, ncols)
        axes[r, c].axis("off")

    # Add a common title
    fig.suptitle("UMAP colored by features", fontsize=14)
    plt.show()     
    
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
                    ax.axis("off")
                    plot_idx += 1
                    continue

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
                    ax.text(0.5, 0.5, f"{feat}\n(plot unavailable)", ha="center", va="center")
                    ax.axis("off")
                    plot_idx += 1
                    continue

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




import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, PathPatch, Polygon
from matplotlib.path import Path
import matplotlib.cm as cm

def plot_directional_chord(
    M,
    labels=None,
    ax=None,
    *,
    sort_labels=True,
    pad=0.02,                 # gap between outer arcs (radians)
    start_angle=np.pi / 2,    # where label 1 starts (12 o'clock)
    r_outer=1.00,
    r_inner=0.78,
    ribbon_alpha=0.65,
    arc_width=0.16,
    min_value=0.0,            # hide ribbons below this
    cmap_name="tab20",
    arrow_size=0.03,          # relative to radius
    arrow_at=0.93,            # where along ribbon centerline the arrow sits (0..1 toward target)
    figsize=(7, 7),
):
    """
    Draw a chord-like circular plot with *directional* transitions.

    Parameters
    ----------
    M : (n,n) array-like OR pandas.DataFrame
        Transition matrix where M[i,j] is flow from i -> j.
    labels : list[str] | None
        Node labels (if M is not a DataFrame).
    ax : matplotlib axis | None
        Provide an axis to draw on; otherwise creates a new figure/axis.
    sort_labels : bool
        If True, sorts labels; if False preserves input order.
    pad, start_angle : float
        Arc layout controls in radians.
    r_outer, r_inner : float
        Outer radius and inner radius where ribbons attach.
    ribbon_alpha : float
        Ribbon transparency.
    arc_width : float
        Thickness of the outer arcs.
    min_value : float
        Minimum flow to draw a ribbon.
    cmap_name : str
        Colormap used for node colors.
    arrow_size : float
        Arrow triangle size (in radius units).
    arrow_at : float
        Position of arrow along centerline (0 source -> 1 target).
    figsize : tuple
        Figure size if ax is None.

    Returns
    -------
    ax : matplotlib axis
    """

    # ---- load / normalize input
    if isinstance(M, pd.DataFrame):
        df = M.copy()
        if labels is None:
            labels = list(df.index)
        M = df.values.astype(float)
    else:
        M = np.asarray(M, dtype=float)
        if labels is None:
            labels = [str(i) for i in range(M.shape[0])]

    n = M.shape[0]
    if M.shape[1] != n:
        raise ValueError("M must be square (n x n).")

    # optional sorting (keeps matrix aligned)
    if sort_labels:
        order = np.argsort(np.array(labels, dtype=str))
        labels = [labels[i] for i in order]
        M = M[np.ix_(order, order)]

    # remove diagonal for chord rendering (optional; comment out if you want self-loops)
    np.fill_diagonal(M, 0.0)

    outflow = M.sum(axis=1)
    inflow = M.sum(axis=0)
    node_weight = outflow + inflow  # used to size outer arcs
    total = node_weight.sum()
    if total <= 0:
        raise ValueError("Matrix has no positive flow to plot.")

    # ---- axis setup
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(aspect="equal"))
    ax.set_xlim(-1.15, 1.15)
    ax.set_ylim(-1.15, 1.15)
    ax.axis("off")

    # ---- color mapping
    cmap = cm.get_cmap(cmap_name, n)
    node_color = [cmap(i) for i in range(n)]

    # ---- allocate outer arc spans for each node
    full = 2 * np.pi
    usable = full - n * pad
    spans = usable * (node_weight / total)

    node_arc = []  # [(a0,a1)]
    ang = start_angle
    for i in range(n):
        a0 = ang
        a1 = ang - spans[i]
        node_arc.append((a0, a1))
        ang = a1 - pad

    # helper: angle to xy on circle
    def pol2xy(r, theta):
        return np.array([r * np.cos(theta), r * np.sin(theta)])

    # ---- draw outer arcs + labels
    for i, (a0, a1) in enumerate(node_arc):
        theta1 = np.degrees(a1)
        theta2 = np.degrees(a0)
        w = Wedge(
            (0, 0),
            r_outer,
            theta1,
            theta2,
            width=arc_width,
            facecolor=node_color[i],
            edgecolor="white",
            lw=1.2,
            alpha=0.95,
        )
        ax.add_patch(w)

        # label
        am = (a0 + a1) / 2
        p = pol2xy(r_outer + 0.08, am)
        rot = np.degrees(am)
        # keep text upright-ish
        if rot < -90 or rot > 90:
            rot += 180
            ha = "right"
        else:
            ha = "left"
        ax.text(
            p[0],
            p[1],
            labels[i],
            rotation=rot,
            rotation_mode="anchor",
            ha=ha,
            va="center",
            fontsize=10,
        )

    # ---- allocate sub-intervals on each node for outgoing ribbons
    # For node i, split its arc by outgoing magnitudes M[i,j]
    out_sub = [[None] * n for _ in range(n)]  # (s0,s1) angles for each i->j on source arc
    for i in range(n):
        a0, a1 = node_arc[i]
        span = a0 - a1
        if outflow[i] <= 0:
            continue
        cursor = a0
        for j in range(n):
            v = M[i, j]
            if v <= 0:
                continue
            d = span * (v / outflow[i])
            s0 = cursor
            s1 = cursor - d
            out_sub[i][j] = (s0, s1)
            cursor = s1

    # For target placement: also split each node by incoming so ribbons land with less overlap
    in_sub = [[None] * n for _ in range(n)]  # (t0,t1) angles for each i->j on target arc
    for j in range(n):
        a0, a1 = node_arc[j]
        span = a0 - a1
        if inflow[j] <= 0:
            continue
        cursor = a0
        for i in range(n):
            v = M[i, j]
            if v <= 0:
                continue
            d = span * (v / inflow[j])
            t0 = cursor
            t1 = cursor - d
            in_sub[i][j] = (t0, t1)
            cursor = t1

    # ---- ribbon builder (a curved "band" using two Bezier curves)
    def add_ribbon(i, j, s0, s1, t0, t1, value):
        # edge points on inner radius
        p_s0 = pol2xy(r_inner, s0)
        p_s1 = pol2xy(r_inner, s1)
        p_t0 = pol2xy(r_inner, t0)
        p_t1 = pol2xy(r_inner, t1)

        # control points pull toward center to create chord-like curvature
        c_scale = 0.60
        c1 = pol2xy(r_inner * c_scale, (s0 + t0) / 2)
        c2 = pol2xy(r_inner * c_scale, (s1 + t1) / 2)

        verts = [
            p_s0,
            c1,
            p_t0,
            p_t1,
            c2,
            p_s1,
            p_s0,
        ]
        codes = [
            Path.MOVETO,
            Path.CURVE3,
            Path.CURVE3,
            Path.LINETO,
            Path.CURVE3,
            Path.CURVE3,
            Path.CLOSEPOLY,
        ]
        patch = PathPatch(
            Path(verts, codes),
            facecolor=node_color[i],     # color by source
            edgecolor="none",
            alpha=ribbon_alpha,
            zorder=2,
        )
        ax.add_patch(patch)

        # ---- direction arrow (triangle) on the ribbon centerline toward target
        sm = (s0 + s1) / 2
        tm = (t0 + t1) / 2

        # sample a quadratic-ish center curve using two control points
        # (simple interpolation in Cartesian space)
        ps = pol2xy(r_inner * 0.99, sm)
        pt = pol2xy(r_inner * 0.99, tm)
        cc = pol2xy(r_inner * 0.55, (sm + tm) / 2)

        def bezier2(a, b, c, u):
            return (1 - u) ** 2 * a + 2 * (1 - u) * u * b + u**2 * c

        u = arrow_at
        p = bezier2(ps, cc, pt, u)
        p2 = bezier2(ps, cc, pt, min(1.0, u + 0.01))
        v = p2 - p
        nv = np.linalg.norm(v)
        if nv == 0:
            return
        v = v / nv
        # perpendicular
        w = np.array([-v[1], v[0]])

        L = arrow_size
        tri = np.vstack(
            [
                p + v * L * 1.4,
                p - v * L * 0.8 + w * L * 0.8,
                p - v * L * 0.8 - w * L * 0.8,
            ]
        )
        ax.add_patch(
            Polygon(tri, closed=True, facecolor="black", edgecolor="none", alpha=0.55, zorder=3)
        )

    # ---- draw ribbons (iterate in descending value so thick ones go below)
    edges = []
    for i in range(n):
        for j in range(n):
            v = M[i, j]
            if v > min_value and out_sub[i][j] and in_sub[i][j]:
                edges.append((v, i, j, out_sub[i][j], in_sub[i][j]))
    edges.sort(reverse=True, key=lambda x: x[0])

    for v, i, j, (s0, s1), (t0, t1) in edges:
        add_ribbon(i, j, s0, s1, t0, t1, v)

    return ax