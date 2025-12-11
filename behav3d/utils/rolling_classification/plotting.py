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
     
def stratified_pick_examples(
    df_windows: pd.DataFrame, 
    X: int, 
    seed: int = 0
    ):
    """
    For each ClusterID in `windows_df`, pick X rows while approximating the
    `sample_name` distribution within that cluster.
    """
    rng = np.random.default_rng(seed)
    selections = []

    for cluster_id, sub in df_windows.groupby("ClusterID"):
        counts = sub["sample_name"].value_counts(normalize=True)
        # Round targets; ensure >=1 per present sample
        targets = {s: max(1, int(round(p * X))) for s, p in counts.items()}
        # Adjust to sum to X
        while sum(targets.values()) > X:
            key = max(targets, key=lambda k: targets[k])
            targets[key] -= 1
        while sum(targets.values()) < X:
            residuals = {s: (counts[s] * X) - targets[s] for s in targets}
            key = max(residuals, key=lambda k: residuals[k])
            targets[key] += 1


        parts = []
        for s, k in targets.items():
            pool = sub[sub["sample_name"] == s]
            if len(pool) <= k:
                parts.append(pool)
            else:
                parts.append(pool.sample(n=k, random_state=seed))
        picked = pd.concat(parts, ignore_index=True)
        if len(picked) > X:
            picked = picked.sample(n=X, random_state=seed)
        selections.append(picked)
    return pd.concat(selections, ignore_index=True)

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

def create_centered_max_projection_cutout(
    df_window,
    df_positions,
    output_folder,
    margin=10,              # half-size in pixels/voxels for X/Y/Z (use ints or (mz,my,mx))
    normalize_per_channel=True,
    pmin=0.0,
    pmax=99.99,
    ):
    """
    For the (sample_name, TrackID) window in df_window, build per-time centered crops
    of size (2*mz+1, 2*my+1, 2*mx+1) around the cell position. Pads at edges so every
    frame is the same size and centered on the cell. Returns RGB max-projection stacks
    for XY, XZ, YZ with shape (T_window, H, W, 3).

    Requires helpers:
      - load_zarr(path) -> array with shape (T, C, Z, Y, X)
      - calc_z_projection(arr, z_axis, projection='max')
      - colorize_channels_to_rgb(arr, percentiles=...)
    """

    # Normalize margins to (mz, my, mx)
    if isinstance(margin, (list, tuple)):
        assert len(margin) == 3, "margin must be int or (mz,my,mx)"
        mz, my, mx = map(int, margin)
    else:
        mz = my = mx = int(margin)

    sample_name = df_window["sample_name"]
    start_t = int(df_window["window_start_position_t"])
    end_t   = int(df_window["window_end_position_t"])

    df_track = df_positions[
        (df_positions["sample_name"] == sample_name) &
        (df_positions["TrackID"] == df_window["TrackID"])
    ]

    zarr_path = Path(output_folder, "images", sample_name, f"{sample_name}.zarr")
    zarr = load_zarr(zarr_path)  # (T, C, Z, Y, X)
    T, C, Z, Y, X = zarr.shape

    if normalize_per_channel:
        # Percentiles per channel from the last timepoint (like your current fn)
        p_img = np.asarray(zarr[-1])
        percentiles = {}
        for c in range(C):
            ch = p_img[c]
            lo = np.percentile(ch, pmin)
            hi = np.percentile(ch, pmax)
            if hi <= lo:
                hi = lo + 1e-6
            percentiles[c] = (float(lo), float(hi))
    else:
        percentiles=None

    # Helper: centered crop with padding (Z,Y,X)
    def _centered_crop_zyx(vol_czyx, cz, cy, cx, mz, my, mx):
        z0 = cz - mz; z1 = cz + mz + 1
        y0 = cy - my; y1 = cy + my + 1
        x0 = cx - mx; x1 = cx + mx + 1

        pad_before = [max(0, -z0), max(0, -y0), max(0, -x0)]
        pad_after  = [max(0,  z1 - vol_czyx.shape[1]),
                    max(0,  y1 - vol_czyx.shape[2]),
                    max(0,  x1 - vol_czyx.shape[3])]

        Z, Y, X = vol_czyx.shape[1], vol_czyx.shape[2], vol_czyx.shape[3]

        sz0 = max(0, z0); sz1 = min(Z, z1)
        sy0 = max(0, y0); sy1 = min(Y, y1)
        sx0 = max(0, x0); sx1 = min(X, x1)

        crop = vol_czyx[:, sz0:sz1, sy0:sy1, sx0:sx1]

        # # choose a padding value that will become bright/white after normalization
        # if np.issubdtype(vol_czyx.dtype, np.integer):
        #     pad_val = np.iinfo(vol_czyx.dtype).max  # e.g. 65535 for uint16
        # else:
        #     pad_val = 1.0  # for float data
        pad_val = 0
        if any(pad_before) or any(pad_after):
            pad_spec = (
                (0, 0),  # C
                (pad_before[0], pad_after[0]),
                (pad_before[1], pad_after[1]),
                (pad_before[2], pad_after[2]),
            )
            crop = np.pad(
                crop,
                pad_spec,
                mode="constant",
                constant_values=pad_val,
            )

        # Ensure exact target shape (C, 2*mz+1, 2*my+1, 2*mx+1)
        target = (vol_czyx.shape[0], 2*mz+1, 2*my+1, 2*mx+1)
        if crop.shape != target:
            pad_z = max(0, target[1] - crop.shape[1])
            pad_y = max(0, target[2] - crop.shape[2])
            pad_x = max(0, target[3] - crop.shape[3])
            crop = np.pad(
                crop,
                ((0, 0), (0, pad_z), (0, pad_y), (0, pad_x)),
                mode="constant",
                constant_values=pad_val,
            )
            crop = crop[:, :target[1], :target[2], :target[3]]

        return crop


    # Build the time stack of centered crops: (T_w, C, Zc, Yc, Xc)
    crops = []
    for t in range(start_t, end_t + 1):
        df_t = df_track[df_track["position_t"] == t]
        if len(df_t) == 0:
            # If no position for this time, repeat last or put zeros
            if len(crops) == 0:
                crops.append(np.zeros((C, 2*mz+1, 2*my+1, 2*mx+1), dtype=zarr.dtype))
            else:
                crops.append(crops[-1])
            continue

        cx = int(df_t["pixel_position_x"].values[0])
        cy = int(df_t["pixel_position_y"].values[0])
        cz = int(df_t["pixel_position_z"].values[0])

        vol_czyx = zarr[t]  # (C, Z, Y, X)
        crop = _centered_crop_zyx(vol_czyx, cz, cy, cx, mz, my, mx)
        crops.append(crop)

    crops = np.stack(crops, axis=0)  # (T_w, C, Zc, Yc, Xc)

    # Projections along the Z/Y/X for XY, XZ, YZ respectively
    # crops axes: (T, C, Z, Y, X)
    xy_proj = calc_z_projection(crops, z_axis=-3, projection='max')  # -> (T, C, Y, X)
    xz_proj = calc_z_projection(crops, z_axis=-2, projection='max')  # -> (T, C, Z, X)
    yz_proj = calc_z_projection(crops, z_axis=-1, projection='max')  # -> (T, C, Z, Y)

    # Colorize per your helper (expects channel_axis=1)
    xy_proj_rgb = colorize_channels_to_rgb(xy_proj, percentiles=percentiles)
    xz_proj_rgb = colorize_channels_to_rgb(xz_proj, percentiles=percentiles)
    yz_proj_rgb = colorize_channels_to_rgb(yz_proj, percentiles=percentiles)

    return xy_proj_rgb, xz_proj_rgb, yz_proj_rgb

def colorize_channels_to_rgb(
    img, 
    channel_axis=1, 
    colors=None, 
    clip=True,
    percentiles=None,
    ):

    if colors is None:
        colors = [
            (0, 255, 255),   # cyan
            (255, 255, 0),   # yellow
            (255, 0, 0),     # red
            (0, 255, 0),     # green
            (255, 105, 180), # pink
        ]

    img = np.asarray(img)
    C = img.shape[channel_axis]
    img_moved = np.moveaxis(img, channel_axis, -1)  # shape (..., C)
    
    if percentiles is not None:
        img_moved = img_moved.astype(np.float32)
        for c in range(C):
            lo, hi = percentiles[c]
            ch = img_moved[..., c]
            img_moved[..., c] = np.clip((ch - lo) / (hi - lo), 0.0, 1.0)
    else:
        dtype_max = np.iinfo(img.dtype).max if np.issubdtype(img.dtype, np.integer) else 1.0
        img_moved = img_moved.astype(np.float32)
        img_moved = img_moved / dtype_max
        
    rgb = np.zeros((*img_moved.shape[:-1], 3), dtype=np.float32)
    for c in range(C):
        color = np.array(colors[c % len(colors)], dtype=np.float32) / 255.0
        rgb += img_moved[..., c, None] * color  # broadcast color

    if clip:
        np.clip(rgb, 0.0, 1.0, out=rgb)

    rgb = (rgb * 255.0).astype(np.uint8)
    return rgb

def _pca_project_xyz(xyz: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    PCA over 3D points (N,3). Returns (proj_2d, mean, components_3x3).
    proj_2d uses the first two principal axes (PC1, PC2).
    """
    assert xyz.ndim == 2 and xyz.shape[1] == 3, "xyz must be (N,3)"
    mean = xyz.mean(axis=0)
    X = xyz - mean
    # SVD for PCA
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    components = Vt  # rows are principal directions, shape (3,3)
    proj2 = X @ components[:2].T  # (N,2)
    return proj2, mean, components

def _prep_traj_limits(traj2: np.ndarray, pad_frac: float) -> Tuple[float, float, float, float]:
    """
    Compute symmetric axes limits for a pleasant view, with a small padding.
    """
    xmin, ymin = traj2.min(axis=0)
    xmax, ymax = traj2.max(axis=0)
    dx = xmax - xmin
    dy = ymax - ymin
    span = max(dx, dy)
    if span == 0:
        span = 1.0
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    pad = span * pad_frac
    return cx - span/2 - pad, cx + span/2 + pad, cy - span/2 - pad, cy + span/2 + pad

def _render_frame_row(
    xy_rgb: np.ndarray, xz_rgb: np.ndarray, yz_rgb: np.ndarray,
    traj2: np.ndarray, t_idx: int, fig_w: float, fig_h: float, dpi: int,
    xy_title: str = "XY", xz_title: str = "XZ", yz_title: str = "YZ",
    pad_frac: float = 0.05
) -> np.ndarray:
    """
    Render a single frame showing three projections and the trajectory (up to t_idx).
    Returns an RGB uint8 image array.
    """
    # Create figure with 1x4 grid
    fig, axes = plt.subplots(1, 4, figsize=(fig_w, fig_h), dpi=dpi,
                             gridspec_kw={"width_ratios": [1,1,1,1]}, constrained_layout=True)
    ax_xy, ax_xz, ax_yz, ax_traj = axes

    # --- Projections (RGB images expected: HxWx3) ---
    for ax, img, title in zip([ax_xy, ax_xz, ax_yz], [xy_rgb, xz_rgb, yz_rgb], [xy_title, xz_title, yz_title]):
        ax.imshow(img)
        h, w = img.shape[:2]
        rect = plt.Rectangle(
            (0, 0), w, h,
            linewidth=2,
            edgecolor='white',
            facecolor='none'
        )
        ax.add_patch(rect)
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])

    # --- Trajectory panel ---
    ax_traj.set_title("Track (PC1–PC2)")
    # draw history up to t_idx (inclusive)
    if traj2.shape[0] > 0:
        t_clamped = max(0, min(t_idx, traj2.shape[0] - 1))

        # Fixed global axes in "pixel units"
        ax_min, ax_max = 0.0, 50.0
        center = 0.5 * (ax_min + ax_max)  # 50

        # PCA output is mean-centered (around 0), so we just shift it
        # so that the mean of the track sits around the center (50, 50).
        traj_shift = traj2 + center

        seg = traj_shift[:t_clamped + 1]
        ax_traj.plot(seg[:, 0], seg[:, 1], linewidth=1.5)
        ax_traj.scatter([seg[-1, 0]], [seg[-1, 1]], s=10)  # current point

        # Use the SAME axes for all tracks, all clusters
        ax_traj.set_xlim(ax_min, ax_max)
        ax_traj.set_ylim(ax_min, ax_max)
        ax_traj.set_aspect("equal", adjustable="box")

        ax_traj.set_xticks([])
        ax_traj.set_yticks([])

    # Draw to array
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    buf = np.asarray(canvas.buffer_rgba())
    plt.close(fig)
    # Convert RGBA -> RGB
    rgb = buf[...,:3].copy()
    return rgb

def _stack_frame_row(xy_t: np.ndarray, xz_t: np.ndarray, yz_t: np.ndarray,
                     traj2: np.ndarray, t_idx: int,
                     fig_w: float, fig_h: float, dpi: int,
                     pad_frac: float,
                     labels: Tuple[str, str, str]) -> np.ndarray:
    """
    Build one row (three projections + trajectory) as a single RGB frame.
    """
    return _render_frame_row(
        xy_rgb=xy_t, xz_rgb=xz_t, yz_rgb=yz_t,
        traj2=traj2, t_idx=t_idx,
        fig_w=fig_w, fig_h=fig_h, dpi=dpi,
        xy_title=labels[0], xz_title=labels[1], yz_title=labels[2],
        pad_frac=pad_frac
    )

def _concat_rows_vertically(rows: List[np.ndarray]) -> np.ndarray:
    """
    Concatenate row images vertically, padding widths as needed.
    """
    widths = [r.shape[1] for r in rows]
    max_w = max(widths)
    out_h = sum(r.shape[0] for r in rows)
    out = np.zeros((out_h, max_w, 3), dtype=np.uint8)
    y = 0
    for r in rows:
        h, w = r.shape[:2]
        out[y:y+h, :w] = r
        y += h
    return out

def _concat_examples_horizontally(examples: List[np.ndarray]) -> np.ndarray:
    """
    Concatenate example images horizontally, padding heights as needed.
    Each example is a full row from _stack_frame_row (projections + track).
    """
    heights = [e.shape[0] for e in examples]
    max_h = max(heights)
    out_w = sum(e.shape[1] for e in examples)
    out = np.zeros((max_h, out_w, 3), dtype=np.uint8)

    x = 0
    for e in examples:
        h, w = e.shape[:2]
        out[:h, x:x+w] = e
        x += w
    return out

def _build_traj2_for_window(df_track: pd.DataFrame, t_start: int, t_end: int) -> np.ndarray:
    """
    Extract XYZ over [t_start, t_end], do PCA once over all timepoints,
    and return the 2D projected coords (N,2) in window order.
    Missing time steps are filled by repeating the last known position.
    """
    # Gather per-time points
    xyz = []
    times = list(range(t_start, t_end+1))
    last = None
    for t in times:
        row = df_track[df_track["position_t"] == t]
        if len(row) == 0:
            if last is None:
                xyz.append([0.0, 0.0, 0.0])
            else:
                xyz.append(last)
        else:
            p = [float(row["pixel_position_x"].values[0]),
                 float(row["pixel_position_y"].values[0]),
                 float(row["pixel_position_z"].values[0])]
            xyz.append(p)
            last = p
    xyz = np.asarray(xyz, dtype=float)
    traj2, mean, comps = _pca_project_xyz(xyz)
    return traj2

def _add_row_title(row_img: np.ndarray, title: str,
                   fontsize: int = 40) -> np.ndarray:
    """
    Add a bold title centered above a row image using matplotlib.
    Returns a new RGB image.
    """
    h, w = row_img.shape[:2]

    # Choose figure size proportional to pixel size so resolution stays reasonable
    dpi = 100
    fig_w = w / dpi
    fig_h = (h + 40) / dpi  # a bit of extra room for the title

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    ax.imshow(row_img)
    ax.set_xticks([])
    ax.set_yticks([])

    ax.set_title(title, fontsize=fontsize, fontweight="bold", pad=10)

    fig.tight_layout(pad=0.1)
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    buf = np.asarray(canvas.buffer_rgba())[..., :3]
    plt.close(fig)
    return buf

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

def _render_fulltrack_frame_with_overlay(
    xy_rgb: np.ndarray,
    xz_rgb: np.ndarray,
    yz_rgb: np.ndarray,
    track_xy: np.ndarray,
    track_xz: np.ndarray,
    track_yz: np.ndarray,
    t_idx: int,
    fig_w: float,
    fig_h: float,
    dpi: int,
    titles: Tuple[str, str, str] = ("XY", "XZ", "YZ"),
    track_color: str = "green",
    ):
    """
    Render a single frame with three projections and an overlaid track that grows
    over time up to t_idx.
    """
    fig, axes = plt.subplots(
        1, 3, figsize=(fig_w, fig_h), dpi=dpi,
        constrained_layout=True
    )
    ax_xy, ax_xz, ax_yz = axes

    # Helper for plotting a growing polyline with a current point
    def _plot_track(ax, track):
        if track is None:
            return
        t_clamped = max(0, min(t_idx, track.shape[0] - 1))
        coords = track[: t_clamped + 1]
        mask = ~np.isnan(coords[:, 0])
        coords = coords[mask]
        if coords.shape[0] == 0:
            return
        ax.plot(coords[:, 0], coords[:, 1], linewidth=1.5, color=track_color)
        ax.scatter(
            coords[-1, 0],
            coords[-1, 1],
            s=10,
            color=track_color,
        )

    # XY
    ax_xy.imshow(xy_rgb)
    _plot_track(ax_xy, track_xy)
    ax_xy.set_title(titles[0])
    ax_xy.set_xticks([])
    ax_xy.set_yticks([])

    # XZ
    ax_xz.imshow(xz_rgb)
    _plot_track(ax_xz, track_xz)
    ax_xz.set_title(titles[1])
    ax_xz.set_xticks([])
    ax_xz.set_yticks([])

    # YZ
    ax_yz.imshow(yz_rgb)
    _plot_track(ax_yz, track_yz)
    ax_yz.set_title(titles[2])
    ax_yz.set_xticks([])
    ax_yz.set_yticks([])

    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    buf = np.asarray(canvas.buffer_rgba())[..., :3]
    plt.close(fig)
    return buf

def create_cluster_videos(
    df_windows: pd.DataFrame,
    df_positions: pd.DataFrame,
    output_folder: str,
    clusters: Optional[List[int]] = None,
    out_dir: str = "cluster_videos",
    normalize_per_channel: bool = False,
    fps: int = 12,
    dpi: int = 200,
    margin=(10, 10, 10),
    pmin: float = 0.0,
    pmax: float = 100.0,
    examples_per_cluster: int = 3,
    seed: int = 0,
    figsize_per_row=(12.0, 4.0),
    traj_pad_frac: float = 0.05,
) -> Dict[int, str]:

    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # Filter windows to selected clusters
    if clusters is not None:
        dfw = df_windows[df_windows["ClusterID"].isin(clusters)].copy()
    else:
        dfw = df_windows.copy()

    # Choose exemplars (X per cluster) using stratified sampler (respects sample_name distribution)
    picks = stratified_pick_examples(dfw, X=examples_per_cluster, seed=seed)

    # Process per cluster
    results: Dict[int, str] = {}
    for cluster_id, sub in picks.groupby("ClusterID"):
        print(f"Processing Cluster {cluster_id} with {len(sub)} exemplars...")
        # Collect rows for this video (each row is one exemplar track)
        rows_info = []
        max_T = 0

        for _, w in sub.iterrows():
            sample_name = w["sample_name"]
            t0 = int(w["window_start_position_t"])
            t1 = int(w["window_end_position_t"])
            track_id = int(w["TrackID"])

            df_track = df_positions[
                (df_positions["sample_name"] == sample_name) &
                (df_positions["TrackID"] == track_id)
            ].copy()

            # Build projections stacks centered per-time
            xy_stack, xz_stack, yz_stack = create_centered_max_projection_cutout(
                w, df_positions, output_folder, normalize_per_channel=normalize_per_channel,
                margin=margin, pmin=pmin, pmax=pmax
            )
            # Shapes: (T_w, H, W, 3)
            T_w = xy_stack.shape[0]
            max_T = max(max_T, T_w)

            # Build projected trajectory (PC1-PC2) for the window
            traj2 = _build_traj2_for_window(df_track, t0, t1)

            rows_info.append({
                "label": f"{sample_name} • Track {track_id}",
                "xy": xy_stack, "xz": xz_stack, "yz": yz_stack,
                "traj2": traj2,
                "T": T_w
            })

        # Build frames by time index
        video_path = out_path / f"cluster_{int(cluster_id)}.mp4"
        
        if video_path.exists():
            video_path.unlink()
        writer = imageio.get_writer(
            video_path, fps=fps, codec="libx264",
            ffmpeg_log_level="warning", quality=9
        )

        try:
            for t in range(max_T):
                row_imgs = []
                for r in rows_info:
                    # Clamp t to row length
                    ti = min(t, r["T"] - 1)
                    frame = _stack_frame_row(
                        xy_t=r["xy"][ti], xz_t=r["xz"][ti], yz_t=r["yz"][ti],
                        traj2=r["traj2"], t_idx=ti,
                        fig_w=figsize_per_row[0], fig_h=figsize_per_row[1], dpi=dpi,
                        pad_frac=traj_pad_frac,
                        labels=('XY', "XZ", "YZ")
                    )
                    row_imgs.append(frame)
                # Vertical concat rows -> full frame
                full_frame = _concat_rows_vertically(row_imgs)
                writer.append_data(full_frame)
        finally:
            writer.close()

        results[int(cluster_id)] = str(video_path)

    return results

def create_cluster_overview_video(
    df_windows: pd.DataFrame,
    df_positions: pd.DataFrame,
    output_folder: str,
    clusters: Optional[List[int]] = None,
    out_dir: str = "cluster_videos",
    normalize_per_channel: bool = False,
    fps: int = 12,
    dpi: int = 200,
    margin=(10, 10, 10),
    pmin: float = 0.0,
    pmax: float = 100.0,
    examples_per_cluster: int = 1,
    seed: int = 0,
    figsize_per_example=(12.0, 4.0),
    traj_pad_frac: float = 0.05,
) -> str:
    
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    video_path = out_path / "clusters_overview.mp4"

    # Filter windows to selected clusters
    if clusters is not None:
        dfw = df_windows[df_windows["ClusterID"].isin(clusters)].copy()
    else:
        dfw = df_windows.copy()

    # Pick exemplars per cluster (same stratified logic as create_cluster_videos)
    picks = stratified_pick_examples(dfw, X=examples_per_cluster, seed=seed)

    # cluster_id -> list of example dicts
    clusters_info: Dict[int, List[Dict]] = {}
    global_max_T = 0

    for cluster_id, sub in picks.groupby("ClusterID"):
        cluster_rows = []
        for _, w in sub.iterrows():
            sample_name = w["sample_name"]
            t0 = int(w["window_start_position_t"])
            t1 = int(w["window_end_position_t"])
            track_id = int(w["TrackID"])

            df_track = df_positions[
                (df_positions["sample_name"] == sample_name) &
                (df_positions["TrackID"] == track_id)
            ].copy()

            xy_stack, xz_stack, yz_stack = create_centered_max_projection_cutout(
                w, df_positions, output_folder,
                normalize_per_channel=normalize_per_channel,
                margin=margin, pmin=pmin, pmax=pmax
            )
            T_w = xy_stack.shape[0]
            global_max_T = max(global_max_T, T_w)

            traj2 = _build_traj2_for_window(df_track, t0, t1)

            cluster_rows.append({
                "cluster_id": int(cluster_id),
                "track_id": track_id,
                "xy": xy_stack, "xz": xz_stack, "yz": yz_stack,
                "traj2": traj2,
                "T": T_w
            })

        clusters_info[int(cluster_id)] = cluster_rows[:examples_per_cluster]

    # Prepare writer for a single overview video
    if video_path.exists():
        video_path.unlink()
    writer = imageio.get_writer(
        video_path, fps=fps, codec="libx264",
        ffmpeg_log_level="warning", quality=9
    )

    try:
        # Loop over time, showing all clusters at once
        for t in range(global_max_T):
            cluster_row_imgs = []

            # Sort clusters by ID for consistent row order
            for cluster_id in sorted(clusters_info.keys()):
                examples = clusters_info[cluster_id]
                example_imgs = []

                for r in examples:
                    ti = min(t, r["T"] - 1)
                    frame = _stack_frame_row(
                        xy_t=r["xy"][ti],
                        xz_t=r["xz"][ti],
                        yz_t=r["yz"][ti],
                        traj2=r["traj2"],
                        t_idx=ti,
                        fig_w=figsize_per_example[0],
                        fig_h=figsize_per_example[1],
                        dpi=dpi,
                        pad_frac=traj_pad_frac,
                        # No sample name here; just track id + projection
                        labels=(f'Track {r["track_id"]} • XY', "XZ", "YZ"),
                    )
                    example_imgs.append(frame)

                # Concatenate examples horizontally to make one row for this cluster
                row_img = _concat_examples_horizontally(example_imgs)

                # Add big bold "Cluster X" above the row
                row_with_title = _add_row_title(row_img, f"Cluster {cluster_id}")
                cluster_row_imgs.append(row_with_title)

            # Stack all cluster rows vertically to make the full frame
            full_frame = _concat_rows_vertically(cluster_row_imgs)
            writer.append_data(full_frame)
    finally:
        writer.close()

    return str(video_path)

def create_fulltrack_cluster_videos(
    df_windows: pd.DataFrame,
    df_positions: pd.DataFrame,
    output_folder: str,
    clusters: Optional[List[int]] = None,
    out_dir: str = "cluster_videos_fulltrack",
    fps: int = 12,
    dpi: int = 200,
    margin: int = 10,
    pmin: float = 0.0,
    pmax: float = 99.99,
    examples_per_cluster: int = 3,
    seed: int = 0,
    figsize_per_row: Tuple[float, float] = (9.0, 3.0),
    normalize_per_channel: bool = True,
    mask_margin: bool = False,
    track_color: str = "green",
    ):
    """
    Create one video per cluster where, for each exemplar:
      * we crop a fixed 3D box that contains the entire track (with `margin`)
      * we compute per-time max projections (XY, XZ, YZ)
      * we overlay the cumulative track of the selected cell center directly
        on top of those max-projection images.

    Returns a non-centered video of selected example T cells
    """
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # Filter windows to selected clusters
    if clusters is not None:
        dfw = df_windows[df_windows["ClusterID"].isin(clusters)].copy()
    else:
        dfw = df_windows.copy()

    # Choose exemplars using same sampler as before
    picks = stratified_pick_examples(dfw, X=examples_per_cluster, seed=seed)

    results: Dict[int, str] = {}

    for cluster_id, sub in picks.groupby("ClusterID"):
        print(f"[fulltrack] Processing Cluster {cluster_id} with {len(sub)} exemplars...")
        rows_info = []
        max_T = 0

        for _, w in sub.iterrows():
            sample_name = w["sample_name"]
            track_id = int(w["TrackID"])
            t0 = int(w["window_start_position_t"])
            t1 = int(w["window_end_position_t"])

            (
                xy_stack,
                xz_stack,
                yz_stack,
                track_xy,
                track_xz,
                track_yz,
            ) = create_fulltrack_max_projection_stacks_with_track(
                w,
                df_positions,
                output_folder,
                margin=margin,
                pmin=pmin,
                pmax=pmax,
                mask_margin=mask_margin,
                normalize_per_channel=normalize_per_channel,
            )

            T_w = xy_stack.shape[0]
            max_T = max(max_T, T_w)

            rows_info.append(
                {
                    "label": f"{sample_name} • Track {track_id}",
                    "xy": xy_stack,
                    "xz": xz_stack,
                    "yz": yz_stack,
                    "track_xy": track_xy,
                    "track_xz": track_xz,
                    "track_yz": track_yz,
                    "T": T_w,
                }
            )

        video_path = out_path / f"cluster_{int(cluster_id)}_fulltrack.mp4"
        if video_path.exists():
            video_path.unlink()

        writer = imageio.get_writer(
            video_path,
            fps=fps,
            codec="libx264",
            ffmpeg_log_level="warning",
            quality=9,
        )

        try:
            for t in range(max_T):
                row_imgs = []
                for r in rows_info:
                    ti = min(t, r["T"] - 1)
                    frame_row = _render_fulltrack_frame_with_overlay(
                        xy_rgb=r["xy"][ti],
                        xz_rgb=r["xz"][ti],
                        yz_rgb=r["yz"][ti],
                        track_xy=r["track_xy"],
                        track_xz=r["track_xz"],
                        track_yz=r["track_yz"],
                        t_idx=ti,
                        fig_w=figsize_per_row[0],
                        fig_h=figsize_per_row[1],
                        dpi=dpi,
                        titles=("XY", "XZ", "YZ"),
                        track_color=track_color,
                    )
                    row_imgs.append(frame_row)

                # stack exemplar rows vertically
                full_frame = _concat_rows_vertically(row_imgs)
                writer.append_data(full_frame)
        finally:
            writer.close()

        results[int(cluster_id)] = str(video_path)

    return results

def test():
    # Example usage (assuming df_analysis is defined)
    import pandas as pd
    import numpy as np
    
    # ssd_dir = r"F:/"
    # ssd_dir = Path(ssd_dir)
    # output_dir = Path(ssd_dir, r"BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE")
    
    downloads_folder = Path("/Users/s.deblank-3/Downloads")
    # downloads_folder = Path(r"C:\Users\Samde\Downloads")
    df_windows_path = downloads_folder / "df_windows.csv"
    df_positions_path = downloads_folder / "df_positions.csv"
    
    # df_analysis.to_csv(df_windows_path)
    # df_tracks_orig.to_csv(df_positions_path)

    df_windows = pd.read_csv(df_windows_path)
    df_positions = pd.read_csv(df_positions_path)

    # output_folder=r"F:/BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE"
    output_folder=r"/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE"

    test = stratified_pick_examples(
        df_windows,
        X=3
    )
    
    chosen_idx=3
    for idx, w in test.iterrows():
        if idx != chosen_idx:
            continue
        xy_stack, xz_stack, yz_stack = create_centered_max_projection_cutout(
            w, df_positions, output_folder, pmax=100
        )
        break
    
    import napari
    viewer = napari.Viewer()
    viewer.add_image(projections[0], rgb=True)
    
    
    create_cluster_videos(
        df_windows,
        df_positions,
        output_folder= output_folder,
        out_dir = r"C:\Users\Samde\Downloads",
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
        df_windows,
        df_positions,
        output_folder=output_folder,
        out_dir=downloads_folder,
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
        df_windows=df_windows,
        df_positions=df_positions,
        output_folder=output_folder,  # folder containing images/<sample>/<sample>.zarr
        out_dir=downloads_folder,
        clusters=None,                # or e.g. [0, 1, 2]
        fps=6,
        margin=20,
        track_color="#63ff33",
        pmin=0.0,
        pmax=99.99,
        examples_per_cluster=4,
        seed=1234,
        normalize_per_channel=True,
        mask_margin=False,
    )