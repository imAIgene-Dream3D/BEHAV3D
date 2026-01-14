import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path
from collections import defaultdict

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
    cluster_col="ClusterID",
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
