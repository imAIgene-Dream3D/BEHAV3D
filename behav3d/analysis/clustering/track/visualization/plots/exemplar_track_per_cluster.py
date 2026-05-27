from pathlib import Path
import textwrap
import time

import imageio.v2 as imageio
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib import colormaps
from matplotlib.collections import LineCollection
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_pdf import PdfPages
from behav3d.analysis.clustering.utils import (
    _mixed_label_sort_key,
    _sanitize_filename_token,
)

from behav3d.analysis.clustering.state.visualization.videos.track_max_projection import (
    create_fulltrack_max_projection_stacks_with_track,
    prepare_fulltrack_max_projection_bundle,
    resolve_tracked_zarr_path,
)
from behav3d.io.images import load_zarr

def _validate_required_columns(df_like, required_cols, name):
    missing = [str(c) for c in required_cols if str(c) not in df_like.columns]
    if len(missing) > 0:
        raise ValueError(f"{name} missing required columns: {missing}")


def _build_state_color_map(adata_full, state_key, cmap_name="tab20"):
    if state_key not in adata_full.obs.columns:
        raise ValueError(f"adata_full.obs missing required state column '{state_key}'.")
    states_series = adata_full.obs[state_key].dropna()
    if len(states_series) == 0:
        raise ValueError(f"No non-null states found in adata_full.obs['{state_key}'].")

    if pd.api.types.is_categorical_dtype(states_series):
        state_values = [str(v) for v in list(states_series.cat.categories)]
    else:
        uniq = [str(v) for v in pd.unique(states_series)]
        state_values = sorted(uniq, key=_mixed_label_sort_key)

    cmap = colormaps.get_cmap(cmap_name).resampled(max(1, len(state_values)))
    color_map = {str(st): cmap(i) for i, st in enumerate(state_values)}
    return state_values, color_map


def select_exemplar_tracks_by_cluster(
    adata_tracks,
    n_per_cluster=5,
    sample_key="sample_name",
    track_key="TrackID",
    cluster_key="ClusterID",
    tmin_key="position_t_min",
    tmax_key="position_t_max",
    seed=0,
    query=None,
):
    """
    Select exemplar tracks per cluster and return:
      chosen_df, cluster_total_counts
    """
    rng = np.random.default_rng(seed)
    needed = {sample_key, track_key, tmin_key, tmax_key, cluster_key}
    missing = [c for c in needed if c not in adata_tracks.obs.columns]
    if missing:
        raise ValueError(f"adata_tracks.obs missing required columns: {missing}")

    full_tracks_df = adata_tracks.obs[[sample_key, track_key, tmin_key, tmax_key, cluster_key]].copy()
    grouped_full = (
        full_tracks_df.groupby([cluster_key, sample_key, track_key], observed=True, as_index=False)
        .agg(**{tmin_key: (tmin_key, "min"), tmax_key: (tmax_key, "max")})
    )
    if len(grouped_full) == 0:
        raise ValueError("No tracks available in adata_tracks.obs.")

    cluster_total_counts = grouped_full.groupby(cluster_key, observed=True).size().to_dict()

    if query is None:
        tracks_df = grouped_full
    else:
        tracks_df = adata_tracks.obs.query(query)[
            [sample_key, track_key, tmin_key, tmax_key, cluster_key]
        ].copy()
        if len(tracks_df) == 0:
            raise ValueError("No tracks left after applying `query` (if any).")
        tracks_df = (
            tracks_df.groupby([cluster_key, sample_key, track_key], observed=True, as_index=False)
            .agg(**{tmin_key: (tmin_key, "min"), tmax_key: (tmax_key, "max")})
        )

    chosen_parts = []
    for _, df_cl in tracks_df.groupby(cluster_key, sort=False, observed=True):
        k = min(int(n_per_cluster), len(df_cl))
        if k <= 0:
            continue
        idx = rng.choice(len(df_cl), size=k, replace=False)
        chosen_parts.append(df_cl.iloc[idx])

    if not chosen_parts:
        raise ValueError("No tracks were selected (n_per_cluster too small or no clusters available).")

    chosen = pd.concat(chosen_parts, axis=0, ignore_index=True)
    return chosen, cluster_total_counts

def plot_tracks_bars_on_ax(
    adata_full,
    chosen_df,
    ax,
    sample_key="sample_name",
    track_key="TrackID",
    time_key="position_t",
    state_key="ClusterID",
    tmin_key="position_t_min",
    tmax_key="position_t_max",
    cmap_name="tab20",
    title=None,
    x_mode="absolute",   # "absolute" or "relative"
    state_color_map=None,
):
    """
    If x_mode="relative": each track is mapped to [0, rel_max] where rel_max is auto-chosen:
      - default: 100
      - if timepoints look discrete and track is short: rel_max = (n_timepoints - 1)
    """
    key_cols = [sample_key, track_key]

    obs = adata_full.obs[[sample_key, track_key, time_key, state_key]].copy()
    chosen_keys = chosen_df[key_cols].drop_duplicates()
    obs = obs.merge(chosen_keys, on=key_cols, how="inner")

    windows = chosen_df[key_cols + [tmin_key, tmax_key]].drop_duplicates()
    obs = obs.merge(windows, on=key_cols, how="left")

    obs["_orig"] = np.arange(len(obs))

    # Auto rel_max helper (per track)
    def _auto_rel_max(tvals: np.ndarray) -> float:
        # If integer-like timepoints and short track, use n-1; else use 100
        if len(tvals) <= 1:
            return 1.0
        # integer-like check (within tolerance)
        int_like = np.all(np.isclose(tvals, np.round(tvals), atol=1e-6))
        if int_like and len(tvals) <= 120:   # heuristic: "short" discrete tracks
            return float(len(tvals) - 1)
        return 100.0

    if x_mode == "relative":
        # compute per-track relative x with per-track rel_max
        obs["_x"] = np.nan

        for (s, t), df in obs.groupby(key_cols, observed=True, sort=False):
            tmin = float(df[tmin_key].iloc[0])
            tmax = float(df[tmax_key].iloc[0])
            denom = (tmax - tmin)

            tvals = df[time_key].astype(float).to_numpy()
            rel_max = _auto_rel_max(tvals)

            if denom == 0:
                rel = np.zeros_like(tvals, dtype=float)
            else:
                rel = (tvals - tmin) / denom
            rel = np.clip(rel, 0.0, 1.0) * rel_max

            obs.loc[df.index, "_x"] = rel

        x_label = "relative time"
        xlim = (0.0, float(np.nanmax(obs["_x"])) if np.isfinite(np.nanmax(obs["_x"])) else 1.0)

    else:
        obs["_x"] = obs[time_key].astype(float)
        x_label = time_key
        xlim = (obs["_x"].min(), obs["_x"].max())

    # sort within each track
    obs = obs.sort_values(key_cols + ["_x", "_orig"])

    if state_color_map is None:
        states = pd.Categorical(obs[state_key]).categories
        cmap = colormaps.get_cmap(cmap_name).resampled(max(1, len(states)))
        color_map = {str(st): cmap(i) for i, st in enumerate(states)}
    else:
        color_map = {str(k): v for k, v in dict(state_color_map).items()}

    tracks = chosen_df[key_cols].drop_duplicates().reset_index(drop=True)

    bar_h = 0.8
    y_gap = 0.25

    for i, row in tracks.iterrows():
        s = row[sample_key]
        t = row[track_key]
        df = obs[(obs[sample_key] == s) & (obs[track_key] == t)]
        if df.empty:
            continue

        x = df["_x"].to_numpy()
        st = df[state_key].to_numpy()

        dx = np.diff(x)
        default_w = np.median(dx[dx > 0]) if np.any(dx > 0) else 1.0
        w = np.r_[dx, default_w]

        y0 = i * (bar_h + y_gap)

        start = x[0]
        cur_state = st[0]
        cur_end = x[0] + w[0]

        for j in range(1, len(x)):
            nxt_state = st[j]
            nxt_start = x[j]
            nxt_end = x[j] + w[j]

            if nxt_state == cur_state and np.isclose(nxt_start, cur_end):
                cur_end = nxt_end
            else:
                ax.broken_barh([(start, cur_end - start)], (y0, bar_h),
                               facecolors=color_map.get(str(cur_state), "grey"))
                start = nxt_start
                cur_state = nxt_state
                cur_end = nxt_end

        ax.broken_barh([(start, cur_end - start)], (y0, bar_h),
                       facecolors=color_map.get(str(cur_state), "grey"))

    ax.set_ylim(-0.2, len(tracks) * (bar_h + y_gap))
    ax.set_xlim(*xlim)
    ax.set_xlabel(x_label)
    ax.set_yticks([])
    if title:
        ax.set_title(title)

def plot_exemplar_tracks_by_cluster(
    adata_full,
    adata_tracks,
    n_per_cluster=5,
    sample_key="sample_name",
    track_key="TrackID",
    time_key="position_t",
    state_key="ClusterID",
    cluster_key="ClusterID",
    tmin_key="position_t_min",
    tmax_key="position_t_max",
    x_mode="relative",
    seed=0,
    query=None,
    cmap_name="tab20",
    max_cols=2,
    legend=True,
    legend_loc="center left",
    legend_bbox_to_anchor=(0.98, 0.5),
    legend_ncol=1,
):
    chosen, cluster_total_counts = select_exemplar_tracks_by_cluster(
        adata_tracks=adata_tracks,
        n_per_cluster=n_per_cluster,
        sample_key=sample_key,
        track_key=track_key,
        cluster_key=cluster_key,
        tmin_key=tmin_key,
        tmax_key=tmax_key,
        seed=seed,
        query=query,
    )

    clusters = list(chosen[cluster_key].dropna().unique())
    n_clusters = len(clusters)
    ncols = min(max_cols, n_clusters) if n_clusters > 0 else 1
    nrows = int(np.ceil(n_clusters / ncols)) if n_clusters > 0 else 1

    fig_h = (
        sum(
            0.45 * min(n_per_cluster, (chosen[chosen[cluster_key] == cl].shape[0])) + 2
            for cl in clusters
        )
        if clusters
        else 4
    )
    fig_w = 14 * ncols
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(fig_w, max(4, fig_h / max(1, ncols))),
        squeeze=False,
    )

    state_values, state_color_map = _build_state_color_map(
        adata_full=adata_full,
        state_key=state_key,
        cmap_name=cmap_name,
    )
    legend_handles = None
    if legend:
        legend_handles = [
            Patch(facecolor=state_color_map.get(str(v), "grey"), edgecolor="k", label=str(v))
            for v in state_values
        ]

    for i, cl in enumerate(clusters):
        r, c = divmod(i, ncols)
        ax = axes[r, c]
        df_cl = chosen.loc[
            chosen[cluster_key] == cl,
            [sample_key, track_key, tmin_key, tmax_key],
        ].copy()

        plot_tracks_bars_on_ax(
            adata_full=adata_full,
            chosen_df=df_cl,
            ax=ax,
            sample_key=sample_key,
            track_key=track_key,
            time_key=time_key,
            state_key=state_key,
            tmin_key=tmin_key,
            tmax_key=tmax_key,
            cmap_name=cmap_name,
            title=f"{cluster_key} = {cl} (N_total={cluster_total_counts.get(cl, 0)})",
            x_mode=x_mode,
            state_color_map=state_color_map,
        )

    for j in range(n_clusters, nrows * ncols):
        r, c = divmod(j, ncols)
        axes[r, c].axis("off")

    if legend and legend_handles:
        fig.legend(
            handles=legend_handles,
            title=state_key,
            loc=legend_loc,
            bbox_to_anchor=legend_bbox_to_anchor,
            frameon=False,
            ncol=legend_ncol,
            fontsize=30,
            title_fontsize=40,
            handlelength=5,
            handleheight=3,
            labelspacing=0.8,
            borderpad=0.3
        )
        fig.tight_layout(rect=(0, 0, 0.92, 1))
    else:
        fig.tight_layout()

    return fig, axes, chosen


def _extract_track_window_obs(
    adata_full,
    *,
    sample_name,
    track_id,
    tmin,
    tmax,
    sample_key="sample_name",
    track_key="TrackID",
    time_key="position_t",
    extra_cols=None,
):
    extra_cols = [] if extra_cols is None else [str(c) for c in list(extra_cols)]
    cols = [sample_key, track_key, time_key]
    cols.extend([c for c in extra_cols if c not in set(cols)])
    obs = adata_full.obs[cols]
    sample_vals = obs[sample_key].astype("string")
    track_vals = obs[track_key].astype("string")
    mask = (sample_vals == str(sample_name)) & (track_vals == str(track_id))
    sub = obs.loc[mask].copy()
    if len(sub) == 0:
        raise ValueError(
            f"No rows found for exemplar sample='{sample_name}', track='{track_id}' in adata_full.obs."
        )
    tvals = pd.to_numeric(sub[time_key], errors="coerce")
    sub = sub.loc[(tvals >= float(tmin)) & (tvals <= float(tmax))].copy()
    if len(sub) == 0:
        raise ValueError(
            "No rows in requested time window for exemplar "
            f"sample='{sample_name}', track='{track_id}', tmin={tmin}, tmax={tmax}."
        )
    sub["_time"] = pd.to_numeric(sub[time_key], errors="coerce")
    sub = sub.dropna(subset=["_time"]).sort_values("_time")
    if len(sub) == 0:
        raise ValueError(
            "No numeric time values for exemplar "
            f"sample='{sample_name}', track='{track_id}', time_key='{time_key}'."
        )
    return sub


def _compute_state_bar_segments(track_df, state_key, time_key, state_color_map):
    times = pd.to_numeric(track_df[time_key], errors="coerce").to_numpy(dtype=float)
    states = track_df[state_key].astype(str).to_numpy()
    if len(times) == 0:
        raise ValueError("Cannot render state bar: no timepoints found.")

    dx = np.diff(times)
    default_w = np.median(dx[dx > 0]) if np.any(dx > 0) else 1.0
    widths = np.r_[dx, default_w]

    segments = []
    start = float(times[0])
    cur_state = str(states[0])
    cur_end = float(times[0] + widths[0])
    for j in range(1, len(times)):
        nxt_state = str(states[j])
        nxt_start = float(times[j])
        nxt_end = float(times[j] + widths[j])
        if (nxt_state == cur_state) and np.isclose(nxt_start, cur_end):
            cur_end = nxt_end
        else:
            segments.append((start, max(0.0, cur_end - start), state_color_map.get(cur_state, "grey")))
            start = nxt_start
            cur_state = nxt_state
            cur_end = nxt_end
    segments.append((start, max(0.0, cur_end - start), state_color_map.get(cur_state, "grey")))

    xmin = float(np.min(times))
    xmax = float(np.max(times + widths))
    if not np.isfinite(xmin) or not np.isfinite(xmax):
        raise ValueError("Invalid state-bar x limits (non-finite values).")
    if np.isclose(xmin, xmax):
        xmax = xmin + 1.0
    return segments, (xmin, xmax)


def _plot_statebar_segments_on_ax(ax, segments, xlim, title=None, cursor_x=None):
    y0 = 0.0
    bar_h = 0.8
    for start, width, color in segments:
        ax.broken_barh([(start, width)], (y0, bar_h), facecolors=color)
    if cursor_x is not None:
        ax.axvline(float(cursor_x), color="black", linewidth=1.5, alpha=0.9)

    ax.set_xlim(*xlim)
    ax.set_ylim(-0.1, 1.1)
    ax.set_yticks([])
    ax.set_xlabel("position_t")
    if title is not None:
        ax.set_title(title, fontsize=8)


def _project_xyz_to_pc12(xyz):
    if xyz.ndim != 2 or xyz.shape[1] != 3:
        raise ValueError(f"xyz must be shaped (N,3), got {xyz.shape}.")
    if xyz.shape[0] < 2:
        raise ValueError("Need at least 2 timepoints to project track trajectory.")
    centered = xyz - xyz.mean(axis=0, keepdims=True)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    proj2 = centered @ vt[:2].T
    return proj2


def _plot_projected_track_state_colored_on_ax(
    ax,
    proj2,
    states,
    *,
    state_color_map,
    title=None,
    fixed_lim=None,
):
    if proj2 is None:
        proj2 = np.zeros((0, 2), dtype=float)
    proj2 = np.asarray(proj2, dtype=float)
    if proj2.ndim != 2 or proj2.shape[1] != 2:
        raise ValueError(f"proj2 must be shaped (N,2), got {proj2.shape}.")

    if proj2.shape[0] >= 2:
        segs = np.stack([proj2[:-1], proj2[1:]], axis=1)
        state_arr = np.asarray(states).astype(str)
        if state_arr.shape[0] < proj2.shape[0]:
            pad = np.repeat(state_arr[-1], proj2.shape[0] - state_arr.shape[0]) if state_arr.shape[0] > 0 else np.repeat("unknown", proj2.shape[0])
            state_arr = np.concatenate([state_arr, pad], axis=0)
        seg_colors = [state_color_map.get(str(v), "grey") for v in state_arr[:-1]]
        lc = LineCollection(segs, colors=seg_colors, linewidths=2.0, alpha=0.95)
        ax.add_collection(lc)
        ax.scatter([proj2[-1, 0]], [proj2[-1, 1]], s=18, color="black", zorder=3)
    elif proj2.shape[0] == 1:
        ax.scatter([proj2[-1, 0]], [proj2[-1, 1]], s=18, color="black", zorder=3)

    lim = None
    if fixed_lim is not None:
        lim = float(fixed_lim)
    elif proj2.shape[0] > 0 and np.isfinite(proj2).all():
        lim = float(np.nanmax(np.abs(proj2)))
    if (lim is None) or (not np.isfinite(lim)) or (lim <= 0):
        lim = 1.0

    pad = 0.08 * lim
    ax.set_xlim(-lim - pad, lim + pad)
    ax.set_ylim(-lim - pad, lim + pad)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    if title is not None:
        ax.set_title(title)


def _plot_track_state_colored_on_ax(
    ax,
    track_df,
    *,
    state_key,
    state_color_map,
    x_col="pixel_position_x",
    y_col="pixel_position_y",
    z_col="pixel_position_z",
    title=None,
    fixed_lim=None,
):
    xyz = track_df[[x_col, y_col, z_col]].to_numpy(dtype=float, copy=False)
    proj2 = _project_xyz_to_pc12(xyz)
    _plot_projected_track_state_colored_on_ax(
        ax=ax,
        proj2=proj2,
        states=track_df[state_key].astype(str).to_numpy(),
        state_color_map=state_color_map,
        title=title,
        fixed_lim=fixed_lim,
    )


def _resolve_trajectory_coordinate_columns(obs):
    pos_triplet = ("position_x", "position_y", "position_z")
    pix_triplet = ("pixel_position_x", "pixel_position_y", "pixel_position_z")
    has_pos = all(c in obs.columns for c in pos_triplet)
    has_pix = all(c in obs.columns for c in pix_triplet)
    if has_pos:
        return pos_triplet
    if has_pix:
        return pix_triplet

    missing_pos = [c for c in pos_triplet if c not in obs.columns]
    missing_pix = [c for c in pix_triplet if c not in obs.columns]
    raise ValueError(
        "Exemplar PDF trajectory panel requires either coordinate triplet "
        "['position_x','position_y','position_z'] or "
        "['pixel_position_x','pixel_position_y','pixel_position_z'] in adata_full.obs. "
        f"Missing position triplet columns: {missing_pos}; missing pixel triplet columns: {missing_pix}."
    )


def _concat_rows_vertically(rows):
    widths = [r.shape[1] for r in rows]
    max_w = max(widths)
    out_h = sum(r.shape[0] for r in rows)
    out = np.zeros((out_h, max_w, 3), dtype=np.uint8)
    y = 0
    for r in rows:
        h, w = r.shape[:2]
        out[y : y + h, :w] = r
        y += h
    return out


def _pad_frame_to_macro_block(frame_rgb, macro_block_size=16):
    frame = np.asarray(frame_rgb)
    if frame.ndim != 3:
        raise ValueError(f"Expected 3D RGB frame, got shape={frame.shape}.")
    h, w = int(frame.shape[0]), int(frame.shape[1])
    block = max(1, int(macro_block_size))
    h_pad = (block - (h % block)) % block
    w_pad = (block - (w % block)) % block
    if h_pad == 0 and w_pad == 0:
        return frame, False, (h, w), (h, w)

    padded = np.pad(
        frame,
        ((0, h_pad), (0, w_pad), (0, 0)),
        mode="constant",
        constant_values=0,
    )
    return padded, True, (h, w), (int(padded.shape[0]), int(padded.shape[1]))


def save_exemplar_statebar_track_pdf_per_cluster(
    adata_full,
    *,
    out_dir,
    chosen_df=None,
    adata_tracks=None,
    n_per_cluster=10,
    sample_key="sample_name",
    track_key="TrackID",
    time_key="position_t",
    state_key="ClusterID",
    cluster_key="ClusterID",
    tmin_key="position_t_min",
    tmax_key="position_t_max",
    rows_per_page=6,
    plot_dpi=300,
    seed=0,
    cmap_name="tab20",
    layout_mode="per_cluster",
    num_example_ranks=5,
):
    mode = _normalize_layout_mode(layout_mode)
    if chosen_df is None:
        if adata_tracks is None:
            raise ValueError("Pass either chosen_df or adata_tracks.")
        chosen_df, _ = select_exemplar_tracks_by_cluster(
            adata_tracks=adata_tracks,
            n_per_cluster=n_per_cluster,
            sample_key=sample_key,
            track_key=track_key,
            cluster_key=cluster_key,
            tmin_key=tmin_key,
            tmax_key=tmax_key,
            seed=seed,
        )

    chosen_df = chosen_df.copy().reset_index(drop=True)
    _validate_required_columns(
        chosen_df,
        [sample_key, track_key, tmin_key, tmax_key, cluster_key],
        name="chosen_df",
    )
    _validate_required_columns(
        adata_full.obs,
        [sample_key, track_key, time_key, state_key],
        name="adata_full.obs",
    )
    traj_x_col, traj_y_col, traj_z_col = _resolve_trajectory_coordinate_columns(adata_full.obs)

    out_dir = Path(out_dir)
    if out_dir.name == "pdf" and out_dir.parent.name in {"per_cluster", "per_example"}:
        out_base = out_dir.parent.parent
    else:
        out_base = out_dir
    out_base.mkdir(parents=True, exist_ok=True)
    per_cluster_dir = out_base / "per_cluster" / "pdf"
    per_example_dir = out_base / "per_example" / "pdf"
    _, state_color_map = _build_state_color_map(
        adata_full=adata_full,
        state_key=state_key,
        cmap_name=cmap_name,
    )

    cluster_vals = list(chosen_df[cluster_key].dropna().unique())
    cluster_vals = sorted(cluster_vals, key=_mixed_label_sort_key)
    paths_by_cluster = {}
    paths_by_example = {}
    prepared_rows = []
    global_lim_candidates = []

    rows_per_page = max(1, int(rows_per_page))

    for _, row in chosen_df.iterrows():
        cluster_val = row[cluster_key]
        if pd.isna(cluster_val):
            continue
        sample_name = row[sample_key]
        track_id = row[track_key]
        tmin = row[tmin_key]
        tmax = row[tmax_key]
        track_df = _extract_track_window_obs(
            adata_full,
            sample_name=sample_name,
            track_id=track_id,
            tmin=tmin,
            tmax=tmax,
            sample_key=sample_key,
            track_key=track_key,
            time_key=time_key,
            extra_cols=[state_key, traj_x_col, traj_y_col, traj_z_col],
        )
        if state_key not in track_df.columns:
            raise ValueError(
                f"Missing '{state_key}' for exemplar sample='{sample_name}', track='{track_id}'."
            )
        segments, xlim = _compute_state_bar_segments(
            track_df=track_df,
            state_key=state_key,
            time_key=time_key,
            state_color_map=state_color_map,
        )
        xyz = track_df[[traj_x_col, traj_y_col, traj_z_col]].to_numpy(dtype=float, copy=False)
        try:
            proj2 = _project_xyz_to_pc12(xyz)
        except Exception:
            if xyz.ndim == 2 and xyz.shape[0] >= 1 and xyz.shape[1] == 3 and np.isfinite(xyz).all():
                proj2 = np.zeros((xyz.shape[0], 2), dtype=float)
            else:
                proj2 = np.zeros((0, 2), dtype=float)

        if proj2.shape[0] > 0 and np.isfinite(proj2).all():
            local_lim = float(np.nanmax(np.abs(proj2)))
            if np.isfinite(local_lim) and local_lim > 0:
                global_lim_candidates.append(local_lim)

        prepared_rows.append(
            {
                "cluster_value": cluster_val,
                "sample_name": sample_name,
                "track_id": track_id,
                "segments": segments,
                "xlim": xlim,
                "proj2": proj2,
                "states": track_df[state_key].astype(str).to_numpy(),
            }
        )

    global_track_panel_limit = (
        float(max(global_lim_candidates)) if len(global_lim_candidates) > 0 else 1.0
    )
    if (not np.isfinite(global_track_panel_limit)) or (global_track_panel_limit <= 0):
        global_track_panel_limit = 1.0

    def _save_group_pdf(pdf_path, group_rows, title_fallback):
        if len(group_rows) == 0:
            fig, ax = plt.subplots(figsize=(10, 2.2))
            ax.axis("off")
            ax.text(
                0.5,
                0.5,
                f"No example tracks available for {title_fallback}.",
                ha="center",
                va="center",
                fontsize=10,
            )
            with PdfPages(pdf_path) as pdf:
                pdf.savefig(fig, dpi=int(plot_dpi), bbox_inches="tight")
            plt.close(fig)
            return

        with PdfPages(pdf_path) as pdf:
            for start in range(0, len(group_rows), rows_per_page):
                page_rows = group_rows[start : start + rows_per_page]
                n_rows = len(page_rows)
                fig, axes = plt.subplots(
                    n_rows,
                    2,
                    figsize=(14, max(3.0, 2.8 * n_rows)),
                    squeeze=False,
                    gridspec_kw={"width_ratios": [1.2, 1.0]},
                    constrained_layout=True,
                )

                for row_i, row_data in enumerate(page_rows):
                    sample_name = row_data["sample_name"]
                    track_id = row_data["track_id"]
                    cluster_value = row_data["cluster_value"]
                    _plot_statebar_segments_on_ax(
                        axes[row_i, 0],
                        segments=row_data["segments"],
                        xlim=row_data["xlim"],
                        title=f"Cluster {cluster_value} | {sample_name} | Track {track_id}",
                    )
                    _plot_projected_track_state_colored_on_ax(
                        axes[row_i, 1],
                        proj2=row_data["proj2"],
                        states=row_data["states"],
                        state_color_map=state_color_map,
                        title="Track trajectory (PC1-PC2), colored by state",
                        fixed_lim=float(global_track_panel_limit),
                    )

                pdf.savefig(fig, dpi=int(plot_dpi), bbox_inches="tight")
                plt.close(fig)

    if mode in {"per_cluster", "both"}:
        per_cluster_dir.mkdir(parents=True, exist_ok=True)
        for cluster_val in cluster_vals:
            cluster_token = _sanitize_filename_token(cluster_val, fallback="cluster")
            pdf_path = per_cluster_dir / f"example_track_cluster_{cluster_token}.pdf"
            sub_rows = [r for r in prepared_rows if r["cluster_value"] == cluster_val]
            if len(sub_rows) == 0:
                continue
            _save_group_pdf(pdf_path, sub_rows, f"cluster '{cluster_val}'")
            paths_by_cluster[str(cluster_val)] = str(pdf_path)

    if mode in {"per_example", "both"}:
        per_example_dir.mkdir(parents=True, exist_ok=True)
        rank_limit = max(1, int(num_example_ranks))
        rank_groups = {int(r): [] for r in range(1, rank_limit + 1)}
        cluster_ranked = {}
        for cluster_val in cluster_vals:
            key = str(cluster_val)
            cluster_ranked[key] = [r for r in prepared_rows if str(r["cluster_value"]) == key]
        for cluster_key_text, rows in cluster_ranked.items():
            for rk, row in enumerate(rows, start=1):
                if rk > rank_limit:
                    break
                rank_groups[int(rk)].append(row)

        for rank in range(1, rank_limit + 1):
            pdf_path = per_example_dir / f"example_track_example_{int(rank):02d}.pdf"
            rank_rows = rank_groups.get(int(rank), [])
            if len(rank_rows) > 0:
                cluster_idx = {str(v): i for i, v in enumerate(cluster_vals)}
                rank_rows = sorted(
                    rank_rows,
                    key=lambda r: cluster_idx.get(str(r["cluster_value"]), 10**9),
                )
            _save_group_pdf(pdf_path, rank_rows, f"example rank {int(rank):02d}")
            paths_by_example[f"{int(rank):02d}"] = str(pdf_path)

    return {
        "pdf_paths_by_cluster": paths_by_cluster,
        "pdf_paths_by_example_rank": paths_by_example,
        "n_clusters": int(len(paths_by_cluster)),
        "n_tracks_selected": int(len(chosen_df)),
        "layout_mode": str(mode),
        "global_track_panel_limit": float(global_track_panel_limit),
    }


def _resolve_sample_zarr_path(output_dir, sample_name):
    expected = Path(output_dir, "images", str(sample_name), f"{sample_name}.zarr")
    if not expected.exists():
        raise ValueError(
            "Missing sample raw zarr for exemplar backprojection video: "
            f"sample='{sample_name}', expected_path='{expected}'."
        )
    return expected


def _normalize_layout_mode(layout_mode):
    mode = str(layout_mode).strip().lower()
    valid = {"per_cluster", "per_example", "both"}
    if mode not in valid:
        raise ValueError(
            "layout_mode must be one of {'per_cluster','per_example','both'}, "
            f"got '{layout_mode}'."
        )
    return mode


def _pad_xy_stack_and_track_to_shape(
    xy_stack,
    track_xy,
    target_h,
    target_w,
    *,
    xy_segment_outline_rgba=None,
):
    if xy_stack.ndim != 4:
        raise ValueError(f"xy_stack must be 4D (T,H,W,3), got shape {xy_stack.shape}.")
    h = int(xy_stack.shape[1])
    w = int(xy_stack.shape[2])
    if h > int(target_h) or w > int(target_w):
        raise ValueError(
            f"Target shape {(target_h, target_w)} smaller than xy shape {(h, w)}."
        )

    pad_top = int((int(target_h) - h) // 2)
    pad_bottom = int(target_h - h - pad_top)
    pad_left = int((int(target_w) - w) // 2)
    pad_right = int(target_w - w - pad_left)

    padded = np.pad(
        xy_stack,
        ((0, 0), (pad_top, pad_bottom), (pad_left, pad_right), (0, 0)),
        mode="constant",
        constant_values=0,
    )

    padded_outline = None
    if xy_segment_outline_rgba is not None:
        xy_segment_outline_rgba = np.asarray(xy_segment_outline_rgba)
        if xy_segment_outline_rgba.ndim != 4 or int(xy_segment_outline_rgba.shape[-1]) != 4:
            raise ValueError(
                "xy_segment_outline_rgba must be 4D (T,H,W,4), "
                f"got shape {xy_segment_outline_rgba.shape}."
            )
        if int(xy_segment_outline_rgba.shape[1]) != h or int(xy_segment_outline_rgba.shape[2]) != w:
            raise ValueError(
                "xy_segment_outline_rgba spatial shape must match xy_stack. "
                f"Got outline={xy_segment_outline_rgba.shape}, xy_stack={xy_stack.shape}."
            )
        padded_outline = np.pad(
            xy_segment_outline_rgba,
            ((0, 0), (pad_top, pad_bottom), (pad_left, pad_right), (0, 0)),
            mode="constant",
            constant_values=0,
        )

    track_shifted = np.array(track_xy, dtype=float, copy=True)
    valid = ~np.isnan(track_shifted[:, 0]) & ~np.isnan(track_shifted[:, 1])
    track_shifted[valid, 0] = track_shifted[valid, 0] + float(pad_left)
    track_shifted[valid, 1] = track_shifted[valid, 1] + float(pad_top)
    return padded, track_shifted, padded_outline


def _render_state_legend_panel(
    *,
    state_values,
    state_color_map,
    state_key,
    height_px,
    dpi=180,
    panel_width_px=280,
):
    panel_h = max(64, int(height_px))
    panel_w = max(120, int(panel_width_px))

    fig = plt.figure(figsize=(panel_w / float(dpi), panel_h / float(dpi)), dpi=int(dpi))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.axis("off")

    handles = [
        Patch(facecolor=state_color_map.get(str(v), "grey"), edgecolor="k", label=str(v))
        for v in list(state_values)
    ]
    if len(handles) > 0:
        ax.legend(
            handles=handles,
            title=str(state_key),
            loc="upper left",
            bbox_to_anchor=(0.02, 0.98),
            borderaxespad=0.0,
            frameon=False,
            fontsize=8,
            title_fontsize=9,
            handlelength=1.8,
            labelspacing=0.45,
        )

    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    panel = np.asarray(canvas.buffer_rgba())[..., :3].copy()
    plt.close(fig)

    if panel.shape[0] == panel_h:
        return panel
    out = np.zeros((panel_h, panel.shape[1], 3), dtype=np.uint8)
    h_copy = min(panel_h, panel.shape[0])
    out[:h_copy, :, :] = panel[:h_copy, :, :]
    return out


def _append_right_panel(base_rgb, panel_rgb, gap_px=8):
    h = int(base_rgb.shape[0])
    panel_h = int(panel_rgb.shape[0])
    if panel_h != h:
        fixed = np.zeros((h, panel_rgb.shape[1], 3), dtype=np.uint8)
        copy_h = min(h, panel_h)
        fixed[:copy_h, :, :] = panel_rgb[:copy_h, :, :]
        panel_rgb = fixed

    gap = np.zeros((h, max(0, int(gap_px)), 3), dtype=np.uint8)
    return np.concatenate([base_rgb, gap, panel_rgb], axis=1)


def _append_left_panel(base_rgb, panel_rgb, gap_px=8):
    h = int(base_rgb.shape[0])
    panel_h = int(panel_rgb.shape[0])
    if panel_h != h:
        fixed = np.zeros((h, panel_rgb.shape[1], 3), dtype=np.uint8)
        copy_h = min(h, panel_h)
        fixed[:copy_h, :, :] = panel_rgb[:copy_h, :, :]
        panel_rgb = fixed

    gap = np.zeros((h, max(0, int(gap_px)), 3), dtype=np.uint8)
    return np.concatenate([panel_rgb, gap, base_rgb], axis=1)


def _render_cluster_label_panel(
    *,
    cluster_value,
    height_px,
    dpi=180,
    panel_width_px=170,
    fontsize=8,
):
    panel_h = max(64, int(height_px))
    panel_w = max(90, int(panel_width_px))
    fig = plt.figure(figsize=(panel_w / float(dpi), panel_h / float(dpi)), dpi=int(dpi))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.axis("off")
    text = textwrap.fill(f"Cluster {cluster_value}", width=14)
    ax.text(0.5, 0.5, text, ha="center", va="center", fontsize=float(fontsize))
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    panel = np.asarray(canvas.buffer_rgba())[..., :3].copy()
    plt.close(fig)
    return panel


def _save_rgb_page_to_pdf(pdf, rgb_img, dpi=220):
    h, w = int(rgb_img.shape[0]), int(rgb_img.shape[1])
    fig = plt.figure(figsize=(w / float(dpi), h / float(dpi)), dpi=int(dpi))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.imshow(rgb_img, interpolation="nearest", aspect="equal")
    ax.axis("off")
    pdf.savefig(fig, dpi=int(dpi))
    plt.close(fig)


def _render_placeholder_frame(message, width_px=1280, height_px=320, dpi=180):
    w = max(480, int(width_px))
    h = max(160, int(height_px))
    fig = plt.figure(figsize=(w / float(dpi), h / float(dpi)), dpi=int(dpi))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.axis("off")
    ax.text(
        0.5,
        0.5,
        str(message),
        ha="center",
        va="center",
        fontsize=14,
        color="white",
        wrap=True,
    )
    fig.patch.set_facecolor("black")
    ax.set_facecolor("black")
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    buf = np.asarray(canvas.buffer_rgba())[..., :3].copy()
    plt.close(fig)
    return buf


def _prepare_ranked_backprojection_rows(
    chosen_df,
    *,
    cluster_key,
    sample_key,
    track_key,
    tmin_key,
    tmax_key,
    examples_per_cluster,
):
    ranked = chosen_df.copy().reset_index(drop=True)
    cluster_vals = list(ranked[cluster_key].dropna().unique())
    cluster_vals = sorted(cluster_vals, key=_mixed_label_sort_key)
    cluster_order = [str(v) for v in cluster_vals]

    ranked["_cluster_text"] = ranked[cluster_key].astype(str)
    ranked["_cluster_order"] = pd.Categorical(
        ranked["_cluster_text"], categories=cluster_order, ordered=True
    )
    ranked = ranked.sort_values(
        ["_cluster_order", str(sample_key), str(track_key), str(tmin_key), str(tmax_key)],
        kind="stable",
    ).reset_index(drop=True)
    ranked["example_rank"] = (
        ranked.groupby(cluster_key, observed=True).cumcount() + 1
    ).astype(int)
    ranked = ranked.loc[ranked["example_rank"] <= int(examples_per_cluster)].reset_index(drop=True)
    return ranked, cluster_order


def _compute_sample_percentiles(zarr_img, *, pmin, pmax):
    _, n_channels, _, _, _ = zarr_img.shape
    p_img = np.asarray(zarr_img[-1])
    percentiles = {}
    for c in range(int(n_channels)):
        ch = p_img[c]
        lo = np.percentile(ch, pmin)
        hi = np.percentile(ch, pmax)
        hi = max(hi, 30000)
        if hi <= lo:
            hi = lo + 1e-6
        percentiles[int(c)] = (float(lo), float(hi))
    return percentiles


def _prepare_exemplar_backprojection_data(
    adata_full,
    *,
    output_dir,
    cell_type=None,
    chosen_df,
    sample_key,
    track_key,
    time_key,
    state_key,
    cluster_key,
    tmin_key,
    tmax_key,
    margin,
    pmin,
    pmax,
    cmap_name,
    coordinate_source_hint,
    examples_per_cluster,
    show_segment_outlines=False,
    segment_style="outline",
    segment_color="#ffffff",
):
    _validate_required_columns(
        chosen_df,
        [sample_key, track_key, tmin_key, tmax_key, cluster_key],
        name="chosen_df",
    )
    _validate_required_columns(
        adata_full.obs,
        [sample_key, track_key, time_key, state_key],
        name="adata_full.obs",
    )
    pixel_triplet = ["pixel_position_x", "pixel_position_y", "pixel_position_z"]
    missing_pixel = [c for c in pixel_triplet if c not in adata_full.obs.columns]
    if len(missing_pixel) > 0:
        src_hint = (
            f" coordinate_source='{coordinate_source_hint}'."
            if (coordinate_source_hint is not None and str(coordinate_source_hint).strip() != "")
            else ""
        )
        raise ValueError(
            "Exemplar backprojection export requires pixel coordinates "
            "['pixel_position_x','pixel_position_y','pixel_position_z'] in adata_full.obs "
            "to align tracks with raw image voxels. "
            f"Missing columns: {missing_pixel}.{src_hint}"
        )

    ranked_df, cluster_order = _prepare_ranked_backprojection_rows(
        chosen_df=chosen_df,
        cluster_key=cluster_key,
        sample_key=sample_key,
        track_key=track_key,
        tmin_key=tmin_key,
        tmax_key=tmax_key,
        examples_per_cluster=examples_per_cluster,
    )
    if len(ranked_df) == 0:
        raise ValueError("No exemplars available after applying per-cluster rank limits.")

    state_values, state_color_map = _build_state_color_map(
        adata_full=adata_full,
        state_key=state_key,
        cmap_name=cmap_name,
    )

    df_positions = adata_full.obs[
        [sample_key, track_key, time_key, "pixel_position_x", "pixel_position_y", "pixel_position_z"]
    ].copy()
    sample_cache = {}
    segment_outline_errors = {}

    rows_info = []
    max_h = 0
    max_w = 0
    for row_id, row in ranked_df.reset_index(drop=True).iterrows():
        sample_name = row[sample_key]
        track_id = row[track_key]
        cluster_value = row[cluster_key]
        rank_value = int(row["example_rank"])
        t0 = int(row[tmin_key])
        t1 = int(row[tmax_key])
        if t1 < t0:
            raise ValueError(
                f"Invalid exemplar window for sample='{sample_name}', track='{track_id}': tmin={t0}, tmax={t1}."
            )

        sample_name_key = str(sample_name)
        if sample_name_key not in sample_cache:
            zarr_path = _resolve_sample_zarr_path(output_dir=output_dir, sample_name=sample_name)
            zarr_img = load_zarr(zarr_path)
            sample_cache[sample_name_key] = {
                "zarr_img": zarr_img,
                "percentiles": _compute_sample_percentiles(
                    zarr_img,
                    pmin=float(pmin),
                    pmax=float(pmax),
                ),
                "tracked_img": None,
                "tracked_img_path": None,
                "tracked_img_error": None,
            }
            if bool(show_segment_outlines):
                try:
                    tracked_img_path = resolve_tracked_zarr_path(
                        output_folder=output_dir,
                        sample_name=sample_name_key,
                        cell_type=cell_type,
                    )
                    if tracked_img_path is None:
                        raise FileNotFoundError(
                            f"Could not resolve tracked zarr for sample='{sample_name_key}'."
                        )
                    sample_cache[sample_name_key]["tracked_img_path"] = str(Path(tracked_img_path))
                    sample_cache[sample_name_key]["tracked_img"] = load_zarr(tracked_img_path)
                except Exception as exc:
                    sample_cache[sample_name_key]["tracked_img_error"] = str(exc)
        cached = sample_cache[sample_name_key]

        track_df = _extract_track_window_obs(
            adata_full,
            sample_name=sample_name,
            track_id=track_id,
            tmin=t0,
            tmax=t1,
            sample_key=sample_key,
            track_key=track_key,
            time_key=time_key,
            extra_cols=[state_key, "pixel_position_x", "pixel_position_y", "pixel_position_z"],
        )
        if state_key not in track_df.columns:
            raise ValueError(
                f"Missing '{state_key}' for exemplar sample='{sample_name}', track='{track_id}'."
            )
        if track_df[pixel_triplet].isna().any().any():
            src_hint = (
                f" coordinate_source='{coordinate_source_hint}'."
                if (coordinate_source_hint is not None and str(coordinate_source_hint).strip() != "")
                else ""
            )
            raise ValueError(
                "Exemplar backprojection export requires non-null pixel coordinates for each exemplar row. "
                f"Found missing pixel coordinates for cluster='{cluster_value}', "
                f"sample='{sample_name}', track='{track_id}'.{src_hint}"
            )

        segments, xlim = _compute_state_bar_segments(
            track_df=track_df,
            state_key=state_key,
            time_key=time_key,
            state_color_map=state_color_map,
        )

        window_row = {
            "sample_name": sample_name,
            "TrackID": track_id,
            "window_start_position_t": t0,
            "window_end_position_t": t1,
        }
        try:
            bundle = prepare_fulltrack_max_projection_bundle(
                df_window=window_row,
                df_positions=df_positions,
                output_folder=output_dir,
                margin=margin,
                pmin=pmin,
                pmax=pmax,
                mask_margin=False,
                normalize_per_channel=True,
                zarr_img=cached["zarr_img"],
                percentiles=cached["percentiles"],
                show_segment_outlines=bool(show_segment_outlines) and (cached["tracked_img_error"] is None),
                segment_style=segment_style,
                segment_color=segment_color,
                tracked_img=cached["tracked_img"],
                tracked_img_path=cached["tracked_img_path"],
                cell_type=cell_type,
            )
        except Exception as exc:
            raise ValueError(
                "Failed to create exemplar backprojection stack for "
                f"cluster='{cluster_value}', sample='{sample_name}', track='{track_id}': {exc}"
            ) from exc

        xy_stack = bundle["xy_proj_rgb"]
        track_xy = bundle["track_xy"]
        xy_segment_outline_rgba = bundle["segment_xy_rgba"]
        row_error = cached.get("tracked_img_error", None) or bundle.get("segment_outline_error", None)
        if row_error is not None:
            row_key = f"{sample_name_key}|track={track_id}|t={t0}-{t1}"
            segment_outline_errors[row_key] = str(row_error)

        if xy_stack.ndim != 4 or int(xy_stack.shape[0]) == 0:
            raise ValueError(
                "Empty XY stack for exemplar "
                f"cluster='{cluster_value}', sample='{sample_name}', track='{track_id}'."
            )

        max_h = max(max_h, int(xy_stack.shape[1]))
        max_w = max(max_w, int(xy_stack.shape[2]))
        rows_info.append(
            {
                "row_id": int(row_id),
                "sample_name": sample_name,
                "track_id": track_id,
                "cluster_value": cluster_value,
                "cluster_text": str(cluster_value),
                "example_rank": int(rank_value),
                "xy_stack": xy_stack,
                "xy_segment_outline_rgba": xy_segment_outline_rgba,
                "track_xy": track_xy,
                "segments": segments,
                "xlim": xlim,
                "t0": int(t0),
                "T": int(xy_stack.shape[0]),
            }
        )

    if len(rows_info) == 0:
        raise ValueError("No exemplar rows were prepared for backprojection export.")

    target_h = int(max_h)
    target_w = int(max_w)
    for row in rows_info:
        xy_pad, track_pad, xy_outline_pad = _pad_xy_stack_and_track_to_shape(
            row["xy_stack"],
            row["track_xy"],
            target_h=target_h,
            target_w=target_w,
            xy_segment_outline_rgba=row.get("xy_segment_outline_rgba", None),
        )
        row["xy_stack"] = xy_pad
        row["track_xy"] = track_pad
        row["xy_segment_outline_rgba"] = xy_outline_pad

    return {
        "rows_info": rows_info,
        "ranked_df": ranked_df,
        "cluster_order": cluster_order,
        "state_values": state_values,
        "state_color_map": state_color_map,
        "global_xy_shape": [int(target_h), int(target_w)],
        "segment_outline_errors": dict(segment_outline_errors),
    }


def _render_statebar_xy_frame(
    *,
    xy_rgb,
    xy_segment_outline_rgba=None,
    track_xy,
    t_idx,
    segments,
    xlim,
    cursor_x,
    sample_name,
    track_id,
    cluster_value,
    track_color="#63ff33",
    fig_w=8.8,
    fig_h=2.6,
    dpi=180,
    show_raw_image=True,
    show_cluster_side_label=True,
    cluster_label_side="left",
    cluster_label_width_px=170,
    cluster_label_fontsize=8,
):
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(fig_w, fig_h),
        dpi=dpi,
        gridspec_kw={"width_ratios": [1.1, 1.4]},
        constrained_layout=True,
    )
    ax_bar, ax_xy = axes

    _plot_statebar_segments_on_ax(
        ax_bar,
        segments=segments,
        xlim=xlim,
        title=f"{sample_name} | Track {track_id}",
        cursor_x=cursor_x,
    )

    base = xy_rgb if bool(show_raw_image) else np.zeros_like(xy_rgb)
    ax_xy.imshow(base, interpolation="nearest", aspect="equal")
    if xy_segment_outline_rgba is not None:
        ax_xy.imshow(xy_segment_outline_rgba, interpolation="nearest", aspect="equal")
    h, w = int(base.shape[0]), int(base.shape[1])
    ax_xy.set_xlim(-0.5, float(w) - 0.5)
    ax_xy.set_ylim(float(h) - 0.5, -0.5)

    t_clamped = max(0, min(int(t_idx), int(track_xy.shape[0]) - 1))
    path = track_xy[: t_clamped + 1]
    valid = ~np.isnan(path[:, 0]) & ~np.isnan(path[:, 1])
    path = path[valid]
    if path.shape[0] >= 2:
        ax_xy.plot(path[:, 0], path[:, 1], linewidth=1.8, color=track_color)
        ax_xy.scatter([path[-1, 0]], [path[-1, 1]], s=14, color=track_color)
    elif path.shape[0] == 1:
        ax_xy.scatter([path[-1, 0]], [path[-1, 1]], s=14, color=track_color)

    ax_xy.set_title(
        "Raw XY + track" if bool(show_raw_image) else "Track overlay",
        fontsize=8,
    )
    ax_xy.set_xticks([])
    ax_xy.set_yticks([])

    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    buf = np.asarray(canvas.buffer_rgba())[..., :3].copy()
    plt.close(fig)
    if bool(show_cluster_side_label):
        panel = _render_cluster_label_panel(
            cluster_value=cluster_value,
            height_px=buf.shape[0],
            dpi=dpi,
            panel_width_px=int(cluster_label_width_px),
            fontsize=float(cluster_label_fontsize),
        )
        side = str(cluster_label_side).strip().lower()
        if side == "right":
            buf = _append_right_panel(buf, panel, gap_px=8)
        else:
            buf = _append_left_panel(buf, panel, gap_px=8)
    return buf


def _init_statebar_xy_row_renderer(
    *,
    xy_shape,
    segments,
    xlim,
    sample_name,
    track_id,
    track_color="#63ff33",
    fig_w=8.8,
    fig_h=2.6,
    dpi=180,
):
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(fig_w, fig_h),
        dpi=dpi,
        gridspec_kw={"width_ratios": [1.1, 1.4]},
        constrained_layout=True,
    )
    ax_bar, ax_xy = axes

    _plot_statebar_segments_on_ax(
        ax_bar,
        segments=segments,
        xlim=xlim,
        title=f"{sample_name} | Track {track_id}",
        cursor_x=None,
    )
    cursor_line = ax_bar.axvline(float(xlim[0]), color="black", linewidth=1.5, alpha=0.9)

    h, w = int(xy_shape[0]), int(xy_shape[1])
    blank = np.zeros((h, w, 3), dtype=np.uint8)
    blank_outline = np.zeros((h, w, 4), dtype=float)
    img_artist = ax_xy.imshow(blank, interpolation="nearest", aspect="equal")
    segment_artist = ax_xy.imshow(blank_outline, interpolation="nearest", aspect="equal")
    ax_xy.set_xlim(-0.5, float(w) - 0.5)
    ax_xy.set_ylim(float(h) - 0.5, -0.5)
    ax_xy.set_title("Raw XY + track", fontsize=8)
    ax_xy.set_xticks([])
    ax_xy.set_yticks([])
    path_line, = ax_xy.plot([], [], linewidth=1.8, color=track_color)
    point_line, = ax_xy.plot([], [], linestyle="None", marker="o", markersize=4, color=track_color)

    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    return {
        "fig": fig,
        "canvas": canvas,
        "img_artist": img_artist,
        "segment_artist": segment_artist,
        "path_line": path_line,
        "point_line": point_line,
        "cursor_line": cursor_line,
        "blank": blank,
        "blank_outline": blank_outline,
    }


def _render_statebar_xy_frame_from_renderer(
    renderer,
    *,
    xy_rgb,
    xy_segment_outline_rgba=None,
    track_xy,
    t_idx,
    cursor_x,
    show_raw_image=True,
    cluster_panel=None,
    cluster_label_side="left",
):
    if bool(show_raw_image):
        renderer["img_artist"].set_data(xy_rgb)
    else:
        renderer["img_artist"].set_data(renderer["blank"])
    if xy_segment_outline_rgba is not None:
        renderer["segment_artist"].set_data(xy_segment_outline_rgba)
    else:
        renderer["segment_artist"].set_data(renderer["blank_outline"])

    t_clamped = max(0, min(int(t_idx), int(track_xy.shape[0]) - 1))
    path = track_xy[: t_clamped + 1]
    valid = ~np.isnan(path[:, 0]) & ~np.isnan(path[:, 1])
    path = path[valid]

    if path.shape[0] >= 2:
        renderer["path_line"].set_data(path[:, 0], path[:, 1])
        renderer["point_line"].set_data([path[-1, 0]], [path[-1, 1]])
    elif path.shape[0] == 1:
        renderer["path_line"].set_data([], [])
        renderer["point_line"].set_data([path[-1, 0]], [path[-1, 1]])
    else:
        renderer["path_line"].set_data([], [])
        renderer["point_line"].set_data([], [])

    renderer["cursor_line"].set_xdata([float(cursor_x), float(cursor_x)])
    renderer["canvas"].draw()
    buf = np.asarray(renderer["canvas"].buffer_rgba())[..., :3].copy()

    if cluster_panel is not None:
        side = str(cluster_label_side).strip().lower()
        if side == "right":
            buf = _append_right_panel(buf, cluster_panel, gap_px=8)
        else:
            buf = _append_left_panel(buf, cluster_panel, gap_px=8)
    return buf


def _close_statebar_xy_row_renderer(renderer):
    fig = renderer.get("fig", None)
    if fig is not None:
        plt.close(fig)


def _group_rows_by_cluster(rows_info, cluster_order):
    cluster_idx = {str(v): i for i, v in enumerate(cluster_order)}
    grouped = {}
    for row in rows_info:
        grouped.setdefault(str(row["cluster_text"]), []).append(row)
    for k in grouped:
        grouped[k] = sorted(grouped[k], key=lambda r: int(r["example_rank"]))
    ordered_keys = sorted(grouped.keys(), key=lambda c: cluster_idx.get(str(c), 10**9))
    return ordered_keys, grouped


def _group_rows_by_rank(rows_info, cluster_order, num_example_ranks):
    cluster_idx = {str(v): i for i, v in enumerate(cluster_order)}
    grouped = {int(r): [] for r in range(1, int(num_example_ranks) + 1)}
    for row in rows_info:
        rk = int(row["example_rank"])
        if rk in grouped:
            grouped[rk].append(row)
    for rk in grouped:
        grouped[rk] = sorted(
            grouped[rk], key=lambda r: cluster_idx.get(str(r["cluster_text"]), 10**9)
        )
    return grouped


def save_exemplar_statebar_backprojection_pdf(
    adata_full,
    *,
    output_dir,
    cell_type=None,
    out_dir,
    chosen_df=None,
    adata_tracks=None,
    n_per_cluster=10,
    sample_key="sample_name",
    track_key="TrackID",
    time_key="position_t",
    state_key="ClusterID",
    cluster_key="ClusterID",
    tmin_key="position_t_min",
    tmax_key="position_t_max",
    rows_per_page=6,
    plot_dpi=220,
    margin=20,
    pmin=0.0,
    pmax=99.99,
    track_color="#63ff33",
    show_segment_outlines=False,
    segment_style="outline",
    segment_color="#ffffff",
    coordinate_source_hint=None,
    seed=0,
    cmap_name="tab20",
    layout_mode="both",
    examples_per_cluster=5,
    num_example_ranks=5,
    show_raw_image=True,
    show_state_legend=True,
    prepared_data=None,
    verbose=False,
):
    t_start = time.perf_counter()
    mode = _normalize_layout_mode(layout_mode)
    if prepared_data is None and chosen_df is None:
        if adata_tracks is None:
            raise ValueError("Pass either chosen_df or adata_tracks.")
        chosen_df, _ = select_exemplar_tracks_by_cluster(
            adata_tracks=adata_tracks,
            n_per_cluster=n_per_cluster,
            sample_key=sample_key,
            track_key=track_key,
            cluster_key=cluster_key,
            tmin_key=tmin_key,
            tmax_key=tmax_key,
            seed=seed,
        )

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    per_cluster_dir = out_dir / "per_cluster" / "pdf"
    per_example_dir = out_dir / "per_example" / "pdf"
    if prepared_data is None:
        t_prep = time.perf_counter()
        prep = _prepare_exemplar_backprojection_data(
            adata_full=adata_full,
            output_dir=output_dir,
            cell_type=cell_type,
            chosen_df=chosen_df.copy().reset_index(drop=True),
            sample_key=sample_key,
            track_key=track_key,
            time_key=time_key,
            state_key=state_key,
            cluster_key=cluster_key,
            tmin_key=tmin_key,
            tmax_key=tmax_key,
            margin=margin,
            pmin=pmin,
            pmax=pmax,
            cmap_name=cmap_name,
            coordinate_source_hint=coordinate_source_hint,
            examples_per_cluster=examples_per_cluster,
            show_segment_outlines=bool(show_segment_outlines),
            segment_style=segment_style,
            segment_color=segment_color,
        )
        if bool(verbose):
            print(
                "[example-track-backprojection][pdf] "
                f"prepared_data in {time.perf_counter() - t_prep:.2f}s"
            )
    else:
        prep = prepared_data
        if bool(verbose):
            print("[example-track-backprojection][pdf] reusing prepared_data")

    rows_info = prep["rows_info"]
    state_values = prep["state_values"]
    state_color_map = prep["state_color_map"]
    cluster_order = prep["cluster_order"]
    rows_per_page = max(1, int(rows_per_page))

    paths_by_cluster = {}
    paths_by_example = {}

    def _save_group_pdf(pdf_path, group_rows, title_fallback):
        with PdfPages(pdf_path) as pdf:
            if len(group_rows) == 0:
                placeholder = _render_placeholder_frame(
                    message=f"No exemplars available for {title_fallback}.",
                    width_px=1200,
                    height_px=320,
                    dpi=plot_dpi,
                )
                if bool(show_state_legend):
                    legend = _render_state_legend_panel(
                        state_values=state_values,
                        state_color_map=state_color_map,
                        state_key=state_key,
                        height_px=placeholder.shape[0],
                        dpi=plot_dpi,
                    )
                    placeholder = _append_right_panel(placeholder, legend)
                _save_rgb_page_to_pdf(pdf, placeholder, dpi=plot_dpi)
                return

            for start in range(0, len(group_rows), rows_per_page):
                chunk = group_rows[start : start + rows_per_page]
                row_frames = []
                for row in chunk:
                    t_last = int(row["T"]) - 1
                    row_frame = _render_statebar_xy_frame(
                        xy_rgb=row["xy_stack"][t_last],
                        xy_segment_outline_rgba=None
                        if row.get("xy_segment_outline_rgba", None) is None
                        else row["xy_segment_outline_rgba"][t_last],
                        track_xy=row["track_xy"],
                        t_idx=t_last,
                        segments=row["segments"],
                        xlim=row["xlim"],
                        cursor_x=float(row["t0"] + t_last),
                        sample_name=row["sample_name"],
                        track_id=row["track_id"],
                        cluster_value=row["cluster_value"],
                        track_color=track_color,
                        dpi=plot_dpi,
                        show_raw_image=show_raw_image,
                    )
                    row_frames.append(row_frame)
                page = _concat_rows_vertically(row_frames)
                if bool(show_state_legend):
                    legend = _render_state_legend_panel(
                        state_values=state_values,
                        state_color_map=state_color_map,
                        state_key=state_key,
                        height_px=page.shape[0],
                        dpi=plot_dpi,
                    )
                    page = _append_right_panel(page, legend)
                _save_rgb_page_to_pdf(pdf, page, dpi=plot_dpi)

    if mode in {"per_cluster", "both"}:
        per_cluster_dir.mkdir(parents=True, exist_ok=True)
        cluster_keys, rows_by_cluster = _group_rows_by_cluster(rows_info, cluster_order)
        for cluster_text in cluster_keys:
            group_rows = rows_by_cluster.get(cluster_text, [])
            cluster_token = _sanitize_filename_token(cluster_text, fallback="cluster")
            pdf_path = per_cluster_dir / f"example_track_backprojection_cluster_{cluster_token}.pdf"
            _save_group_pdf(pdf_path, group_rows, f"cluster '{cluster_text}'")
            paths_by_cluster[str(cluster_text)] = str(pdf_path)

    if mode in {"per_example", "both"}:
        per_example_dir.mkdir(parents=True, exist_ok=True)
        rows_by_rank = _group_rows_by_rank(
            rows_info=rows_info,
            cluster_order=cluster_order,
            num_example_ranks=int(num_example_ranks),
        )
        for rank in range(1, int(num_example_ranks) + 1):
            rank_rows = rows_by_rank.get(rank, [])
            pdf_path = per_example_dir / f"example_track_backprojection_example_{int(rank):02d}.pdf"
            _save_group_pdf(pdf_path, rank_rows, f"example rank {int(rank):02d}")
            paths_by_example[f"{int(rank):02d}"] = str(pdf_path)

    out = {
        "pdf_paths_by_cluster": paths_by_cluster,
        "pdf_paths_by_example_rank": paths_by_example,
        "n_tracks_selected": int(len(prep["ranked_df"])),
        "layout_mode": str(mode),
        "global_xy_shape": list(prep["global_xy_shape"]),
        "segment_outline_errors": dict(prep.get("segment_outline_errors", {})),
    }
    if bool(verbose):
        elapsed = time.perf_counter() - t_start
        print(
            "[example-track-backprojection][pdf] "
            f"completed in {elapsed:.2f}s | per_cluster={len(paths_by_cluster)} | "
            f"per_example={len(paths_by_example)}"
        )
    return out


def save_exemplar_statebar_backprojection_video_per_cluster(
    adata_full,
    *,
    output_dir,
    cell_type=None,
    out_dir,
    chosen_df=None,
    adata_tracks=None,
    n_per_cluster=10,
    sample_key="sample_name",
    track_key="TrackID",
    time_key="position_t",
    state_key="ClusterID",
    cluster_key="ClusterID",
    tmin_key="position_t_min",
    tmax_key="position_t_max",
    fps=6,
    dpi=180,
    margin=20,
    pmin=0.0,
    pmax=99.99,
    track_color="#63ff33",
    show_segment_outlines=False,
    segment_style="outline",
    segment_color="#ffffff",
    coordinate_source_hint=None,
    seed=0,
    cmap_name="tab20",
    layout_mode="both",
    examples_per_cluster=5,
    num_example_ranks=5,
    show_state_legend=True,
    ffmpeg_preset="veryfast",
    ffmpeg_crf=21,
    prepared_data=None,
    verbose=False,
):
    t_start = time.perf_counter()
    mode = _normalize_layout_mode(layout_mode)
    if prepared_data is None and chosen_df is None:
        if adata_tracks is None:
            raise ValueError("Pass either chosen_df or adata_tracks.")
        chosen_df, _ = select_exemplar_tracks_by_cluster(
            adata_tracks=adata_tracks,
            n_per_cluster=n_per_cluster,
            sample_key=sample_key,
            track_key=track_key,
            cluster_key=cluster_key,
            tmin_key=tmin_key,
            tmax_key=tmax_key,
            seed=seed,
        )

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    per_cluster_dir = out_dir / "per_cluster" / "mp4"
    per_example_dir = out_dir / "per_example" / "mp4"

    if prepared_data is None:
        t_prep = time.perf_counter()
        prep = _prepare_exemplar_backprojection_data(
            adata_full=adata_full,
            output_dir=output_dir,
            cell_type=cell_type,
            chosen_df=chosen_df.copy().reset_index(drop=True),
            sample_key=sample_key,
            track_key=track_key,
            time_key=time_key,
            state_key=state_key,
            cluster_key=cluster_key,
            tmin_key=tmin_key,
            tmax_key=tmax_key,
            margin=margin,
            pmin=pmin,
            pmax=pmax,
            cmap_name=cmap_name,
            coordinate_source_hint=coordinate_source_hint,
            examples_per_cluster=examples_per_cluster,
            show_segment_outlines=bool(show_segment_outlines),
            segment_style=segment_style,
            segment_color=segment_color,
        )
    else:
        prep = prepared_data
    if bool(verbose):
        print(
            "[example-track-backprojection][mp4] "
            f"Prepared data in {time.perf_counter() - t_prep:.2f}s"
            if prepared_data is None
            else "[example-track-backprojection][mp4] Reusing prepared data"
        )
    rows_info = prep["rows_info"]
    state_values = prep["state_values"]
    state_color_map = prep["state_color_map"]
    cluster_order = prep["cluster_order"]
    paths_by_cluster = {}
    paths_by_example = {}

    legend_cache = {}

    def _legend_for_height(height_px):
        if not bool(show_state_legend):
            return None
        key = int(height_px)
        if key not in legend_cache:
            legend_cache[key] = _render_state_legend_panel(
                state_values=state_values,
                state_color_map=state_color_map,
                state_key=state_key,
                height_px=key,
                dpi=dpi,
            )
        return legend_cache[key]

    def _write_group_video(video_path, group_rows, title_fallback):
        t_group = time.perf_counter()
        if video_path.exists():
            video_path.unlink()
        if bool(verbose):
            print(f"[example-track-backprojection][mp4] Backprojecting {title_fallback}")

        writer = imageio.get_writer(
            video_path,
            fps=int(fps),
            codec="libx264",
            ffmpeg_log_level="warning",
            macro_block_size=16,
            pixelformat="yuv420p",
            output_params=[
                "-preset",
                str(ffmpeg_preset),
                "-crf",
                str(int(ffmpeg_crf)),
            ],
        )
        try:
            if len(group_rows) == 0:
                frame = _render_placeholder_frame(
                    message=f"No exemplars available for {title_fallback}.",
                    width_px=1200,
                    height_px=320,
                    dpi=dpi,
                )
                legend = _legend_for_height(frame.shape[0])
                if legend is not None:
                    frame = _append_right_panel(frame, legend)
                frame, was_padded, shape_in, shape_out = _pad_frame_to_macro_block(frame, 16)
                writer.append_data(frame)
                return

            renderers = []
            try:
                t_build = time.perf_counter()
                for row in group_rows:
                    renderer = _init_statebar_xy_row_renderer(
                        xy_shape=row["xy_stack"][0].shape,
                        segments=row["segments"],
                        xlim=row["xlim"],
                        sample_name=row["sample_name"],
                        track_id=row["track_id"],
                        track_color=track_color,
                        dpi=dpi,
                    )
                    cluster_panel = _render_cluster_label_panel(
                        cluster_value=row["cluster_value"],
                        height_px=renderer["blank"].shape[0],
                        dpi=dpi,
                        panel_width_px=170,
                        fontsize=8,
                    )
                    renderers.append(
                        {
                            "row": row,
                            "renderer": renderer,
                            "cluster_panel": cluster_panel,
                        }
                    )
                build_elapsed = time.perf_counter() - t_build

                max_t = max(int(r["T"]) for r in group_rows)
                t_frames = time.perf_counter()
                for t_idx in range(max_t):
                    row_frames = []
                    for item in renderers:
                        row = item["row"]
                        ti = min(int(t_idx), int(row["T"]) - 1)
                        frame = _render_statebar_xy_frame_from_renderer(
                            item["renderer"],
                            xy_rgb=row["xy_stack"][ti],
                            xy_segment_outline_rgba=None
                            if row.get("xy_segment_outline_rgba", None) is None
                            else row["xy_segment_outline_rgba"][ti],
                            track_xy=row["track_xy"],
                            t_idx=ti,
                            cursor_x=float(row["t0"] + ti),
                            show_raw_image=True,
                            cluster_panel=item["cluster_panel"],
                            cluster_label_side="left",
                        )
                        row_frames.append(frame)

                    full_frame = _concat_rows_vertically(row_frames)
                    legend = _legend_for_height(full_frame.shape[0])
                    if legend is not None:
                        full_frame = _append_right_panel(full_frame, legend)
                    full_frame, was_padded, shape_in, shape_out = _pad_frame_to_macro_block(
                        full_frame, 16
                    )
                    writer.append_data(full_frame)
            finally:
                for item in renderers:
                    _close_statebar_xy_row_renderer(item["renderer"])
        finally:
            writer.close()
        if bool(verbose):
            elapsed_group = time.perf_counter() - t_group
            print(
                "[example-track-backprojection][mp4] "
                f"Finished {title_fallback} in {elapsed_group:.2f}s"
            )

    if mode in {"per_cluster", "both"}:
        per_cluster_dir.mkdir(parents=True, exist_ok=True)
        cluster_keys, rows_by_cluster = _group_rows_by_cluster(rows_info, cluster_order)
        for cluster_text in cluster_keys:
            cluster_token = _sanitize_filename_token(cluster_text, fallback="cluster")
            video_path = per_cluster_dir / f"example_track_backprojection_cluster_{cluster_token}.mp4"
            _write_group_video(video_path, rows_by_cluster.get(cluster_text, []), f"cluster '{cluster_text}'")
            paths_by_cluster[str(cluster_text)] = str(video_path)

    if mode in {"per_example", "both"}:
        per_example_dir.mkdir(parents=True, exist_ok=True)
        rows_by_rank = _group_rows_by_rank(
            rows_info=rows_info,
            cluster_order=cluster_order,
            num_example_ranks=int(num_example_ranks),
        )
        for rank in range(1, int(num_example_ranks) + 1):
            video_path = per_example_dir / f"example_track_backprojection_example_{int(rank):02d}.mp4"
            _write_group_video(video_path, rows_by_rank.get(rank, []), f"example rank {int(rank):02d}")
            paths_by_example[f"{int(rank):02d}"] = str(video_path)

    out = {
        "video_paths_by_cluster": paths_by_cluster,
        "video_paths_by_example_rank": paths_by_example,
        "n_clusters": int(len(paths_by_cluster)),
        "n_example_ranks": int(len(paths_by_example)),
        "n_tracks_selected": int(len(prep["ranked_df"])),
        "layout_mode": str(mode),
        "global_xy_shape": list(prep["global_xy_shape"]),
        "segment_outline_errors": dict(prep.get("segment_outline_errors", {})),
    }
    if bool(verbose):
        elapsed = time.perf_counter() - t_start
        print(
            "[example-track-backprojection][mp4] "
            f"completed in {elapsed:.2f}s | per_cluster={len(paths_by_cluster)} | "
            f"per_example={len(paths_by_example)}"
        )
    return out
