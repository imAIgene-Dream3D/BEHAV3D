from pathlib import Path

import imageio.v2 as imageio
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib import colormaps
from matplotlib.collections import LineCollection
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_pdf import PdfPages

from behav3d.analysis.clustering.state.visualization.videos.track_max_projection import (
    create_fulltrack_max_projection_stacks_with_track,
)


def _mixed_label_sort_key(value):
    text = str(value)
    return (0, int(text)) if text.isdigit() else (1, text)


def _sanitize_filename_token(value, fallback="plot"):
    token = str(value).strip()
    if token == "":
        token = str(fallback)
    token = "".join(ch if (ch.isalnum() or ch in "._-") else "_" for ch in token)
    token = token.strip("._-")
    return token if token != "" else str(fallback)


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
        ax.set_title(title)


def _project_xyz_to_pc12(xyz):
    if xyz.ndim != 2 or xyz.shape[1] != 3:
        raise ValueError(f"xyz must be shaped (N,3), got {xyz.shape}.")
    if xyz.shape[0] < 2:
        raise ValueError("Need at least 2 timepoints to project track trajectory.")
    centered = xyz - xyz.mean(axis=0, keepdims=True)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    proj2 = centered @ vt[:2].T
    return proj2


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
):
    xyz = track_df[[x_col, y_col, z_col]].to_numpy(dtype=float, copy=False)
    proj2 = _project_xyz_to_pc12(xyz)

    if proj2.shape[0] < 2:
        raise ValueError("Track projection needs at least 2 points.")

    segs = np.stack([proj2[:-1], proj2[1:]], axis=1)
    st = track_df[state_key].astype(str).to_numpy()
    seg_colors = [state_color_map.get(str(v), "grey") for v in st[:-1]]
    lc = LineCollection(segs, colors=seg_colors, linewidths=2.0, alpha=0.95)
    ax.add_collection(lc)
    ax.scatter([proj2[-1, 0]], [proj2[-1, 1]], s=18, color="black", zorder=3)

    xmin, ymin = np.min(proj2, axis=0)
    xmax, ymax = np.max(proj2, axis=0)
    lim = max(abs(xmin), abs(xmax), abs(ymin), abs(ymax))
    lim = 1.0 if lim <= 0 else lim
    pad = 0.08 * lim
    ax.set_xlim(-lim - pad, lim + pad)
    ax.set_ylim(-lim - pad, lim + pad)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    if title is not None:
        ax.set_title(title)


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
):
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
    out_dir.mkdir(parents=True, exist_ok=True)
    _, state_color_map = _build_state_color_map(
        adata_full=adata_full,
        state_key=state_key,
        cmap_name=cmap_name,
    )

    cluster_vals = list(chosen_df[cluster_key].dropna().unique())
    cluster_vals = sorted(cluster_vals, key=_mixed_label_sort_key)
    paths = {}

    cluster_key_token = _sanitize_filename_token(cluster_key, fallback="cluster")
    state_key_token = _sanitize_filename_token(state_key, fallback="state")
    rows_per_page = max(1, int(rows_per_page))

    for cluster_val in cluster_vals:
        cluster_token = _sanitize_filename_token(cluster_val, fallback="cluster")
        pdf_path = out_dir / (
            f"exemplar_statebar_track_cluster_{cluster_token}_"
            f"cluster_{cluster_key_token}_state_{state_key_token}.pdf"
        )

        sub = chosen_df.loc[chosen_df[cluster_key] == cluster_val].reset_index(drop=True)
        if len(sub) == 0:
            continue

        with PdfPages(pdf_path) as pdf:
            for start in range(0, len(sub), rows_per_page):
                page = sub.iloc[start : start + rows_per_page].reset_index(drop=True)
                n_rows = len(page)
                fig, axes = plt.subplots(
                    n_rows,
                    2,
                    figsize=(14, max(3.0, 2.8 * n_rows)),
                    squeeze=False,
                    gridspec_kw={"width_ratios": [1.2, 1.0]},
                    constrained_layout=True,
                )

                for row_i, row in page.iterrows():
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
                    _plot_statebar_segments_on_ax(
                        axes[row_i, 0],
                        segments=segments,
                        xlim=xlim,
                        title=f"Cluster {cluster_val} | {sample_name} | Track {track_id}",
                    )
                    _plot_track_state_colored_on_ax(
                        axes[row_i, 1],
                        track_df=track_df,
                        state_key=state_key,
                        state_color_map=state_color_map,
                        x_col=traj_x_col,
                        y_col=traj_y_col,
                        z_col=traj_z_col,
                        title="Track trajectory (PC1-PC2), colored by state",
                    )

                pdf.savefig(fig, dpi=int(plot_dpi), bbox_inches="tight")
                plt.close(fig)

        paths[str(cluster_val)] = str(pdf_path)

    return {
        "pdf_paths_by_cluster": paths,
        "n_clusters": int(len(paths)),
        "n_tracks_selected": int(len(chosen_df)),
    }


def _resolve_sample_zarr_path(output_dir, sample_name):
    expected = Path(output_dir, "images", str(sample_name), f"{sample_name}.zarr")
    if not expected.exists():
        raise ValueError(
            "Missing sample raw zarr for exemplar backprojection video: "
            f"sample='{sample_name}', expected_path='{expected}'."
        )
    return expected


def _render_statebar_xy_frame(
    *,
    xy_rgb,
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
        title=f"Cluster {cluster_value} | {sample_name} | Track {track_id}",
        cursor_x=cursor_x,
    )

    ax_xy.imshow(xy_rgb)
    t_clamped = max(0, min(int(t_idx), int(track_xy.shape[0]) - 1))
    path = track_xy[: t_clamped + 1]
    valid = ~np.isnan(path[:, 0])
    path = path[valid]
    if path.shape[0] >= 2:
        ax_xy.plot(path[:, 0], path[:, 1], linewidth=1.8, color=track_color)
        ax_xy.scatter([path[-1, 0]], [path[-1, 1]], s=14, color=track_color)
    elif path.shape[0] == 1:
        ax_xy.scatter([path[-1, 0]], [path[-1, 1]], s=14, color=track_color)

    ax_xy.set_title("Raw XY + track overlay")
    ax_xy.set_xticks([])
    ax_xy.set_yticks([])

    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    buf = np.asarray(canvas.buffer_rgba())[..., :3].copy()
    plt.close(fig)
    return buf


def save_exemplar_statebar_backprojection_video_per_cluster(
    adata_full,
    *,
    output_dir,
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
    coordinate_source_hint=None,
    seed=0,
    cmap_name="tab20",
):
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
    pixel_triplet = ["pixel_position_x", "pixel_position_y", "pixel_position_z"]
    missing_pixel = [c for c in pixel_triplet if c not in adata_full.obs.columns]
    if len(missing_pixel) > 0:
        src_hint = (
            f" coordinate_source='{coordinate_source_hint}'."
            if (coordinate_source_hint is not None and str(coordinate_source_hint).strip() != "")
            else ""
        )
        raise ValueError(
            "Exemplar backprojection video requires pixel coordinates "
            "['pixel_position_x','pixel_position_y','pixel_position_z'] in adata_full.obs "
            "to align tracks with raw image voxels. "
            f"Missing columns: {missing_pixel}.{src_hint}"
        )

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    _, state_color_map = _build_state_color_map(
        adata_full=adata_full,
        state_key=state_key,
        cmap_name=cmap_name,
    )

    df_positions = adata_full.obs[
        [sample_key, track_key, time_key, "pixel_position_x", "pixel_position_y", "pixel_position_z"]
    ].copy()

    cluster_vals = list(chosen_df[cluster_key].dropna().unique())
    cluster_vals = sorted(cluster_vals, key=_mixed_label_sort_key)
    paths = {}

    cluster_key_token = _sanitize_filename_token(cluster_key, fallback="cluster")
    state_key_token = _sanitize_filename_token(state_key, fallback="state")

    for cluster_val in cluster_vals:
        sub = chosen_df.loc[chosen_df[cluster_key] == cluster_val].reset_index(drop=True)
        if len(sub) == 0:
            continue

        rows_info = []
        max_t = 0
        for _, row in sub.iterrows():
            sample_name = row[sample_key]
            track_id = row[track_key]
            t0 = int(row[tmin_key])
            t1 = int(row[tmax_key])
            if t1 < t0:
                raise ValueError(
                    f"Invalid exemplar window for sample='{sample_name}', track='{track_id}': tmin={t0}, tmax={t1}."
                )

            _resolve_sample_zarr_path(output_dir=output_dir, sample_name=sample_name)

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
                    "Exemplar backprojection video requires non-null pixel coordinates for each exemplar row. "
                    f"Found missing pixel coordinates for cluster='{cluster_val}', "
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
                xy_stack, _, _, track_xy, _, _ = create_fulltrack_max_projection_stacks_with_track(
                    df_window=window_row,
                    df_positions=df_positions,
                    output_folder=output_dir,
                    margin=margin,
                    pmin=pmin,
                    pmax=pmax,
                    mask_margin=False,
                    normalize_per_channel=True,
                )
            except Exception as exc:
                raise ValueError(
                    "Failed to create exemplar backprojection stack for "
                    f"cluster='{cluster_val}', sample='{sample_name}', track='{track_id}': {exc}"
                ) from exc

            if xy_stack.ndim != 4 or int(xy_stack.shape[0]) == 0:
                raise ValueError(
                    "Empty XY stack for exemplar "
                    f"cluster='{cluster_val}', sample='{sample_name}', track='{track_id}'."
                )

            t_len = int(xy_stack.shape[0])
            max_t = max(max_t, t_len)
            rows_info.append(
                {
                    "sample_name": sample_name,
                    "track_id": track_id,
                    "xy_stack": xy_stack,
                    "track_xy": track_xy,
                    "segments": segments,
                    "xlim": xlim,
                    "t0": int(t0),
                    "T": t_len,
                }
            )

        cluster_token = _sanitize_filename_token(cluster_val, fallback="cluster")
        video_path = out_dir / (
            f"exemplar_statebar_backprojection_cluster_{cluster_token}_"
            f"cluster_{cluster_key_token}_state_{state_key_token}.mp4"
        )
        if video_path.exists():
            video_path.unlink()

        writer = imageio.get_writer(
            video_path,
            fps=int(fps),
            codec="libx264",
            ffmpeg_log_level="warning",
            quality=9,
        )
        try:
            for t_idx in range(max_t):
                row_frames = []
                for row_info in rows_info:
                    ti = min(int(t_idx), int(row_info["T"]) - 1)
                    frame = _render_statebar_xy_frame(
                        xy_rgb=row_info["xy_stack"][ti],
                        track_xy=row_info["track_xy"],
                        t_idx=ti,
                        segments=row_info["segments"],
                        xlim=row_info["xlim"],
                        cursor_x=float(row_info["t0"] + ti),
                        sample_name=row_info["sample_name"],
                        track_id=row_info["track_id"],
                        cluster_value=cluster_val,
                        track_color=track_color,
                        dpi=dpi,
                    )
                    row_frames.append(frame)
                writer.append_data(_concat_rows_vertically(row_frames))
        finally:
            writer.close()

        paths[str(cluster_val)] = str(video_path)

    return {
        "video_paths_by_cluster": paths,
        "n_clusters": int(len(paths)),
        "n_tracks_selected": int(len(chosen_df)),
    }
