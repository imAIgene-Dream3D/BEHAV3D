from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg

import imageio

from behav3d.io.images import load_zarr
from behav3d.preprocessing import calc_z_projection

def stratified_pick_examples(df_windows, X, seed=0):
    """
    For each ClusterID in `windows_df`, pick X rows while approximating the
    `sample_name` distribution within that cluster.
    """
    rng = np.random.default_rng(seed)
    selections = []

    for cluster_id, sub in df_windows.groupby("ClusterID"):
        counts = sub["sample_name"].value_counts(normalize=True)
        targets = {s: max(1, int(round(p * X))) for s, p in counts.items()}

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
    margin=10,
    pmin=0,
    pmax=99.99,
    mask_margin=False,
):
    sample_name = df_window["sample_name"]
    start_t = int(df_window["window_start_position_t"])
    end_t = int(df_window["window_end_position_t"])

    df_track = df_positions[
        (df_positions["sample_name"] == sample_name)
        & (df_positions["TrackID"] == df_window["TrackID"])
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
            hi = lo + 1e-6
        percentiles[c] = (float(lo), float(hi))

    if mask_margin:
        masked = np.zeros_like(zarr)
        for t in range(start_t, end_t + 1):
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
        x_min:x_max,
    ]

    cropped_img = np.asarray(cropped_img)
    xy_proj = calc_z_projection(cropped_img, z_axis=-3, projection="max")
    xz_proj = calc_z_projection(cropped_img, z_axis=-2, projection="max")
    yz_proj = calc_z_projection(cropped_img, z_axis=-1, projection="max")

    xy_proj_rgb = colorize_channels_to_rgb(xy_proj, percentiles=percentiles)
    xz_proj_rgb = colorize_channels_to_rgb(xz_proj, percentiles=percentiles)
    yz_proj_rgb = colorize_channels_to_rgb(yz_proj, percentiles=percentiles)

    return xy_proj_rgb, xz_proj_rgb, yz_proj_rgb


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
    sample_name = df_window["sample_name"]
    start_t = int(df_window["window_start_position_t"])
    end_t = int(df_window["window_end_position_t"])

    df_track = df_positions[
        (df_positions["sample_name"] == sample_name)
        & (df_positions["TrackID"].astype(int) == int(df_window["TrackID"]))
    ]

    zarr_path = Path(output_folder, "images", sample_name, f"{sample_name}.zarr")
    zarr = load_zarr(zarr_path)
    T, C, Z, Y, X = zarr.shape

    if normalize_per_channel:
        p_img = np.asarray(zarr[-1])
        percentiles = {}
        for c in range(C):
            ch = p_img[c]
            lo = np.percentile(ch, pmin)
            hi = np.percentile(ch, pmax)
            hi_floor = 30000
            hi = max(hi, hi_floor)
            if hi <= lo:
                hi = lo + 1e-6
            percentiles[c] = (float(lo), float(hi))
    else:
        percentiles = None

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

    df_track_window = df_track[
        (df_track["position_t"] >= start_t) & (df_track["position_t"] <= end_t)
    ]
    df_bbox = df_track_window if len(df_track_window) else df_track

    z_min = int(df_bbox["pixel_position_z"].min())
    z_max = int(df_bbox["pixel_position_z"].max())
    y_min = int(df_bbox["pixel_position_y"].min())
    y_max = int(df_bbox["pixel_position_y"].max())
    x_min = int(df_bbox["pixel_position_x"].min())
    x_max = int(df_bbox["pixel_position_x"].max())

    z_min = max(0, int(z_min - margin))
    y_min = max(0, int(y_min - margin))
    x_min = max(0, int(x_min - margin))

    z_max = min(Z, int(z_max + margin + 1))
    y_max = min(Y, int(y_max + margin + 1))
    x_max = min(X, int(x_max + margin + 1))

    cropped_img = masked[
        start_t : end_t + 1,
        :,
        z_min:z_max,
        y_min:y_max,
        x_min:x_max,
    ]
    cropped_img = np.asarray(cropped_img)

    xy_proj = calc_z_projection(cropped_img, z_axis=-3, projection="max")
    xz_proj = calc_z_projection(cropped_img, z_axis=-2, projection="max")
    yz_proj = calc_z_projection(cropped_img, z_axis=-1, projection="max")

    xy_proj_rgb = colorize_channels_to_rgb(xy_proj, percentiles=percentiles)
    xz_proj_rgb = colorize_channels_to_rgb(xz_proj, percentiles=percentiles)
    yz_proj_rgb = colorize_channels_to_rgb(yz_proj, percentiles=percentiles)

    T_window = cropped_img.shape[0]
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

        x_local = pos_x - x_min
        y_local = pos_y - y_min
        z_local = pos_z - z_min

        track_xy[idx, 0] = x_local
        track_xy[idx, 1] = y_local

        track_xz[idx, 0] = x_local
        track_xz[idx, 1] = z_local

        track_yz[idx, 0] = y_local
        track_yz[idx, 1] = z_local

    return xy_proj_rgb, xz_proj_rgb, yz_proj_rgb, track_xy, track_xz, track_yz


def _pca_project_xyz(xyz):
    assert xyz.ndim == 2 and xyz.shape[1] == 3, "xyz must be (N,3)"
    mean = xyz.mean(axis=0)
    X = xyz - mean
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    components = Vt
    proj2 = X @ components[:2].T
    return proj2, mean, components


def _render_frame_row(
    xy_rgb,
    xz_rgb,
    yz_rgb,
    traj2,
    t_idx,
    fig_w,
    fig_h,
    dpi,
    xy_title="XY",
    xz_title="XZ",
    yz_title="YZ",
    pad_frac=0.05,
):
    fig, axes = plt.subplots(
        1, 4, figsize=(fig_w, fig_h), dpi=dpi,
        gridspec_kw={"width_ratios": [1, 1, 1, 1]},
        constrained_layout=True,
    )
    ax_xy, ax_xz, ax_yz, ax_traj = axes

    for ax, img, title in zip(
        [ax_xy, ax_xz, ax_yz],
        [xy_rgb, xz_rgb, yz_rgb],
        [xy_title, xz_title, yz_title],
    ):
        ax.imshow(img)
        h, w = img.shape[:2]
        rect = plt.Rectangle((0, 0), w, h, linewidth=2, edgecolor="white", facecolor="none")
        ax.add_patch(rect)
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])

    ax_traj.set_title("Track (PC1–PC2)")
    if traj2.shape[0] > 0:
        t_clamped = max(0, min(t_idx, traj2.shape[0] - 1))
        ax_min, ax_max = 0.0, 50.0
        center = 0.5 * (ax_min + ax_max)
        traj_shift = traj2 + center
        seg = traj_shift[: t_clamped + 1]

        ax_traj.plot(seg[:, 0], seg[:, 1], linewidth=1.5)
        ax_traj.scatter([seg[-1, 0]], [seg[-1, 1]], s=10)

        ax_traj.set_xlim(ax_min, ax_max)
        ax_traj.set_ylim(ax_min, ax_max)
        ax_traj.set_aspect("equal", adjustable="box")
        ax_traj.set_xticks([])
        ax_traj.set_yticks([])

    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    buf = np.asarray(canvas.buffer_rgba())
    plt.close(fig)
    return buf[..., :3].copy()


def _stack_frame_row(xy_t, xz_t, yz_t, traj2, t_idx, fig_w, fig_h, dpi, pad_frac, labels):
    return _render_frame_row(
        xy_rgb=xy_t,
        xz_rgb=xz_t,
        yz_rgb=yz_t,
        traj2=traj2,
        t_idx=t_idx,
        fig_w=fig_w,
        fig_h=fig_h,
        dpi=dpi,
        xy_title=labels[0],
        xz_title=labels[1],
        yz_title=labels[2],
        pad_frac=pad_frac,
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


def _concat_examples_horizontally(examples):
    heights = [e.shape[0] for e in examples]
    max_h = max(heights)
    out_w = sum(e.shape[1] for e in examples)
    out = np.zeros((max_h, out_w, 3), dtype=np.uint8)

    x = 0
    for e in examples:
        h, w = e.shape[:2]
        out[:h, x : x + w] = e
        x += w
    return out


def _build_traj2_for_window(df_track, t_start, t_end):
    xyz = []
    times = list(range(t_start, t_end + 1))
    last = None
    for t in times:
        row = df_track[df_track["position_t"] == t]
        if len(row) == 0:
            xyz.append([0.0, 0.0, 0.0] if last is None else last)
        else:
            p = [
                float(row["pixel_position_x"].values[0]),
                float(row["pixel_position_y"].values[0]),
                float(row["pixel_position_z"].values[0]),
            ]
            xyz.append(p)
            last = p

    xyz = np.asarray(xyz, dtype=float)
    traj2, mean, comps = _pca_project_xyz(xyz)
    return traj2


def _add_row_title(row_img, title, fontsize=40):
    h, w = row_img.shape[:2]
    dpi = 100
    fig_w = w / dpi
    fig_h = (h + 40) / dpi

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


def _render_fulltrack_frame_with_overlay(
    xy_rgb,
    xz_rgb,
    yz_rgb,
    track_xy,
    track_xz,
    track_yz,
    t_idx,
    fig_w,
    fig_h,
    dpi,
    titles=("XY", "XZ", "YZ"),
    track_color="green",
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
        ax.scatter(coords[-1, 0], coords[-1, 1], s=10, color=track_color)

    ax_xy.imshow(xy_rgb)
    _plot_track(ax_xy, track_xy)
    ax_xy.set_title(titles[0])
    ax_xy.set_xticks([])
    ax_xy.set_yticks([])

    ax_xz.imshow(xz_rgb)
    _plot_track(ax_xz, track_xz)
    ax_xz.set_title(titles[1])
    ax_xz.set_xticks([])
    ax_xz.set_yticks([])

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
    df_windows,
    df_positions,
    output_folder,
    clusters=None,
    out_dir="cluster_videos",
    normalize_per_channel=False,
    fps=12,
    dpi=200,
    margin=(10, 10, 10),
    pmin=0.0,
    pmax=100.0,
    examples_per_cluster=3,
    seed=0,
    figsize_per_row=(12.0, 4.0),
    traj_pad_frac=0.05,
):
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    if clusters is not None:
        dfw = df_windows[df_windows["ClusterID"].isin(clusters)].copy()
    else:
        dfw = df_windows.copy()

    picks = stratified_pick_examples(dfw, X=examples_per_cluster, seed=seed)

    results = {}
    for cluster_id, sub in picks.groupby("ClusterID"):
        print(f"Processing Cluster {cluster_id} with {len(sub)} exemplars...")
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

            xy_stack, xz_stack, yz_stack = create_centered_max_projection_cutout(
                w,
                df_positions,
                output_folder,
                normalize_per_channel=normalize_per_channel,
                margin=margin,
                pmin=pmin,
                pmax=pmax,
            )

            T_w = xy_stack.shape[0]
            max_T = max(max_T, T_w)

            traj2 = _build_traj2_for_window(df_track, t0, t1)

            rows_info.append({
                "label": f"{sample_name} • Track {track_id}",
                "xy": xy_stack,
                "xz": xz_stack,
                "yz": yz_stack,
                "traj2": traj2,
                "T": T_w,
            })

        video_path = out_path / f"cluster_{int(cluster_id)}.mp4"
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
                    frame = _stack_frame_row(
                        xy_t=r["xy"][ti],
                        xz_t=r["xz"][ti],
                        yz_t=r["yz"][ti],
                        traj2=r["traj2"],
                        t_idx=ti,
                        fig_w=figsize_per_row[0],
                        fig_h=figsize_per_row[1],
                        dpi=dpi,
                        pad_frac=traj_pad_frac,
                        labels=("XY", "XZ", "YZ"),
                    )
                    row_imgs.append(frame)

                full_frame = _concat_rows_vertically(row_imgs)
                writer.append_data(full_frame)
        finally:
            writer.close()

        results[int(cluster_id)] = str(video_path)

    return results


def create_cluster_overview_video(
    df_windows,
    df_positions,
    output_folder,
    clusters=None,
    out_dir="cluster_videos",
    normalize_per_channel=False,
    fps=12,
    dpi=200,
    margin=(10, 10, 10),
    pmin=0.0,
    pmax=100.0,
    examples_per_cluster=1,
    seed=0,
    figsize_per_example=(12.0, 4.0),
    traj_pad_frac=0.05,
):
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    video_path = out_path / "clusters_overview.mp4"

    if clusters is not None:
        dfw = df_windows[df_windows["ClusterID"].isin(clusters)].copy()
    else:
        dfw = df_windows.copy()

    picks = stratified_pick_examples(dfw, X=examples_per_cluster, seed=seed)

    clusters_info = {}
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
                w,
                df_positions,
                output_folder,
                normalize_per_channel=normalize_per_channel,
                margin=margin,
                pmin=pmin,
                pmax=pmax,
            )

            T_w = xy_stack.shape[0]
            global_max_T = max(global_max_T, T_w)

            traj2 = _build_traj2_for_window(df_track, t0, t1)

            cluster_rows.append({
                "cluster_id": int(cluster_id),
                "track_id": track_id,
                "xy": xy_stack,
                "xz": xz_stack,
                "yz": yz_stack,
                "traj2": traj2,
                "T": T_w,
            })

        clusters_info[int(cluster_id)] = cluster_rows[:examples_per_cluster]

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
        for t in range(global_max_T):
            cluster_row_imgs = []

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
                        labels=(f'Track {r["track_id"]} • XY', "XZ", "YZ"),
                    )
                    example_imgs.append(frame)

                row_img = _concat_examples_horizontally(example_imgs)
                row_with_title = _add_row_title(row_img, f"Cluster {cluster_id}")
                cluster_row_imgs.append(row_with_title)

            full_frame = _concat_rows_vertically(cluster_row_imgs)
            writer.append_data(full_frame)
    finally:
        writer.close()

    return str(video_path)


def create_fulltrack_cluster_videos(
    df_windows,
    df_positions,
    output_folder,
    clusters=None,
    out_dir="cluster_videos_fulltrack",
    fps=12,
    dpi=200,
    margin=10,
    pmin=0.0,
    pmax=99.99,
    examples_per_cluster=3,
    seed=0,
    figsize_per_row=(9.0, 3.0),
    normalize_per_channel=True,
    mask_margin=False,
    track_color="green",
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

    if clusters is not None:
        dfw = df_windows[df_windows["ClusterID"].isin(clusters)].copy()
    else:
        dfw = df_windows.copy()

    picks = stratified_pick_examples(dfw, X=examples_per_cluster, seed=seed)

    results = {}

    for cluster_id, sub in picks.groupby("ClusterID"):
        print(f"[fulltrack] Processing Cluster {cluster_id} with {len(sub)} exemplars...")
        rows_info = []
        max_T = 0

        for _, w in sub.iterrows():
            sample_name = w["sample_name"]
            track_id = int(w["TrackID"])

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

            rows_info.append({
                "label": f"{sample_name} • Track {track_id}",
                "xy": xy_stack,
                "xz": xz_stack,
                "yz": yz_stack,
                "track_xy": track_xy,
                "track_xz": track_xz,
                "track_yz": track_yz,
                "T": T_w,
            })

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
    