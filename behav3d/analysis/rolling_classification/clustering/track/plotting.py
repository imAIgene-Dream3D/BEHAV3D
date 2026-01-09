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

    # state colors
    states = pd.Categorical(obs[state_key]).categories
    # cmap = cm.get_cmap(cmap_name, max(1, len(states)))
    cmap = colormaps.get_cmap(cmap_name).resampled(max(1, len(states)))

    color_map = {st: cmap(i) for i, st in enumerate(states)}

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
                               facecolors=color_map.get(cur_state, "grey"))
                start = nxt_start
                cur_state = nxt_state
                cur_end = nxt_end

        ax.broken_barh([(start, cur_end - start)], (y0, bar_h),
                       facecolors=color_map.get(cur_state, "grey"))

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
    rng = np.random.default_rng(seed)

    needed = {sample_key, track_key, tmin_key, tmax_key, cluster_key}
    missing = [c for c in needed if c not in adata_tracks.obs.columns]
    if missing:
        raise ValueError(f"adata_tracks.obs missing required columns: {missing}")

    tracks_df = adata_tracks.obs[[sample_key, track_key, tmin_key, tmax_key, cluster_key]].copy()
    if query is not None:
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
    for cl, df_cl in tracks_df.groupby(cluster_key, sort=False, observed=True):
        k = min(int(n_per_cluster), len(df_cl))
        if k <= 0:
            continue
        idx = rng.choice(len(df_cl), size=k, replace=False)
        chosen_parts.append(df_cl.iloc[idx])

    if not chosen_parts:
        raise ValueError("No tracks were selected (n_per_cluster too small or no clusters available).")

    chosen = pd.concat(chosen_parts, axis=0, ignore_index=True)

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

    legend_handles = None
    if legend:
        if state_key not in adata_full.obs.columns:
            raise ValueError(
                f"adata_full.obs missing required column for legend/state coloring: {state_key}"
            )

        states_series = adata_full.obs[state_key].dropna()

        if pd.api.types.is_categorical_dtype(states_series):
            state_values = list(states_series.cat.categories)
        else:
            uniq = pd.unique(states_series)
            try:
                state_values = sorted(uniq)
            except Exception:
                state_values = list(uniq)

        cmap = colormaps.get_cmap(cmap_name).resampled(max(1, len(state_values)))

        legend_handles = [
            Patch(facecolor=cmap(i), edgecolor="k", label=str(v))
            for i, v in enumerate(state_values)
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
            title=f"{cluster_key} = {cl} (n={len(df_cl)})",
            x_mode=x_mode,
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
