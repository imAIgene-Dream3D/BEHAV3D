import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

def plot_exemplar_track_bars(
    adata,
    n_tracks=10,
    sample_key="sample_name",
    track_key="TrackID",
    time_key="position_t",
    state_key="ClusterID",
    seed=0,
    selection="random",   # "random" or "first_time"
    cmap_name="tab20",
    ax=None,
):
    """
    Plot exemplar tracks as horizontal bars colored by state over time.

    Uniqueness / selection is by (sample_key, track_key).
    """

    obs = adata.obs[[sample_key, track_key, time_key, state_key]].copy()

    # Stable categorical mapping for states -> colors
    obs[state_key] = obs[state_key].astype("category")
    states = list(obs[state_key].cat.categories)

    # Unique track identifiers are (sample, track)
    keys_df = obs[[sample_key, track_key]].drop_duplicates().reset_index(drop=True)

    rng = np.random.default_rng(seed)

    if selection == "random":
        idx = rng.choice(len(keys_df), size=min(n_tracks, len(keys_df)), replace=False)
        chosen = keys_df.iloc[idx].copy()
    elif selection == "first_time":
        subset_idx = rng.choice(len(keys_df), size=min(n_tracks, len(keys_df)), replace=False)
        subset = keys_df.iloc[subset_idx].copy()

        # sort chosen subset by each track's first timepoint
        tmp = obs.merge(subset, on=[sample_key, track_key], how="inner")
        starts = (
            tmp.groupby([sample_key, track_key])[time_key]
            .min()
            .sort_values()
            .reset_index()
        )
        chosen = starts[[sample_key, track_key]].copy()
    else:
        raise ValueError("selection must be 'random' or 'first_time'")

    # Color mapping
    cmap = plt.get_cmap(cmap_name, max(len(states), 1))
    state_to_color = {s: cmap(i) for i, s in enumerate(states)}

    # Create axis
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 0.45 * len(chosen) + 2))
    else:
        fig = ax.figure

    # Plot each (sample, track) as a row
    ylabels = []
    for row_i, (sname, tid) in enumerate(zip(chosen[sample_key], chosen[track_key])):
        df = obs.loc[
            (obs[sample_key] == sname) & (obs[track_key] == tid),
            [time_key, state_key],
        ].sort_values(time_key)

        t = df[time_key].to_numpy()
        st = df[state_key].to_numpy()

        if len(t) == 0:
            continue

        # boundaries where state changes
        if len(st) == 1:
            # single point -> tiny segment
            start = float(t[0])
            width = 1.0
            ax.broken_barh([(start, width)], (row_i - 0.4, 0.8),
                           facecolors=state_to_color[df[state_key].iloc[0]],
                           edgecolors="none")
        else:
            change_idx = np.flatnonzero(st[1:] != st[:-1]) + 1
            boundaries = np.r_[0, change_idx, len(st)]

            for a, b in zip(boundaries[:-1], boundaries[1:]):
                state = st[a]
                start = float(t[a])

                if b < len(t):
                    end = float(t[b])
                else:
                    dt = np.median(np.diff(t))
                    end = float(t[-1] + (dt if np.isfinite(dt) and dt > 0 else 1.0))

                ax.broken_barh([(start, max(0.0, end - start))],
                               (row_i - 0.4, 0.8),
                               facecolors=state_to_color[state],
                               edgecolors="none")

        ylabels.append(f"{sname} | {tid}")

    ax.set_yticks(range(len(ylabels)))
    ax.set_yticklabels(ylabels)
    ax.set_xlabel(time_key)
    ax.set_ylabel(f"{sample_key} | {track_key}")
    ax.set_title(f"Exemplar tracks (seed={seed})")

    # Legend
    handles = [mpl.patches.Patch(color=state_to_color[s], label=str(s)) for s in states]
    ax.legend(handles=handles, title=state_key, bbox_to_anchor=(1.02, 1), loc="upper left")

    fig.tight_layout()
    return ax

def compare_cluster_distribution(df, col_a, col_b):
    counts = pd.crosstab(df[col_a], df[col_b])
    props = counts.div(counts.sum(axis=1), axis=0)

    # Majority mapping (for each cluster in A, which B is most common?)
    maj_target = counts.idxmax(axis=1)
    maj_count = counts.max(axis=1)
    purity = (maj_count / counts.sum(axis=1)).fillna(0.0)
    mapping_summary = pd.DataFrame({
        f'{col_a}': counts.index,
        f'major_{col_b}': maj_target.values,
        'major_count': maj_count.values,
        'purity': purity.values
    }).sort_values('purity', ascending=False).reset_index(drop=True)

    # Plot heatmap (matplotlib, no external styling)
    fig, ax = plt.subplots(figsize=(max(6, props.shape[1]*0.8), max(4, props.shape[0]*0.6)))
    im = ax.imshow(props.values, aspect='auto')
    ax.set_xticks(np.arange(props.shape[1]))
    ax.set_xticklabels(props.columns, rotation=45, ha='right')
    ax.set_yticks(np.arange(props.shape[0]))
    ax.set_yticklabels(props.index)
    ax.set_xlabel(col_b)
    ax.set_title(f'Proportions of {col_b} within each {col_a} (row-normalized)')

    cbar = plt.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('proportion', rotation=270, labelpad=15)

    # Optional: annotate cells
    for i in range(props.shape[0]):
        for j in range(props.shape[1]):
            val = props.values[i, j]
            text = f'{val:.2f}'
            ax.text(j, i, text, ha='center', va='center', fontsize=8)

    plt.tight_layout()
    plt.show()
    
