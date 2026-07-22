import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

from behav3d.analysis.behavior.general.visualization.plots.proportion_bars import (
    hash_stable_label_color_map,
)
from behav3d.analysis.behavior.state.utils import _apply_state_order

A4_LANDSCAPE = (11.69, 8.27)


def plot_feature_violin_by_group(
    df,
    feature_col,
    group_col,
    *,
    group_order=None,
    colors=None,
    ylabel=None,
    title=None,
    show_points=False,
    point_alpha=0.4,
    point_size=6,
    mean_marker_size=60,
    figsize=A4_LANDSCAPE,
):
    """One-panel violin of a numeric per-row feature split by a categorical group column.

    Colors follow the same hash-stable convention used for proportion-bar plots
    (``hash_stable_label_color_map``), so a violin of e.g. mean contact fraction per
    ``ClusterID`` uses the same per-cluster colors as the cluster-proportion bars.
    """
    sub = df[[group_col, feature_col]].copy()
    sub[group_col] = sub[group_col].astype(str)
    sub[feature_col] = pd.to_numeric(sub[feature_col], errors="coerce")
    sub = sub.dropna(subset=[group_col, feature_col])

    if group_order is not None:
        resolved_order = [str(g) for g in group_order if str(g) in set(sub[group_col])]
    else:
        resolved_order = sorted(sub[group_col].unique().tolist())
        resolved_order = _apply_state_order(resolved_order, [])

    resolved_colors = hash_stable_label_color_map(resolved_order, colors=colors)

    fig, ax = plt.subplots(figsize=figsize)
    if len(sub) == 0 or len(resolved_order) == 0:
        ax.text(0.5, 0.5, "No data available", ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
        return fig

    palette = [resolved_colors[g] for g in resolved_order]
    sns.violinplot(
        data=sub, x=group_col, y=feature_col, hue=group_col, legend=False,
        order=resolved_order, hue_order=resolved_order, palette=palette, inner=None, cut=0, ax=ax,
    )

    if show_points:
        sns.stripplot(
            data=sub, x=group_col, y=feature_col,
            order=resolved_order, ax=ax,
            color="black", dodge=False, jitter=0.2, alpha=point_alpha, size=point_size,
        )

    means = sub.groupby(group_col, observed=True)[feature_col].mean().reindex(resolved_order)
    ax.scatter(
        np.arange(len(resolved_order)), means.to_numpy(),
        s=mean_marker_size, color="white", edgecolor="black", linewidths=1.0, zorder=3,
    )

    ax.set_xlabel(str(group_col))
    ax.set_ylabel(str(ylabel) if ylabel is not None else str(feature_col))
    if title is not None:
        ax.set_title(str(title))
    fig.tight_layout()
    return fig
