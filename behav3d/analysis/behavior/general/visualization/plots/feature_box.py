import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from scipy import stats

from behav3d.analysis.behavior.general.visualization.plots.proportion_bars import (
    hash_stable_label_color_map,
    welch_ttest_stars,
    SIGNIFICANCE_LEGEND_TEXT,
)
from behav3d.analysis.behavior.state.utils import _apply_state_order

A4_LANDSCAPE = (11.69, 8.27)

_MIN_GROUP_N_FOR_TEST = 3
_ALPHA = 0.05


def _pairwise_mannwhitney(sub, feature_col, group_col, order):
    """Two-sided Mann-Whitney U test for every pair of ``order`` groups with enough data.

    Returns a list of ``(group_a, group_b, p_value, stars)``. Pairs where either side has
    fewer than ``_MIN_GROUP_N_FOR_TEST`` values get ``p_value=nan`` (not tested, not reported).
    """
    values = {g: pd.to_numeric(sub.loc[sub[group_col] == g, feature_col], errors="coerce").dropna().to_numpy() for g in order}
    rows = []
    for i in range(len(order)):
        for j in range(i + 1, len(order)):
            group_a, group_b = order[i], order[j]
            vals_a, vals_b = values[group_a], values[group_b]
            if len(vals_a) >= _MIN_GROUP_N_FOR_TEST and len(vals_b) >= _MIN_GROUP_N_FOR_TEST:
                _, p_value = stats.mannwhitneyu(vals_a, vals_b, alternative="two-sided")
                p_value = float(p_value)
            else:
                p_value = float("nan")
            rows.append((group_a, group_b, p_value, welch_ttest_stars(p_value)))
    return rows


def _draw_significance_brackets(ax, order, sig_pairs, y_max, data_range):
    """Draw one stacked bracket + star label per significant pair, closest pairs lowest."""
    positions = {g: i for i, g in enumerate(order)}
    step = max(data_range * 0.07, 1e-9)
    ordered_pairs = sorted(
        sig_pairs,
        key=lambda row: (abs(positions[row[1]] - positions[row[0]]), min(positions[row[0]], positions[row[1]])),
    )
    placed = []  # (x1, x2, y) of brackets already drawn
    top = y_max
    for group_a, group_b, _p_value, stars in ordered_pairs:
        x1, x2 = sorted((positions[group_a], positions[group_b]))
        y = y_max + step
        for px1, px2, py in placed:
            if x1 <= px2 and x2 >= px1:
                y = max(y, py + step)
        tick = step * 0.15
        ax.plot([x1, x1, x2, x2], [y, y + tick, y + tick, y], lw=1.0, color="black")
        ax.text((x1 + x2) / 2, y + tick, stars, ha="center", va="bottom", fontsize=9)
        placed.append((x1, x2, y + tick))
        top = max(top, y + tick)
    ax.set_ylim(top=top + step)


def _format_ns_footnote(ns_pairs):
    if not ns_pairs:
        return None
    listed = ", ".join(f"{group_a} vs {group_b}" for group_a, group_b, _p, _s in ns_pairs)
    return f"Mann-Whitney U, not significant (p≥0.05): {listed}"


def plot_feature_box_by_group(
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
    """One-panel boxplot of a numeric per-row feature split by a categorical group column, with
    pairwise Mann-Whitney U significance testing between every pair of groups.

    Colors follow the same hash-stable convention used for proportion-bar plots
    (``hash_stable_label_color_map``), so a box of e.g. mean contact fraction per
    ``ClusterID`` uses the same per-cluster colors as the cluster-proportion bars.

    Only pairs that test significant (p<0.05) get a bracket drawn above the boxes, to keep the
    plot readable when there are many groups. Non-significant pairs (and skipped brackets) are
    listed in a footnote below the plot instead of being drawn.
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
    sns.boxplot(
        data=sub, x=group_col, y=feature_col, hue=group_col, legend=False,
        order=resolved_order, hue_order=resolved_order, palette=palette, ax=ax,
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

    below_lines = []
    if len(resolved_order) >= 2:
        pairwise = _pairwise_mannwhitney(sub, feature_col, group_col, resolved_order)
        tested_pairs = [row for row in pairwise if np.isfinite(row[2])]
        sig_pairs = [row for row in tested_pairs if row[2] < _ALPHA]
        ns_pairs = [row for row in tested_pairs if row[2] >= _ALPHA]
        if sig_pairs:
            y_max = float(np.nanmax(sub[feature_col].to_numpy()))
            y_min = float(np.nanmin(sub[feature_col].to_numpy()))
            _draw_significance_brackets(ax, resolved_order, sig_pairs, y_max, max(y_max - y_min, 1e-9))
        if tested_pairs:
            below_lines.append(SIGNIFICANCE_LEGEND_TEXT)
        footnote = _format_ns_footnote(ns_pairs)
        if footnote:
            below_lines.append(footnote)

    if below_lines:
        fig.tight_layout(rect=(0.0, 0.03 + 0.03 * len(below_lines), 1.0, 1.0))
        fig.text(0.5, 0.005, "\n".join(below_lines), ha="center", va="bottom", fontsize=7, wrap=True)
    else:
        fig.tight_layout()
    return fig
