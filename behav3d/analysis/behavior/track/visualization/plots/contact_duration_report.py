"""Standalone contact-duration comparison.

Pulled out of the bundled ``save_track_contact_group_analysis`` report (see ``reports.py``) so it
can be run on its own: for each track, how long it stayed in contact with each touched class of a
chosen target cell type (e.g. a T-cell/organoid track's contacts broken down by which macrophage
morphology class — round / elongated / plastic — it touched), compared:

- pairwise between every two classes, and
- each class vs. every other class pooled together ("rest"),

using either an unpaired Welch's t-test (the default, matching the rest of the app) or a paired
t-test where the pairs are formed by averaging within a user-chosen column (typically
``sample_name``) before pairing — mirroring a manual R workflow of
``group_by(sample_name, class) |> summarise(mean(...)) |> pivot_wider() |> t.test(paired=TRUE)``.
"""
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from behav3d.analysis.behavior.utils import _mixed_label_sort_key
from behav3d.analysis.behavior.general.visualization.plots.proportion_bars import (
    hash_stable_label_color_map,
    welch_ttest_stars,
    SIGNIFICANCE_LEGEND_TEXT,
    _chunk_list,
)
from behav3d.analysis.behavior.track.contact_grouping import (
    compute_track_contact_features,
    merge_track_contact_features_into_obs,
    compute_track_contact_target_class_features,
    _contact_group_col_name,
    _contact_class_max_bout_col_name,
)

_MIN_GROUP_N_FOR_TEST = 2


def _rest_label(other_classes):
    if len(other_classes) <= 3:
        return f"rest ({'+'.join(other_classes)})"
    return f"rest ({len(other_classes)} classes)"


def _build_comparisons(class_order):
    """``[(label_a, classes_a, label_b, classes_b), ...]`` — every pairwise combo between
    individual classes, plus (only once there are >= 3 classes, since with 2 it would just
    duplicate the pairwise entry) one "class vs. every other class pooled" entry per class."""
    comparisons = [(a, [a], b, [b]) for a, b in combinations(class_order, 2)]
    if len(class_order) >= 3:
        for cls in class_order:
            others = [c for c in class_order if c != cls]
            comparisons.append((cls, [cls], _rest_label(others), others))
    return comparisons


def _welch_stats(values_a, values_b):
    values_a = np.asarray(values_a, dtype=float)
    values_b = np.asarray(values_b, dtype=float)
    n_a, n_b = int(len(values_a)), int(len(values_b))
    mean_a = float(np.mean(values_a)) if n_a else float("nan")
    mean_b = float(np.mean(values_b)) if n_b else float("nan")
    if n_a >= _MIN_GROUP_N_FOR_TEST and n_b >= _MIN_GROUP_N_FOR_TEST:
        t_stat, p_value = stats.ttest_ind(values_b, values_a, equal_var=False, nan_policy="omit")
        t_stat, p_value = float(t_stat), float(p_value)
    else:
        t_stat, p_value = float("nan"), float("nan")
    return dict(n_a=n_a, n_b=n_b, mean_a=mean_a, mean_b=mean_b, t_stat=t_stat, p_value=p_value)


def _paired_stats(values_a, values_b):
    """``values_a``/``values_b`` are already the aligned per-pairing-unit mean durations (equal
    length, same pairing-unit order) — the caller does the inner join."""
    n = int(len(values_a))
    mean_a = float(np.mean(values_a)) if n else float("nan")
    mean_b = float(np.mean(values_b)) if n else float("nan")
    if n >= _MIN_GROUP_N_FOR_TEST:
        t_stat, p_value = stats.ttest_rel(values_b, values_a, nan_policy="omit")
        t_stat, p_value = float(t_stat), float(p_value)
    else:
        t_stat, p_value = float("nan"), float("nan")
    return dict(n_a=n, n_b=n, mean_a=mean_a, mean_b=mean_b, t_stat=t_stat, p_value=p_value)


def _comparison_row(duration_df, *, label_a, classes_a, label_b, classes_b, test_mode, pairing_col):
    """Compute stats plus the actual (timepoints) values tested, for one comparison."""
    if test_mode == "paired":
        per_unit = duration_df.groupby([pairing_col, "target_class"])["duration_timepoints"].mean()
        target_level = per_unit.index.get_level_values("target_class")
        side_a = per_unit[target_level.isin(classes_a)].groupby(level=pairing_col).mean()
        side_b = per_unit[target_level.isin(classes_b)].groupby(level=pairing_col).mean()
        joined = pd.concat({"a": side_a, "b": side_b}, axis=1).dropna()
        values_a = joined["a"].to_numpy()
        values_b = joined["b"].to_numpy()
        stats_row = _paired_stats(values_a, values_b)
    else:
        values_a = duration_df.loc[duration_df["target_class"].isin(classes_a), "duration_timepoints"].to_numpy()
        values_b = duration_df.loc[duration_df["target_class"].isin(classes_b), "duration_timepoints"].to_numpy()
        stats_row = _welch_stats(values_a, values_b)

    diff = (
        stats_row["mean_b"] - stats_row["mean_a"]
        if np.isfinite(stats_row["mean_a"]) and np.isfinite(stats_row["mean_b"])
        else float("nan")
    )
    stats_row.update(
        group_a=label_a,
        group_b=label_b,
        diff=diff,
        stars=welch_ttest_stars(stats_row["p_value"]),
        test_mode=test_mode,
        pairing_col=pairing_col if test_mode == "paired" else None,
        values_a=values_a,
        values_b=values_b,
    )
    return stats_row


def _draw_pair_box(ax, values_a, values_b, *, label_a, label_b, ylabel, stars, colors):
    data = [np.asarray(values_a, dtype=float), np.asarray(values_b, dtype=float)]
    bp = ax.boxplot(data, tick_labels=[label_a, label_b], patch_artist=True, widths=0.6, showfliers=False)
    for patch, label in zip(bp["boxes"], (label_a, label_b)):
        patch.set_facecolor(colors.get(label, "#808080"))
        patch.set_alpha(0.55)
    rng = np.random.default_rng(0)
    for i, vals in enumerate(data, start=1):
        if len(vals) == 0:
            continue
        jitter = (rng.random(len(vals)) - 0.5) * 0.15
        ax.scatter(np.full(len(vals), i) + jitter, vals, s=8, color="black", alpha=0.4, zorder=3)
    ax.set_ylabel(ylabel, fontsize=7)
    ax.tick_params(axis="x", labelsize=6.5)
    ax.tick_params(axis="y", labelsize=7)

    finite_vals = [v for arr in data if len(arr) for v in arr[np.isfinite(arr)]]
    if not finite_vals:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes, fontsize=7)
        return
    y_max, y_min = max(finite_vals), min(finite_vals)
    span = max(y_max - y_min, 1e-9)
    if stars and stars != "n.s.":
        y = y_max + span * 0.08
        tick = span * 0.03
        ax.plot([1, 1, 2, 2], [y, y + tick, y + tick, y], lw=0.8, color="black")
        ax.text(1.5, y + tick, stars, ha="center", va="bottom", fontsize=8)
        ax.set_ylim(top=y + tick + span * 0.15)
    elif stars == "n.s.":
        ax.text(0.5, 1.02, "n.s.", transform=ax.transAxes, ha="center", va="bottom", fontsize=6.5, color="#666")


def _plot_duration_comparison_page(
    page_rows, *, contact_col, target_cell_type_label, test_mode, pairing_col, minutes_per_frame,
    ncols, label_colors,
):
    n = len(page_rows)
    nrows = int(np.ceil(n / ncols))
    show_minutes = minutes_per_frame is not None
    fig = plt.figure(figsize=(11.69, 8.27))
    outer = fig.add_gridspec(nrows=max(1, nrows), ncols=max(1, ncols), hspace=1.0, wspace=0.6, top=0.86, bottom=0.10)

    for i, row in enumerate(page_rows):
        r, c = divmod(i, ncols)
        inner = outer[r, c].subgridspec(1, 2 if show_minutes else 1, wspace=0.7)
        ax_tp = fig.add_subplot(inner[0, 0])
        _draw_pair_box(
            ax_tp, row["values_a"], row["values_b"],
            label_a=row["group_a"], label_b=row["group_b"],
            ylabel="Duration (timepoints)", stars=row["stars"], colors=label_colors,
        )
        ax_tp.set_title(f"{row['group_a']} vs {row['group_b']}", fontsize=7)
        if show_minutes:
            ax_min = fig.add_subplot(inner[0, 1])
            _draw_pair_box(
                ax_min,
                np.asarray(row["values_a"], dtype=float) * minutes_per_frame,
                np.asarray(row["values_b"], dtype=float) * minutes_per_frame,
                label_a=row["group_a"], label_b=row["group_b"],
                ylabel="Duration (minutes)", stars=row["stars"], colors=label_colors,
            )

    subtitle = f"contact_col={contact_col}  target={target_cell_type_label}  test={test_mode}"
    if test_mode == "paired":
        subtitle += f"  pairing_col={pairing_col}"
    if not show_minutes:
        subtitle += "  (minutes unavailable — no time metadata)"
    fig.suptitle(
        f"Contact duration comparison — max sustained contact-bout length by {target_cell_type_label} class\n{subtitle}",
        fontsize=10, fontweight="bold",
    )
    fig.text(0.5, 0.02, SIGNIFICANCE_LEGEND_TEXT, ha="center", va="bottom", fontsize=7)
    return fig


def save_track_contact_duration_comparison(
    adata_tracks,
    df_timepoints,
    out_dir,
    *,
    contact_col,
    min_bout_length,
    target_class_lookup,
    touching_col,
    time_varying,
    target_cell_type_label,
    class_order=None,
    class_colors=None,
    test_mode="welch",
    pairing_col=None,
    minutes_per_frame=None,
    sample_col="sample_name",
    groupby_cols=("sample_name", "TrackID"),
    comparisons_per_page=12,
    verbose=False,
):
    """For every class of ``target_cell_type_label`` actually touched, compare how long tracks
    stay in contact with that class — pairwise between classes, and each class vs. every other
    class pooled ("rest") — using an unpaired Welch's t-test or a paired t-test (pairs formed by
    averaging within ``pairing_col``, e.g. ``"sample_name"``, before pairing).

    Requires the per-cell contact-attribution path (``target_class_lookup``/``touching_col`` —
    see ``contact_grouping.build_target_class_lookup_from_state_adata`` /
    ``build_target_class_lookup_from_track_adata``); there is nothing to compare "by class"
    without it.

    Writes one combined PDF (``contact_duration_comparison.pdf``, small boxplot pairs — timepoints
    and minutes side by side when ``minutes_per_frame`` is given — paginated
    ``comparisons_per_page`` per page) plus a CSV with one row per comparison, into the same
    ``{out_dir}/contact_analysis/{contact_col}/`` folder used by ``save_track_contact_group_analysis``.

    Returns a dict of artifact paths plus ``n_comparisons``/``n_pages``/``class_order``.
    """
    test_mode = str(test_mode).strip().lower()
    if test_mode not in ("welch", "paired"):
        raise ValueError(f"test_mode must be 'welch' or 'paired', got {test_mode!r}.")
    if test_mode == "paired" and not pairing_col:
        raise ValueError("pairing_col is required when test_mode='paired'.")

    groupby_cols = [str(c) for c in list(groupby_cols)]
    group_col = _contact_group_col_name(contact_col)
    max_bout_col = _contact_class_max_bout_col_name(contact_col)

    contact_features = compute_track_contact_features(
        df_timepoints, adata_tracks, contact_col=contact_col, min_bout_length=min_bout_length,
        groupby_cols=groupby_cols, verbose=verbose,
    )
    merge_track_contact_features_into_obs(
        adata_tracks, contact_features, contact_col=contact_col, min_bout_length=min_bout_length,
        groupby_cols=groupby_cols,
    )
    long_target_df, _group_df = compute_track_contact_target_class_features(
        df_timepoints, adata_tracks, target_class_lookup,
        contact_col=contact_col, touching_col=touching_col, time_varying=bool(time_varying),
        contact_group_col=group_col, groupby_cols=groupby_cols, verbose=verbose,
    )

    duration_df = long_target_df.reset_index().rename(columns={max_bout_col: "duration_timepoints"})
    duration_df["target_class"] = duration_df["target_class"].astype(str)

    if test_mode == "paired" and pairing_col not in duration_df.columns:
        if pairing_col not in adata_tracks.obs.columns:
            raise KeyError(
                f"pairing_col '{pairing_col}' not found in adata_tracks.obs — merge it in before "
                f"calling this function (e.g. via core.metadata.merge_condition_columns_into_obs)."
            )
        pairing_lookup = (
            adata_tracks.obs[groupby_cols + [pairing_col]].drop_duplicates(subset=groupby_cols)
        )
        duration_df = duration_df.merge(pairing_lookup, on=groupby_cols, how="left")

    touched_classes = sorted(
        duration_df["target_class"].dropna().unique().tolist(), key=_mixed_label_sort_key,
    )
    resolved_class_order = [str(c) for c in class_order] if class_order is not None else touched_classes
    resolved_class_order = [c for c in resolved_class_order if c in touched_classes]

    out_dir = Path(out_dir) / "contact_analysis" / str(contact_col)
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_dir = out_dir / "csv"
    csv_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / "contact_duration_comparison.pdf"
    csv_path = csv_dir / "contact_duration_comparison.csv"

    csv_columns = [
        "group_a", "group_b", "n_a", "n_b", "mean_a_timepoints", "mean_b_timepoints",
        "diff_timepoints", "mean_a_minutes", "mean_b_minutes", "t_stat", "p_value", "stars",
        "test_mode", "pairing_col",
    ]

    if len(resolved_class_order) < 2:
        pd.DataFrame(columns=csv_columns).to_csv(csv_path, index=False)
        with PdfPages(pdf_path) as pdf:
            fig, ax = plt.subplots(figsize=(8, 4))
            ax.text(
                0.5, 0.5,
                f"Fewer than 2 touched '{target_cell_type_label}' classes found — nothing to compare.",
                ha="center", va="center", wrap=True,
            )
            ax.axis("off")
            pdf.savefig(fig)
            plt.close(fig)
        if verbose:
            print(f"Contact duration comparison: fewer than 2 touched classes, nothing to compare ({pdf_path}).")
        return {
            "contact_col": str(contact_col), "pdf_path": str(pdf_path), "csv_path": str(csv_path),
            "csv_dir": str(csv_dir), "n_comparisons": 0, "n_pages": 1, "test_mode": test_mode,
            "pairing_col": pairing_col, "class_order": resolved_class_order,
            "minutes_per_frame": minutes_per_frame,
        }

    comparisons = _build_comparisons(resolved_class_order)
    rows = [
        _comparison_row(
            duration_df, label_a=label_a, classes_a=classes_a, label_b=label_b, classes_b=classes_b,
            test_mode=test_mode, pairing_col=pairing_col,
        )
        for label_a, classes_a, label_b, classes_b in comparisons
    ]

    all_labels = sorted({row["group_a"] for row in rows} | {row["group_b"] for row in rows})
    label_colors = hash_stable_label_color_map(all_labels, colors=class_colors)

    comparisons_per_page = max(1, int(comparisons_per_page))
    ncols = min(4, max(1, int(np.ceil(np.sqrt(comparisons_per_page)))))
    pages = _chunk_list(rows, comparisons_per_page)
    with PdfPages(pdf_path) as pdf:
        for page_rows in pages:
            fig = _plot_duration_comparison_page(
                page_rows, contact_col=contact_col, target_cell_type_label=target_cell_type_label,
                test_mode=test_mode, pairing_col=pairing_col, minutes_per_frame=minutes_per_frame,
                ncols=ncols, label_colors=label_colors,
            )
            pdf.savefig(fig)
            plt.close(fig)

    csv_rows = []
    for row in rows:
        m_a = row["mean_a"] * minutes_per_frame if minutes_per_frame is not None and np.isfinite(row["mean_a"]) else float("nan")
        m_b = row["mean_b"] * minutes_per_frame if minutes_per_frame is not None and np.isfinite(row["mean_b"]) else float("nan")
        csv_rows.append({
            "group_a": row["group_a"], "group_b": row["group_b"],
            "n_a": row["n_a"], "n_b": row["n_b"],
            "mean_a_timepoints": row["mean_a"], "mean_b_timepoints": row["mean_b"],
            "diff_timepoints": row["diff"], "mean_a_minutes": m_a, "mean_b_minutes": m_b,
            "t_stat": row["t_stat"], "p_value": row["p_value"], "stars": row["stars"],
            "test_mode": row["test_mode"], "pairing_col": row["pairing_col"],
        })
    pd.DataFrame(csv_rows, columns=csv_columns).to_csv(csv_path, index=False)

    if verbose:
        print(
            f"Saved contact duration comparison ({len(rows)} comparisons, {len(pages)} page(s), "
            f"test_mode={test_mode}): {pdf_path}"
        )

    return {
        "contact_col": str(contact_col),
        "min_bout_length": int(min_bout_length),
        "target_cell_type_label": str(target_cell_type_label),
        "pdf_path": str(pdf_path),
        "csv_path": str(csv_path),
        "csv_dir": str(csv_dir),
        "n_comparisons": len(rows),
        "n_pages": len(pages),
        "test_mode": test_mode,
        "pairing_col": pairing_col,
        "class_order": resolved_class_order,
        "minutes_per_frame": minutes_per_frame,
    }
