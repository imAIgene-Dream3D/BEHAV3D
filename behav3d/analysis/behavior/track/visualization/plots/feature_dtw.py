from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler

from behav3d.analysis.filtering import round_legend_ticks


def plot_cluster_percentage_bars(
    df_clust_perc: pd.DataFrame,
    outprefix,
    group_cols=None,
    plot_results=True,
    grid_figsize=(20, 10),
    ncols_samples=3,
    sample_panel_width=6.0,
    sample_row_height=2.8,
    title_fontsize=24,
    label_fontsize=12,
    info_fontsize=10,
):
    """
    Produce a PDF with cluster percentage visualizations.

    Page 1: 'Combined' view (samples pooled). Grid layout: 1 or 2 user-selected
            `group_cols` (e.g. organoid_line_condition, tcell_line_condition)
            arranged as a true 2D grid on one page — the column with more unique
            values becomes rows, the other becomes columns. If only 1 column is
            given, panels are laid out as a single column of rows. If more than
            2 columns are given, only the first 2 are used. If group_cols is
            empty/None, this page is skipped.

    Last Page: 'Per-sample' view. Panels laid out in 3 columns (configurable).
               Each panel is a single horizontal stacked bar of ClusterID percentages
               for one sample_name.
    """
    outprefix = Path(outprefix)
    outpdf = outprefix.with_suffix(".pdf")

    if group_cols is None:
        group_cols = []
    elif isinstance(group_cols, str):
        group_cols = [group_cols]

    valid_group_cols = [c for c in group_cols if c in df_clust_perc.columns]
    if len(valid_group_cols) > 2:
        print(f"  Note: {len(valid_group_cols)} group_cols given — using only the first 2 for the grid.")
        valid_group_cols = valid_group_cols[:2]

    grid_config = None
    if valid_group_cols:
        pooled_group_cols = valid_group_cols + ["ClusterID"]

        pooled = (
            df_clust_perc
            .groupby(pooled_group_cols, as_index=False, observed=True)["count"]
            .sum()
        )
        pooled_combo_totals = (
            pooled.groupby(valid_group_cols, observed=True)["count"]
            .sum()
            .reset_index(name="combo_total_pooled")
        )
        pooled = pooled.merge(pooled_combo_totals, on=valid_group_cols, how="left")
        pooled["percentage_pooled"] = pooled["count"] / pooled["combo_total_pooled"]

        if len(valid_group_cols) == 2:
            row_col, col_col = sorted(valid_group_cols, key=lambda c: -df_clust_perc[c].nunique())
            row_values = df_clust_perc[row_col].drop_duplicates().sort_values().tolist()
            col_values = df_clust_perc[col_col].drop_duplicates().sort_values().tolist()
        else:
            row_col = valid_group_cols[0]
            col_col = None
            row_values = df_clust_perc[row_col].drop_duplicates().sort_values().tolist()
            col_values = [None]
        grid_config = {
            "row_col": row_col,
            "row_values": row_values,
            "col_col": col_col,
            "col_values": col_values,
            "nrows": len(row_values),
            "ncols": len(col_values),
        }

    per_sample = (
        df_clust_perc
        .groupby(["sample_name", "ClusterID"], as_index=False, observed=True)["count"]
        .sum()
    )
    per_sample_totals = (
        per_sample.groupby("sample_name", observed=True)["count"]
        .sum()
        .reset_index(name="sample_total_all")
    )
    per_sample = per_sample.merge(per_sample_totals, on="sample_name", how="left")
    per_sample["percentage_overall"] = per_sample["count"] / per_sample["sample_total_all"]
    samples = per_sample["sample_name"].drop_duplicates().tolist()

    cluster_order = (
        df_clust_perc["ClusterID"]
        .drop_duplicates()
        .sort_values()
        .tolist()
    )

    with PdfPages(outpdf) as pdf:
        if grid_config is not None:
            row_col = grid_config["row_col"]
            row_values = grid_config["row_values"]
            col_col = grid_config["col_col"]
            col_values = grid_config["col_values"]
            nrows = grid_config["nrows"]
            ncols = grid_config["ncols"]

            fig1, axes1 = plt.subplots(
                nrows, ncols,
                figsize=grid_figsize,
                sharex=True, sharey=True,
                squeeze=False,
            )

            legend_handles_1, legend_labels_1 = None, None

            for i, row_val in enumerate(row_values):
                for j, col_val in enumerate(col_values):
                    ax = axes1[i, j]
                    if col_col is not None:
                        subset = pooled[
                            (pooled[row_col] == row_val) &
                            (pooled[col_col] == col_val)
                        ]
                    else:
                        subset = pooled[pooled[row_col] == row_val]

                    if i == 0 and col_col is not None:
                        ax.set_title(f"{col_val}", fontsize=title_fontsize)
                    if j == 0:
                        ax.text(
                            -0.05, 0.5, f"{row_val}",
                            ha="right", va="center",
                            rotation=0, transform=ax.transAxes,
                            fontsize=title_fontsize,
                        )

                    if subset.empty:
                        ax.axis("off")
                        continue

                    subset_unique = subset.drop_duplicates(subset=["ClusterID"], keep="first")
                    row = (
                        subset_unique
                        .set_index("ClusterID")
                        .reindex(cluster_order)["percentage_pooled"]
                        .fillna(0.0)
                    )
                    pivot = row.to_frame().T
                    pivot.index = ["combined"]

                    bar_ax = pivot.plot(kind="barh", stacked=True, ax=ax, legend=False)

                    if legend_handles_1 is None:
                        legend_handles_1, legend_labels_1 = bar_ax.get_legend_handles_labels()

                    num_cells = int(subset["combo_total_pooled"].iloc[0]) if not subset.empty else 0
                    ax.text(
                        0.5, 0.08, f"# Cells (pooled): {num_cells}",
                        ha="center", va="center",
                        transform=ax.transAxes, fontsize=info_fontsize,
                    )
                    ax.axis("off")

            row_col_name = row_col.replace("_line_condition", "").replace("_", " ").title()
            if col_col is not None:
                col_col_name = col_col.replace("_line_condition", "").replace("_", " ").title()
                fig1.suptitle(f"Cluster percentages: {row_col_name} × {col_col_name}", fontsize=title_fontsize, y=0.995)
            else:
                fig1.suptitle(f"Cluster percentages: {row_col_name}", fontsize=title_fontsize, y=0.995)

            if legend_handles_1 and legend_labels_1:
                fig1.legend(
                    legend_handles_1, legend_labels_1,
                    fontsize=label_fontsize,
                    title="ClusterID", title_fontsize=label_fontsize,
                    bbox_to_anchor=(0.92, 0.5), loc="center left",
                )
                fig1.tight_layout(rect=[0, 0, 0.88, 0.96])
            else:
                fig1.tight_layout(rect=[0, 0, 1, 0.96])

            if plot_results:
                plt.show()
            pdf.savefig(fig1, bbox_inches="tight")
            plt.close(fig1)

        n_panels = len(samples)
        ncols = ncols_samples
        nrows = int(np.ceil(n_panels / ncols))
        fig2, axes2 = plt.subplots(
            nrows, ncols,
            figsize=(sample_panel_width * ncols, sample_row_height * nrows),
            squeeze=False,
        )

        legend_handles_2, legend_labels_2 = None, None

        for k, sample in enumerate(samples):
            i, j = divmod(k, ncols)
            ax = axes2[i, j]
            sub = per_sample[per_sample["sample_name"] == sample]
            if sub.empty:
                ax.axis("off")
                continue

            sub_unique = sub.drop_duplicates(subset=["ClusterID"], keep="first")
            row = (
                sub_unique
                .set_index("ClusterID")
                .reindex(cluster_order)["percentage_overall"]
                .fillna(0.0)
            )
            pivot = row.to_frame().T
            pivot.index = [sample]

            bar_ax = pivot.plot(kind="barh", stacked=True, ax=ax, legend=False)

            if legend_handles_2 is None:
                legend_handles_2, legend_labels_2 = bar_ax.get_legend_handles_labels()

            num_cells = int(sub["sample_total_all"].iloc[0]) if not sub.empty else 0
            ax.set_title(sample, fontsize=label_fontsize)
            ax.text(
                0.5, 0.08, f"# Cells: {num_cells}",
                ha="center", va="center",
                transform=ax.transAxes, fontsize=info_fontsize,
            )
            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.set_yticks([])
            ax.set_xticks([])
            ax.spines[["top", "right", "bottom", "left"]].set_visible(False)

        for k in range(n_panels, nrows * ncols):
            i, j = divmod(k, ncols)
            axes2[i, j].axis("off")

        fig2.suptitle("Per-sample cluster composition", fontsize=title_fontsize, y=0.995)

        if legend_handles_2 and legend_labels_2:
            fig2.legend(
                legend_handles_2, legend_labels_2,
                fontsize=label_fontsize,
                title="ClusterID", title_fontsize=label_fontsize,
                bbox_to_anchor=(0.92, 0.5), loc="center left",
            )
            fig2.tight_layout(rect=[0, 0, 0.88, 0.96])
        else:
            fig2.tight_layout(rect=[0, 0, 1, 0.96])

        if plot_results:
            plt.show()
        pdf.savefig(fig2, bbox_inches="tight")
        plt.close(fig2)


def plot_feature_umap(
    df_umap,
    info_cols,
    sample_cols,
    outpath,
    rows_per_page=4,
    nr_cols=2,
    rows_first_img=2,
    figsize=(8.27, 11.69),
    plot_results=True,
):
    n_plots = len(info_cols)
    n_rows = (n_plots // nr_cols) + (1 if n_plots % nr_cols != 0 else 0) + rows_first_img
    nr_pages = (n_rows // rows_per_page) + (1 if n_rows % rows_per_page != 0 else 0)

    with PdfPages(outpath) as pdf:
        plot_idx = 0

        for page in range(nr_pages):
            fig = plt.figure(figsize=figsize)
            gs = GridSpec(rows_per_page, nr_cols, figure=fig, wspace=0.3)

            if page == 0:
                ax = fig.add_subplot(gs[:rows_first_img, :])
                sns.scatterplot(
                    data=df_umap, x="UMAP1", y="UMAP2",
                    hue="ClusterID", s=20, alpha=0.5, palette="Set1", ax=ax,
                )
                ax.legend(
                    loc="upper left", prop={"size": 8},
                    bbox_to_anchor=(1, 1), borderpad=0.3,
                    labelspacing=0.4, columnspacing=0.1, frameon=False,
                )
                ax.set_title("ClusterID", fontsize=10, loc="center")
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                ax.set_xlabel("")
                ax.set_ylabel("")
                if plot_results:
                    plt.show()

            remaining_axes = [
                fig.add_subplot(gs[i, j])
                for i in range(rows_first_img if page == 0 else 0, rows_per_page)
                for j in range(nr_cols)
            ]

            for ax in remaining_axes:
                if plot_idx >= len(info_cols):
                    ax.remove()
                    continue
                colorcol = info_cols[plot_idx]
                if colorcol in sample_cols or df_umap.dtypes[colorcol] == bool:
                    sns.scatterplot(
                        data=df_umap, x="UMAP1", y="UMAP2",
                        s=20, hue=colorcol, alpha=0.5, palette="Set2", ax=ax,
                    )
                else:
                    sns.scatterplot(
                        data=df_umap, x="UMAP1", y="UMAP2",
                        s=20, hue=colorcol, palette="viridis", alpha=0.5, ax=ax,
                    )

                ax.legend(
                    loc="upper left", prop={"size": 8},
                    bbox_to_anchor=(1, 1), borderpad=0.3,
                    labelspacing=0.4, columnspacing=0.1, frameon=False,
                )
                ax.set_title(colorcol, fontsize=10, loc="center")
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_xticklabels([])
                ax.set_yticklabels([])
                ax.set_xlabel("")
                ax.set_ylabel("")
                plot_idx += 1

            fig.subplots_adjust(left=0.05, right=0.85, top=0.95, bottom=0.05)
            pdf.savefig(fig, dpi=600)
            plt.close(fig)


def plot_clustering_feature_heatmap(
    df_umap,
    info_cols,
    sample_cols,
    outpath,
    rows_per_page=7,
    nr_cols=2,
    figsize=(8.27, 11.69),
    plot_results=True,
    show_points=False,
    point_alpha=0.5,
    point_size=8,
    mean_marker_size=60,
):
    """
    Produce a PDF with:
      - Page 1: full-page min-max scaled heatmap of cluster means.
      - Subsequent pages: per-feature violin plots tiled across pages.
    """
    info_cols = list(info_cols) if info_cols is not None else []
    sample_cols = list(sample_cols) if sample_cols is not None else []

    def _round_legend_ticks(max_val):
        try:
            return round_legend_ticks(max_val)
        except Exception:
            if not np.isfinite(max_val) or max_val <= 0:
                return 1.0
            magnitude = 10.0 ** np.floor(np.log10(max_val))
            return float(np.ceil(max_val / magnitude) * magnitude)

    _cols_for_means = list(dict.fromkeys(list(info_cols) + ["ClusterID"]))
    df_for_means = df_umap[_cols_for_means].drop(columns=sample_cols, errors="ignore")
    cluster_means = (
        df_for_means
        .groupby("ClusterID", observed=False)
        .mean(numeric_only=True)
        .reset_index()
    )

    cluster_means_scaled = cluster_means.copy()
    scale_columns = [c for c in cluster_means.columns if c != "ClusterID"]

    X = cluster_means_scaled[scale_columns].apply(pd.to_numeric, errors="coerce")
    X = X.replace([np.inf, -np.inf], np.nan)
    all_nan_cols = X.columns[X.isna().all()].tolist()
    if all_nan_cols:
        X = X.drop(columns=all_nan_cols)
        scale_columns = [c for c in scale_columns if c not in all_nan_cols]

    if len(scale_columns) > 0:
        X_filled = X.copy()
        med = X_filled.median(numeric_only=True)
        X_filled = X_filled.fillna(med)
        cluster_means_scaled[scale_columns] = MinMaxScaler().fit_transform(X_filled[scale_columns])
        df_heatmap_scaled = cluster_means_scaled.melt(id_vars="ClusterID", var_name="var", value_name="AU")
        overall_heatmap_data = df_heatmap_scaled.pivot(index="var", columns="ClusterID", values="AU")
    else:
        overall_heatmap_data = pd.DataFrame()

    value_cols = [c for c in info_cols if c not in set(sample_cols) and c != "ClusterID"]
    df_values = df_umap[["ClusterID"] + value_cols].copy()
    df_values[value_cols] = df_values[value_cols].apply(pd.to_numeric, errors="coerce")
    df_values.replace([np.inf, -np.inf], np.nan, inplace=True)
    df_long = df_values.melt(id_vars="ClusterID", var_name="var", value_name="value")

    cluster_order = sorted(df_values["ClusterID"].dropna().unique().tolist())
    feat_names = [c for c in value_cols if c in df_long["var"].unique()]
    n_plots = len(feat_names)
    rows_for_plots = (n_plots + nr_cols - 1) // nr_cols
    nr_pages = max(1, (rows_for_plots + rows_per_page - 1) // rows_per_page)

    with PdfPages(outpath) as pdf:
        fig, ax = plt.subplots(figsize=figsize)
        if not overall_heatmap_data.empty:
            try:
                heatmap = sns.heatmap(
                    overall_heatmap_data, ax=ax, cmap="viridis", cbar=True, yticklabels=True,
                )
                ax.set_title("Min–Max Scaled Cluster Means", fontsize=14, pad=14)
                ax.set_xlabel("ClusterID", fontsize=10)
                ax.set_ylabel("", fontsize=10)
                ax.tick_params(axis="y", labelsize=6)
                ax.tick_params(axis="x", labelsize=8, rotation=0)
                cbar = heatmap.collections[0].colorbar
                cbar.ax.tick_params(labelsize=8)
                fig.tight_layout(pad=2.0)
            except Exception:
                ax.text(0.5, 0.5, "Overview heatmap unavailable", ha="center", va="center")
                ax.axis("off")
        else:
            ax.text(0.5, 0.5, "No features available for overview scaling", ha="center", va="center")
            ax.axis("off")

        pdf.savefig(fig, dpi=600)
        plt.close(fig)

        plot_idx = 0
        for page in range(nr_pages):
            fig = plt.figure(figsize=figsize)
            gs = GridSpec(rows_per_page, nr_cols, figure=fig, hspace=1.5, wspace=0.3)
            remaining_axes = [
                fig.add_subplot(gs[r, c])
                for r in range(rows_per_page)
                for c in range(nr_cols)
            ]

            for ax in remaining_axes:
                if plot_idx >= n_plots:
                    ax.remove()
                    continue

                feat = feat_names[plot_idx]
                sub = df_long.loc[
                    df_long["var"] == feat, ["ClusterID", "value"]
                ].dropna(subset=["ClusterID", "value"])
                if sub.empty:
                    ax.text(0.5, 0.5, f"{feat}\n(no finite data)", ha="center", va="center")
                    ax.axis("off")
                    plot_idx += 1
                    continue

                try:
                    sns.violinplot(
                        data=sub, x="ClusterID", y="value",
                        order=cluster_order, inner=None, ax=ax, cut=0,
                    )
                except Exception:
                    ax.text(0.5, 0.5, f"{feat}\n(plot unavailable)", ha="center", va="center")
                    ax.axis("off")
                    plot_idx += 1
                    continue

                if show_points:
                    sns.stripplot(
                        data=sub, x="ClusterID", y="value",
                        order=cluster_order, ax=ax,
                        dodge=False, jitter=0.2, alpha=point_alpha, size=point_size,
                    )

                means = sub.groupby("ClusterID", observed=False)["value"].mean().reindex(cluster_order)
                ax.scatter(
                    np.arange(len(cluster_order)), means.values,
                    s=mean_marker_size, edgecolor="black", linewidths=0.8, zorder=3,
                )

                ax.set_title(feat, fontsize=9)
                ax.set_xlabel("ClusterID", fontsize=8)
                ax.set_ylabel("Value", fontsize=8)
                ax.tick_params(axis="x", rotation=0, labelsize=7)
                ax.tick_params(axis="y", labelsize=7)
                plot_idx += 1

            fig.subplots_adjust(left=0.20, right=0.98, top=0.95, bottom=0.08)
            pdf.savefig(fig, dpi=600)
            plt.close(fig)

    if plot_results:
        print(f"Saved PDF to: {outpath}")
