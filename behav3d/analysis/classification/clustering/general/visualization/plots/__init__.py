import math
from pathlib import Path  # only needed if you pass Path objects as outpath

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec

import seaborn as sns
import scanpy as sc

from pandas.api.types import is_numeric_dtype
from sklearn.preprocessing import MinMaxScaler


def plot_per_cluster_proportions(
    df,
    groupby = ["ClusterID", "sample_name"],
    show=True
    ):
    prop_df = (
        df
        .groupby(groupby)
        .size()
        .groupby(level=0)
        .apply(lambda x: x / x.sum())  # normalize within each cluster
        .unstack(fill_value=0)         # make sample_name columns
    )
    # Plot stacked bar chart
    prop_df.index = prop_df.index.get_level_values(0)
    # Create stacked bar plot
    ax = prop_df.plot(kind='bar', stacked=True, figsize=(10, 6))

    plt.title("Proportion of sample_name per ClusterID")
    plt.xlabel("ClusterID")
    plt.ylabel("Proportion")
    plt.legend(title="Sample Name", bbox_to_anchor=(1.05, 1), loc='upper left')

    # Ensure upright tick labels
    plt.xticks(ticks=range(len(prop_df.index)), labels=prop_df.index.astype(str), rotation=0, ha='center')

    plt.tight_layout()

    if show:
        plt.show()
    else:
        return ax
 
def plot_number_per_clusters(
    df,
    cluster_col="ClusterID",
    show=True
    ):
    plt.figure(figsize=(6, 3))
    h_counts = df[cluster_col].value_counts().sort_index()
    ax = sns.barplot(x=h_counts.index.astype(str), y=h_counts.values, color="tab:green")
    ax.bar_label(ax.containers[0], padding=3)
    plt.title("Cluster sizes")
    plt.xlabel("Cluster")
    plt.ylabel("Count")
    ax.margins(y=0.15)
    plt.tight_layout()
    
    if show:
        plt.show()
    else:
        return ax
     
def plot_top_ranking_features(
    adata,
    groupby="ClusterID",
    n_features=15,
    method="wilcoxon",
):
    # Run differential expression per group
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
    )

    key = "rank_genes_groups"
    # Get group names from the rank_genes_groups result
    groups = adata.uns[key]["names"].dtype.names
    n_groups = len(groups)

    # Define a roughly square grid layout
    n_cols = math.ceil(np.sqrt(n_groups))
    n_rows = math.ceil(n_groups / n_cols)

    # Create figure and axes grid
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(4.5 * n_cols, 5 * n_rows),
        squeeze=False,
    )

    # Plot top-ranked genes for each group
    for ax, group in zip(axes.flat, groups):
        df = sc.get.rank_genes_groups_df(adata, group=group, key=key).copy()

        # Rank genes by absolute score
        df["abs_score"] = df["scores"].abs()
        top = df.sort_values("abs_score", ascending=False).head(n_features).iloc[::-1]

        # Horizontal bar plot of scores
        ax.barh(top["names"], top["scores"])
        ax.axvline(0, linewidth=0.8)

        ax.set_title(f"Cluster {group}")
        ax.set_xlabel("score")

        # Color gene labels by score direction
        for label, score in zip(ax.get_yticklabels(), top["scores"]):
            label.set_color("blue" if score > 0 else "red")

    # Disable unused axes
    for ax in axes.flat[len(groups):]:
        ax.axis("off")

    # Add shared legend for score direction
    legend_elements = [
        Line2D([0], [0], color="blue", lw=2, label="Present / high"),
        Line2D([0], [0], color="red", lw=2, label="Absent / low"),
    ]

    fig.legend(
        handles=legend_elements,
        loc="lower center",
        ncol=2,
        frameon=False,
    )

    # Adjust layout to make room for legend
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.show()

def plot_umap_feature_grid(
    df: pd.DataFrame,
    feature_cols: list[str],
    x_col: str = "UMAP1",
    y_col: str = "UMAP2",
    ncols: int = 4,
    max_plots: int | None = None,
    point_size: int = 5,
    alpha: float = 0.5,
    numeric_cmap: str = "viridis",
    categorical_palette: str = "tab20",
    add_colorbar: bool = True,
    page: int = 0,   # for pagination: 0-based page index
    ):
    """
    Creates a multi-row, multi-column grid of UMAP scatterplots colored by each feature in feature_cols.
    Filters out non-scalar or missing features automatically. Supports pagination via `page`.
    """
    # Filter valid features (exist, scalar, not all NaN)
    valid = []
    for c in feature_cols:
        if c in df.columns and df[c].notna().any():
            valid.append(c)
    if max_plots is not None:
        valid = valid[:max_plots]

    if len(valid) == 0:
        raise ValueError("No valid features to plot.")

    n = len(valid)
    nrows = math.ceil(n / ncols)

    # Pagination support: choose a slice of features per page
    per_page = nrows * ncols
    start = page * per_page
    end = min(start + per_page, len(valid))
    feats = valid[start:end]
    if len(feats) == 0:
        raise ValueError(f"No features to plot on page {page} (only {math.ceil(len(valid)/per_page)} page(s) available).")

    # Axes limits shared across panels
    x_min, x_max = df[x_col].min(), df[x_col].max()
    y_min, y_max = df[y_col].min(), df[y_col].max()

    # Build grid
    fig, axes = plt.subplots(
        nrows=math.ceil(len(feats)/ncols),
        ncols=ncols,
        figsize=(4*ncols, 3.5*math.ceil(len(feats)/ncols)),
        squeeze=False,
        constrained_layout=True
    )

    for i, feat in enumerate(feats):
        r, c = divmod(i, ncols)
        ax = axes[r, c]

        s = df[feat]
        # Numeric vs categorical handling
        if is_numeric_dtype(s):
            # Numeric: use matplotlib scatter for easy colorbar handling
            scanpy = ax.scatter(
                df[x_col], df[y_col],
                s=point_size, alpha=alpha,
                c=s, cmap=numeric_cmap, edgecolors="none"
            )
            if add_colorbar:
                cb = plt.colorbar(scanpy, ax=ax, fraction=0.046, pad=0.04)
                cb.ax.tick_params(labelsize=8)
        else:
            # Categorical: enforce category dtype and use seaborn palette
            s_cat = s.astype("category")
            tmp = df.copy()
            tmp[feat] = s_cat
            sns.scatterplot(
                data=tmp, x=x_col, y=y_col, hue=feat,
                palette=categorical_palette, s=point_size, alpha=alpha,
                legend=False, ax=ax
            )

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_title(feat, fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("")

    # Hide any unused axes (if grid not full)
    total_cells = axes.size
    for j in range(len(feats), total_cells):
        r, c = divmod(j, ncols)
        axes[r, c].axis("off")

    # Add a common title
    fig.suptitle("UMAP colored by features", fontsize=14)
    plt.show()     
    
def plot_clustering_feature_heatmap(
    df_umap,
    info_cols,
    sample_cols,
    outpath,
    rows_per_page = 7,
    nr_cols = 2,
    figsize = (8.27, 11.69),
    plot_results=True,
    show_points=False,
    point_alpha=0.5,
    point_size=8,
    mean_marker_size=60,
):
    """
    Produce a PDF with:
      • Page 1: full-page min–max scaled heatmap of cluster means.
      • Subsequent pages: per-feature violin plots tiled across pages.
    """

    info_cols   = list(info_cols) if info_cols is not None else []
    sample_cols = list(sample_cols) if sample_cols is not None else []

    # ---- Helper ----
    def _round_legend_ticks(max_val):
        try:
            return round_legend_ticks(max_val)
        except Exception:
            if not np.isfinite(max_val) or max_val <= 0:
                return 1.0
            magnitude = 10.0 ** np.floor(np.log10(max_val))
            return float(np.ceil(max_val / magnitude) * magnitude)

    # ---- Cluster means ----
    df_for_means = (
        df_umap[list(info_cols) + ["ClusterID"]]
        .drop(columns=sample_cols, errors="ignore")
    )
    cluster_means = (
        df_for_means
        .groupby("ClusterID", observed=False)
        .mean(numeric_only=True)
        .reset_index()
    )

    # ---- Min-max scaling ----
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

    # ---- Prepare violin plot data ----
    value_cols = [c for c in info_cols if c not in set(sample_cols)]
    df_values = df_umap[["ClusterID"] + value_cols].copy()
    for c in value_cols:
        df_values[c] = pd.to_numeric(df_values[c], errors="coerce")
    df_values.replace([np.inf, -np.inf], np.nan, inplace=True)
    df_long = df_values.melt(id_vars="ClusterID", var_name="var", value_name="value")

    cluster_order = sorted(df_values["ClusterID"].dropna().unique().tolist())
    feat_names = [c for c in value_cols if c in df_long["var"].unique()]
    n_plots = len(feat_names)
    rows_for_plots = (n_plots + nr_cols - 1) // nr_cols
    nr_pages = max(1, (rows_for_plots + rows_per_page - 1) // rows_per_page)

    with PdfPages(outpath) as pdf:
        # ---- Page 1: Full-page heatmap ----
        fig, ax = plt.subplots(figsize=figsize)
        if not overall_heatmap_data.empty:
            try:
                heatmap = sns.heatmap(
                    overall_heatmap_data,
                    ax=ax,
                    cmap="viridis",
                    cbar=True,
                    yticklabels=True
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

        # ---- Remaining pages: Violin plots ----
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
                sub = df_long.loc[df_long["var"] == feat, ["ClusterID", "value"]].dropna(subset=["ClusterID", "value"])
                if sub.empty:
                    ax.text(0.5, 0.5, f"{feat}\n(no finite data)", ha="center", va="center")
                    ax.axis("off")
                    plot_idx += 1
                    continue

                try:
                    sns.violinplot(
                        data=sub,
                        x="ClusterID",
                        y="value",
                        order=cluster_order,
                        inner=None,
                        ax=ax,
                        cut=0
                    )
                except Exception:
                    ax.text(0.5, 0.5, f"{feat}\n(plot unavailable)", ha="center", va="center")
                    ax.axis("off")
                    plot_idx += 1
                    continue

                if show_points:
                    sns.stripplot(
                        data=sub,
                        x="ClusterID",
                        y="value",
                        order=cluster_order,
                        ax=ax,
                        dodge=False,
                        jitter=0.2,
                        alpha=point_alpha,
                        size=point_size
                    )

                means = sub.groupby("ClusterID", observed=False)["value"].mean().reindex(cluster_order)
                ax.scatter(
                    np.arange(len(cluster_order)),
                    means.values,
                    s=mean_marker_size,
                    edgecolor="black",
                    linewidths=0.8,
                    zorder=3
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