def plot_cluster_feature_heatmap(
    df: pd.DataFrame,
    cluster_col: str = "ClusterID",
    feature_cols=None,
    agg: str = "mean",              # "mean", "median", or callable
    sort_clusters: str = "size",    # "size", "index", "hierarchical", or None
    linkage_method: str = "ward",   # used if sort_clusters="hierarchical"
    metric: str = "euclidean",      # used if sort_clusters="hierarchical"
    figsize_scale=(0.35, 0.4),      # (per_cluster_width, per_feature_height)
    title: str | None = None
):
    """
    Heatmap of aggregated feature values per cluster.

    Output axes:
      X = clusters
      Y = features

    Colors:
      blue  = values < 0
      white = 0
      red   = values > 0
    """
    if feature_cols is None:
        raise ValueError("Pass feature_cols (e.g., descriptive_feature_cols).")

    sub = df[[cluster_col] + feature_cols].dropna(subset=[cluster_col]).copy()
    grouped = sub.groupby(cluster_col)[feature_cols].agg(agg)

    # ---- cluster sorting options ----
    if sort_clusters == "size":
        order = sub[cluster_col].value_counts().index
        grouped = grouped.loc[order]

    elif sort_clusters == "index":
        grouped = grouped.sort_index()

    elif sort_clusters == "hierarchical":
        X = grouped.values
        if linkage_method.lower() == "ward":
            Z = linkage(X, method="ward")
        else:
            d = pdist(X, metric=metric)
            Z = linkage(d, method=linkage_method)
        grouped = grouped.iloc[leaves_list(Z)]

    # enforce feature order
    grouped = grouped[feature_cols]

    # TRANSPOSE so rows=features (Y), cols=clusters (X)
    plot_data = grouped.T  # shape: (features, clusters)

    per_clust_w, per_feat_h = figsize_scale

    # --- tweak row height here ---
    per_feat_h = min(per_feat_h, 0.14)   # cap default row height
    # you can try 0.15–0.22 depending on taste

    fig_w = per_clust_w * plot_data.shape[1] + 4
    fig_h = max(2.0, per_feat_h * plot_data.shape[0] + 2)

    # ----- diverging cmap centered at 0 -----
    # symmetric scaling around 0 (best for z-scored inputs)
    max_abs = np.nanmax(np.abs(plot_data.values))
    norm = TwoSlopeNorm(vmin=-max_abs, vcenter=0.0, vmax=max_abs)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(
        plot_data.values,
        aspect="auto",
        interpolation="nearest",
        cmap="RdBu_r",   # blue for <0, red for >0
        norm=norm
    )

    cbar = plt.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label(f"{agg} (input already z-normalized)")

    # X axis = clusters
    ax.set_xticks(np.arange(plot_data.shape[1]))
    ax.set_xticklabels(plot_data.columns, rotation=45, ha="right", fontsize=9)

    # Y axis = features
    ax.set_yticks(np.arange(plot_data.shape[0]))
    ax.set_yticklabels(plot_data.index, fontsize=8)

    ax.set_xlabel(cluster_col)
    ax.set_ylabel("Features")
    ax.set_title(title or f"{agg} feature values per {cluster_col} ({sort_clusters} sorted)")

    plt.tight_layout()
    plt.show()

    return grouped

def plot_feature_cluster_heatmap(
    df_analysis: pd.DataFrame,
    feature_cols,
    cluster_col: str = "ClusterID",
    drop_noise_label: int | None = -1,            # set None to keep all clusters
    feature_category_map: dict[str, str] | None = None,
    zscore_rows: bool = True,                     # z-score features across clusters
    row_distance: str = "abs_correlation",        # {"abs_correlation","correlation","euclidean"}
    col_distance: str = "correlation",            # {"correlation","euclidean"}
    row_linkage: str = "average",
    col_linkage: str = "average",
    cmap: str = "vlag",
    figsize=(12, 14),
    savepath: str | None = None,
):
    """
    RNA-style clustermap:
      - rows = features clustered by (absolute) correlation
      - cols = clusters ordered by similarity
      - cells = mean(feature) within each cluster
    """
    # --- safety checks
    assert cluster_col in df_analysis.columns, f"{cluster_col} not found in df_analysis"
    for c in feature_cols:
        if c not in df_analysis.columns:
            raise ValueError(f"Feature '{c}' not found in df_analysis")

    # --- optionally remove HDBSCAN noise
    df = df_analysis.copy()
    if drop_noise_label is not None:
        df = df[df[cluster_col] != drop_noise_label]

    # --- compute feature x cluster matrix of means
    cluster_means = df.groupby(cluster_col)[list(feature_cols)].mean().T
    cluster_means = cluster_means.replace([np.inf, -np.inf], np.nan).dropna(how="all")

    # remove zero-variance rows (cannot compute correlation)
    nonconst = cluster_means.loc[cluster_means.std(axis=1) > 0]
    if nonconst.empty:
        raise ValueError("All features have zero variance across clusters.")
    M = nonconst.copy()

    # optional row z-scoring (emphasize patterns)
    if zscore_rows:
        M = (M - M.mean(axis=1).to_numpy()[:, None]) / (M.std(axis=1, ddof=0).to_numpy()[:, None] + 1e-9)

    # --- row (feature) distances
    if row_distance == "euclidean":
        d_rows = pdist(M.values, metric="euclidean")
    elif row_distance == "correlation":
        d_rows = pdist(M.values, metric="correlation")  # 1 - corr
    elif row_distance == "abs_correlation":
        C = np.corrcoef(M.values)                       # (n_features x n_features)
        D = 1.0 - np.abs(C)
        np.fill_diagonal(D, 0.0)
        d_rows = squareform(D, checks=False)
    else:
        raise ValueError(f"Unsupported row_distance='{row_distance}'")
    row_link = linkage(d_rows, method=row_linkage)
    row_order = leaves_list(row_link)

    # --- column (cluster) distances
    if col_distance == "euclidean":
        d_cols = pdist(M.T.values, metric="euclidean")
    elif col_distance == "correlation":
        d_cols = pdist(M.T.values, metric="correlation")
    else:
        raise ValueError(f"Unsupported col_distance='{col_distance}'")
    col_link = linkage(d_cols, method=col_linkage)
    col_order = leaves_list(col_link)

    # --- annotation bars
    row_colors = None
    if feature_category_map is not None:
        cats = pd.Series([feature_category_map.get(f, "other") for f in M.index], index=M.index)
        palette = sns.color_palette("tab20", n_colors=cats.nunique())
        lut = dict(zip(cats.unique(), palette))
        row_colors = cats.map(lut)

    cluster_counts = df[cluster_col].value_counts().reindex(M.columns).fillna(0).astype(int)
    norm = Normalize(vmin=cluster_counts.min(), vmax=max(cluster_counts.max(), 1))
    col_colors = [plt.cm.Blues(norm(v)) for v in cluster_counts.values]

    n_rows = M.shape[0]
    row_height = 0.25  # inches per feature row (adjust to taste)
    fig_height = max(6, n_rows * row_height)
    fig_width = figsize[0]
    # --- plot
    g = sns.clustermap(
        M,
        row_linkage=row_link,
        col_linkage=col_link,
        row_colors=row_colors,
        col_colors=col_colors,
        cmap=cmap,
        figsize=(fig_width, fig_height),
        xticklabels=True,
        yticklabels=True,
        cbar_kws={"label": "mean (row z-score)" if zscore_rows else "mean"},
        dendrogram_ratio=(0.12, 0.12),
        colors_ratio=(0.02, 0.04),
    )
    g.ax_heatmap.set_xlabel("Cluster")
    g.ax_heatmap.set_ylabel("Feature")

    # row legend (feature categories), if provided
    if feature_category_map is not None:
        for cat, color in lut.items():
            g.ax_col_dendrogram.bar(0, 0, color=color, label=cat, linewidth=0)
        g.ax_col_dendrogram.legend(title="Feature category", ncols=min(3, len(lut)), loc="center",
                                   bbox_to_anchor=(0.5, 1.2), frameon=False)

    # column legend (cluster sizes)
    import matplotlib.patches as mpatches
    ticks = np.unique(np.linspace(cluster_counts.min(),
                                  max(cluster_counts.max(), 1), 3).astype(int))
    if len(ticks):
        handles = [mpatches.Patch(color=plt.cm.Blues(norm(v)), label=f"N={v}") for v in ticks]
        g.ax_row_dendrogram.legend(handles=handles, title="Cluster size",
                                   loc="center", bbox_to_anchor=(0.5, 1.15), frameon=False)

    plt.tight_layout()
    if savepath:
        g.savefig(savepath, dpi=300, bbox_inches="tight")

    return g, {
        "matrix": M,
        "row_order": row_order,
        "col_order": col_order,
        "row_linkage": row_link,
        "col_linkage": col_link,
        "cluster_counts": cluster_counts,
    }

def plot_hmm_transition_matrix(hmm_model):
    K = hmm_model.n_components
    state_names = [f"S{i}" for i in range(K)]
    T_df = pd.DataFrame(hmm_model.transmat_, index=state_names, columns=state_names)
    T_df

    plt.figure(figsize=(5,4))
    plt.imshow(hmm_model.transmat_, aspect="auto")
    plt.colorbar(label="P(next state)")
    plt.xticks(range(K), state_names)
    plt.yticks(range(K), state_names)
    plt.xlabel("to state")
    plt.ylabel("from state")
    plt.title("HMM transition matrix")
    plt.show()