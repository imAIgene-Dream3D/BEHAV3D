def run_pca(
    adata, 
    ncomps=50, 
    pca_var_selection=0.95, 
    svd_solver="full",
    zero_center=True,
    random_state=None
    ):
    scanpy.pp.pca(
        adata,
        zero_center=zero_center,
        n_comps=ncomps,
        svd_solver=svd_solver,
        random_state=random_state,
    )
    var_ratio = adata.uns["pca"]["variance_ratio"]
    n_pcs = int(np.searchsorted(np.cumsum(var_ratio), pca_var_selection) + 1)
    adata.obsm["X_pca"] = adata.obsm["X_pca"][:, :n_pcs]
    return adata

def run_leiden_clustering(
    X,
    n_neighbors=30,
    flavor="igraph",
    n_iterations=-1,
    directed=False,
    metric="euclidean",
    method="umap",
    use_rep="X_pca",
    resolution=1.0,
    random_state=0,
    stability_resolutions=(0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0),
    n_stability_repeats=8,
    key_added="clusters_leiden",
):
    """
    If X is AnnData:
        - modifies X in-place
        - returns only the updated AnnData

    If X is pandas DataFrame or numpy/array-like:
        - creates temporary AnnData
        - returns labels only
    """

    inplace = isinstance(X, scanpy.AnnData)

    if inplace:
        adata = X  # edit in place
    else:
        if hasattr(X, "values"):
            adata = scanpy.AnnData(X.values)
            try:
                adata.obs_names = X.index.astype(str)
            except Exception:
                pass
        else:
            adata = scanpy.AnnData(X)

    scanpy.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        metric=metric,
        method=method,
        knn=True,
        use_rep=use_rep,
        random_state=random_state
    )

    if resolution in ("auto", None):
        labels, best_res, summary = _leiden_stability_search(
            adata,
            resolutions=stability_resolutions,
            n_repeats=n_stability_repeats,
        )
        adata.uns["leiden_stability_summary"] = summary
        adata.uns["leiden_stability_best_res"] = best_res
        print(f"[auto] Selected Leiden resolution = {best_res:.3f}")

        adata.obs[key_added] = pd.Categorical(labels)

    else:
        scanpy.tl.leiden(
            adata,
            flavor=flavor,
            n_iterations=n_iterations,
            directed=directed,
            resolution=float(resolution),
            random_state=random_state,
            key_added=key_added,
        )
        labels = adata.obs[key_added].astype("category").cat.codes.to_numpy()
        adata.uns["leiden_stability_best_res"] = float(resolution)

    # Return labels based from 1 isntead of 0
    adata.obs[key_added] = (
        adata.obs[key_added].astype(int) + 1
    ).astype(str).astype("category")
    
    if inplace:
        return adata
    else:
        return labels

def _leiden_stability_search(
    adata,
    resolutions=(0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0),
    n_repeats=8,
    random_states=None,
    restrict_k_range=(3, 80),
    
):
    """
    Internal helper: find the most stable Leiden resolution via mean pairwise ARI.
    Returns (best_labels, best_res, summary_df)
    """
    if "neighbors" not in adata.uns:
        raise ValueError("adata must have neighbors computed before stability search.")

    if random_states is None:
        random_states = list(range(1, n_repeats + 1))

    results, per_res_labels = [], {}

    for res in resolutions:
        print("Testing resolution:", res)
        runs = []
        for rs in random_states:
            key = f"leiden_tmp_{res:.3f}_{rs}"
            scanpy.tl.leiden(adata, resolution=res, random_state=rs, key_added=key)
            runs.append(adata.obs[key].to_numpy())

        # pairwise ARIs
        aris = [
            adjusted_rand_score(runs[i], runs[j])
            for i, j in combinations(range(len(runs)), 2)
        ]
        k_counts = [len(np.unique(r)) for r in runs]

        results.append(
            dict(
                resolution=res,
                mean_ari=np.mean(aris),
                std_ari=np.std(aris),
                median_k=np.median(k_counts),
            )
        )
        per_res_labels[res] = runs

    summary = pd.DataFrame(results)
    # restrict range of cluster counts
    if restrict_k_range is not None:
        lo, hi = restrict_k_range
        mask = (summary["median_k"] >= lo) & (summary["median_k"] <= hi)
        if mask.any():
            summary = summary.loc[mask]

    # pick best resolution
    best_res = summary.loc[summary["mean_ari"].idxmax(), "resolution"]

    # pick representative labeling at best_res
    best_runs = per_res_labels[best_res]
    avg_aris = [
        np.mean(
            [adjusted_rand_score(a, b) for j, b in enumerate(best_runs) if j != i]
        )
        for i, a in enumerate(best_runs)
    ]
    best_labels = best_runs[int(np.argmax(avg_aris))]

    return best_labels, float(best_res), summary

def merge_small_clusters(
    adata,
    key="ClusterID",
    min_size=200,
    use_rep="X_pca",  # or "X" or "X_umap" depending on what you trust
):
    labels = adata.obs[key].astype("category")
    counts = labels.value_counts()

    small_clusters = counts[counts < min_size].index.tolist()
    if len(small_clusters) == 0:
        print("No small clusters to merge.")
        return adata

    X = adata.obsm[use_rep] if use_rep in adata.obsm else adata.X

    # compute centroids for all clusters
    centroids = {}
    for cl in labels.cat.categories:
        idx = (labels == cl).to_numpy()
        centroids[cl] = np.asarray(X[idx].mean(axis=0)).ravel()

    centroids = pd.DataFrame(centroids).T  # (n_clusters, n_features)

    # for each small cluster, find the closest *larger* cluster
    new_labels = labels.copy()
    for cl in small_clusters:
        cl_size = counts[cl]
        cl_centroid = centroids.loc[cl].to_numpy()

        # candidate target clusters: not the same and >= min_size
        big_clusters = [
            c for c, s in counts.items()
            if (c != cl) and (s >= min_size)
        ]
        if len(big_clusters) == 0:
            print(f"Cluster {cl} is small but there is no larger cluster to merge into.")
            continue

        big_centroids = centroids.loc[big_clusters].to_numpy()
        dists = np.linalg.norm(big_centroids - cl_centroid, axis=1)
        target = big_clusters[int(np.argmin(dists))]

        print(f"Merging small cluster {cl} (n={cl_size}) into {target} (n={counts[target]}).")
        new_labels[new_labels == cl] = target

    # re-categorize to keep it tidy
    new_labels = new_labels.astype("category").cat.remove_unused_categories()
    adata.obs[key] = new_labels
    
    return adata