import numpy as np
import pandas as pd

import scanpy  # used as scanpy.pp.*, scanpy.tl.*, scanpy.AnnData

from itertools import combinations
from sklearn.metrics import adjusted_rand_score

def run_pca(
    adata, 
    ncomps=50, 
    pca_var_selection=0.95, 
    svd_solver="full",
    zero_center=True,
    random_state=None,
    verbose=False,
    ):
    n_obs = int(adata.n_obs)
    n_vars = int(adata.n_vars)
    if n_obs <= 0 or n_vars <= 0:
        raise ValueError(
            f"PCA requires non-empty matrix; received n_obs={n_obs}, n_vars={n_vars}."
        )
    max_allowed = int(min(n_obs, n_vars))
    requested = max_allowed if ncomps is None else int(ncomps)
    safe_ncomps = max(1, min(requested, max_allowed))
    if bool(verbose) and requested != safe_ncomps:
        print(
            "[track-clustering] INFO "
            f"PCA ncomps clamped from {requested} to {safe_ncomps} "
            f"(n_obs={n_obs}, n_vars={n_vars})."
        )

    scanpy.pp.pca(
        adata,
        zero_center=zero_center,
        n_comps=safe_ncomps,
        svd_solver=svd_solver,
        random_state=random_state,
    )
    var_ratio = adata.uns["pca"]["variance_ratio"]
    n_pcs = int(np.searchsorted(np.cumsum(var_ratio), pca_var_selection) + 1)
    available_pcs = int(adata.obsm["X_pca"].shape[1])
    n_pcs = max(1, min(n_pcs, available_pcs))
    
    # truncate X_pca
    adata.obsm["X_pca"] = adata.obsm["X_pca"][:, :n_pcs]
    
    # truncate loadings so projection/ingest is consistent
    if "PCs" in adata.varm:
        adata.varm["PCs"] = adata.varm["PCs"][:, :n_pcs]

    # truncate PCA metadata arrays to match
    if "variance" in adata.uns["pca"]:
        adata.uns["pca"]["variance"] = adata.uns["pca"]["variance"][:n_pcs]
    adata.uns["pca"]["variance_ratio"] = adata.uns["pca"]["variance_ratio"][:n_pcs]

    # keep params consistent (optional but nice)
    if "params" in adata.uns["pca"]:
        adata.uns["pca"]["params"]["n_comps"] = n_pcs
        
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
    stability_resolutions=(0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5),
    n_stability_repeats=10,
    n_subsample_repeats=10,
    subsample_fraction=0.8,
    subsample_random_state=0,
    min_mean_ari_seed=0.7,
    min_mean_ari_subsample=0.7,
    max_std_ari_seed=0.1,
    max_std_ari_subsample=0.1,
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

    has_precomputed_neighbors = (
        ("neighbors" in adata.uns)
        and ("connectivities" in adata.obsp)
    )

    resolved_n_neighbors = n_neighbors

    if n_neighbors is None:
        if not has_precomputed_neighbors:
            raise ValueError(
                "n_neighbors=None but no precomputed neighbors graph found in AnnData. "
                "Precompute neighbors first (scanpy.pp.neighbors) or provide n_neighbors."
            )
        params = adata.uns.get("neighbors", {}).get("params", {})
        inferred = params.get("n_neighbors", None)
        if inferred is not None:
            resolved_n_neighbors = int(inferred)
        else:
            raise ValueError(
                "n_neighbors=None and precomputed neighbors graph lacks "
                "`adata.uns['neighbors']['params']['n_neighbors']`. "
                "Provide n_neighbors explicitly or recompute neighbors with Scanpy."
            )
        # print(f"Using precomputed neighbors graph (n_neighbors={resolved_n_neighbors}).")
    else:
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
        labels, best_res, summary = leiden_stability_search(
            adata,
            resolutions=stability_resolutions,
            n_repeats=n_stability_repeats,
            n_subsample_repeats=n_subsample_repeats,
            subsample_fraction=subsample_fraction,
            subsample_random_state=subsample_random_state,
            n_neighbors=resolved_n_neighbors,
            metric=metric,
            method=method,
            use_rep=use_rep,
            flavor=flavor,
            n_iterations=n_iterations,
            directed=directed,
            min_mean_ari_seed=min_mean_ari_seed,
            min_mean_ari_subsample=min_mean_ari_subsample,
            max_std_ari_seed=max_std_ari_seed,
            max_std_ari_subsample=max_std_ari_subsample,
        )
        adata.uns["leiden_stability_summary"] = summary
        adata.uns["leiden_stability_best_res"] = best_res
        adata.uns["leiden_stability_labels_by_resolution"] = summary.attrs.get("labels_by_resolution", {})
        adata.uns["leiden_stability_plot_results"] = summary.attrs.get("plot_results", False)
        best_row = summary.loc[summary["resolution"] == best_res].iloc[0]
        sub_txt = ""
        if pd.notna(best_row.get("mean_ari_subsample", np.nan)):
            sub_txt = (
                f" | subsample_mean_ARI={best_row['mean_ari_subsample']:.4f}"
                f" | subsample_std_ARI={best_row['std_ari_subsample']:.4f}"
                f" | subsample_median_k={int(best_row['median_k_subsample'])}"
            )
        print(
            "[auto] Selected Leiden resolution = "
            f"{best_res:.3f} | seed_mean_ARI={best_row['mean_ari']:.4f} | seed_std_ARI={best_row['std_ari']:.4f} | seed_median_k={int(best_row['median_k'])}"
            f"{sub_txt}"
        )

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

def leiden_stability_search(
    adata,
    resolutions=(0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5),
    n_repeats=10,
    n_subsample_repeats=10,
    subsample_fraction=0.8,
    subsample_random_state=0,
    n_neighbors=30,
    metric="euclidean",
    method="umap",
    use_rep="X_pca",
    flavor="igraph",
    n_iterations=-1,
    directed=False,
    random_states=None,
    restrict_k_range=(1, 10),
    min_mean_ari_seed=0.7,
    min_mean_ari_subsample=0.7,
    max_std_ari_seed=0.1,
    max_std_ari_subsample=0.1,
    plot_results=False,
    prefer="highest_resolution", # options: "highest_resolution" or "highest_k"
    cleanup_tmp_keys=True,
):
    """
    Find a Leiden resolution that is stable across repeated runs.
    Stability is evaluated in two ways:
      1) Seed stability: mean pairwise ARI across different random seeds on one graph.
      2) Subsampling stability: repeatedly subsample cells, recompute neighbors,
         rerun Leiden, and compute ARI on overlapping cells.

    Selection strategy (recommended):
      - Filter by restrict_k_range (median_k within [lo, hi])
      - Keep resolutions with seed mean_ari >= min_mean_ari_seed
      - Keep resolutions with subsample mean_ari >= min_mean_ari_subsample
      - Keep resolutions with seed std_ari <= max_std_ari_seed
      - Keep resolutions with subsample std_ari <= max_std_ari_subsample
      - Choose the *highest* resolution among those (or highest median_k)
      - If none pass both thresholds, fall back to argmax(seed mean_ari)

    Returns: (best_labels, best_res, summary_df)
    """
    if "neighbors" not in adata.uns:
        raise ValueError("adata must have neighbors computed before stability search.")

    if random_states is None:
        random_states = list(range(1, n_repeats + 1))
    else:
        if len(random_states) != n_repeats:
            # not fatal, but it's clearer if these match
            n_repeats = len(random_states)

    results = []
    per_res_keys = {}     # store obs keys (so we can cleanup)
    per_res_runs = {}     # store label arrays per resolution
    per_res_repr_labels = {}  # representative labels per tested resolution
    obs_names = adata.obs_names.to_numpy()
    rng = np.random.default_rng(subsample_random_state)

    for res in resolutions:
        runs = []
        keys = []
        for rs in random_states:
            key = f"leiden_tmp_{res:.3f}_{rs}"
            scanpy.tl.leiden(
                adata,
                resolution=float(res),
                random_state=int(rs),
                key_added=key,
                flavor=flavor,
                n_iterations=n_iterations,
                directed=directed,
            )
            runs.append(adata.obs[key].to_numpy())
            keys.append(key)

        # Seed stability on a fixed graph
        aris = [
            adjusted_rand_score(runs[i], runs[j])
            for i, j in combinations(range(len(runs)), 2)
        ]
        k_counts = [len(np.unique(r)) for r in runs]

        mean_ari = float(np.mean(aris)) if len(aris) else 1.0
        std_ari = float(np.std(aris)) if len(aris) else 0.0
        median_k = float(np.median(k_counts))
        # Representative labeling for this resolution: run with max avg ARI to others
        if len(runs) == 1:
            rep_labels = runs[0]
        else:
            avg_aris_res = []
            for i_run, a_run in enumerate(runs):
                others = [runs[j_run] for j_run in range(len(runs)) if j_run != i_run]
                avg_aris_res.append(float(np.mean([adjusted_rand_score(a_run, b_run) for b_run in others])))
            rep_labels = runs[int(np.argmax(avg_aris_res))]
        per_res_repr_labels[float(res)] = np.asarray(rep_labels)

        # Biological/subsampling stability with neighbor recomputation
        aris_subsample = []
        k_counts_subsample = []
        if n_subsample_repeats is not None and n_subsample_repeats >= 2 and 0 < subsample_fraction < 1:
            n_total = adata.n_obs
            n_sub = max(2, int(np.floor(n_total * subsample_fraction)))
            subsample_runs = []

            for i in range(int(n_subsample_repeats)):
                idx = np.sort(rng.choice(n_total, size=n_sub, replace=False))
                sub_names = obs_names[idx]
                ad_sub = adata[idx].copy()

                scanpy.pp.neighbors(
                    ad_sub,
                    n_neighbors=min(n_neighbors, max(2, ad_sub.n_obs - 1)),
                    metric=metric,
                    method=method,
                    knn=True,
                    use_rep=use_rep,
                    random_state=subsample_random_state + i,
                )
                scanpy.tl.leiden(
                    ad_sub,
                    resolution=float(res),
                    random_state=subsample_random_state,
                    key_added="leiden_subsample_tmp",
                    flavor=flavor,
                    n_iterations=n_iterations,
                    directed=directed,
                )
                sub_labels = ad_sub.obs["leiden_subsample_tmp"].to_numpy()
                k_counts_subsample.append(len(np.unique(sub_labels)))
                subsample_runs.append({"names": sub_names, "labels": sub_labels})

            for i, j in combinations(range(len(subsample_runs)), 2):
                common, idx_i, idx_j = np.intersect1d(
                    subsample_runs[i]["names"],
                    subsample_runs[j]["names"],
                    return_indices=True,
                )
                if len(common) < 2:
                    continue
                aris_subsample.append(
                    adjusted_rand_score(
                        subsample_runs[i]["labels"][idx_i],
                        subsample_runs[j]["labels"][idx_j],
                    )
                )

        mean_ari_subsample = float(np.mean(aris_subsample)) if len(aris_subsample) else np.nan
        std_ari_subsample = float(np.std(aris_subsample)) if len(aris_subsample) else np.nan
        median_k_subsample = float(np.median(k_counts_subsample)) if len(k_counts_subsample) else np.nan

        print(
            f"Testing resolution={res}: "
            f"seed(mean_ARI={mean_ari:.4f}, std_ARI={std_ari:.4f}, median_k={int(median_k)}, run_k={k_counts}); "
            f"subsample(mean_ARI={mean_ari_subsample:.4f}, std_ARI={std_ari_subsample:.4f}, "
            f"median_k={int(median_k_subsample) if not np.isnan(median_k_subsample) else 'NA'}, run_k={k_counts_subsample})"
        )

        results.append(
            dict(
                resolution=float(res),
                # Backward-compatible: these are seed-stability metrics
                mean_ari=mean_ari,
                std_ari=std_ari,
                median_k=median_k,
                # Additional biological/subsampling metrics
                mean_ari_seed=mean_ari,
                std_ari_seed=std_ari,
                median_k_seed=median_k,
                mean_ari_subsample=mean_ari_subsample,
                std_ari_subsample=std_ari_subsample,
                median_k_subsample=median_k_subsample,
            )
        )
        per_res_runs[float(res)] = runs
        per_res_keys[float(res)] = keys

    summary_all = pd.DataFrame(results).sort_values("resolution").reset_index(drop=True)
    summary_for_selection = summary_all.copy()

    # Restrict by cluster count range (using median_k) for selection only
    if restrict_k_range is not None:
        lo, hi = restrict_k_range
        summary_for_selection = summary_for_selection[
            (summary_for_selection["median_k"] >= lo) & (summary_for_selection["median_k"] <= hi)
        ].copy()

    if summary_for_selection.empty:
        raise ValueError("No resolutions left after applying restrict_k_range filter.")

    # ---- Selection: highest stable by both seed + subsample ARI thresholds ----
    stable = summary_for_selection[
        (summary_for_selection["mean_ari"] >= float(min_mean_ari_seed))
        & (summary_for_selection["std_ari"] <= float(max_std_ari_seed))
    ].copy()
    if "mean_ari_subsample" in stable.columns and not stable["mean_ari_subsample"].isna().all():
        stable = stable[
            (stable["mean_ari_subsample"] >= float(min_mean_ari_subsample))
            & (stable["std_ari_subsample"] <= float(max_std_ari_subsample))
        ].copy()
    print(
        f"Selection thresholds: min_seed_ARI={float(min_mean_ari_seed):.3f}, "
        f"min_subsample_ARI={float(min_mean_ari_subsample):.3f}, "
        f"max_seed_std_ARI={float(max_std_ari_seed):.3f}, "
        f"max_subsample_std_ARI={float(max_std_ari_subsample):.3f} "
        f"-> passing_resolutions={stable['resolution'].tolist() if not stable.empty else []}"
    )

    if not stable.empty:
        if prefer == "highest_k":
            # choose the most detailed stable partition
            best_row = stable.sort_values(["median_k", "resolution"], ascending=[False, False]).iloc[0]
        else:
            # default: choose the highest resolution that is stable enough
            best_row = stable.sort_values("resolution", ascending=False).iloc[0]
    else:
        # fallback: if nothing meets threshold, pick the best ARI
        best_row = summary_for_selection.loc[summary_for_selection["mean_ari"].idxmax()]

    best_res = float(best_row["resolution"])
    best_runs = per_res_runs[best_res]

    # Representative labeling: pick run most similar to others (max avg ARI)
    avg_aris = []
    for i, a in enumerate(best_runs):
        others = [best_runs[j] for j in range(len(best_runs)) if j != i]
        avg_aris.append(float(np.mean([adjusted_rand_score(a, b) for b in others])))

    best_labels = best_runs[int(np.argmax(avg_aris))]

    # Cleanup temporary columns
    if cleanup_tmp_keys:
        for res, keys in per_res_keys.items():
            for k in keys:
                if k in adata.obs.columns:
                    del adata.obs[k]

    # Expose tested-resolution representative labels for downstream visualization
    # (kept in DataFrame attrs for backward-compatible return signature)
    summary_all.attrs["labels_by_resolution"] = {
        f"{float(res):.10g}": labels
        for res, labels in per_res_repr_labels.items()
    }
    summary_all.attrs["plot_results"] = bool(plot_results)

    return best_labels, best_res, summary_all


