import pandas as pd
import numpy as np

import scanpy
from sklearn.metrics import adjusted_rand_score
from itertools import combinations

from hmmlearn.hmm import GaussianHMM

import matplotlib.pyplot as plt
from behav3d.utils.rolling_classification import *

def run_pca(
    adata, 
    ncomps=50, 
    pca_var=0.95, 
    svd_solver="full",
    zero_center=True,
    random_state=None
    ):
    sc.pp.pca(
        adata,
        zero_center=zero_center,
        n_comps=ncomps,
        svd_solver=svd_solver,
        random_state=random_state,
    )
    var_ratio = adata.uns["pca"]["variance_ratio"]
    n_pcs = int(np.searchsorted(np.cumsum(var_ratio), pca_var) + 1)
    adata.obsm["X_pca"] = adata.obsm["X_pca"][:, :n_pcs]
    return adata

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
    
def run_leiden_clustering(
    X,
    n_neighbors=30,
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
    key="clusters_leiden",
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

def run_hmm_state_classification(
    df_features: pd.DataFrame,
    feature_cols: list[str],
    model=None,
    id_cols: list[str] = ["sample_name", "TrackID"],
    time_col: str = "position_t",

    # n_states can be int OR "auto"
    n_states: int | str = "auto",
    k_min: int = 1,
    k_max: int = 8,

    covariance_type: str = "full",
    n_iter: int = 200,
    tol: float = 1e-3,
    random_state: int = 0,
    verbose: bool = False,
    # If you *really* want to allow NaNs through, set this False
    error_on_nans: bool = True,
    out_col_name: str = "hmm_state",
):
    """
    Fit an HMM to track time series and assign a hidden state to each timepoint.
    If n_states="auto", choose number of states via BIC sweep from k_min..k_max.

    Returns
    -------
    (df_with_states, fitted_model, selection_df)
        df_with_states is df_features plus a column "hmm_state".
        fitted_model is the chosen global model.
        selection_df is a DataFrame with columns ["k","bic","logL"] if auto,
        otherwise None.
    """

    def _prep_X_and_lengths(df: pd.DataFrame):
        req = set(id_cols + [time_col] + feature_cols)
        missing = req - set(df.columns)
        if missing:
            raise ValueError(f"df_features missing columns: {missing}")

        df_sorted = df.sort_values(id_cols + [time_col], kind="mergesort").copy()

        # Build concatenated observation matrix + lengths for HMM
        sequences = []
        lengths = []

        for _, g in df_sorted.groupby(id_cols, sort=False):
            X = g[feature_cols].to_numpy(dtype=float)

            if error_on_nans and np.isnan(X).any():
                raise ValueError(
                    "NaNs found in feature matrix. "
                    "Please interpolate/impute before calling run_hmm_state_classification."
                )

            sequences.append(X)
            lengths.append(len(g))

        X_all = np.vstack(sequences)
        return X_all, lengths, df_sorted

    def _num_params(K: int, D: int, cov_type: str):
        if cov_type == "full":
            cov_params = K * (D * (D + 1) / 2)
        elif cov_type == "diag":
            cov_params = K * D
        elif cov_type == "spherical":
            cov_params = K
        elif cov_type == "tied":
            cov_params = (D * (D + 1) / 2)
        else:
            raise ValueError(f"Unknown covariance_type: {cov_type}")
        return (K - 1) + K * (K - 1) + K * D + cov_params

    def _bic(model: GaussianHMM, X: np.ndarray, lengths: list[int]):
        logL = model.score(X, lengths=lengths)
        N, D = X.shape
        K = model.n_components
        p = _num_params(K, D, model.covariance_type)
        return -2 * logL + p * np.log(N)

    X_all, lengths, df_sorted = _prep_X_and_lengths(df_features)

    selection_df = None

    if model is None:
        # ---------- choose K ----------
        if isinstance(n_states, str):
            if n_states.lower() != "auto":
                raise ValueError("n_states must be int or 'auto'")
            if k_max < k_min:
                raise ValueError("k_max must be >= k_min")

            rows = []
            best_k, best_bic, best_model = None, np.inf, None

            init_params = "stmc"  # exclude 't'
            params = "stmc"
        
            for k in range(k_min, k_max + 1):
                m = GaussianHMM(
                    n_components=k,
                    covariance_type=covariance_type,
                    n_iter=n_iter,
                    tol=tol,
                    random_state=random_state,
                    verbose=verbose,
                    params=params,
                    init_params=init_params,
                )

                m.fit(X_all, lengths=lengths)
                bic_k = _bic(m, X_all, lengths)
                logL_k = m.score(X_all, lengths=lengths)
                rows.append({"k": k, "bic": bic_k, "logL": logL_k})

                if bic_k < best_bic:
                    best_k, best_bic, best_model = k, bic_k, m

            model = best_model
            selection_df = pd.DataFrame(rows)

        else:
            model = GaussianHMM(
                n_components=int(n_states),
                covariance_type=covariance_type,
                n_iter=n_iter,
                tol=tol,
                random_state=random_state,
                verbose=verbose,
            )
            model.fit(X_all, lengths=lengths)
            
    # ---------- decode ----------
    states_all = model.predict(X_all, lengths=lengths)
    # Make HMM state start from 1, instead of 0
    states_all = states_all + 1
    out = df_sorted.copy()
    out[out_col_name] = states_all

    state_map = out[id_cols + [time_col, out_col_name]]

    df_features_clean = df_features.drop(columns=[out_col_name], errors="ignore")
    df_out = df_features_clean.merge(
        state_map,
        on=id_cols + [time_col],
        how="left",
        sort=False,
        validate="many_to_one",
    )
    return df_out, model, selection_df

def run_sticky_hmm_state_classification(
    df_features: pd.DataFrame,
    feature_cols: list[str],
    model=None,
    id_cols: list[str] = ["sample_name", "TrackID"],
    time_col: str = "position_t",

    # n_states can be int OR "auto"
    n_states: int | str = "auto",
    k_min: int = 1,
    k_max: int = 8,

    covariance_type: str = "full",
    n_iter: int = 200,
    tol: float = 1e-3,
    random_state: int = 0,
    verbose: bool = False,
    error_on_nans: bool = True,
    out_col_name: str = "hmm_state",

    # --- sticky HMM knobs ---
    stickiness_kappa: float = 8,   # extra prior mass on self transitions (diagonal)
    transmat_alpha: float = 1.0,       # base prior mass everywhere
):
    """
    Fit a *sticky* Gaussian HMM to track time series and assign a hidden state
    to each timepoint.

    Stickiness is implemented by a Dirichlet prior on each transition row:
      prior_ij = alpha  (i != j)
      prior_ii = alpha + kappa

    If n_states="auto", choose number of states via BIC sweep from k_min..k_max.

    Returns
    -------
    (df_with_states, fitted_model, selection_df)
    """

    def _prep_X_and_lengths(df: pd.DataFrame):
        req = set(id_cols + [time_col] + feature_cols)
        missing = req - set(df.columns)
        if missing:
            raise ValueError(f"df_features missing columns: {missing}")

        df_sorted = df.sort_values(id_cols + [time_col], kind="mergesort").copy()

        sequences = []
        lengths = []
        for _, g in df_sorted.groupby(id_cols, sort=False):
            X = g[feature_cols].to_numpy(dtype=float)

            if error_on_nans and np.isnan(X).any():
                raise ValueError(
                    "NaNs found in feature matrix. "
                    "Please interpolate/impute before calling run_sticky_hmm_state_classification."
                )

            sequences.append(X)
            lengths.append(len(g))

        X_all = np.vstack(sequences)
        return X_all, lengths, df_sorted

    def _num_params(K: int, D: int, cov_type: str):
        if cov_type == "full":
            cov_params = K * (D * (D + 1) / 2)
        elif cov_type == "diag":
            cov_params = K * D
        elif cov_type == "spherical":
            cov_params = K
        elif cov_type == "tied":
            cov_params = (D * (D + 1) / 2)
        else:
            raise ValueError(f"Unknown covariance_type: {cov_type}")
        return (K - 1) + K * (K - 1) + K * D + cov_params

    def _bic(model: GaussianHMM, X: np.ndarray, lengths: list[int]):
        logL = model.score(X, lengths=lengths)
        N, D = X.shape
        K = model.n_components
        p = _num_params(K, D, model.covariance_type)
        return -2 * logL + p * np.log(N)

    def _make_sticky_prior(K: int, alpha: float, kappa: float) -> np.ndarray:
        if alpha <= 0:
            raise ValueError("transmat_alpha must be > 0 for a valid Dirichlet prior.")
        if kappa < 0:
            raise ValueError("stickiness_kappa must be >= 0.")
        prior = np.full((K, K), float(alpha), dtype=float)
        np.fill_diagonal(prior, float(alpha) + float(kappa))
        return prior

    def _make_sticky_init_transmat(K: int, alpha: float, kappa: float) -> np.ndarray:
        prior = _make_sticky_prior(K, alpha, kappa)
        tm = prior / prior.sum(axis=1, keepdims=True)
        return tm

    def _build_model(K: int):
        """
        Build a GaussianHMM with sticky transition prior.
        Uses transmat_prior if available; otherwise, falls back to post-fit smoothing.
        """
        sticky_prior = _make_sticky_prior(K, transmat_alpha, stickiness_kappa)
        sticky_init = _make_sticky_init_transmat(K, transmat_alpha, stickiness_kappa)

        # In hmmlearn, init_params controls what gets (re)initialized before fit.
        # We set transmat_ ourselves, so exclude 't' from init_params.
        init_params = "smc"  # exclude 't'
        params = "stmc"

        m = GaussianHMM(
            n_components=K,
            covariance_type=covariance_type,
            n_iter=n_iter,
            tol=tol,
            random_state=random_state,
            verbose=verbose,
            params=params,
            init_params=init_params,
            transmat_prior=sticky_prior,
            min_covar=1e-3# matrix-valued sticky Dirichlet prior
        )
        # Provide an explicitly sticky starting point
        m.startprob_ = np.full(K, 1.0 / K)
        m.transmat_ = sticky_init
        return m, True  # True => has transmat_prior support

    def _apply_postfit_sticky_smoothing(model: GaussianHMM, K: int):
        """
        Fallback if transmat_prior isn't supported:
        After fitting, bias the learned transmat_ toward self transitions and renormalize.
        This is not as principled as a true sticky prior during EM, but is often useful.
        """
        tm = model.transmat_.copy()
        tm = np.maximum(tm, 1e-15)

        # Add pseudo-count style mass then renormalize row-wise
        add = np.full((K, K), float(transmat_alpha))
        np.fill_diagonal(add, float(transmat_alpha) + float(stickiness_kappa))
        tm = tm + add
        tm = tm / tm.sum(axis=1, keepdims=True)
        model.transmat_ = tm

    X_all, lengths, df_sorted = _prep_X_and_lengths(df_features)
    selection_df = None

    if model is None:
        if isinstance(n_states, str):
            if n_states.lower() != "auto":
                raise ValueError("n_states must be int or 'auto'")
            if k_max < k_min:
                raise ValueError("k_max must be >= k_min")

            rows = []
            best_k, best_bic, best_model = None, np.inf, None

            for k in range(k_min, k_max + 1):
                m, has_prior = _build_model(k)
                m.fit(X_all, lengths=lengths)

                if not has_prior:
                    _apply_postfit_sticky_smoothing(m, k)

                bic_k = _bic(m, X_all, lengths)
                logL_k = m.score(X_all, lengths=lengths)
                rows.append({"k": k, "bic": bic_k, "logL": logL_k})

                if bic_k < best_bic:
                    best_k, best_bic, best_model = k, bic_k, m

            model = best_model
            selection_df = pd.DataFrame(rows)

        else:
            K = int(n_states)
            model, has_prior = _build_model(K)
            model.fit(X_all, lengths=lengths)
            if not has_prior:
                _apply_postfit_sticky_smoothing(model, K)

    # ---------- decode ----------
    states_all = model.predict(X_all, lengths=lengths) + 1  # 1-based states

    out = df_sorted.copy()
    out[out_col_name] = states_all

    state_map = out[id_cols + [time_col, out_col_name]]

    df_features_clean = df_features.drop(columns=[out_col_name], errors="ignore")
    df_out = df_features_clean.merge(
        state_map,
        on=id_cols + [time_col],
        how="left",
        sort=False,
        validate="many_to_one",
    )
    return df_out, model, selection_df

def compute_cluster_transition_matrix(
    adata,
    cluster_key: str,
    id_cols=("sample_name", "TrackID"),
    time_key="position_t",
    normalize: bool = True,
    plot: bool = False,
    ax: plt.Axes = None,
    only_transitions: bool = False,
):
    """
    Compute a transition matrix between cluster states from tracked objects over time.

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object containing tracking and clustering info in .obs.
    cluster_key : str
        Column in adata.obs with cluster labels (e.g. 'ClusterID', 'leiden').
    id_cols : sequence of str, default ("sample_name", "TrackID")
        Columns in adata.obs that together identify each track/object.
    time_key : str, default "position_t"
        Column in adata.obs giving time or frame index (must be sortable).
    normalize : bool, default True
        If True, return row-normalized probabilities (HMM-style).
        If False, returns raw transition counts.
    plot : bool, default False
        If True, plot the transition matrix as a heatmap.
    ax : matplotlib.axes.Axes, optional
        Axis to plot on. If None and plot=True, a new figure/axis is created.
    only_transitions : bool, default False
        If True, remove self-transitions (diagonal) from the *returned* matrices
        by setting diagonal counts to 0 and re-normalizing across off-diagonal
        entries (so rows sum to 1 when there is at least one off-diagonal transition).
        Also makes the diagonal appear empty in the plot.

    Returns
    -------
    transition_counts : pandas.DataFrame
        Matrix of transition counts, shape (n_states, n_states).
        Rows = current state, columns = next state.
    transition_probs : pandas.DataFrame
        Row-normalized transition matrix.
        If only_transitions=True, this is P(next | current, next != current).
    """
    df = adata.obs[list(id_cols) + [cluster_key, time_key]].copy()
    df = df.sort_values(list(id_cols) + [time_key])

    # Next state within each track
    df["next_state"] = df.groupby(list(id_cols))[cluster_key].shift(-1)

    # Drop last timepoint of each track (no "next" state)
    transitions = df.dropna(subset=["next_state"]).copy()

    # Build transition counts
    transition_counts = pd.crosstab(
        transitions[cluster_key],
        transitions["next_state"]
    )

    # Ensure all states present on both axes
    states = sorted(df[cluster_key].unique())
    transition_counts = transition_counts.reindex(
        index=states, columns=states, fill_value=0
    )

    # Optionally remove self-transitions in the returned matrices
    if only_transitions:
        np.fill_diagonal(transition_counts.values, 0)

    # Row-normalize to probabilities
    row_sums = transition_counts.sum(axis=1)
    transition_probs = transition_counts.div(row_sums.replace(0, np.nan), axis=0)

    if plot:
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 5))

        data_to_plot = transition_probs if normalize else transition_counts
        plot_arr = data_to_plot.to_numpy(dtype=float, copy=True)

        # Make diagonal appear empty ONLY when only_transitions=True
        if only_transitions:
            diag_mask = np.eye(plot_arr.shape[0], dtype=bool)
            plot_arr = np.ma.array(plot_arr, mask=diag_mask)

        imshow_kwargs = {}
        if normalize:  # probability heatmap
            imshow_kwargs.update(dict(vmin=0.0, vmax=1.0))
            
        im = ax.imshow(plot_arr, aspect="auto", **imshow_kwargs)

        # Make masked values transparent (so diagonal looks empty)
        if only_transitions:
            cmap = im.get_cmap()
            cmap.set_bad(alpha=0)

        ax.set_xticks(np.arange(len(states)))
        ax.set_yticks(np.arange(len(states)))
        ax.set_xticklabels(states, rotation=90)
        ax.set_yticklabels(states)

        plt.colorbar(im, ax=ax, label="P(next | current)" if normalize else "Count")
        ax.set_xlabel("Next state")
        ax.set_ylabel("Current state")
        ax.set_title("Cluster transition matrix")
        plt.tight_layout()

    return transition_counts, transition_probs

def filter_short_state_runs(
    adata,
    cluster_key: str,
    id_cols=("sample_name", "TrackID"),
    time_key: str = "position_t",
    length_removed: int = 1,
    new_key: str | None = None,
    inplace: bool = False,
):
    """
    Progressive + ordered smoothing:
      - For each track, for k = 1..length_removed:
          scan left-to-right;
          whenever you see a run of length == k, replace it immediately,
          then step back and re-check (so merges change what comes next).
    """
    df = adata.obs[list(id_cols) + [cluster_key, time_key]].copy()
    df = df.sort_values(list(id_cols) + [time_key])

    out = df[cluster_key].copy()

    def _smooth_exact_k_in_order(states: np.ndarray, k: int) -> np.ndarray:
        """Scan left-to-right; replace runs of exact length k immediately; re-check after each change."""
        n = len(states)
        if n == 0:
            return states

        i = 0
        while i < n:
            # find run [i:j)
            j = i + 1
            while j < n and states[j] == states[i]:
                j += 1
            run_len = j - i

            if run_len == k:
                prev_state = states[i - 1] if i > 0 else None
                next_state = states[j] if j < n else None

                if prev_state is None and next_state is None:
                    replacement = states[i]          # only run
                elif prev_state is None:
                    replacement = next_state         # start
                elif next_state is None:
                    replacement = prev_state         # end
                else:
                    # your rule: if neighbors differ, take previous
                    replacement = prev_state if prev_state != next_state else prev_state

                if replacement != states[i]:
                    states[i:j] = replacement
                    # step back to re-check merges around the edit
                    i = max(0, i - 1)
                    continue  # re-evaluate from (possibly) merged region

            # move to next run
            i = j

        return states

    for _, sub in df.groupby(list(id_cols), sort=False):
        idx = sub.index.to_numpy()
        states = sub[cluster_key].to_numpy().copy()

        for k in range(1, length_removed + 1):
            states = _smooth_exact_k_in_order(states, k)

        out.loc[idx] = states
        
    if new_key is not None:
        adata.obs[new_key] = out
    else:
        adata.obs[cluster_key] = out
    return adata

def condense_state_runs(
    adata,
    cluster_key: str,
    id_cols=["sample_name", "TrackID"],
    time_key: str = "position_t",
):
    """
    Collapse consecutive identical states within each track into runs.

    Example:
        For a given track:
            t:   0  1  2  3  4  5
            s:   1  1  1  2  2  1
        becomes runs:
            run 1: state=1, start_t=0, end_t=2, length=3
            run 2: state=2, start_t=3, end_t=4, length=2
            run 3: state=1, start_t=5, end_t=5, length=1

    Parameters
    ----------
    adata : anndata.AnnData
        AnnData with tracking and clustering info in .obs.
    cluster_key : str
        Column in adata.obs with state/cluster labels.
    id_cols : list of str, optional
        Columns that together define a track (e.g. ['sample_name', 'TrackID']).
        If None, defaults to ['TrackID'].
    time_key : str, default 'position_t'
        Column giving time/frame index (must be sortable).

    Returns
    -------
    runs_df : pandas.DataFrame
        One row per run, with columns:
            - all id_cols
            - cluster_key (state of the run)
            - 'run_start_time', 'run_end_time', 'run_length'
            - 'run_index' (order of the run within each track)
    """

    df = adata.obs[id_cols + [cluster_key, time_key]].copy()
    df = df.sort_values(id_cols + [time_key])

    runs_list = []

    for keys, sub in df.groupby(id_cols, sort=False):
        sub = sub.sort_values(time_key)

        # Identify where state changes
        # run_id increases whenever cluster changes
        run_id = (sub[cluster_key] != sub[cluster_key].shift()).cumsum()
        sub = sub.assign(_run_id=run_id)

        # Aggregate per run
        agg = sub.groupby("_run_id").agg(
            {
                cluster_key: "first",
                time_key: ["min", "max", "count"],
            }
        )

        # Flatten columns
        agg.columns = [
            "state",
            "run_start_time",
            "run_end_time",
            "run_length",
        ]

        # Add id_cols back & run index
        for col, val in zip(id_cols, keys if isinstance(keys, tuple) else (keys,)):
            agg[col] = val

        agg["run_index"] = np.arange(len(agg))

        runs_list.append(agg.reset_index(drop=True))

    runs_df = pd.concat(runs_list, axis=0, ignore_index=True)

    # Rename 'state' back to cluster_key so we can reuse easily
    runs_df = runs_df.rename(columns={"state": cluster_key})

    return runs_df


def rle_encode(states):
    """Return list of (state, run_length)."""
    if len(states) == 0:
        return []
    runs = []
    cur = states[0]
    length = 1
    for s in states[1:]:
        if s == cur:
            length += 1
        else:
            runs.append((cur, length))
            cur = s
            length = 1
    runs.append((cur, length))
    return runs

def fractions_from_runs(runs, states):
    total = sum(l for _, l in runs) if runs else 0
    cols = [f"overall_fraction_{s}" for s in states]
    out = {c: 0.0 for c in cols}
    if total == 0:
        return out, cols
    for s, l in runs:
        out[f"overall_fraction_{s}"] += l / total
    return out, cols

def bout_stats_from_runs(runs, states):
    lengths = {s: [] for s in states}
    for s, l in runs:
        lengths[s].append(l)

    total_bouts = len(runs)
    track_length = sum(l for _, l in runs)

    cols = []
    for s in states:
        cols.extend([f"bouts_nr_{s}", f"bouts_mean_length_{s}", f"bouts_max_length_{s}"])

    if total_bouts == 0 or track_length == 0:
        return {c: 0.0 for c in cols}, cols

    out = {}
    for s in states:
        arr = np.asarray(lengths[s], dtype=float)
        out[f"bouts_nr_{s}"] = float(arr.size) / float(total_bouts)
        out[f"bouts_mean_length_{s}"] = float(arr.mean()) / float(track_length) if arr.size else 0.0
        out[f"bouts_max_length_{s}"] = float(arr.max()) / float(track_length) if arr.size else 0.0

    return out, cols

def transition_probs_from_runs(runs, state_to_idx, states):
    K = len(states)
    M = np.zeros((K, K), dtype=float)
    labels = [s for s, _ in runs]
    if len(labels) < 2:
        return M  # all zeros
    for a, b in zip(labels[:-1], labels[1:]):
        if a in state_to_idx and b in state_to_idx:
            M[state_to_idx[a], state_to_idx[b]] += 1.0
    row_sums = M.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    return M / row_sums

def ngram_counts_from_runs(runs, n=3, weight="count"):
    labels = [s for s, _ in runs]
    lens = [l for _, l in runs]
    out = {}
    if len(labels) < n:
        return out
    for i in range(len(labels) - n + 1):
        g = tuple(labels[i : i + n])
        w = float(lens[i]) if weight == "duration" else 1.0
        out[g] = out.get(g, 0.0) + w
    return out

def l2_normalize_block(adata, cols, layer=None):
    """Row-wise L2 normalize selected columns. Leaves all-zero rows unchanged."""
    adata_norm = adata.copy()
    if layer is not None:
        X = adata_norm[:, cols].layers[layer]
    else:
        X = adata_norm[:, cols].X
    norms = np.linalg.norm(X, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    adata_norm[:, cols] = X / norms
    return adata_norm

def l2_normalize_all(adata):
    """Row-wise L2 normalize all columns. Leaves all-zero rows unchanged."""
   
    X = adata.X.copy()
    norms = np.linalg.norm(X, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    adata.X = X / norms
    return(adata)

def l2_normalize_features_blocks(
    adata,
    blocks,
    layer=None
    ):
    for block in blocks:
        adata = l2_normalize_block(adata, block, layer=layer)
    return adata
 
def extract_descibing_track_state_features(
    adata,
    group_col=("sample_name", "TrackID"),
    time_col="position_t",
    state_col="ClusterID",
    use_fractions=True,
    use_bout_stats=True,
    use_transitions=True,
    use_bigrams=True,
    use_trigrams=True,
    ngram_weight="count",
):
    """
    Track-level feature extraction -> returns (track_adata, blocks)

    track_adata.X        = features (n_tracks x n_features)
    track_adata.obs      = per-track metadata (sample_name, TrackID, position_t_min/max, n_timepoints)
    track_adata.var_names= feature names
    """

    group_col = list(group_col) if isinstance(group_col, (list, tuple)) else [group_col]

    obs = adata.obs[group_col + [time_col, state_col]].copy()
    obs = obs.sort_values(group_col + [time_col])

    # Stable state universe
    states = pd.Index(obs[state_col].astype("category").cat.categories).tolist()
    state_to_idx = {s: i for i, s in enumerate(states)}

    # Use observed=True to avoid unobserved categorical groups exploding results
    gb = obs.groupby(group_col, sort=False, observed=True)

    # Collect runs + n-gram vocab
    runs_by_track = []
    bigram_set, trigram_set = set(), set()

    for tid, df_t in gb:
        seq = df_t[state_col].tolist()
        runs = rle_encode(seq)  # <-- your RLE function

        runs_by_track.append((tid, runs))

        if use_bigrams:
            bigram_set |= set(ngram_counts_from_runs(runs, n=2, weight=ngram_weight).keys())
        if use_trigrams:
            trigram_set |= set(ngram_counts_from_runs(runs, n=3, weight=ngram_weight).keys())

    bigram_list = sorted(bigram_set)
    trigram_list = sorted(trigram_set)

    def bg_col(g): return f"bigram_{g[0]}>{g[1]}"
    def tg_col(g): return f"trigram_{g[0]}>{g[1]}>{g[2]}"

    # Track block columns explicitly while building features
    fraction_cols = []
    bout_cols = []
    transition_cols = []
    bigram_cols = [bg_col(g) for g in bigram_list] if use_bigrams else []
    trigram_cols = [tg_col(g) for g in trigram_list] if use_trigrams else []

    rows = []
    id_rows = []

    for tid, runs in runs_by_track:
        # tid is scalar if len(group_col)==1, else tuple
        if len(group_col) == 1:
            id_rows.append([tid])
        else:
            id_rows.append(list(tid))

        feats = {}

        if use_fractions:
            f, fcols = fractions_from_runs(runs, states)
            feats.update(f)
            if not fraction_cols:
                fraction_cols = fcols

        if use_bout_stats:
            b, bcols = bout_stats_from_runs(runs, states)
            feats.update(b)
            if not bout_cols:
                bout_cols = bcols

        if use_transitions:
            T = transition_probs_from_runs(runs, state_to_idx, states)
            if not transition_cols:
                transition_cols = [f"transitions_{a}>{b}" for a in states for b in states]
            for a in states:
                ia = state_to_idx[a]
                for b in states:
                    ib = state_to_idx[b]
                    feats[f"transitions_{a}>{b}"] = float(T[ia, ib])

        if use_bigrams:
            bg = ngram_counts_from_runs(runs, n=2, weight=ngram_weight)
            total = sum(bg.values()) or 1.0
            for g in bigram_list:
                feats[bg_col(g)] = bg.get(g, 0.0) / total

        if use_trigrams:
            tg = ngram_counts_from_runs(runs, n=3, weight=ngram_weight)
            total = sum(tg.values()) or 1.0
            for g in trigram_list:
                feats[tg_col(g)] = tg.get(g, 0.0) / total

        rows.append(feats)

    # Features
    df_feat = pd.DataFrame(rows).fillna(0.0)

    # IDs in the exact same order as df_feat rows
    df_ids = pd.DataFrame(id_rows, columns=group_col)

    # Time summaries (also observed=True to avoid unobserved category combos)
    df_time = (
        obs.groupby(group_col, sort=False, observed=True)[time_col]
           .agg(position_t_min="min", position_t_max="max", n_timepoints="size")
           .reset_index()
    )

    # Build obs table in the same row order as features
    df_meta = df_ids.merge(df_time, on=group_col, how="left")

    # Combine for conversion using your helper
    df_out = pd.concat([df_meta.reset_index(drop=True), df_feat.reset_index(drop=True)], axis=1)

    blocks = [fraction_cols, bout_cols, transition_cols, bigram_cols, trigram_cols]

    feature_cols = df_feat.columns.tolist()
    obs_cols = df_meta.columns.tolist()

    track_adata = df_to_adata(df_out, feature_cols=feature_cols, obs_cols=obs_cols)

    # Optional bookkeeping
    track_adata.uns["feature_blocks"] = {
        name: cols for name, cols in zip(
            ["fractions", "bout_stats", "transitions", "bigrams", "trigrams"], blocks
        ) if cols
    }
    track_adata.uns["feature_params"] = {
        "group_col": group_col,
        "time_col": time_col,
        "state_col": state_col,
        "ngram_weight": ngram_weight,
        "use_fractions": use_fractions,
        "use_bout_stats": use_bout_stats,
        "use_transitions": use_transitions,
        "use_bigrams": use_bigrams,
        "use_trigrams": use_trigrams,
    }

    return track_adata, blocks

def subset_full_adata_from_summary(
    adata_full,
    adata_tracks,
    n_tracks: int | None = None,
    cluster_key: str = "ClusterID",   
    n_tracks_per_cluster = None,
    sample_key: str = "sample_name",
    track_key: str = "TrackID",
    time_key: str = "position_t",
    tmin_key: str = "position_t_min",
    tmax_key: str = "position_t_max",
    random_state: int = 0,
    query: str | None = None, 
):
    """
    Subset a per-timepoint AnnData (adata_full) using a per-track summary AnnData (adata_tracks).

    Keeps rows in adata_full.obs that:
      - match (sample_key, track_key) present in chosen summary tracks
      - and have time_key within [tmin_key, tmax_key] from adata_tracks.obs for that track.

    Track selection
    ---------------
    - If n_tracks_per_cluster is provided:
        * sample up to N tracks *within each cluster* from adata_tracks.obs[cluster_key]
        * if n_tracks_per_cluster is an int: use the same N for every cluster
        * if it's a dict: use per-cluster values (clusters not in dict -> 0)
      In this mode, `n_tracks` is ignored.
    - If n_tracks is provided: 
        * sample n_tracks total across all tracks.
    - Else: keep all tracks (after optional query).


    Returns
    -------
    adata_sub : AnnData
    chosen : pd.DataFrame
        The chosen track windows (+ cluster column if present/used).
    """

    # --- tracks table (optionally filtered)
    base_cols = [sample_key, track_key, tmin_key, tmax_key]
    need_cluster = n_tracks_per_cluster is not None
    if need_cluster and cluster_key not in adata_tracks.obs.columns:
        raise ValueError(f"cluster_key='{cluster_key}' not found in adata_tracks.obs")

    cols = base_cols + ([cluster_key] if need_cluster and cluster_key not in base_cols else [])
    tracks_df = adata_tracks.obs[cols].copy()

    if query is not None:
        tracks_df = adata_tracks.obs.query(query)[cols].copy()

    tracks_df = tracks_df.reset_index(drop=True)

    if len(tracks_df) == 0:
        raise ValueError("No tracks available after applying query/filtering.")

    rng = np.random.default_rng(random_state)

    # --- choose which tracks
    if n_tracks_per_cluster is not None:
        # normalize to dict: cluster -> n
        if isinstance(n_tracks_per_cluster, dict):
            n_map = dict(n_tracks_per_cluster)
            default_n = 0
        else:
            default_n = int(n_tracks_per_cluster)
            n_map = {}

        chosen_parts = []
        # observed=True avoids unobserved categorical groups exploding
        for cl, df_cl in tracks_df.groupby(cluster_key, sort=False, observed=True):
            n_cl = n_map.get(cl, default_n)
            if n_cl <= 0 or len(df_cl) == 0:
                continue
            k = min(n_cl, len(df_cl))
            idx = rng.choice(len(df_cl), size=k, replace=False)
            chosen_parts.append(df_cl.iloc[idx])

        chosen = pd.concat(chosen_parts, axis=0, ignore_index=True) if chosen_parts else tracks_df.iloc[0:0].copy()

        if len(chosen) == 0:
            raise ValueError("No tracks selected with n_tracks_per_cluster (check cluster values / mapping).")
    elif n_tracks is not None:
        # Sample n_tracks total
        n = min(int(n_tracks), len(tracks_df))
        chosen = tracks_df.iloc[rng.choice(len(tracks_df), size=n, replace=False)].copy()
    else:
        chosen = tracks_df
      
    chosen = chosen.reset_index(drop=True)

    # --- build mask by merging timepoints with chosen windows
    obs_full = adata_full.obs[[sample_key, track_key, time_key]].copy()

    tmp = obs_full.reset_index(drop=False).merge(
        chosen[base_cols],  # only need keys + windows for filtering
        on=[sample_key, track_key],
        how="inner",
    )

    in_win = (tmp[time_key] >= tmp[tmin_key]) & (tmp[time_key] <= tmp[tmax_key])

    # IMPORTANT: subset by integer row positions to be robust even if obs index has duplicates
    kept_rowpos = tmp.loc[in_win].index.to_numpy()
    adata_sub = adata_full[kept_rowpos, :].copy()

    return adata_sub, chosen