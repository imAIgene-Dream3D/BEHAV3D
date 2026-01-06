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
    pca_var_selection=0.95, 
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
    n_pcs = int(np.searchsorted(np.cumsum(var_ratio), pca_var_selection) + 1)
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

def relabel_cluster_ids(
    adata,
    mapping,
    cluster_key="ClusterID",
    original_key="ClusterID_original",
    new_key=None,
    keep_unmapped=True,
    unmapped_label="unlabeled",
    overwrite_original=False,
):
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"{cluster_key} not found in adata.obs")

    if (original_key in adata.obs.columns) and (not overwrite_original):
        raise ValueError(f"{original_key} already exists")

    adata.obs[original_key] = adata.obs[cluster_key].astype(str).copy()
    current = adata.obs[cluster_key].astype(str)

    if isinstance(mapping, dict):
        map_dict = {str(k): v for k, v in mapping.items()}
        mapped = current.map(map_dict)
    else:
        uniq = np.array(sorted(current.unique(), key=lambda x: (int(x) if x.isdigit() else x)))
        labels = list(mapping)
        if len(labels) < len(uniq):
            raise ValueError("Not enough labels for number of clusters")
        map_dict = {uniq[i]: labels[i] for i in range(len(uniq))}
        mapped = current.map(map_dict)

    if keep_unmapped:
        out = mapped.where(mapped.notna(), current)
    else:
        out = mapped.fillna(unmapped_label)

    out_col = new_key if new_key is not None else cluster_key
    adata.obs[out_col] = pd.Categorical(out)

    return adata