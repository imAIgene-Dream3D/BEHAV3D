import pandas as pd
import numpy as np

import scanpy
from sklearn.metrics import adjusted_rand_score
from itertools import combinations

from hmmlearn.hmm import GaussianHMM

import matplotlib.pyplot as plt

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
        method="umap",
        knn=True,
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
            transmat_prior=sticky_prior,  # matrix-valued sticky Dirichlet prior
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