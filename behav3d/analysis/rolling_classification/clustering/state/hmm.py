import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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
    stickiness_kappa: float = 8.0,     # extra prior mass on self transitions (diagonal)
    transmat_alpha: float = 1.0,       # base prior mass everywhere
    min_covar: float = 1e-3,           # covariance floor for numerical stability
):
    """
    Fit a *sticky* Gaussian HMM to track time series and assign a hidden state
    to each timepoint.

    Stickiness is implemented via a Dirichlet prior on each transition row:
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
        # startprob_ is FIXED => do not count (K-1)
        # transmat_ is learned => count K*(K-1) free params (row-normalized)
        # emissions: means + covars
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
        trans_params = K * (K - 1)
        mean_params = K * D
        return trans_params + mean_params + cov_params

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
        sticky_prior = _make_sticky_prior(K, transmat_alpha, stickiness_kappa)
        sticky_init = _make_sticky_init_transmat(K, transmat_alpha, stickiness_kappa)

        # We set startprob_ and transmat_ ourselves -> exclude 's' and 't' from init_params.
        # We do NOT want to train startprob_ -> exclude 's' from params.
        init_params = "mc"
        params = "tmc"  # learn transmat + emissions, keep startprob_ fixed

        m = GaussianHMM(
            n_components=K,
            covariance_type=covariance_type,
            n_iter=n_iter,
            tol=tol,
            random_state=random_state,
            verbose=verbose,
            params=params,
            init_params=init_params,
            transmat_prior=sticky_prior,  # matrix-valued Dirichlet prior (sticky)
            min_covar=min_covar,
        )

        # Fixed uniform initial state distribution (not trained, not overwritten)
        m.startprob_ = np.full(K, 1.0 / K, dtype=float)

        # Sticky starting point for transitions (not overwritten)
        m.transmat_ = sticky_init
        return m

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
                m = _build_model(k)
                m.fit(X_all, lengths=lengths)

                bic_k = _bic(m, X_all, lengths)
                logL_k = m.score(X_all, lengths=lengths)
                rows.append({"k": k, "bic": bic_k, "logL": logL_k})

                if bic_k < best_bic:
                    best_k, best_bic, best_model = k, bic_k, m

            model = best_model
            selection_df = pd.DataFrame(rows)

        else:
            K = int(n_states)
            model = _build_model(K)
            model.fit(X_all, lengths=lengths)

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

def hmm_emission_distance_symmetric_kl(model, eps=1e-6):
    K, D = model.means_.shape
    mus = np.asarray(model.means_)

    if model.covariance_type == "full":
        covs = [np.asarray(model.covars_[k]) for k in range(K)]
    elif model.covariance_type == "diag":
        covs = [np.diag(np.asarray(model.covars_[k])) for k in range(K)]
    else:
        raise ValueError("Supported: covariance_type 'full' or 'diag'")

    covs = [C + np.eye(D) * eps for C in covs]

    def kl_gauss(mu0, S0, mu1, S1):
        invS1 = np.linalg.inv(S1)
        diff = (mu1 - mu0).reshape(-1, 1)
        term_trace = np.trace(invS1 @ S0)
        term_quad = float(diff.T @ invS1 @ diff)
        sign0, logdet0 = np.linalg.slogdet(S0)
        sign1, logdet1 = np.linalg.slogdet(S1)
        if sign0 <= 0 or sign1 <= 0:
            raise ValueError("Non-positive definite covariance encountered. Increase min_covar/eps or use diag.")
        term_logdet = (logdet1 - logdet0)
        return 0.5 * (term_trace + term_quad - D + term_logdet)

    dist = np.zeros((K, K), float)
    for i in range(K):
        for j in range(i + 1, K):
            kl_ij = kl_gauss(mus[i], covs[i], mus[j], covs[j])
            kl_ji = kl_gauss(mus[j], covs[j], mus[i], covs[i])
            d = 0.5 * (kl_ij + kl_ji)
            dist[i, j] = d
            dist[j, i] = d

    return pd.DataFrame(dist, index=np.arange(1, K + 1), columns=np.arange(1, K + 1))

def plot_state_distance_heatmap(dist_df, title="State distance (symmetric KL)", cmap="viridis"):
    M = dist_df.to_numpy(float)
    fig, ax = plt.subplots(figsize=(0.65 * len(dist_df) + 3, 0.6 * len(dist_df) + 2.5))
    im = ax.imshow(M, aspect="auto", interpolation="nearest", cmap=cmap)

    ax.set_title(title)
    ax.set_xlabel("State")
    ax.set_ylabel("State")

    ax.set_xticks(np.arange(len(dist_df)))
    ax.set_yticks(np.arange(len(dist_df)))
    ax.set_xticklabels(dist_df.columns.tolist())
    ax.set_yticklabels(dist_df.index.tolist())

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("distance")

    plt.tight_layout()
    plt.show()

def merge_close_states_by_distance(
    df_out,
    dist_df,
    threshold,
    state_col="hmm_state",
    new_col="hmm_state_merged",
    relabel_compact=True,
):
    states = dist_df.index.to_list()
    K = len(states)
    idx = {s: i for i, s in enumerate(states)}

    parent = list(range(K))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    M = dist_df.to_numpy(float)
    for i in range(K):
        for j in range(i + 1, K):
            if M[i, j] <= threshold:
                union(i, j)

    groups = {}
    for s in states:
        r = find(idx[s])
        groups.setdefault(r, []).append(s)

    groups_list = list(groups.values())
    groups_list = [sorted(g) for g in groups_list]
    groups_list.sort(key=lambda g: (len(g), g[0]), reverse=True)

    if relabel_compact:
        new_labels = {tuple(g): i + 1 for i, g in enumerate(groups_list)}
        mapping = {}
        for g in groups_list:
            lab = new_labels[tuple(g)]
            for s in g:
                mapping[s] = lab
    else:
        mapping = {}
        for g in groups_list:
            lab = min(g)
            for s in g:
                mapping[s] = lab

    out = df_out.copy()
    out[new_col] = out[state_col].map(mapping)

    merge_df = pd.DataFrame(
        [{"merged_state": mapping[g[0]], "members": g} for g in groups_list]
    ).sort_values("merged_state")

    return out, mapping, merge_df