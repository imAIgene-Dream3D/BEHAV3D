import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from hmmlearn.hmm import GaussianHMM


def run_hmm_state_classification(
    df_features,
    feature_cols,
    model=None,
    id_cols=("sample_name", "TrackID"),
    time_col="position_t",
    n_states="auto",
    k_min=1,
    k_max=8,
    covariance_type="full",
    n_iter=200,
    tol=1e-3,
    random_state=0,
    verbose=False,
    error_on_nans=True,
    out_col_name="hmm_state",
    # --- sticky knobs (optional) ---
    sticky=False,
    stickiness_kappa=8.0,
    transmat_alpha=1.0,
    min_covar=1e-3,
):
    """
    Fit a (Gaussian) HMM to track time series and assign a hidden state to each timepoint.

    If n_states="auto", choose number of states via BIC sweep from k_min..k_max.

    If sticky=True, apply a "sticky" Dirichlet prior to the transition matrix:
      prior_ij = alpha         (i != j)
      prior_ii = alpha + kappa (diagonal)

    Returns
    -------
    (df_with_states, fitted_model, selection_df)
        df_with_states is df_features plus a column out_col_name.
        fitted_model is the chosen global model.
        selection_df is a DataFrame with columns ["k","bic","logL"] if auto, otherwise None.
    """

    def _prep_X_and_lengths(df):
        req = set(list(id_cols) + [time_col] + list(feature_cols))
        missing = req - set(df.columns)
        if missing:
            raise ValueError(f"df_features missing columns: {missing}")

        df_sorted = df.sort_values(list(id_cols) + [time_col], kind="mergesort").copy()

        sequences = []
        lengths = []
        for _, g in df_sorted.groupby(list(id_cols), sort=False):
            X = g[list(feature_cols)].to_numpy(dtype=float)

            if error_on_nans and np.isnan(X).any():
                raise ValueError(
                    "NaNs found in feature matrix. "
                    "Please interpolate/impute before calling run_hmm_state_classification."
                )

            sequences.append(X)
            lengths.append(len(g))

        X_all = np.vstack(sequences) if sequences else np.empty((0, len(feature_cols)), dtype=float)
        return X_all, lengths, df_sorted

    def _num_params(K, D, cov_type, count_startprob=True):
        # startprob_ params: (K-1) if learned, else 0
        # transmat_ free params: K*(K-1) (row-normalized)
        # means: K*D
        # covars: depends on type
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

        start_params = (K - 1) if count_startprob else 0
        trans_params = K * (K - 1)
        mean_params = K * D
        return start_params + trans_params + mean_params + cov_params

    def _bic(m, X, lengths_, count_startprob=True):
        logL = m.score(X, lengths=lengths_)
        N, D = X.shape
        K = m.n_components
        p = _num_params(K, D, m.covariance_type, count_startprob=count_startprob)
        return -2 * logL + p * np.log(N)

    def _make_sticky_prior(K, alpha, kappa):
        if alpha <= 0:
            raise ValueError("transmat_alpha must be > 0 for a valid Dirichlet prior.")
        if kappa < 0:
            raise ValueError("stickiness_kappa must be >= 0.")
        prior = np.full((K, K), float(alpha), dtype=float)
        np.fill_diagonal(prior, float(alpha) + float(kappa))
        return prior

    def _make_sticky_init_transmat(K, alpha, kappa):
        prior = _make_sticky_prior(K, alpha, kappa)
        return prior / prior.sum(axis=1, keepdims=True)

    def _build_model(K):
        if sticky:
            sticky_prior = _make_sticky_prior(K, transmat_alpha, stickiness_kappa)
            sticky_init = _make_sticky_init_transmat(K, transmat_alpha, stickiness_kappa)

            # We set startprob_ and transmat_ ourselves -> exclude 's' and 't' from init_params.
            # We do NOT want to train startprob_ -> exclude 's' from params.
            init_params = "mc"
            params = "tmc"

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
                min_covar=min_covar,
            )

            # Fixed uniform initial state distribution
            m.startprob_ = np.full(K, 1.0 / K, dtype=float)

            # Sticky starting point for transitions
            m.transmat_ = sticky_init
            return m

        # non-sticky
        return GaussianHMM(
            n_components=K,
            covariance_type=covariance_type,
            n_iter=n_iter,
            tol=tol,
            random_state=random_state,
            verbose=verbose,
            min_covar=min_covar,
        )

    X_all, lengths, df_sorted = _prep_X_and_lengths(df_features)

    if X_all.shape[0] == 0:
        out = df_sorted.copy()
        out[out_col_name] = np.nan
        selection_df = None
        df_features_clean = df_features.drop(columns=[out_col_name], errors="ignore")
        df_out = df_features_clean.merge(
            out[list(id_cols) + [time_col, out_col_name]],
            on=list(id_cols) + [time_col],
            how="left",
            sort=False,
            validate="many_to_one",
        )
        return df_out, model, selection_df

    selection_df = None

    if model is None:
        if isinstance(n_states, str):
            if n_states.lower() != "auto":
                raise ValueError("n_states must be int or 'auto'")
            if k_max < k_min:
                raise ValueError("k_max must be >= k_min")

            rows = []
            best_bic = np.inf
            best_model = None

            for k in range(int(k_min), int(k_max) + 1):
                m = _build_model(k)
                m.fit(X_all, lengths=lengths)
                if hasattr(m, "monitor_") and not m.monitor_.converged:
                    warnings.warn(
                        f"HMM (n_states={k}) stopped after reaching n_iter={n_iter} without converging. "
                        "The result may be suboptimal — consider increasing n_iter or relaxing tol.",
                        UserWarning,
                        stacklevel=2,
                    )

                # if sticky=True and startprob fixed, don't count startprob params in BIC
                count_startprob = not sticky
                bic_k = _bic(m, X_all, lengths, count_startprob=count_startprob)
                logL_k = m.score(X_all, lengths=lengths)

                rows.append({"k": k, "bic": bic_k, "logL": logL_k})

                if bic_k < best_bic:
                    best_bic = bic_k
                    best_model = m

            model = best_model
            selection_df = pd.DataFrame(rows)
        else:
            model = _build_model(int(n_states))
            model.fit(X_all, lengths=lengths)
            if hasattr(model, "monitor_") and not model.monitor_.converged:
                warnings.warn(
                    f"HMM (n_states={n_states}) stopped after reaching n_iter={n_iter} without converging. "
                    "The result may be suboptimal — consider increasing n_iter or relaxing tol.",
                    UserWarning,
                    stacklevel=2,
                )

    # ---------- decode ----------
    states_all = model.predict(X_all, lengths=lengths) + 1  # 1-based
    out = df_sorted.copy()
    out[out_col_name] = states_all

    state_map = out[list(id_cols) + [time_col, out_col_name]]
    df_features_clean = df_features.drop(columns=[out_col_name], errors="ignore")

    df_out = df_features_clean.merge(
        state_map,
        on=list(id_cols) + [time_col],
        how="left",
        sort=False,
        validate="many_to_one",
    )
    return df_out, model, selection_df


def run_sticky_hmm_state_classification(
    df_features,
    feature_cols,
    model=None,
    id_cols=("sample_name", "TrackID"),
    time_col="position_t",
    n_states="auto",
    k_min=1,
    k_max=8,
    covariance_type="full",
    n_iter=200,
    tol=1e-3,
    random_state=0,
    verbose=False,
    error_on_nans=True,
    out_col_name="hmm_state",
    # sticky defaults
    stickiness_kappa=8.0,
    transmat_alpha=1.0,
    min_covar=1e-3,
):
    """
    Convenience wrapper around run_hmm_state_classification with sticky=True and
    default sticky hyperparameters.
    """
    return run_hmm_state_classification(
        df_features=df_features,
        feature_cols=feature_cols,
        model=model,
        id_cols=id_cols,
        time_col=time_col,
        n_states=n_states,
        k_min=k_min,
        k_max=k_max,
        covariance_type=covariance_type,
        n_iter=n_iter,
        tol=tol,
        random_state=random_state,
        verbose=verbose,
        error_on_nans=error_on_nans,
        out_col_name=out_col_name,
        sticky=True,
        stickiness_kappa=stickiness_kappa,
        transmat_alpha=transmat_alpha,
        min_covar=min_covar,
    )


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
            raise ValueError(
                "Non-positive definite covariance encountered. Increase min_covar/eps or use diag."
            )
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
