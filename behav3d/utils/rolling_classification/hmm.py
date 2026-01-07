import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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