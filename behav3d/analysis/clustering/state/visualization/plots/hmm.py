import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

def hmm_state_feature_zscores(df_out, feature_cols, model, out_col="hmm_state"):
    """
    Returns:
      z_df: (K x D) z-scores per state/feature
      mu_df: (K x D) raw state means
    """
    # hmm_state is 1-based in your pipeline; hmmlearn is 0-based internally.
    K, D = model.means_.shape

    mu = model.means_  # (K, D)

    X = df_out[feature_cols].to_numpy(float)
    mu_global = np.nanmean(X, axis=0)
    sd_global = np.nanstd(X, axis=0, ddof=1)
    sd_global = np.where(sd_global == 0, 1.0, sd_global)

    z = (mu - mu_global) / sd_global

    states = np.arange(1, K+1)
    z_df  = pd.DataFrame(z,  index=states, columns=feature_cols)
    mu_df = pd.DataFrame(mu, index=states, columns=feature_cols)
    return z_df, mu_df

def top_hmm_features_per_state(z_df, top_n=5):
    rows = []
    for s in z_df.index:
        ser = z_df.loc[s]
        top = ser.reindex(ser.abs().sort_values(ascending=False).index).head(top_n)
        for feat, val in top.items():
            rows.append({"state": s, "feature": feat, "z": float(val)})
    return pd.DataFrame(rows)

def plot_hmm_top_ranking_features(
    df_out: pd.DataFrame,
    feature_cols: list[str],
    model,
    out_col: str = "hmm_state",
    n_features: int = 15,
    rank_by: str = "abs",   # "abs" (default) or "pos" or "neg"
):
    """
    Like Scanpy rank_genes_groups plot, but for hmmlearn GaussianHMM emissions.

    We rank features per state using z-scored emission means:
        z_{k,j} = (mu_{k,j} - mu_global_j) / sd_global_j

    Blue labels = high vs global, Red labels = low vs global.
    """

    # ---- sanity checks ----
    if not hasattr(model, "means_"):
        raise ValueError("Model must have means_ (fitted hmmlearn GaussianHMM).")
    K, D = model.means_.shape
    if D != len(feature_cols):
        raise ValueError(
            f"model.means_ has D={D} features, but feature_cols has {len(feature_cols)}. "
            "Ensure the same order/columns were used for fitting."
        )
    for c in feature_cols:
        if c not in df_out.columns:
            raise ValueError(f"Feature column missing from df_out: {c}")

    # ---- compute z-scores of emission means ----
    X = df_out[feature_cols].to_numpy(float)
    mu_global = np.nanmean(X, axis=0)
    sd_global = np.nanstd(X, axis=0, ddof=1)
    sd_global = np.where(sd_global == 0, 1.0, sd_global)

    mu_states = np.asarray(model.means_)  # (K, D)
    z = (mu_states - mu_global) / sd_global
    z_df = pd.DataFrame(z, index=[str(i) for i in range(1, K + 1)], columns=feature_cols)

    # ---- grid layout ----
    n_groups = K
    n_cols = math.ceil(np.sqrt(n_groups))
    n_rows = math.ceil(n_groups / n_cols)

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(4.5 * n_cols, 5 * n_rows),
        squeeze=False,
    )

    # ---- plot top-ranked features per state ----
    for ax, state in zip(axes.flat, z_df.index):
        ser = z_df.loc[state].copy()

        if rank_by == "abs":
            feats = ser.abs().sort_values(ascending=False).head(n_features).index
            top = ser.loc[feats].reindex(ser.loc[feats].abs().sort_values().index)  # reverse for barh
        elif rank_by == "pos":
            top = ser.sort_values(ascending=False).head(n_features).sort_values()
        elif rank_by == "neg":
            top = ser.sort_values(ascending=True).head(n_features).sort_values()
        else:
            raise ValueError("rank_by must be one of: 'abs', 'pos', 'neg'")

        ax.barh(top.index, top.values)
        ax.axvline(0, linewidth=0.8)
        ax.set_title(f"State {state}")
        ax.set_xlabel("emission mean z-score")

        # Color labels by direction
        for label, score in zip(ax.get_yticklabels(), top.values):
            label.set_color("blue" if score > 0 else "red")

    # ---- disable unused axes ----
    for ax in axes.flat[len(z_df.index):]:
        ax.axis("off")

    # ---- shared legend ----
    legend_elements = [
        Line2D([0], [0], color="blue", lw=2, label="High vs global"),
        Line2D([0], [0], color="red", lw=2, label="Low vs global"),
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=2, frameon=False)

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.show()

def plot_hmm_top_ranking_features(
    df_out,
    feature_cols,
    model,
    out_col="hmm_state",
    n_features=15,
    rank_by="abs",  # "abs", "pos", "neg"
):
    z_df, _ = hmm_state_feature_zscores(
        df_out=df_out,
        feature_cols=feature_cols,
        model=model,
        out_col=out_col,
    )

    states = z_df.index.tolist()
    n_groups = len(states)

    n_cols = math.ceil(np.sqrt(n_groups))
    n_rows = math.ceil(n_groups / n_cols)

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(4.5 * n_cols, 5 * n_rows),
        squeeze=False,
    )

    # ---- global x-axis limit (shared) ----
    if rank_by == "abs":
        max_val = z_df.abs().values.max()
    elif rank_by == "pos":
        max_val = z_df.values.max()
    elif rank_by == "neg":
        max_val = abs(z_df.values.min())
    else:
        raise ValueError("rank_by must be 'abs', 'pos', or 'neg'")

    xlim = (-max_val * 1.05, max_val * 1.05)

    for ax, state in zip(axes.flat, states):
        ser = z_df.loc[state].copy()

        if rank_by == "abs":
            feats = ser.abs().sort_values(ascending=False).head(n_features).index
            top = ser.loc[feats].reindex(
                ser.loc[feats].abs().sort_values().index
            )
        elif rank_by == "pos":
            top = ser.sort_values(ascending=False).head(n_features).sort_values()
        elif rank_by == "neg":
            top = ser.sort_values(ascending=True).head(n_features).sort_values()

        ax.barh(top.index, top.values)
        ax.axvline(0, linewidth=0.8)
        ax.set_xlim(xlim)

        ax.set_title(f"State {state}")
        ax.set_xlabel("emission mean z-score")

        for label, score in zip(ax.get_yticklabels(), top.values):
            label.set_color("blue" if score > 0 else "red")

    for ax in axes.flat[len(states):]:
        ax.axis("off")

    legend_elements = [
        Line2D([0], [0], color="blue", lw=2, label="High vs global"),
        Line2D([0], [0], color="red", lw=2, label="Low vs global"),
    ]

    fig.legend(
        handles=legend_elements,
        loc="lower center",
        ncol=2,
        frameon=False,
    )

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.show()
 
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
