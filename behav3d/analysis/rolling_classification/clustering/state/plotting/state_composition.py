import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_state_composition_over_time(
    adata,
    *,
    time_col: str = "position_t",
    state_col: str = "hmm_state",
    sample_col: str = "sample_name",
    relative: bool = True,
    group_by_sample: bool = False,
    state_order=None,
    figsize=(8, 4),
    show: bool = True,
    ax=None,
    legend: bool = True,
):
    """
    Compute (and optionally plot) continuous stacked state composition over time.

    Returns
    -------
    If show == False:
        dict of DataFrames:
            {sample_name -> DataFrame(time x state)}
            or {"all": DataFrame} if not grouped
    If show == True:
        (data_dict, fig, ax/axes)
    """

    obs = adata.obs

    # --- validate ---
    for c in [time_col, state_col]:
        if c not in obs.columns:
            raise KeyError(f"{c} not found in adata.obs")
    if group_by_sample and sample_col not in obs.columns:
        raise KeyError(f"{sample_col} not found in adata.obs")

    # --- clean ---
    cols = [time_col, state_col] + ([sample_col] if group_by_sample else [])
    df = obs[cols].dropna().copy()
    df[time_col] = pd.to_numeric(df[time_col], errors="coerce")
    df = df.dropna(subset=[time_col])
    df[time_col] = df[time_col].astype(int)

    if state_order is None:
        state_order = (
            df[state_col]
            .astype(str)
            .value_counts()
            .index
            .tolist()
        )

    # --- compute per-sample matrices ---
    def _compute(panel_df):
        panel_df = panel_df.copy()
        panel_df[state_col] = panel_df[state_col].astype(str)

        mat = (
            panel_df
            .groupby([time_col, state_col], observed=True)
            .size()
            .unstack(fill_value=0)
            .sort_index()
        )

        for s in state_order:
            if s not in mat.columns:
                mat[s] = 0
        mat = mat[state_order]

        if relative:
            mat = mat.div(
                mat.sum(axis=1).replace(0, np.nan),
                axis=0
            ).fillna(0.0)

        return mat

    if group_by_sample:
        data = {
            str(s): _compute(df[df[sample_col].astype(str) == str(s)])
            for s in df[sample_col].astype(str).unique()
        }
    else:
        data = {"all": _compute(df)}

    # --- early exit if no plotting ---
    if not show:
        return data

    # --- plotting ---
    def _plot_one(mat, ax, title=None):
        x = mat.index.to_numpy()
        bottom = np.zeros(len(x))

        for s in mat.columns:
            v = mat[s].to_numpy()
            ax.bar(
                x,
                v,
                bottom=bottom,
                width=1.0,
                align="edge",
                linewidth=0,
                label=str(s),
            )
            bottom += v

        ax.set_xlim(x.min(), x.max() + 1)
        ax.margins(x=0)
        ax.set_xlabel(time_col)
        ax.set_ylabel("Proportion" if relative else "Count")
        if title:
            ax.set_title(title)

    if group_by_sample:
        n = len(data)
        fig, axes = plt.subplots(
            nrows=n,
            ncols=1,
            sharex=True,
            figsize=(figsize[0], figsize[1] * n),
        )
        if n == 1:
            axes = [axes]

        for ax_i, (samp, mat) in zip(axes, data.items()):
            _plot_one(mat, ax_i, title=f"{sample_col}: {samp}")

        if legend:
            handles, labels = axes[0].get_legend_handles_labels()
            fig.legend(
                handles,
                labels,
                loc="lower center",
                ncol=min(len(labels), 6),
                frameon=False,
            )
            fig.tight_layout(rect=(0, 0.06, 1, 1))
        else:
            fig.tight_layout()

        return data, fig, axes

    # single axis
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    _plot_one(data["all"], ax)

    if legend:
        ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.02, 1),
            frameon=False,
        )
        fig.tight_layout()

    return data, fig, ax


df_fig, fig, ax = plot_state_composition_over_time(
    adata_full, 
    time_col="position_t", 
    state_col="ClusterID", 
    relative=False
    )

df_fig, fig, axes= plot_state_composition_over_time(
    adata_full, 
    time_col="position_t", 
    state_col="ClusterID", 
    relative=True,
    group_by_sample=True
    )

fig, axes = plot_state_composition_over_time(adata, group_by_sample=True, relative=True)
plt.show()