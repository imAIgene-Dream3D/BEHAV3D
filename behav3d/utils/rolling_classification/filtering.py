def filter_and_truncate_tracks(
    adata,
    groupby_cols = ["sample_name", "TrackID"],  # <-- composite track identity
    time_col: str = "position_t",                                  # <-- time/frame column
    min_length: int = 10,                                   # x
    max_length: int = 200                                   # y
):
    missing = [k for k in groupby_cols if k not in adata.obs.columns]
    if missing:
        raise KeyError(f"Missing groupby_cols in adata.obs: {missing}")

    if time_col not in adata.obs.columns:
        raise KeyError(f"'{time_col}' not found in adata.obs. Set time_col to your time column.")

    obs = adata.obs.copy()

    # --- 1) Filter groups with too few positions ---
    # size per group (composite track key)
    group_sizes = obs.groupby(list(groupby_cols), observed=True).size()
    keep_groups = group_sizes[group_sizes >= min_length].index

    # fast membership test: join on multi-index
    obs["_keep"] = obs.set_index(list(groupby_cols)).index.isin(keep_groups)
    obs = obs.loc[obs["_keep"]].copy()

    # --- 2) Keep last max_length per group by time ordering ---
    # stable tie-breaker (in case time_col has ties)
    obs["_orig_idx"] = np.arange(len(obs))
    obs = obs.sort_values(list(groupby_cols) + [time_col, "_orig_idx"])

    # rank within group
    obs["_rank"] = obs.groupby(list(groupby_cols), observed=True).cumcount()
    sizes = obs.groupby(list(groupby_cols), observed=True)["_rank"].transform("max") + 1
    obs = obs.loc[obs["_rank"] >= (sizes - max_length)]

    # subset AnnData with rows in this exact order
    idx = obs.index
    adata_out = adata[idx].copy()

    # cleanup helper cols if they persisted into adata_out.obs
    for col in ["_keep", "_orig_idx", "_rank"]:
        if col in adata_out.obs.columns:
            adata_out.obs.drop(columns=[col], inplace=True)

    return adata_out