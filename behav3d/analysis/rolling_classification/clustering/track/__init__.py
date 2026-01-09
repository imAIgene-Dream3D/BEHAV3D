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


