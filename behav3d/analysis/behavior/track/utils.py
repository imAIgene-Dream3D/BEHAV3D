from pathlib import Path

import numpy as np
import pandas as pd

from behav3d.analysis.behavior.utils import _sanitize_filename_token


def _winfo(prefix, message):
    print(f"[{prefix}] INFO {message}")


def _ordered_unique(values):
    out = []
    seen = set()
    for value in list(values or []):
        text = str(value).strip()
        if text == "" or text in seen:
            continue
        out.append(text)
        seen.add(text)
    return out


def _resolve_dtaidistance_paths(output_dir, cell_type, output_subdir_name="behavorial_trajectories"):
    root = Path(output_dir).expanduser()
    analysis_outdir = root / "analysis" / str(cell_type)
    state_outdir = analysis_outdir / "behavioral_states"
    outfolder = analysis_outdir / str(output_subdir_name)
    clustering_outfolder = outfolder / "clustering"
    quality_control_outfolder = outfolder / "quality_control"
    outfolder.mkdir(parents=True, exist_ok=True)
    clustering_outfolder.mkdir(parents=True, exist_ok=True)
    quality_control_outfolder.mkdir(parents=True, exist_ok=True)
    return {
        "root": root,
        "analysis_outdir": analysis_outdir,
        "state_outdir": state_outdir,
        "outfolder": outfolder,
        "clustering_outfolder": clustering_outfolder,
        "quality_control_outfolder": quality_control_outfolder,
    }


def get_dtaidistance_track_trajectories_filename(cell_type):
    cell_token = _sanitize_filename_token(cell_type, fallback="cell")
    return f"BEHAV3D_{cell_token}_behavioral_trajectories_dtaidistance.h5ad"


def _default_behavioral_states_path(output_dir, cell_type):
    return (
        Path(output_dir).expanduser()
        / "analysis"
        / str(cell_type)
        / "behavioral_states"
        / f"BEHAV3D_{cell_type}_behavioral_states.h5ad"
    )


def _resolve_optional_int(value):
    text = str(value).strip()
    return None if text == "" else int(text)


def _resolve_optional_float(value):
    text = str(value).strip()
    return None if text == "" else float(text)


def _filter_tracks_for_dtaidistance(
    adata,
    *,
    groupby_cols=("sample_name", "TrackID"),
    time_col="position_t",
    trajectory_size=None,
    min_length=None,
    trim_mode="last",
):
    groupby_cols = [str(c) for c in list(groupby_cols)]
    missing = [col for col in groupby_cols if col not in adata.obs.columns]
    if missing:
        raise KeyError(f"Missing groupby_cols in adata.obs: {missing}")
    if str(time_col) not in adata.obs.columns:
        raise KeyError(f"'{time_col}' not found in adata.obs.")

    obs = adata.obs.copy()
    obs["_orig_idx"] = np.arange(len(obs))
    obs = obs.sort_values(groupby_cols + [str(time_col), "_orig_idx"])

    if min_length is not None:
        group_sizes = obs.groupby(groupby_cols, observed=True).size()
        keep_groups = group_sizes[group_sizes >= int(min_length)].index
        obs["_keep"] = obs.set_index(groupby_cols).index.isin(keep_groups)
        obs = obs.loc[obs["_keep"]].copy()

    if trajectory_size is not None:
        trajectory_size = int(trajectory_size)
        if trajectory_size <= 0:
            raise ValueError("trajectory_size must be positive when supplied.")
        trim_mode = str(trim_mode).strip().lower()
        if trim_mode not in {"first", "last"}:
            raise ValueError("trim_mode must be 'first' or 'last'.")
        obs["_rank"] = obs.groupby(groupby_cols, observed=True).cumcount()
        sizes = obs.groupby(groupby_cols, observed=True)["_rank"].transform("max") + 1
        if trim_mode == "first":
            obs = obs.loc[obs["_rank"] < trajectory_size]
        else:
            obs = obs.loc[obs["_rank"] >= (sizes - trajectory_size)]

    idx = obs.index
    adata_out = adata[idx].copy()
    for col in ["_orig_idx", "_keep", "_rank"]:
        if col in adata_out.obs.columns:
            adata_out.obs.drop(columns=[col], inplace=True)
    return adata_out


def _build_identity_cluster_mapping_from_obs(obs, cluster_col="ClusterID"):
    if cluster_col not in obs.columns:
        return {}
    values = pd.Series(obs[cluster_col]).dropna().astype(str).unique().tolist()
    values = sorted(values, key=lambda x: (0, int(x)) if str(x).isdigit() else (1, str(x)))
    return {str(value): str(value) for value in values}
