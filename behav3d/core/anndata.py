import math
from typing import Optional, Tuple, Dict, Literal, List, Iterable
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
from scipy import sparse
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, leaves_list

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

import scanpy as sc
import umap
import igraph as ig
import leidenalg as la
from tqdm import tqdm
import imageio_ffmpeg as iioff

from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler
from sklearn.cluster import KMeans, HDBSCAN, AgglomerativeClustering
from sklearn.feature_selection import VarianceThreshold
from sklearn.impute import SimpleImputer

from behav3d.io.images import load_zarr
from behav3d.preprocessing import calc_z_projection


def relabel_img_from_adata(
    label_image: np.ndarray,
    adata,
    sample_name,
    obs_col: str,
    trackid_col: str = "TrackID",
    time_col: str = "position_t",
    sample_col = "sample_name",
    default_value=0,
    numeric_from_category: bool = True,
):
    cols = [trackid_col, time_col, obs_col]
    if sample_col is not None:
        cols.append(sample_col)

    df = adata.obs[cols].copy()

    df[obs_col] = df[obs_col].astype(np.int64)
    df[trackid_col] = df[trackid_col].astype(np.int64)
    df[time_col] = df[time_col].astype(np.int64)
    
    if sample_col is not None:
        if sample_name is None:
            raise ValueError(
                "sample_col is provided but sample_value is None. "
                "Please specify which sample to select."
            )
        df = df[df[sample_col] == sample_name]

    # If nothing left after filtering, return array filled with default_value
    if df.empty:
        value_dtype = np.array([default_value]).dtype
        return np.full(label_image.shape, default_value, dtype=value_dtype)
   
    value_dtype = np.result_type(df[obs_col].dtype, np.array([default_value]).dtype)

    relabeled = np.empty(label_image.shape, dtype=value_dtype)

    T = label_image.shape[0]
    grouped = df.groupby(time_col)

    for t in range(T):
        labels_t = label_image[t].astype(np.int64, copy=False)
        max_label_t = int(labels_t.max())

        lut = np.full(max_label_t + 1, default_value, dtype=value_dtype)

        if t in grouped.groups:
            df_t = grouped.get_group(t)
            track_ids = df_t[trackid_col].to_numpy()
            values = df_t[obs_col].to_numpy()

            for tid, val in zip(track_ids, values):
                tid_int = int(tid)
                if 0 <= tid_int <= max_label_t:
                    lut[tid_int] = val
        relabeled[t] = lut[labels_t]

    return relabeled

def df_to_adata(df, feature_cols, obs_cols=None):
    """Create AnnData from df, store metadata in .obs."""
    X = df[feature_cols].to_numpy()
    adata = sc.AnnData(X)
    adata.var_names = feature_cols
    if obs_cols is not None:
        adata.obs = df[obs_cols].copy()
    return adata

def adata_add_back_to_df(df, adata, cols_from_obs, prefix=None):
    """Copy selected .obs columns back into df."""
    for c in cols_from_obs:
        newc = f"{prefix}{c}" if prefix else c
        df[newc] = adata.obs[c].to_numpy()
    return df

def adata_to_df(adata, layer=None, dense=True, obs_cols=None, prefix=None):
    X = adata.X if layer is None else adata.layers[layer]

    feature_cols = adata.var_names.to_list()
    if prefix:
        feature_cols = [f"{prefix}{c}" for c in feature_cols]

    out = pd.DataFrame(X, columns=feature_cols, index=adata.obs_names)

    if obs_cols is not None:
        out = pd.concat([adata.obs[obs_cols].copy(), out], axis=1)
    return out

def merge_pandas_cols_into_obs_anndata(
    cols,
    adata,
    df_analysis,
    on=("sample_name", "TrackID", "position_t"),
    obs_index_col="_obs_index",
    how="left",
    make_category=True,
    astype_str=True,
    inplace=True,
    ):
    if isinstance(cols, str):
        cols = [cols]
    else:
        cols = list(cols)

    missing_in_obs = [c for c in on if c not in adata.obs.columns]
    missing_in_df  = [c for c in on if c not in df_analysis.columns]
    if missing_in_obs:
        raise KeyError("Keys missing in adata.obs: %s" % missing_in_obs)
    if missing_in_df:
        raise KeyError("Keys missing in df_analysis: %s" % missing_in_df)

    missing_cols = [c for c in cols if c not in df_analysis.columns]
    if missing_cols:
        raise KeyError("Columns missing in df_analysis: %s" % missing_cols)

    target = adata if inplace else adata.copy()

    lookup = df_analysis[list(on) + cols].copy()
    obs_df = target.obs.reset_index().rename(columns={"index": obs_index_col})
    obs_df = obs_df.drop(columns=[c for c in cols if c in obs_df.columns], errors="ignore")
    merged = obs_df.merge(lookup, on=list(on), how=how)
    merged = merged.set_index(obs_index_col).loc[target.obs.index, cols]

    for c in cols:
        s = merged[c]
        if astype_str:
            s = s.astype(str)
        if make_category:
            s = s.astype("category")
        target.obs[c] = s

    if not inplace:
        return target
  