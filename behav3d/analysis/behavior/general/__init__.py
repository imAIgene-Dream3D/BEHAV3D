import numpy as np
import pandas as pd

from behav3d.analysis.behavior.utils import _mixed_label_sort_key


def relabel_cluster_ids(
    adata,
    mapping,
    cluster_key="ClusterID",
    new_key=None,
    keep_unmapped=True,
    unmapped_label="unlabeled",
    overwrite_original=False,
    categories=None,
):
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"{cluster_key} not found in adata.obs")

    original_key=f"{cluster_key}_original"
    if (original_key in adata.obs.columns) and (not overwrite_original):
        raise ValueError(f"{original_key} already exists")

    adata.obs[original_key] = adata.obs[cluster_key].astype(str).copy()
    current = adata.obs[cluster_key].astype(str)

    if isinstance(mapping, dict):
        map_dict = {str(k): v for k, v in mapping.items()}
        mapped = current.map(map_dict)
    else:
        uniq = np.array(sorted(current.unique(), key=lambda x: (int(x) if x.isdigit() else x)))
        labels = list(mapping)
        if len(labels) < len(uniq):
            raise ValueError("Not enough labels for number of clusters")
        map_dict = {uniq[i]: labels[i] for i in range(len(uniq))}
        mapped = current.map(map_dict)

    if keep_unmapped:
        out = mapped.where(mapped.notna(), current)
    else:
        out = mapped.fillna(unmapped_label)

    out_col = new_key if new_key is not None else cluster_key
    present = sorted({str(x) for x in out.dropna().unique()}, key=_mixed_label_sort_key)
    if categories is not None:
        requested = [str(c) for c in categories]
        present_set = set(present)
        cats = [c for c in requested if c in present_set] + [
            c for c in present if c not in requested
        ]
    else:
        cats = present
    adata.obs[out_col] = pd.Categorical(out, categories=cats)

    return adata