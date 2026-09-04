"""Value-based column classification for track/state feature CSVs.

Deliberately a leaf module: it depends on nothing but ``math`` and ``pandas``.
These helpers used to live in :mod:`behav3d.widgets.base_state_classification`,
whose imports pull in the whole clustering stack (scanpy -> umap -> pynndescent),
costing ~12 s of numba JIT compilation on first import. The napari Single Cell
tab calls them on every metadata load, so that cost landed on the Qt main thread
and froze the GUI. They are re-exported from the old module for backwards
compatibility.
"""
import math

import pandas as pd


def normalize_binary_value(value, tol=1e-9):
    """Return 0/1 if ``value`` normalizes to a boolean, else ``None``.

    Accepts native bools, the strings true/false/t/f, and numeric 0/1.
    Shared by the notebook and napari binary-column detectors so both use
    identical semantics.
    """
    if isinstance(value, (bool, pd.BooleanDtype)):
        return 1 if bool(value) else 0
    if isinstance(value, str):
        sval = value.strip().lower()
        if sval in {"true", "t"}:
            return 1
        if sval in {"false", "f"}:
            return 0
        try:
            value = float(sval)
        except Exception:
            return None
    try:
        fval = float(value)
    except Exception:
        return None
    if not math.isfinite(fval):
        return None
    if abs(fval - 0.0) <= tol:
        return 0
    if abs(fval - 1.0) <= tol:
        return 1
    return None


def detect_binary_columns_from_csv(csv_path, cols, chunksize=50000):
    """Value-based binary-column detection over the *full* CSV.

    A column is binary only if every non-NA value normalizes (via
    :func:`normalize_binary_value`) to ``{0, 1}``. This replaces fragile
    dtype sniffing over a small row sample, which mis-classifies numeric
    columns as binary whenever the sampled rows happen to be NaN/blank.
    """
    if csv_path is None or len(cols) == 0:
        return []

    states = {str(c): {"seen": set(), "invalid": False} for c in cols}
    try:
        for chunk in pd.read_csv(csv_path, usecols=cols, chunksize=chunksize, low_memory=False):
            for col in cols:
                st = states[str(col)]
                if st["invalid"]:
                    continue
                series = chunk[col].dropna()
                if len(series) == 0:
                    continue
                unique_vals = pd.unique(series)
                for raw in unique_vals:
                    norm = normalize_binary_value(raw)
                    if norm is None:
                        st["invalid"] = True
                        break
                    st["seen"].add(int(norm))
                    if len(st["seen"]) > 2:
                        st["invalid"] = True
                        break
    except Exception:
        return []

    return sorted(
        [
            c
            for c in cols
            if (not states[str(c)]["invalid"]) and (len(states[str(c)]["seen"]) > 0)
        ]
    )


def detect_non_numeric_columns_from_csv(csv_path, cols, chunksize=50000):
    """Value-based detection of columns unsuitable as continuous HMM features.

    A column is flagged if any non-NA value fails to parse as a finite float --
    e.g. free-text/categorical labels ("um", "27t") or comma-separated contact-ID
    lists ("45,46"). Such object-dtype columns silently poison the HMM
    observation matrix's dtype when selected as a feature: pandas keeps the
    column as ``object`` even after numeric coercion, so ``adata.X`` ends up
    object-dtype and anndata's h5ad writer -- assuming an object array must be
    strings -- crashes with "Can't implicitly convert non-string objects to
    strings" the moment it hits the actual float values. These columns must be
    excluded before they are ever offered as selectable timepoint features.
    """
    if csv_path is None or len(cols) == 0:
        return []

    invalid = {str(c): False for c in cols}
    try:
        for chunk in pd.read_csv(csv_path, usecols=cols, chunksize=chunksize, low_memory=False):
            for col in cols:
                key = str(col)
                if invalid[key]:
                    continue
                series = chunk[col].dropna()
                if len(series) == 0:
                    continue
                if pd.to_numeric(series, errors="coerce").isna().any():
                    invalid[key] = True
    except Exception:
        return []

    return sorted([c for c in cols if invalid[str(c)]])
