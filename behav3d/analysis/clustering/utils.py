import re
import time
from pathlib import Path

import numpy as np


def _vinfo(verbose, prefix, message):
    if not bool(verbose):
        return
    print(f"[{prefix}] INFO {str(message)}")


def _vstart(verbose, prefix, step_name):
    if bool(verbose):
        print(f"[{prefix}] START {step_name}")
    return time.perf_counter()


def _vdone(verbose, prefix, step_name, t_start):
    if bool(verbose):
        dt = time.perf_counter() - float(t_start)
        print(f"[{prefix}] DONE {step_name} | took {dt:.2f}s")


def _vsave(verbose, prefix, label, path):
    if not bool(verbose):
        return
    print(f"[{prefix}] SAVED {str(label)}: {path}")


def _resolve_output_dir(output_dir):
    if output_dir is None:
        raise ValueError("output_dir is required.")
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    return output_dir_path


def _mixed_label_sort_key(value):
    text = str(value)
    if re.fullmatch(r"-?\d+", text):
        return (0, int(text))
    return (1, text)


def _sanitize_filename_token(value, fallback="value"):
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip())
    token = token.strip("._-")
    return token if len(token) > 0 else str(fallback)


def _to_numpy_2d(X):
    if hasattr(X, "toarray"):
        X = X.toarray()
    arr = np.asarray(X)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    if arr.ndim != 2:
        raise ValueError(f"Expected 2D matrix, got shape={arr.shape}.")
    return arr
