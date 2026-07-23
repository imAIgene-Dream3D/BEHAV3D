import os
import warnings
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/behav3d-mpl-cache")

import anndata as ad
import matplotlib
import pandas as pd

matplotlib.use("Agg", force=True)

from behav3d.analysis.behavior.state.utils import (
    _get_classification_state_order,
    _set_classification_state_order,
)
from behav3d.analysis.behavior.state.visualization.backprojection import (
    _build_code_map,
    export_behavioral_state_backprojection_zarrs,
)
from behav3d.analysis.behavior.state.visualization.plots.state_transitions import (
    compute_cluster_transition_matrix,
    save_state_transition_report,
)

_STATE_COL = "behavioral_state"
_CUSTOM_ORDER = ["mu", "zeta", "alpha"]


def _build_obs():
    return pd.DataFrame(
        {
            "sample_name": ["s1"] * 6,
            "TrackID": [1, 1, 1, 2, 2, 2],
            "position_t": [0, 1, 2, 0, 1, 2],
            _STATE_COL: ["zeta", "alpha", "mu", "zeta", "alpha", "mu"],
        }
    )


def test_build_code_map_defaults_to_alphabetical_order():
    obs = _build_obs()
    code_map = _build_code_map(obs, state_col=_STATE_COL)
    assert list(code_map.keys()) == ["alpha", "mu", "zeta"]


def test_build_code_map_honors_saved_state_order():
    obs = _build_obs()
    code_map = _build_code_map(obs, state_col=_STATE_COL, state_order=_CUSTOM_ORDER)
    assert list(code_map.keys()) == _CUSTOM_ORDER
    assert code_map == {"mu": 1, "zeta": 2, "alpha": 3}


def test_export_behavioral_state_backprojection_zarrs_resolves_saved_order(tmp_path):
    obs = _build_obs()
    adata = ad.AnnData(obs=obs)
    _set_classification_state_order(adata, _STATE_COL, _CUSTOM_ORDER)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        manifest = export_behavioral_state_backprojection_zarrs(
            adata=adata,
            output_dir=tmp_path,
            cell_type="tcell",
            state_col=_STATE_COL,
            n_workers=1,
            verbose=False,
        )

    # No raw/tracked images exist under tmp_path, so every sample is skipped, but the
    # label->code mapping is still computed up front and must follow the saved order.
    assert manifest["label_to_code"] == {"mu": 1, "zeta": 2, "alpha": 3}
    assert manifest["label_map"] == {"1": "mu", "2": "zeta", "3": "alpha"}


def test_compute_cluster_transition_matrix_honors_state_order():
    obs = _build_obs()
    adata = ad.AnnData(obs=obs)
    counts, probs = compute_cluster_transition_matrix(
        adata, cluster_key=_STATE_COL, state_order=_CUSTOM_ORDER,
    )
    assert counts.index.tolist() == _CUSTOM_ORDER
    assert counts.columns.tolist() == _CUSTOM_ORDER
    assert probs.index.tolist() == _CUSTOM_ORDER
    assert probs.columns.tolist() == _CUSTOM_ORDER


def test_compute_cluster_transition_matrix_defaults_to_alphabetical():
    obs = _build_obs()
    adata = ad.AnnData(obs=obs)
    counts, _ = compute_cluster_transition_matrix(adata, cluster_key=_STATE_COL)
    assert counts.index.tolist() == ["alpha", "mu", "zeta"]


def test_save_state_transition_report_honors_saved_state_order(tmp_path):
    obs = _build_obs()
    adata = ad.AnnData(obs=obs)
    _set_classification_state_order(adata, _STATE_COL, _CUSTOM_ORDER)

    save_state_transition_report(
        adata=adata,
        output_dir=tmp_path,
        state_col=_STATE_COL,
        id_cols=("sample_name", "TrackID"),
        time_col="position_t",
        include_ngram_rankings=False,
        state_order=_get_classification_state_order(adata, _STATE_COL),
        verbose=False,
    )

    counts_csv = tmp_path / "transition_matrix" / "data" / "transition_matrix_counts.csv"
    assert counts_csv.exists()
    counts = pd.read_csv(counts_csv, index_col=0)
    assert counts.index.tolist() == _CUSTOM_ORDER
    assert counts.columns.tolist() == _CUSTOM_ORDER
