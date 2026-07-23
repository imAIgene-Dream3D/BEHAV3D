import os
import warnings
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/behav3d-mpl-cache")

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg", force=True)

from behav3d.analysis.behavior.state.utils import (
    _get_classification_state_order,
    _set_classification_state_colors,
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
from behav3d.analysis.behavior.state.visualization.plots import state_composition as state_composition_mod
from behav3d.analysis.behavior.state.visualization.plots.state_composition import (
    save_state_composition_report,
    save_state_condition_comparison_report,
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


def test_get_classification_state_order_survives_h5ad_roundtrip(tmp_path):
    """`AnnData.write`/`read_h5ad` silently turns a plain Python list stored in `.uns` into a
    numpy array. `_get_classification_state_order` must still recognize it as a valid saved
    order (previously it only checked `isinstance(order, (list, tuple))`, which excludes
    `numpy.ndarray` and caused every state-order lookup on a freshly-loaded-from-disk adata to
    silently return `[]`, even though the order was saved correctly)."""
    obs = _build_obs()
    adata = ad.AnnData(obs=obs)
    _set_classification_state_order(adata, _STATE_COL, _CUSTOM_ORDER)

    h5ad_path = tmp_path / "roundtrip.h5ad"
    adata.write(h5ad_path)
    reloaded = ad.read_h5ad(h5ad_path)

    assert isinstance(reloaded.uns["classification"]["state_order"][_STATE_COL], np.ndarray)
    assert _get_classification_state_order(reloaded, _STATE_COL) == _CUSTOM_ORDER


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


_CUSTOM_COLORS = {"mu": "#111111", "zeta": "#222222", "alpha": "#333333"}


def test_save_state_composition_report_self_resolves_saved_order_and_colors(tmp_path, monkeypatch):
    """`state_order`/`state_colors` are optional — when omitted, the report must resolve them
    from `adata.uns["classification"]` instead of falling back to alphabetical order / hash
    colors."""
    obs = _build_obs()
    adata = ad.AnnData(obs=obs)
    _set_classification_state_order(adata, _STATE_COL, _CUSTOM_ORDER)
    _set_classification_state_colors(adata, _STATE_COL, _CUSTOM_COLORS)

    captured_colors = {}
    real_normalize = state_composition_mod._normalize_label_color_map

    def _spy_normalize(labels, colors=None, cmap_name="tab20"):
        captured_colors["colors"] = colors
        return real_normalize(labels, colors=colors, cmap_name=cmap_name)

    monkeypatch.setattr(state_composition_mod, "_normalize_label_color_map", _spy_normalize)

    report_out = save_state_composition_report(
        adata=adata,
        output_pdf_path=tmp_path / "state_composition_report.pdf",
        output_csv_path=tmp_path / "state_composition_report_auc.csv",
        state_col=_STATE_COL,
        sample_col="sample_name",
        time_col="position_t",
        include_pooled_summary=False,
        verbose=False,
    )

    assert captured_colors["colors"] == _CUSTOM_COLORS

    auc = pd.read_csv(report_out["csv_path"])
    first_sample_states = auc[auc["sample_name"] == "s1"]["state_id"].tolist()
    assert first_sample_states == _CUSTOM_ORDER


def test_save_state_condition_comparison_report_self_resolves_saved_order_and_colors(tmp_path, monkeypatch):
    obs = pd.DataFrame(
        {
            "sample_name": ["s1"] * 3 + ["s2"] * 3 + ["s3"] * 3 + ["s4"] * 3,
            "TrackID": [1] * 12,
            "condition": ["A"] * 6 + ["B"] * 6,
            _STATE_COL: (["zeta", "alpha", "mu"]) * 4,
        }
    )
    adata = ad.AnnData(obs=obs)
    _set_classification_state_order(adata, _STATE_COL, _CUSTOM_ORDER)
    _set_classification_state_colors(adata, _STATE_COL, _CUSTOM_COLORS)

    captured_colors = {}
    real_normalize = state_composition_mod._normalize_label_color_map

    def _spy_normalize(labels, colors=None, cmap_name="tab20"):
        captured_colors["colors"] = colors
        return real_normalize(labels, colors=colors, cmap_name=cmap_name)

    monkeypatch.setattr(state_composition_mod, "_normalize_label_color_map", _spy_normalize)

    save_state_condition_comparison_report(
        adata=adata,
        output_pdf_path=tmp_path / "condition_comparison.pdf",
        output_csv_path=tmp_path / "condition_comparison.csv",
        state_col=_STATE_COL,
        sample_col="sample_name",
        condition_col="condition",
        verbose=False,
    )

    assert captured_colors["colors"] == _CUSTOM_COLORS


def test_save_state_transition_report_self_serves_saved_state_order(tmp_path):
    """`state_order`/`state_colors` are optional — when omitted, the report must resolve
    them itself from `adata.uns["classification"]` instead of silently falling back to
    alphabetical order and default colors."""
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
        verbose=False,
    )

    counts_csv = tmp_path / "transition_matrix" / "data" / "transition_matrix_counts.csv"
    assert counts_csv.exists()
    counts = pd.read_csv(counts_csv, index_col=0)
    assert counts.index.tolist() == _CUSTOM_ORDER
    assert counts.columns.tolist() == _CUSTOM_ORDER
