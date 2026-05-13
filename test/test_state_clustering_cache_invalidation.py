from pathlib import Path
from types import SimpleNamespace
import sys

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import yaml

import behav3d.analysis.clustering.state.classification as state_classification
import behav3d.widgets.state_classification as state_classification_widget
from behav3d.analysis.clustering.state.visualization import backprojection as state_backprojection
from behav3d.features import rolling_window_features
from behav3d.analysis.clustering.state.classification import (
    _resolve_state_paths,
    prepare_state_classification_dataset,
    run_state_clustering,
)
from behav3d.analysis.clustering.state.classification import run_hmm_state_clustering


def _load_hmm_widget_module():
    pytest.importorskip("ipywidgets")
    pytest.importorskip("dask.array")
    return state_classification_widget


def _make_positions_df(n_tracks=6, track_len=15):
    rows = []
    for track_id in range(n_tracks):
        group_id = track_id % 3
        sample_name = f"sample_{track_id // 3}"
        for t in range(track_len):
            rows.append(
                {
                    "sample_name": sample_name,
                    "TrackID": track_id,
                    "position_t": t,
                    "position_x": float(t + 0.15 * track_id),
                    "position_y": float(group_id * 2.5 + np.sin(t / 3.0) + 0.05 * track_id),
                    "position_z": float((t % 4) * 0.4 + group_id * 0.2),
                    "speed": float(group_id + 0.25 * t + 0.03 * ((-1) ** t)),
                    "elongation": float(1.0 + group_id * 0.4 + 0.08 * np.cos(t / 2.0) + 0.01 * track_id),
                }
            )
    return pd.DataFrame(rows)


def _make_positions_df_variable_lengths(track_lengths):
    rows = []
    for track_id, track_len in enumerate(list(track_lengths)):
        group_id = track_id % 3
        sample_name = f"sample_{track_id // 3}"
        for t in range(int(track_len)):
            rows.append(
                {
                    "sample_name": sample_name,
                    "TrackID": track_id,
                    "position_t": t,
                    "position_x": float(t + 0.15 * track_id),
                    "position_y": float(group_id * 2.5 + np.sin(t / 3.0) + 0.05 * track_id),
                    "position_z": float((t % 4) * 0.4 + group_id * 0.2),
                    "speed": float(group_id + 0.25 * t + 0.03 * ((-1) ** t)),
                    "elongation": float(1.0 + group_id * 0.4 + 0.08 * np.cos(t / 2.0) + 0.01 * track_id),
                }
            )
    return pd.DataFrame(rows)


def _with_binary_contact_flag(df):
    out = df.copy()
    out["contact_flag"] = (
        (pd.to_numeric(out["TrackID"], errors="coerce").fillna(0).astype(int) % 2) == 0
    ).astype(int)
    return out


def _with_multi_binary_flags(df):
    out = df.copy()
    patterns = {
        0: (1, 0, 0),
        1: (0, 1, 0),
        2: (0, 0, 1),
        3: (1, 1, 0),
    }
    values = pd.to_numeric(out["TrackID"], errors="coerce").fillna(0).astype(int).map(
        lambda tid: patterns[int(tid) % 4]
    )
    out["contact_a"] = values.map(lambda x: int(x[0]))
    out["contact_b"] = values.map(lambda x: int(x[1]))
    out["contact_c"] = values.map(lambda x: int(x[2]))
    return out


def _std_cols(cols):
    return sorted([str(col) for col in list(cols) if str(col).endswith("_standard_deviation")])


def _force_serial_window_generation(monkeypatch):
    def _serial_create_descriptive_track_dataset(*args, **kwargs):
        kwargs.setdefault("n_jobs", 1)
        return rolling_window_features.create_descriptive_track_dataset(*args, **kwargs)

    monkeypatch.setattr(
        state_classification,
        "create_descriptive_track_dataset",
        _serial_create_descriptive_track_dataset,
    )


def test_prepare_cache_rebuilds_when_std_is_unticked(tmp_path, monkeypatch):
    _force_serial_window_generation(monkeypatch)
    df_positions = _make_positions_df()
    prepared_path = tmp_path / "prepared_full.h5ad"

    prepare_state_classification_dataset(
        df_positions=df_positions,
        features=["speed", "elongation"],
        binary_features_to_group=[],
        window_size=4,
        descriptive_features=["mean", "median", "std"],
        incomplete_window_policy="partial",
        scale_features=True,
        reuse_prepared_dataset=False,
        save_prepared_dataset=True,
        prepared_dataset_path=prepared_path,
        verbose=False,
    )

    adata_prepared = prepare_state_classification_dataset(
        df_positions=df_positions,
        features=["speed", "elongation"],
        binary_features_to_group=[],
        window_size=4,
        descriptive_features=["mean", "median"],
        incomplete_window_policy="partial",
        scale_features=True,
        reuse_prepared_dataset=True,
        save_prepared_dataset=True,
        prepared_dataset_path=prepared_path,
        verbose=False,
    )

    assert _std_cols(adata_prepared.var_names) == []
    assert adata_prepared.uns["preprocessing"]["windowing"]["descriptive_features"] == ["mean", "median"]

    saved = sc.read_h5ad(prepared_path)
    assert _std_cols(saved.var_names) == []
    assert list(saved.uns["preprocessing"]["windowing"]["descriptive_features"]) == ["mean", "median"]


def test_prepare_rebuilds_older_cache_missing_window_metadata(tmp_path, monkeypatch):
    _force_serial_window_generation(monkeypatch)
    df_positions = _make_positions_df()
    prepared_path = tmp_path / "prepared_full_missing_meta.h5ad"

    prepare_state_classification_dataset(
        df_positions=df_positions,
        features=["speed", "elongation"],
        binary_features_to_group=[],
        window_size=4,
        descriptive_features=["mean", "median", "std"],
        incomplete_window_policy="partial",
        scale_features=True,
        reuse_prepared_dataset=False,
        save_prepared_dataset=True,
        prepared_dataset_path=prepared_path,
        verbose=False,
    )

    stale = sc.read_h5ad(prepared_path)
    stale.uns["preprocessing"] = {
        key: value
        for key, value in stale.uns["preprocessing"].items()
        if key != "windowing"
    }
    stale.write(prepared_path, compression="gzip")

    adata_prepared = prepare_state_classification_dataset(
        df_positions=df_positions,
        features=["speed", "elongation"],
        binary_features_to_group=[],
        window_size=4,
        descriptive_features=["mean", "median"],
        incomplete_window_policy="partial",
        scale_features=True,
        reuse_prepared_dataset=True,
        save_prepared_dataset=True,
        prepared_dataset_path=prepared_path,
        verbose=False,
    )

    assert _std_cols(adata_prepared.var_names) == []
    assert adata_prepared.uns["preprocessing"]["windowing"]["descriptive_features"] == ["mean", "median"]


def test_run_state_clustering_rebuilds_model_cache_when_std_is_removed(tmp_path, monkeypatch):
    _force_serial_window_generation(monkeypatch)
    df_positions = _make_positions_df()
    output_dir = Path(tmp_path) / "state_cache_case"

    run_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=[],
        output_dir=output_dir,
        cell_type="tcell",
        window_size=4,
        min_spacing=1,
        n_neighbors=10,
        min_dist=0.1,
        resolution=3,
        descriptive_features=["mean", "median", "std"],
        clustering_method="kmeans",
        incomplete_window_policy="partial",
        random_state=7,
        reuse_prepared_dataset=True,
        plot_exemplar_videos=False,
        df_positions=df_positions,
        verbose=False,
    )

    model_adata = run_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=[],
        output_dir=output_dir,
        cell_type="tcell",
        window_size=4,
        min_spacing=1,
        n_neighbors=10,
        min_dist=0.1,
        resolution=3,
        descriptive_features=["mean", "median"],
        clustering_method="kmeans",
        incomplete_window_policy="partial",
        random_state=7,
        reuse_prepared_dataset=True,
        plot_exemplar_videos=False,
        df_positions=df_positions,
        verbose=False,
    )

    assert _std_cols(model_adata.var_names) == []

    state_paths = _resolve_state_paths(output_dir, "tcell")
    saved_model = sc.read_h5ad(state_paths.model_adata_path)
    assert _std_cols(saved_model.var_names) == []

    heatmap_csv = state_paths.state_clustering_outdir / "behavioral_clustering_diagnostics_heatmap_matrix.csv"
    heatmap_df = pd.read_csv(heatmap_csv, index_col=0)
    assert _std_cols(heatmap_df.columns) == []

    feature_distribution_pdf = state_paths.state_clustering_outdir / "behavioral_clustering_feature_distributions.pdf"
    assert feature_distribution_pdf.exists()
    assert model_adata.uns["clustering"]["feature_distribution_pdf"] == str(feature_distribution_pdf)


def test_run_state_clustering_supports_pca_and_no_pca_modes(tmp_path, monkeypatch):
    _force_serial_window_generation(monkeypatch)
    df_positions = _make_positions_df()

    model_with_pca = run_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=[],
        output_dir=Path(tmp_path) / "with_pca",
        cell_type="tcell",
        window_size=4,
        min_spacing=1,
        n_neighbors=10,
        min_dist=0.1,
        resolution=0.5,
        descriptive_features=["mean", "median", "std"],
        pca_var_selection=0.95,
        use_pca=True,
        clustering_method="leiden",
        incomplete_window_policy="partial",
        random_state=11,
        reuse_prepared_dataset=True,
        plot_exemplar_videos=False,
        df_positions=df_positions,
        verbose=False,
    )
    assert "X_pca" in model_with_pca.obsm
    assert "X_umap" in model_with_pca.obsm
    assert "neighbors" in model_with_pca.uns

    model_no_pca_leiden = run_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=[],
        output_dir=Path(tmp_path) / "no_pca_leiden",
        cell_type="tcell",
        window_size=4,
        min_spacing=1,
        n_neighbors=10,
        min_dist=0.1,
        resolution=0.5,
        descriptive_features=["mean", "median", "std"],
        pca_var_selection=0.95,
        use_pca=False,
        clustering_method="leiden",
        incomplete_window_policy="partial",
        random_state=11,
        reuse_prepared_dataset=True,
        plot_exemplar_videos=False,
        df_positions=df_positions,
        verbose=False,
    )
    assert "X_pca" not in model_no_pca_leiden.obsm
    assert "X_umap" in model_no_pca_leiden.obsm
    assert "neighbors" in model_no_pca_leiden.uns
    assert model_no_pca_leiden.uns["preprocessing"]["model_cache"]["neighbors"]["use_rep"] == "X"

    model_no_pca_kmeans = run_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=[],
        output_dir=Path(tmp_path) / "no_pca_kmeans",
        cell_type="tcell",
        window_size=4,
        min_spacing=1,
        n_neighbors=10,
        min_dist=0.1,
        resolution=3,
        descriptive_features=["mean", "median", "std"],
        pca_var_selection=0.95,
        use_pca=False,
        clustering_method="kmeans",
        incomplete_window_policy="partial",
        random_state=11,
        reuse_prepared_dataset=True,
        plot_exemplar_videos=False,
        df_positions=df_positions,
        verbose=False,
    )
    assert "X_pca" not in model_no_pca_kmeans.obsm
    assert "intrinsic_behavioral_cluster" in model_no_pca_kmeans.obs.columns


def test_run_state_clustering_rebuilds_model_cache_when_toggling_pca(tmp_path, monkeypatch):
    _force_serial_window_generation(monkeypatch)
    df_positions = _make_positions_df()
    output_dir = Path(tmp_path) / "toggle_pca_case"

    run_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=[],
        output_dir=output_dir,
        cell_type="tcell",
        window_size=4,
        min_spacing=1,
        n_neighbors=10,
        min_dist=0.1,
        resolution=3,
        descriptive_features=["mean", "median", "std"],
        pca_var_selection=0.95,
        use_pca=True,
        clustering_method="kmeans",
        incomplete_window_policy="partial",
        random_state=13,
        reuse_prepared_dataset=True,
        plot_exemplar_videos=False,
        df_positions=df_positions,
        verbose=False,
    )

    model_no_pca = run_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=[],
        output_dir=output_dir,
        cell_type="tcell",
        window_size=4,
        min_spacing=1,
        n_neighbors=10,
        min_dist=0.1,
        resolution=3,
        descriptive_features=["mean", "median", "std"],
        pca_var_selection=0.95,
        use_pca=False,
        clustering_method="kmeans",
        incomplete_window_policy="partial",
        random_state=13,
        reuse_prepared_dataset=True,
        plot_exemplar_videos=False,
        df_positions=df_positions,
        verbose=False,
    )
    assert "X_pca" not in model_no_pca.obsm
    assert model_no_pca.uns["preprocessing"]["model_cache"]["pca"]["enabled"] is False
    assert model_no_pca.uns["preprocessing"]["model_cache"]["neighbors"]["use_rep"] == "X"
    assert model_no_pca.uns["preprocessing"]["model_cache"]["umap"]["use_rep"] == "X"
    assert model_no_pca.uns["preprocessing"]["model_cache"]["clustering"]["use_rep"] == "X"
    assert model_no_pca.uns["clustering"]["use_pca"] is False
    assert model_no_pca.uns["clustering"]["use_rep"] == "X"

    model_with_pca_again = run_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=[],
        output_dir=output_dir,
        cell_type="tcell",
        window_size=4,
        min_spacing=1,
        n_neighbors=10,
        min_dist=0.1,
        resolution=3,
        descriptive_features=["mean", "median", "std"],
        pca_var_selection=0.95,
        use_pca=True,
        clustering_method="kmeans",
        incomplete_window_policy="partial",
        random_state=13,
        reuse_prepared_dataset=True,
        plot_exemplar_videos=False,
        df_positions=df_positions,
        verbose=False,
    )
    assert "X_pca" in model_with_pca_again.obsm
    assert model_with_pca_again.uns["preprocessing"]["model_cache"]["pca"]["enabled"] is True
    assert model_with_pca_again.uns["preprocessing"]["model_cache"]["neighbors"]["use_rep"] == "X_pca"
    assert model_with_pca_again.uns["clustering"]["use_pca"] is True
    assert model_with_pca_again.uns["clustering"]["use_rep"] == "X_pca"


def test_run_state_clustering_rebuilds_older_model_cache_missing_pca_metadata(tmp_path, monkeypatch):
    _force_serial_window_generation(monkeypatch)
    df_positions = _make_positions_df()
    output_dir = Path(tmp_path) / "older_model_cache_missing_pca_meta"

    run_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=[],
        output_dir=output_dir,
        cell_type="tcell",
        window_size=4,
        min_spacing=1,
        n_neighbors=10,
        min_dist=0.1,
        resolution=3,
        descriptive_features=["mean", "median", "std"],
        pca_var_selection=0.95,
        use_pca=True,
        clustering_method="kmeans",
        incomplete_window_policy="partial",
        random_state=17,
        reuse_prepared_dataset=True,
        plot_exemplar_videos=False,
        df_positions=df_positions,
        verbose=False,
    )

    state_paths = _resolve_state_paths(output_dir, "tcell")
    stale_model = sc.read_h5ad(state_paths.model_adata_path)
    del stale_model.uns["preprocessing"]["model_cache"]["pca"]["enabled"]
    stale_model.write(state_paths.model_adata_path, compression="gzip")

    rebuilt_model = run_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=[],
        output_dir=output_dir,
        cell_type="tcell",
        window_size=4,
        min_spacing=1,
        n_neighbors=10,
        min_dist=0.1,
        resolution=3,
        descriptive_features=["mean", "median", "std"],
        pca_var_selection=0.95,
        use_pca=True,
        clustering_method="kmeans",
        incomplete_window_policy="partial",
        random_state=17,
        reuse_prepared_dataset=True,
        plot_exemplar_videos=False,
        df_positions=df_positions,
        verbose=False,
    )

    assert "X_pca" in rebuilt_model.obsm
    assert rebuilt_model.uns["preprocessing"]["model_cache"]["pca"]["enabled"] is True


def test_run_hmm_state_clustering_fixed_k_outputs(tmp_path):
    df_positions = _make_positions_df()
    output_dir = Path(tmp_path) / "hmm_fixed_k"

    cluster_out = run_hmm_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=[],
        output_dir=output_dir,
        cell_type="tcell",
        n_states=3,
        random_state=19,
        df_positions=df_positions,
        return_details=True,
        verbose=False,
    )
    model_adata = cluster_out["model_adata"]
    artifact = cluster_out["hmm_deployment_artifact"]

    assert list(model_adata.var_names) == ["speed", "elongation"]
    assert state_classification.INTRINSIC_STATE_COL in model_adata.obs.columns
    assert state_classification.HMM_INTRINSIC_RAW_STATE_COL in model_adata.obs.columns
    assert "binary_group" in model_adata.obs.columns
    assert state_classification.FULL_STATE_COL in model_adata.obs.columns
    assert "full_behavioral_cluster" in model_adata.obs.columns
    assert model_adata.uns["clustering"]["clustering_method"] == "hmm"
    assert model_adata.uns["clustering"]["hmm"]["selected_k"] == 3
    assert artifact["artifact_kind"] == "hmm_state_deployment"
    assert "schema_version" not in artifact
    assert "schema_version" not in artifact["pipeline_metadata"]
    assert len(artifact["full_label_mapping"]) > 0

    state_paths = _resolve_state_paths(output_dir, "tcell")
    saved_model = sc.read_h5ad(state_paths.model_adata_path)
    assert saved_model.n_obs == model_adata.n_obs
    assert (state_paths.state_clustering_outdir / "behavioral_clustering_diagnostics.pdf").exists()
    assert (state_paths.state_clustering_outdir / "behavioral_clustering_feature_distributions.pdf").exists()
    assert (state_paths.state_clustering_outdir / "behavioral_clustering_hmm_state_counts.csv").exists()


def test_run_hmm_state_clustering_writes_grouped_binary_constraints_and_enforces(tmp_path):
    df_positions = _with_multi_binary_flags(_make_positions_df(n_tracks=8, track_len=10))
    output_dir = Path(tmp_path) / "hmm_grouped_binary_constraints"

    cluster_out = run_hmm_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=["contact_a", "contact_b", "contact_c"],
        output_dir=output_dir,
        cell_type="tcell",
        n_states=3,
        random_state=31,
        df_positions=df_positions,
        return_details=True,
        verbose=False,
    )
    model_adata = cluster_out["model_adata"]
    clustering_meta = model_adata.uns["clustering"]
    constraints = clustering_meta["binary_group_constraints"]
    forbidden = constraints["forbidden_binary_combinations"]

    assert isinstance(forbidden, dict)
    assert "2" in forbidden
    assert any(str(k) == "3" for k in forbidden.keys())

    state_paths = _resolve_state_paths(output_dir, "tcell")
    saved_model = sc.read_h5ad(state_paths.model_adata_path)
    saved_constraints = saved_model.uns["clustering"]["binary_group_constraints"]
    assert isinstance(saved_constraints["forbidden_binary_combinations"], dict)

    state_classification._rebuild_hmm_full_behavioral_state_from_intrinsic(
        adata=saved_model,
        binary_cols_to_merge=["contact_a", "contact_b", "contact_c"],
        intrinsic_col=state_classification.INTRINSIC_STATE_COL,
        binary_group_constraints=saved_constraints,
        enforce_binary_group_constraints=True,
    )
    assert state_classification.FULL_STATE_COL in saved_model.obs.columns


def test_classification_module_exposes_hmm_without_legacy_classifier_schema():
    assert "behav3d.analysis.clustering.state.classification" in sys.modules
    hmm_module = sys.modules["behav3d.analysis.clustering.state.classification"]
    assert not hasattr(hmm_module, "STATE_CLASSIFIER_PIPELINE_SCHEMA_VERSION")


def test_hmm_widget_refresh_enablement_uses_hmm_intrinsic_column(monkeypatch):
    widget_mod = _load_hmm_widget_module()
    widgets = pytest.importorskip("ipywidgets")
    monkeypatch.setattr(widget_mod._LegacyStateClassificationPanel, "_refresh_enablement", lambda self: None)

    panel = object.__new__(widget_mod.StateClassificationHMMPanel)
    panel.model_adata = SimpleNamespace(
        obs=pd.DataFrame(
            {
                state_classification.INTRINSIC_STATE_COL: ["1", "2"],
                state_classification.FULL_STATE_COL: ["no_contact_1", "no_contact_2"],
            }
        )
    )
    panel._current_cell_type = lambda: "tcell"
    panel._selected_hmm_feature_columns = lambda: ["speed"]
    panel.backproj_sample_dd = SimpleNamespace(value="sample_0")
    panel.apply_hmm_artifact_picker = SimpleNamespace(value="")
    panel._refresh_analysis_plots_status = lambda: None
    panel.btn_cluster = widgets.Button()
    panel.btn_rename_intrinsic = widgets.Button(disabled=True)
    panel.btn_combine_intrinsic = widgets.Button(disabled=True)
    panel.intrinsic_combine_name = widgets.Text(disabled=True)
    panel.btn_rename_full = widgets.Button(disabled=True)
    panel.btn_combine_full = widgets.Button(disabled=True)
    panel.full_combine_name = widgets.Text(disabled=True)
    panel.btn_open_intrinsic_backprojection = widgets.Button()
    panel.btn_apply_hmm_artifact = widgets.Button()
    panel.btn_create_state_composition_plots = widgets.Button()
    panel.btn_create_state_transition_plots = widgets.Button()

    widget_mod.StateClassificationHMMPanel._refresh_enablement(panel)

    assert panel.btn_rename_intrinsic.disabled is False
    assert panel.btn_combine_intrinsic.disabled is False
    assert panel.btn_rename_full.disabled is False


def test_hmm_widget_load_normalization_and_mapping_yaml(tmp_path):
    widgets = pytest.importorskip("ipywidgets")
    widget_mod = _load_hmm_widget_module()
    obs = pd.DataFrame(
        {
            "intrinsic_behavioral_cluster": ["1", "2", "2"],
            "full_behavioral_cluster": ["contact_1", "contact_2", "contact_2"],
            state_classification.BINARY_GROUP_COL: ["contact", "contact", "contact"],
        }
    )
    panel = object.__new__(widget_mod.StateClassificationHMMPanel)
    panel.model_adata = SimpleNamespace(
        obs=obs,
        uns={"clustering": {"clustering_method": "hmm"}},
    )
    normalized = widget_mod.StateClassificationHMMPanel._normalize_loaded_hmm_model_adata(panel)
    assert normalized["normalized"] is True
    assert state_classification.INTRINSIC_STATE_COL in panel.model_adata.obs.columns
    assert state_classification.FULL_STATE_COL in panel.model_adata.obs.columns
    assert state_classification.HMM_INTRINSIC_RAW_STATE_COL in panel.model_adata.obs.columns

    panel.output_dir = str(tmp_path)
    panel._current_cell_type = lambda: "tcell"
    panel._model_adata_path = lambda: tmp_path / "model.h5ad"
    panel._rename_mapping_yaml_path = lambda: tmp_path / "mapping.yml"
    panel._full_color_pickers = {}
    yaml_path = widget_mod.StateClassificationHMMPanel._write_cluster_name_mappings_yaml(panel)
    payload = yaml.safe_load(Path(yaml_path).read_text())
    current = payload["current_mappings"]

    assert current[
        f"{state_classification.HMM_INTRINSIC_RAW_STATE_COL}_to_{state_classification.INTRINSIC_STATE_COL}"
    ] == {"1": "1", "2": "2"}
    assert current[f"generated_{state_classification.FULL_STATE_COL}_to_{state_classification.FULL_STATE_COL}"] == {
        "contact_1": "contact_1",
        "contact_2": "contact_2",
    }
    assert set(payload["current_colors"][state_classification.FULL_STATE_COL]) == {"contact_1", "contact_2"}

    panel.rename_full_rows = widgets.VBox([])
    panel.rename_full_status = widgets.HTML()
    panel.btn_rename_full = widgets.Button()
    panel.btn_combine_full = widgets.Button()
    panel.full_combine_name = widgets.Text()
    widget_mod.StateClassificationHMMPanel._rebuild_full_rename_rows(panel)
    assert set(panel._full_color_pickers) == {"contact_1", "contact_2"}
    assert all(isinstance(picker, widgets.ColorPicker) for picker in panel._full_color_pickers.values())


def test_hmm_widget_rename_mapping_yaml_path_is_available_on_deployment_panel(tmp_path):
    widget_mod = _load_hmm_widget_module()
    panel = object.__new__(widget_mod.StateClassificationHMMDeploymentPanel)
    panel.output_dir = str(tmp_path)
    panel._current_cell_type = lambda: "organoid"

    yaml_path = widget_mod.StateClassificationHMMPanel._rename_mapping_yaml_path(panel)
    assert yaml_path.name == "hmm_cluster_name_mappings_organoid.yml"
    assert "hmm_behavioral_classification" in str(yaml_path)


def test_hmm_widget_mapping_dict_from_obs_handles_missing_and_multi_mapping():
    widget_mod = _load_hmm_widget_module()
    panel = object.__new__(widget_mod.StateClassificationHMMPanel)

    panel.model_adata = SimpleNamespace(obs=pd.DataFrame({"foo": ["1"], "bar": ["a"]}))
    assert widget_mod.StateClassificationHMMPanel._mapping_dict_from_obs(panel, "missing", "bar") == {}
    assert widget_mod.StateClassificationHMMPanel._mapping_dict_from_obs(panel, "foo", "missing") == {}

    panel.model_adata = SimpleNamespace(
        obs=pd.DataFrame(
            {
                "src": ["1", " 1 ", "2", "3", "", "4"],
                "dst": ["alpha", "alpha", "beta", "beta", "ignored", " beta "],
            }
        )
    )
    assert widget_mod.StateClassificationHMMPanel._mapping_dict_from_obs(panel, "src", "dst") == {
        "1": "alpha",
        "2": "beta",
        "3": "beta",
        "4": "beta",
    }

    panel.model_adata = SimpleNamespace(
        obs=pd.DataFrame(
            {
                "src": ["5", "5", "5", "6"],
                "dst": ["gamma_2", "gamma_1", "gamma_1", "delta"],
            }
        )
    )
    assert widget_mod.StateClassificationHMMPanel._mapping_dict_from_obs(panel, "src", "dst") == {
        "5": ["gamma_1", "gamma_2"],
        "6": "delta",
    }


def test_hmm_widget_load_normalization_ignores_non_hmm_model():
    widget_mod = _load_hmm_widget_module()
    panel = object.__new__(widget_mod.StateClassificationHMMPanel)
    panel.model_adata = SimpleNamespace(
        obs=pd.DataFrame(
            {
                "intrinsic_behavioral_cluster": ["1"],
                "full_behavioral_cluster": ["no_contact_1"],
            }
        ),
        uns={"clustering": {"clustering_method": "kmeans"}},
    )
    normalized = widget_mod.StateClassificationHMMPanel._normalize_loaded_hmm_model_adata(panel)
    assert normalized == {"normalized": False, "reason": "not_hmm"}
    assert state_classification.INTRINSIC_STATE_COL not in panel.model_adata.obs.columns
    assert state_classification.FULL_STATE_COL not in panel.model_adata.obs.columns


def test_state_backprojection_code_color_map_uses_state_labels():
    code_colors = state_backprojection._build_state_code_color_map(
        {"rest": 1, "move": 2},
        state_colors={"rest": "#ff0000", "move": "#00ff00"},
    )
    assert code_colors == {"1": "#ff0000", "2": "#00ff00"}


def test_run_hmm_state_clustering_auto_selection_and_smoothing(tmp_path):
    df_positions = _make_positions_df()

    unsmoothed = run_hmm_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=[],
        output_dir=Path(tmp_path) / "hmm_auto_unsmoothed",
        cell_type="tcell",
        n_states="auto",
        k_min=2,
        k_max=3,
        feature_smoothing_window=1,
        random_state=23,
        df_positions=df_positions,
        verbose=False,
    )

    smoothed = run_hmm_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=[],
        output_dir=Path(tmp_path) / "hmm_auto_smoothed",
        cell_type="tcell",
        n_states="auto",
        k_min=2,
        k_max=3,
        feature_smoothing_window=3,
        smoothing_min_periods=1,
        random_state=23,
        df_positions=df_positions,
        verbose=False,
    )

    assert smoothed.uns["preprocessing"]["timepoint_smoothing"]["feature_smoothing_window"] == 3
    assert smoothed.uns["clustering"]["model_selection_csv"] is not None
    assert Path(smoothed.uns["clustering"]["model_selection_csv"]).exists()
    assert not np.allclose(np.asarray(unsmoothed.X), np.asarray(smoothed.X))


def test_run_hmm_state_clustering_validates_inputs(tmp_path):
    df_positions = _make_positions_df()

    with pytest.raises(ValueError, match="n_states is required"):
        run_hmm_state_clustering(
            features=["speed", "elongation"],
            binary_features_to_group=[],
            output_dir=Path(tmp_path) / "hmm_missing_n_states",
            cell_type="tcell",
            df_positions=df_positions,
            verbose=False,
        )

    with pytest.raises(ValueError, match="Missing required columns"):
        run_hmm_state_clustering(
            features=["speed", "missing_feature"],
            binary_features_to_group=[],
            output_dir=Path(tmp_path) / "hmm_missing_feature",
            cell_type="tcell",
            n_states=2,
            df_positions=df_positions,
            verbose=False,
        )


def test_run_hmm_state_clustering_dispatches_sticky_helper(tmp_path, monkeypatch):
    df_positions = _make_positions_df()
    calls = []

    def _fake_plain(**kwargs):
        calls.append("plain")
        out = kwargs["df_features"].copy()
        out[kwargs["out_col_name"]] = ((np.arange(len(out)) % 2) + 1).astype(int)
        return out, SimpleNamespace(n_components=2), None

    def _fake_sticky(**kwargs):
        calls.append("sticky")
        out = kwargs["df_features"].copy()
        out[kwargs["out_col_name"]] = ((np.arange(len(out)) % 2) + 1).astype(int)
        return out, SimpleNamespace(n_components=2), None

    monkeypatch.setattr(state_classification, "run_hmm_state_classification", _fake_plain)
    monkeypatch.setattr(state_classification, "run_sticky_hmm_state_classification", _fake_sticky)

    plain = run_hmm_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=[],
        output_dir=Path(tmp_path) / "hmm_dispatch_plain",
        cell_type="tcell",
        n_states=2,
        sticky=False,
        random_state=29,
        df_positions=df_positions,
        verbose=False,
    )
    sticky = run_hmm_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=[],
        output_dir=Path(tmp_path) / "hmm_dispatch_sticky",
        cell_type="tcell",
        n_states=2,
        sticky=True,
        random_state=29,
        df_positions=df_positions,
        verbose=False,
    )

    assert calls == ["plain", "sticky"]
    assert plain.uns["clustering"]["hmm"]["sticky"] is False
    assert sticky.uns["clustering"]["hmm"]["sticky"] is True


def test_hmm_deployment_artifact_roundtrip_and_apply(tmp_path):
    df_positions = _with_binary_contact_flag(_make_positions_df(n_tracks=6, track_len=10))
    train_output_dir = Path(tmp_path) / "hmm_deployment_train"

    cluster_out = state_classification.run_hmm_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=["contact_flag"],
        output_dir=train_output_dir,
        cell_type="tcell",
        n_states=3,
        random_state=31,
        df_positions=df_positions,
        return_details=True,
        verbose=False,
    )
    model_adata = cluster_out["model_adata"]
    hmm_model = cluster_out["hmm_model"]
    state_paths = cluster_out["state_paths"]

    intrinsic_labels = sorted(
        model_adata.obs[state_classification.INTRINSIC_STATE_COL].astype(str).unique().tolist(),
        key=state_classification._mixed_label_sort_key,
    )
    mapping = {label: f"state_{label}" for label in intrinsic_labels}
    state_classification.relabel_cluster_ids(
        adata=model_adata,
        mapping=mapping,
        cluster_key=state_classification.INTRINSIC_STATE_COL,
        new_key=state_classification.INTRINSIC_STATE_COL,
        keep_unmapped=True,
        overwrite_original=True,
    )
    state_classification._rebuild_hmm_full_behavioral_state_from_intrinsic(
        adata=model_adata,
        binary_cols_to_merge=["contact_flag"],
        intrinsic_col=state_classification.INTRINSIC_STATE_COL,
    )
    full_labels = sorted(
        model_adata.obs[state_classification.FULL_STATE_COL].astype(str).unique().tolist(),
        key=state_classification._mixed_label_sort_key,
    )
    full_mapping = {label: f"curated_{idx}_{label}" for idx, label in enumerate(full_labels, start=1)}
    state_classification.relabel_cluster_ids(
        adata=model_adata,
        mapping=full_mapping,
        cluster_key=state_classification.FULL_STATE_COL,
        overwrite_original=True,
        keep_unmapped=True,
    )
    expected_state_colors = {
        new_label: color
        for new_label, color in zip(
            sorted(full_mapping.values(), key=state_classification._mixed_label_sort_key),
            ["#ff0000", "#00aa00", "#0000ff", "#ffaa00", "#00aaff", "#aa00ff"],
        )
    }
    model_adata.uns.setdefault("classification", {})["state_colors"] = {
        state_classification.FULL_STATE_COL: dict(expected_state_colors)
    }

    artifact_path = Path(tmp_path) / "saved_hmm_deployment.pkl"
    saved_artifact = state_classification.save_hmm_deployment_artifact(
        output_path=artifact_path,
        model_adata=model_adata,
        hmm_model=hmm_model,
        state_paths=state_paths,
        source_model_adata_path=state_paths.model_adata_path,
        verbose=False,
    )
    loaded_artifact = state_classification.load_hmm_deployment_artifact(artifact_path)

    assert saved_artifact["artifact_kind"] == "hmm_state_deployment"
    assert "schema_version" not in saved_artifact
    assert "schema_version" not in saved_artifact["pipeline_metadata"]
    assert loaded_artifact["artifact_kind"] == "hmm_state_deployment"
    assert "schema_version" not in loaded_artifact
    assert "schema_version" not in loaded_artifact["pipeline_metadata"]
    assert loaded_artifact["full_label_strategy"] == "binary_group_plus_intrinsic"
    assert loaded_artifact["observation_feature_cols"] == ["speed", "elongation"]
    assert loaded_artifact["binary_cols_to_merge"] == ["contact_flag"]
    assert all(str(v).startswith("state_") for v in loaded_artifact["intrinsic_label_mapping"].values())
    assert set(loaded_artifact["full_label_mapping"].values()) == set(full_mapping.values())
    assert loaded_artifact["state_colors"][state_classification.FULL_STATE_COL] == expected_state_colors
    expected_raw_full_keys = set(
        (
            model_adata.obs["binary_group"].astype(str)
            + "_"
            + model_adata.obs[state_classification.HMM_INTRINSIC_RAW_STATE_COL].astype(str)
        ).tolist()
    )
    assert set(loaded_artifact["full_label_mapping"].keys()) == expected_raw_full_keys

    apply_output_dir = Path(tmp_path) / "hmm_deployment_apply"
    feature_outdir = apply_output_dir / "analysis" / "tcell" / "track_features"
    feature_outdir.mkdir(parents=True, exist_ok=True)
    df_positions.to_csv(feature_outdir / "BEHAV3D_tcell_combined_track_features_filtered.csv", index=False)

    adata_applied = state_classification.apply_hmm_deployment_artifact_to_full_dataset(
        output_dir=apply_output_dir,
        cell_type="tcell",
        hmm_deployment_artifact=artifact_path,
        verbose=False,
    )

    assert state_classification.INTRINSIC_STATE_COL in adata_applied.obs.columns
    assert state_classification.FULL_STATE_COL in adata_applied.obs.columns
    assert "full_behavioral_cluster" in adata_applied.obs.columns
    assert "intrinsic_behavioral_cluster_confidence" in adata_applied.obs.columns
    assert "full_behavioral_cluster_confidence" in adata_applied.obs.columns
    assert all(
        str(v).startswith("state_")
        for v in adata_applied.obs[state_classification.INTRINSIC_STATE_COL].astype(str).unique().tolist()
    )
    raw_expected_full = (
        adata_applied.obs["binary_group"].astype(str)
        + "_" + adata_applied.obs[state_classification.HMM_INTRINSIC_RAW_STATE_COL].astype(str)
    )
    expected_full = raw_expected_full.map(loaded_artifact["full_label_mapping"])
    assert (adata_applied.obs[state_classification.FULL_STATE_COL].astype(str) == expected_full.astype(str)).all()
    assert (
        adata_applied.uns["classification"]["state_colors"][state_classification.FULL_STATE_COL]
        == expected_state_colors
    )
    assert np.allclose(
        adata_applied.obs["full_behavioral_cluster_confidence"].to_numpy(dtype=float),
        adata_applied.obs["intrinsic_behavioral_cluster_confidence"].to_numpy(dtype=float),
    )


def test_hmm_deployment_artifact_updates_after_intrinsic_rename(tmp_path):
    df_positions = _with_binary_contact_flag(_make_positions_df(n_tracks=6, track_len=8))
    output_dir = Path(tmp_path) / "hmm_deployment_rename"

    cluster_out = state_classification.run_hmm_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=["contact_flag"],
        output_dir=output_dir,
        cell_type="tcell",
        n_states=3,
        random_state=37,
        df_positions=df_positions,
        return_details=True,
        verbose=False,
    )
    model_adata = cluster_out["model_adata"]
    artifact_path = Path(tmp_path) / "hmm_deployment_rename.pkl"

    state_classification.save_hmm_deployment_artifact(
        output_path=artifact_path,
        model_adata=model_adata,
        hmm_model=cluster_out["hmm_model"],
        state_paths=cluster_out["state_paths"],
        verbose=False,
    )

    first_mapping = {}
    for idx, label in enumerate(
        sorted(
            model_adata.obs[state_classification.INTRINSIC_STATE_COL].astype(str).unique().tolist(),
            key=state_classification._mixed_label_sort_key,
        )
    ):
        first_mapping[str(label)] = "merged_state" if idx < 2 else f"kept_{label}"
    state_classification.relabel_cluster_ids(
        adata=model_adata,
        mapping=first_mapping,
        cluster_key=state_classification.INTRINSIC_STATE_COL,
        new_key=state_classification.INTRINSIC_STATE_COL,
        keep_unmapped=True,
        overwrite_original=True,
    )
    state_classification._rebuild_hmm_full_behavioral_state_from_intrinsic(
        adata=model_adata,
        binary_cols_to_merge=["contact_flag"],
        intrinsic_col=state_classification.INTRINSIC_STATE_COL,
    )
    second_mapping = {}
    for label in sorted(
        model_adata.obs[state_classification.INTRINSIC_STATE_COL].astype(str).unique().tolist(),
        key=state_classification._mixed_label_sort_key,
    ):
        second_mapping[str(label)] = "final_merged" if str(label) == "merged_state" else f"final_{label}"
    state_classification.relabel_cluster_ids(
        adata=model_adata,
        mapping=second_mapping,
        cluster_key=state_classification.INTRINSIC_STATE_COL,
        new_key=state_classification.INTRINSIC_STATE_COL,
        keep_unmapped=True,
        overwrite_original=True,
    )
    state_classification._rebuild_hmm_full_behavioral_state_from_intrinsic(
        adata=model_adata,
        binary_cols_to_merge=["contact_flag"],
        intrinsic_col=state_classification.INTRINSIC_STATE_COL,
    )

    state_classification.save_hmm_deployment_artifact(
        output_path=artifact_path,
        model_adata=model_adata,
        hmm_model=cluster_out["hmm_model"],
        state_paths=cluster_out["state_paths"],
        verbose=False,
    )
    updated_artifact = state_classification.load_hmm_deployment_artifact(artifact_path)

    assert "final_merged" in set(updated_artifact["intrinsic_label_mapping"].values())
    assert all(
        label in {"final_merged"} or str(label).startswith("final_kept_")
        for label in updated_artifact["intrinsic_label_mapping"].values()
    )
    assert isinstance(updated_artifact["full_label_mapping"], dict)
    assert len(updated_artifact["full_label_mapping"]) > 0
    expected_raw_full_keys = set(
        (
            model_adata.obs["binary_group"].astype(str)
            + "_"
            + model_adata.obs[state_classification.HMM_INTRINSIC_RAW_STATE_COL].astype(str)
        ).tolist()
    )
    assert set(updated_artifact["full_label_mapping"].keys()) == expected_raw_full_keys


def test_hmm_deployment_artifact_does_not_require_schema_version(tmp_path):
    df_positions = _with_binary_contact_flag(_make_positions_df(n_tracks=4, track_len=6))
    output_dir = Path(tmp_path) / "hmm_no_schema"

    cluster_out = state_classification.run_hmm_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=["contact_flag"],
        output_dir=output_dir,
        cell_type="tcell",
        n_states=2,
        random_state=53,
        df_positions=df_positions,
        return_details=True,
        verbose=False,
    )
    artifact = state_classification._build_hmm_deployment_artifact(
        model_adata=cluster_out["model_adata"],
        hmm_model=cluster_out["hmm_model"],
        state_paths=cluster_out["state_paths"],
    )
    assert "schema_version" not in artifact
    assert "schema_version" not in artifact["pipeline_metadata"]
    assert "state_colors" in artifact

    old_artifact = dict(artifact)
    old_artifact.pop("state_colors", None)
    state_classification._validate_hmm_deployment_artifact(old_artifact)


def test_hmm_deployment_apply_fails_when_canonical_full_mapping_is_missing(tmp_path):
    df_positions = _with_binary_contact_flag(_make_positions_df(n_tracks=4, track_len=6))
    output_dir = Path(tmp_path) / "hmm_missing_full_mapping"

    cluster_out = state_classification.run_hmm_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=["contact_flag"],
        output_dir=output_dir,
        cell_type="tcell",
        n_states=2,
        random_state=59,
        df_positions=df_positions,
        return_details=True,
        verbose=False,
    )
    artifact = state_classification._build_hmm_deployment_artifact(
        model_adata=cluster_out["model_adata"],
        hmm_model=cluster_out["hmm_model"],
        state_paths=cluster_out["state_paths"],
    )
    missing_key = next(iter(artifact["full_label_mapping"].keys()))
    artifact["full_label_mapping"].pop(missing_key)

    feature_outdir = output_dir / "analysis" / "tcell" / "track_features"
    feature_outdir.mkdir(parents=True, exist_ok=True)
    df_positions.to_csv(feature_outdir / "BEHAV3D_tcell_combined_track_features_filtered.csv", index=False)

    with pytest.raises(ValueError, match="not applicable to this data"):
        state_classification.apply_hmm_deployment_artifact_to_full_dataset(
            output_dir=output_dir,
            cell_type="tcell",
            hmm_deployment_artifact=artifact,
            verbose=False,
        )


def test_run_hmm_state_clustering_start_offset_backfills_initial_rows(tmp_path):
    df_positions = _with_binary_contact_flag(_make_positions_df(n_tracks=4, track_len=6))
    output_dir = Path(tmp_path) / "hmm_start_offset_backfill"

    model_adata = state_classification.run_hmm_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=["contact_flag"],
        output_dir=output_dir,
        cell_type="tcell",
        n_states=2,
        start_offset=1,
        start_offset_fill_mode="backfill",
        random_state=43,
        df_positions=df_positions,
        verbose=False,
    )

    assert model_adata.n_obs == len(df_positions)
    assert model_adata.uns["preprocessing"]["start_offset"] == 1
    assert model_adata.uns["preprocessing"]["start_offset_fill_mode"] == "backfill"
    for _, group in model_adata.obs.sort_values(["sample_name", "TrackID", "position_t"]).groupby(
        ["sample_name", "TrackID"],
        sort=False,
        observed=False,
    ):
        labels = pd.Series(
            group[state_classification.INTRINSIC_STATE_COL],
            dtype="string",
        ).reset_index(drop=True)
        assert labels.notna().all()
        assert labels.iloc[0] == labels.iloc[1]


def test_hmm_deployment_apply_start_offset_leave_unassigned(tmp_path):
    df_positions = _with_binary_contact_flag(_make_positions_df_variable_lengths([1, 5, 5, 5]))
    train_output_dir = Path(tmp_path) / "hmm_start_offset_leave_unassigned_train"

    cluster_out = state_classification.run_hmm_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=["contact_flag"],
        output_dir=train_output_dir,
        cell_type="tcell",
        n_states=2,
        start_offset=1,
        start_offset_fill_mode="leave_unassigned",
        random_state=47,
        df_positions=df_positions,
        return_details=True,
        verbose=False,
    )
    artifact_path = Path(tmp_path) / "hmm_start_offset_leave_unassigned.pkl"
    state_classification.save_hmm_deployment_artifact(
        output_path=artifact_path,
        model_adata=cluster_out["model_adata"],
        hmm_model=cluster_out["hmm_model"],
        state_paths=cluster_out["state_paths"],
        verbose=False,
    )

    apply_output_dir = Path(tmp_path) / "hmm_start_offset_leave_unassigned_apply"
    feature_outdir = apply_output_dir / "analysis" / "tcell" / "track_features"
    feature_outdir.mkdir(parents=True, exist_ok=True)
    df_positions.to_csv(feature_outdir / "BEHAV3D_tcell_combined_track_features_filtered.csv", index=False)

    adata_applied = state_classification.apply_hmm_deployment_artifact_to_full_dataset(
        output_dir=apply_output_dir,
        cell_type="tcell",
        hmm_deployment_artifact=artifact_path,
        verbose=False,
    )

    obs = adata_applied.obs.sort_values(["sample_name", "TrackID", "position_t"]).copy()
    assert adata_applied.uns["preprocessing"]["start_offset"] == 1
    assert adata_applied.uns["preprocessing"]["start_offset_fill_mode"] == "leave_unassigned"

    first_rows = obs.groupby(["sample_name", "TrackID"], sort=False, observed=False).head(1)
    multi_frame_first_rows = first_rows[first_rows["TrackID"].astype(int) != 0]
    assert pd.Series(multi_frame_first_rows[state_classification.INTRINSIC_STATE_COL], dtype="string").isna().all()
    assert pd.Series(multi_frame_first_rows[state_classification.FULL_STATE_COL], dtype="string").isna().all()
    assert pd.Series(multi_frame_first_rows["full_behavioral_cluster"], dtype="string").isna().all()
    assert pd.to_numeric(multi_frame_first_rows["intrinsic_behavioral_cluster_confidence"], errors="coerce").isna().all()

    scored_rows = obs.groupby(["sample_name", "TrackID"], sort=False, observed=False).nth(1)
    scored_rows = scored_rows[scored_rows.index.get_level_values("TrackID").astype(int) != 0]
    assert pd.Series(scored_rows[state_classification.INTRINSIC_STATE_COL], dtype="string").notna().all()
    assert pd.Series(scored_rows[state_classification.FULL_STATE_COL], dtype="string").notna().all()

    short_track = obs[obs["TrackID"].astype(int) == 0]
    assert len(short_track) == 1
    assert pd.Series(short_track[state_classification.INTRINSIC_STATE_COL], dtype="string").isna().all()
    assert pd.Series(short_track[state_classification.FULL_STATE_COL], dtype="string").isna().all()
    assert pd.Series(short_track["full_behavioral_cluster"], dtype="string").isna().all()


def test_hmm_deployment_apply_uses_track_lengths(tmp_path, monkeypatch):
    df_positions = _with_binary_contact_flag(_make_positions_df(n_tracks=4, track_len=6))
    output_dir = Path(tmp_path) / "hmm_deployment_lengths"

    cluster_out = state_classification.run_hmm_state_clustering(
        features=["speed", "elongation"],
        binary_features_to_group=["contact_flag"],
        output_dir=output_dir,
        cell_type="tcell",
        n_states=2,
        random_state=41,
        df_positions=df_positions,
        return_details=True,
        verbose=False,
    )
    artifact = state_classification._build_hmm_deployment_artifact(
        model_adata=cluster_out["model_adata"],
        hmm_model=cluster_out["hmm_model"],
        state_paths=cluster_out["state_paths"],
    )

    feature_outdir = output_dir / "analysis" / "tcell" / "track_features"
    feature_outdir.mkdir(parents=True, exist_ok=True)
    df_positions.to_csv(feature_outdir / "BEHAV3D_tcell_combined_track_features_filtered.csv", index=False)

    captured = {"predict": None, "predict_proba": None}
    model = artifact["model"]
    orig_predict = model.predict
    orig_predict_proba = model.predict_proba

    def _predict_spy(X, lengths=None):
        captured["predict"] = None if lengths is None else tuple(int(v) for v in lengths)
        return orig_predict(X, lengths=lengths)

    def _predict_proba_spy(X, lengths=None):
        captured["predict_proba"] = None if lengths is None else tuple(int(v) for v in lengths)
        return orig_predict_proba(X, lengths=lengths)

    monkeypatch.setattr(model, "predict", _predict_spy)
    monkeypatch.setattr(model, "predict_proba", _predict_proba_spy)

    adata_applied = state_classification.apply_hmm_deployment_artifact_to_full_dataset(
        output_dir=output_dir,
        cell_type="tcell",
        hmm_deployment_artifact=artifact,
        verbose=False,
    )

    expected_lengths = tuple(
        int(v)
        for v in (
            adata_applied.obs[["sample_name", "TrackID"]]
            .astype({"sample_name": "string"})
            .groupby(["sample_name", "TrackID"], sort=False, observed=False)
            .size()
            .tolist()
        )
    )
    assert captured["predict"] == expected_lengths
    assert captured["predict_proba"] == expected_lengths
