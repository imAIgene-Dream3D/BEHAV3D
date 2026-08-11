import pandas as pd
import pytest

from behav3d.core.combine_experiments import combine_behav3d_experiments

_METADATA_COLUMNS = [
    "sample_name", "exp_nr", "well", "raw_image_path",
    "pixel_distance_xy", "pixel_distance_z", "distance_unit",
    "time_interval", "time_unit", "or_organoid_line_condition",
]


def _make_experiment(tmp_path, name, sample_names, cell_type=None, feature_columns=None):
    exp_dir = tmp_path / name
    exp_dir.mkdir()

    rows = []
    for sample_name in sample_names:
        rows.append({
            "sample_name": sample_name,
            "exp_nr": 1,
            "well": "A1",
            "raw_image_path": f"/data/{name}/{sample_name}.czi",
            "pixel_distance_xy": 0.5,
            "pixel_distance_z": 2.0,
            "distance_unit": "um",
            "time_interval": 1.0,
            "time_unit": "m",
            "or_organoid_line_condition": "10T_unstim",
        })
    metadata = pd.DataFrame(rows, columns=_METADATA_COLUMNS)
    metadata.to_csv(exp_dir / "metadata.csv", index=False)

    if cell_type is not None:
        feature_columns = feature_columns or ["speed", "volume"]
        feature_dir = exp_dir / "analysis" / cell_type / "track_features"
        feature_dir.mkdir(parents=True)
        feature_rows = []
        for sample_name in sample_names:
            row = {"sample_name": sample_name, "TrackID": 1, "position_t": 0}
            for col in feature_columns:
                row[col] = 1.0
            feature_rows.append(row)
        feature_df = pd.DataFrame(feature_rows)
        feature_df.to_csv(
            feature_dir / f"BEHAV3D_{cell_type}_combined_track_features.csv", index=False
        )

    return exp_dir


def test_happy_path_merges_metadata_and_features(tmp_path):
    exp1 = _make_experiment(tmp_path, "exp1", ["S1", "S2"], cell_type="tcell")
    exp2 = _make_experiment(tmp_path, "exp2", ["S3", "S4"], cell_type="tcell")
    output_dir = tmp_path / "combined"

    result = combine_behav3d_experiments([exp1, exp2], output_dir)

    assert sorted(result["metadata"]["sample_name"]) == ["S1", "S2", "S3", "S4"]
    assert (output_dir / "metadata.csv").exists()
    assert result["warnings"] == []
    assert "tcell" in result["combined_features"]

    merged_features = pd.read_csv(result["combined_features"]["tcell"])
    assert sorted(merged_features["sample_name"]) == ["S1", "S2", "S3", "S4"]


def test_duplicate_sample_names_raise(tmp_path):
    exp1 = _make_experiment(tmp_path, "exp1", ["S1", "S2"])
    exp2 = _make_experiment(tmp_path, "exp2", ["S2", "S3"])
    output_dir = tmp_path / "combined"

    with pytest.raises(ValueError, match="S2"):
        combine_behav3d_experiments([exp1, exp2], output_dir)

    assert not output_dir.exists()


def test_mismatched_feature_columns_warns_but_merges(tmp_path):
    exp1 = _make_experiment(tmp_path, "exp1", ["S1"], cell_type="tcell", feature_columns=["speed", "volume"])
    exp2 = _make_experiment(tmp_path, "exp2", ["S2"], cell_type="tcell", feature_columns=["speed", "contact"])
    output_dir = tmp_path / "combined"

    result = combine_behav3d_experiments([exp1, exp2], output_dir)

    assert any("do not have identical columns" in w for w in result["warnings"])
    merged_features = pd.read_csv(result["combined_features"]["tcell"])
    assert set(merged_features.columns) >= {"sample_name", "speed", "volume", "contact"}
    assert merged_features.loc[merged_features["sample_name"] == "S1", "contact"].isna().all()
    assert merged_features.loc[merged_features["sample_name"] == "S2", "volume"].isna().all()


def test_cell_type_missing_from_one_source_is_skipped(tmp_path):
    exp1 = _make_experiment(tmp_path, "exp1", ["S1"], cell_type="tcell")
    exp2 = _make_experiment(tmp_path, "exp2", ["S2"])  # no feature extraction run
    output_dir = tmp_path / "combined"

    result = combine_behav3d_experiments([exp1, exp2], output_dir)

    assert "tcell" not in result["combined_features"]
    assert result["skipped_cell_types"]["tcell"] == ["exp2"]
    assert any("tcell" in w for w in result["warnings"])
    # metadata itself still merges fine
    assert sorted(result["metadata"]["sample_name"]) == ["S1", "S2"]
