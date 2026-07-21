import importlib.util
import os
from pathlib import Path


os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/behav3d-mpl-cache")

import matplotlib
import pandas as pd

matplotlib.use("Agg", force=True)

from matplotlib import pyplot as plt

MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "behav3d/analysis/behavior/general/visualization/plots/proportion_bars.py"
)
SPEC = importlib.util.spec_from_file_location("proportion_bars", MODULE_PATH)
proportion_bars = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(proportion_bars)


def test_stacked_proportion_barh_grid_accepts_series_values():
    fig = proportion_bars.plot_page_stacked_proportion_barh_grid(
        {
            "sample_a": pd.Series({"rest": 0.25, "move": 0.75}),
            "sample_b": pd.Series({"rest": 0.5, "move": 0.5}),
        },
        row_order=["sample_a", "sample_b"],
        class_order=["rest", "move"],
        colors={"rest": "#4477AA", "move": "#EE6677"},
        title="Series-backed proportions",
    )

    assert fig is not None
    plt.close(fig)


def _build_condition_diff_inputs():
    samples = ["s1", "s2", "s3", "s4", "s5", "s6"]
    per_sample_class_props = pd.DataFrame(
        {
            "rest": [0.2, 0.3, 0.5, 0.6, 0.7, 0.8],
            "move": [0.8, 0.7, 0.5, 0.4, 0.3, 0.2],
        },
        index=samples,
    )
    sample_metadata = pd.DataFrame(
        {
            "condition": ["None", "None", "M21", "M21", "M23", "M23"],
            "batch": ["b1", "b2", "b1", "b2", "b1", "b2"],
        },
        index=samples,
    )
    return per_sample_class_props, sample_metadata


def test_compute_condition_diff_stats_pairwise_orders_pairs_by_first_appearance():
    per_sample_class_props, sample_metadata = _build_condition_diff_inputs()

    out = proportion_bars.compute_condition_diff_stats_pairwise(
        per_sample_class_props,
        sample_metadata,
        class_order=["rest", "move"],
        condition_col="condition",
    )

    assert list(out.keys()) == ["(all)"]
    pairs = list(out["(all)"].keys())
    assert pairs == [("None", "M21"), ("None", "M23"), ("M21", "M23")]
    for diff_df in out["(all)"].values():
        assert list(diff_df["class"]) == ["rest", "move"]


def test_compute_condition_diff_stats_pairwise_splits_by_group_cols():
    per_sample_class_props, sample_metadata = _build_condition_diff_inputs()

    out = proportion_bars.compute_condition_diff_stats_pairwise(
        per_sample_class_props,
        sample_metadata,
        class_order=["rest", "move"],
        condition_col="condition",
        group_cols=["batch"],
    )

    assert set(out.keys()) == {"b1", "b2"}
    for pairs in out.values():
        assert list(pairs.keys()) == [("None", "M21"), ("None", "M23"), ("M21", "M23")]


def test_plot_condition_diff_grid_single_column_and_multi_group(tmp_path):
    per_sample_class_props, sample_metadata = _build_condition_diff_inputs()
    colors = {"rest": "#4477AA", "move": "#EE6677"}

    diff_stats = proportion_bars.compute_condition_diff_stats_pairwise(
        per_sample_class_props,
        sample_metadata,
        class_order=["rest", "move"],
        condition_col="condition",
    )
    out_pdf = tmp_path / "single_row.pdf"
    result = proportion_bars.plot_condition_diff_grid(
        diff_stats,
        class_order=["rest", "move"],
        colors=colors,
        title="Single-row comparison",
        out_pdf=out_pdf,
        out_csv=out_pdf.with_suffix(".csv"),
    )
    assert Path(result["pdf_path"]).exists()
    csv_out = pd.read_csv(result["csv_path"], keep_default_na=False)
    assert set(csv_out["group"].unique()) == {"(all)"}
    assert set(zip(csv_out["level_a"], csv_out["level_b"])) == {
        ("None", "M21"), ("None", "M23"), ("M21", "M23"),
    }

    diff_stats_grouped = proportion_bars.compute_condition_diff_stats_pairwise(
        per_sample_class_props,
        sample_metadata,
        class_order=["rest", "move"],
        condition_col="condition",
        group_cols=["batch"],
    )
    out_pdf_grouped = tmp_path / "grouped.pdf"
    result_grouped = proportion_bars.plot_condition_diff_grid(
        diff_stats_grouped,
        class_order=["rest", "move"],
        colors=colors,
        title="Grouped comparison",
        out_pdf=out_pdf_grouped,
        out_csv=out_pdf_grouped.with_suffix(".csv"),
    )
    assert Path(result_grouped["pdf_path"]).exists()
    csv_out_grouped = pd.read_csv(result_grouped["csv_path"])
    assert set(csv_out_grouped["group"].unique()) == {"b1", "b2"}
