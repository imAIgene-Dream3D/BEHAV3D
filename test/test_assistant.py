"""
Tests for the BEHAV3D napari co-pilot assistant (pure-logic layer).

Runs under pytest *or* standalone:  python test/test_assistant.py
Requires the `behav3d` conda env (numpy, the behav3d package). No GPU / Modal /
Qt needed — only the dependency-free logic is exercised here.
"""
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "chatbot"))

from behav3d.napari._assistant_schema import (
    flatten_config_to_cards, cards_for_step, dump_cards_json,
)
from behav3d.napari._assistant_actions import (
    get_by_dotted, set_by_dotted, validate_value, build_actions, apply_set_ui_value,
)
from behav3d.napari._assistant_context import (
    summarize_metadata, _diff_from_defaults, validate_metadata_records,
    _metadata_builder_state, _experiment_reference_context, build_context,
)
from behav3d.napari._assistant_controls import (
    CONTROL_CONTRACT_VERSION, active_cell_type, control_registry,
)
from behav3d.napari._assistant_recommendations import (
    calculate_edt_recommendations, format_edt_recommendations,
)
from behav3d.napari._assistant import researcher_facing_text
from behav3d.analysis.track_counts import (
    calculate_track_count_table, format_track_count_summary,
    generate_track_count_summary,
)
from behav3d.analysis.behavior.track.utils import _filter_tracks_for_dtaidistance
from behav3d.analysis.behavior.track.visualization.backprojection import (
    _build_track_cluster_lookup,
    _resolve_track_cluster_label,
)
from behav3d.napari._results_catalog import scan_outputs


# --------------------------------------------------------------------------
def _make_dtaidistance_adata(track_lengths):
    rows = []
    for track_id, track_len in enumerate(track_lengths):
        for t in range(int(track_len)):
            rows.append(
                {
                    "sample_name": "sample_0",
                    "TrackID": int(track_id),
                    "position_t": int(t),
                    "full_behavioral_cluster": "state_a" if t % 2 == 0 else "state_b",
                }
            )
    obs = pd.DataFrame(rows)
    return ad.AnnData(X=np.zeros((len(obs), 1)), obs=obs)


def test_dtaidistance_split_long_tracks_first_mode_preserves_original_trackid():
    adata = _make_dtaidistance_adata([350])
    out = _filter_tracks_for_dtaidistance(
        adata,
        trajectory_size=100,
        min_length=100,
        trim_mode="first",
        split_long_tracks=True,
    )

    obs = out.obs.sort_values(["trajectory_window_id", "position_t"])
    assert obs["TrackID"].nunique() == 1
    assert obs["TrackID"].iloc[0] == 0
    assert obs.groupby("trajectory_window_id", observed=True).size().tolist() == [100, 100, 100]
    assert obs["position_t"].min() == 0
    assert obs["position_t"].max() == 299


def test_dtaidistance_split_long_tracks_last_mode_drops_leading_leftover():
    adata = _make_dtaidistance_adata([350, 100, 80])
    out = _filter_tracks_for_dtaidistance(
        adata,
        trajectory_size=100,
        min_length=100,
        trim_mode="last",
        split_long_tracks=True,
    )

    obs = out.obs.sort_values(["TrackID", "trajectory_window_id", "position_t"])
    long_track = obs[obs["TrackID"] == 0]
    exact_track = obs[obs["TrackID"] == 1]
    short_track = obs[obs["TrackID"] == 2]
    assert len(long_track) == 300
    assert long_track["position_t"].min() == 50
    assert long_track["position_t"].max() == 349
    assert long_track.groupby("trajectory_window_id", observed=True).size().tolist() == [100, 100, 100]
    assert len(exact_track) == 100
    assert exact_track["trajectory_window_id"].nunique() == 1
    assert len(short_track) == 0


def test_dtaidistance_split_track_cluster_lookup_uses_clicked_window():
    obs = pd.DataFrame(
        {
            "sample_name": ["sample_0", "sample_0"],
            "TrackID": [7, 7],
            "trajectory_window_id": [0, 1],
            "position_t_min": [0, 100],
            "position_t_max": [99, 199],
            "ClusterID": ["early", "late"],
        }
    )
    adata_tracks = ad.AnnData(X=np.zeros((2, 1)), obs=obs)
    lookup = _build_track_cluster_lookup(adata_tracks)

    assert _resolve_track_cluster_label(lookup, "sample_0", "7", cursor_time=25) == "early"
    assert _resolve_track_cluster_label(lookup, "sample_0", "7", cursor_time=125) == "late"
    assert _resolve_track_cluster_label(lookup, "sample_0", "7", cursor_time=None) == "multiple"


def test_schema_flatten_synthetic():
    cfg = {"tracking": {"immune": {"method": "trackpy",
                                   "lap": {"track_cost_px": 60}}},
           "cellpose": {"channel_labels": {}}, "seed": 42}
    cards = flatten_config_to_cards(config=cfg, calculated_features={})
    keys = {c["key"] for c in cards}
    assert "tracking.immune.method" in keys           # recurses param groups
    assert "tracking.immune.lap.track_cost_px" in keys
    assert "cellpose.channel_labels" in keys          # int-keyed map = leaf
    assert "seed" in keys
    method = next(c for c in cards if c["key"] == "tracking.immune.method")
    assert method["choices"] and "btrack" in method["choices"]
    assert method["category"] == "immune" and method["step"] == "tracking"


def test_schema_real_config_coverage():
    cards = flatten_config_to_cards()
    assert len(cards) > 150
    steps = {c["step"] for c in cards}
    assert {"tracking", "segmentation", "feature_extraction", "filtering"} <= steps
    assert all("default" in c and "description" in c for c in cards)


def test_dotted_get_set():
    d = {}
    set_by_dotted(d, "a.b.c", 7)
    assert get_by_dotted(d, "a.b.c") == 7
    assert get_by_dotted(d, "a.x", "DEF") == "DEF"


def test_validation_types_choices_bounds():
    cards = {c["key"]: c for c in flatten_config_to_cards()}
    c = cards["tracking.immune.trackpy.search_range_px"]
    ok, v, _ = validate_value(c, "50"); assert ok and v == 50
    m = cards["tracking.immune.method"]
    assert validate_value(m, "btrack")[0]
    assert not validate_value(m, "nope")[0]
    b = cards["tracking.immune.btrack.use_optimize"]
    assert validate_value(b, "true") == (True, True, "")
    # synthetic min/max
    c2 = {"key": "x", "type": "int", "min": 1, "max": 10}
    assert not validate_value(c2, 0)[0]
    assert not validate_value(c2, 99)[0]
    assert validate_value(c2, 5)[0]


def test_build_actions_preview_and_rejection():
    cards = flatten_config_to_cards()
    params = {"tracking": {"immune": {"method": "trackpy"}}}
    acts = build_actions([
        {"name": "set_parameter", "arguments": {"key": "tracking.immune.method", "value": "btrack"}},
        {"name": "set_parameter", "arguments": {"key": "bogus.key", "value": 1}},
        {"name": "navigate_to_step", "arguments": {"step": "tracking"}},
    ], cards, params)
    assert acts[0].ok and "immune cells: method" in acts[0].preview
    assert not acts[1].ok                     # unknown key rejected
    assert acts[2].kind == "navigate_to_step"


def test_build_actions_rejects_noops_and_empty():
    cards = flatten_config_to_cards()
    # proposing the value it already has (e.g. '' -> '') is a no-op, no Fill button
    acts = build_actions(
        [{"name": "set_parameter", "arguments": {"key": "paths.metadata_csv", "value": ""}}],
        cards, {"paths": {"metadata_csv": ""}})
    assert not acts[0].ok and "no change" in acts[0].message.lower()
    assert acts[0].data.get("no_op") is True
    # empty string for an unset value is rejected too
    acts2 = build_actions(
        [{"name": "set_parameter", "arguments": {"key": "paths.output_dir", "value": "   "}}],
        cards, {})
    assert not acts2[0].ok
    # a genuine change is still actionable
    acts3 = build_actions(
        [{"name": "set_parameter", "arguments": {"key": "paths.metadata_csv", "value": "/data/m.csv"}}],
        cards, {"paths": {"metadata_csv": ""}})
    assert acts3[0].ok and "/data/m.csv" in acts3[0].preview


def test_build_actions_metadata_builder_dynamic_fields():
    acts = build_actions([
        {"name": "fill_metadata_builder",
         "arguments": {"field": "cell_line", "value": "Jurkat", "index": 0}},
        {"name": "fill_metadata_builder",
         "arguments": {"field": "cell_line", "value": "Jurkat", "index": 0,
                       "cell_type": "tcell"}},
    ], [], {})
    assert not acts[0].ok
    assert acts[1].ok and acts[1].data["cell_type"] == "tcell"
    assert "Sample 1 tcell line" in acts[1].preview


def test_metadata_summary_is_defensive():
    import pandas as pd
    assert summarize_metadata(None)["loaded"] is False
    assert summarize_metadata(pd.DataFrame())["loaded"] is False
    s = summarize_metadata(pd.DataFrame({"x": [1, 2]}))   # no sample_name col
    assert s["loaded"] and s["n_samples"] == 2 and s["sample_names"] == []


def test_diff_from_defaults():
    diffs = _diff_from_defaults({"tracking": {"immune": {"method": "btrack"}}})
    assert diffs["tracking.immune.method"] == {"current": "btrack", "default": "trackpy"}


def test_dump_cards_json(tmp_path=None):
    import tempfile, json
    d = tmp_path or tempfile.mkdtemp()
    path = os.path.join(str(d), "cards.json")
    n = dump_cards_json(path)
    assert n > 150 and json.load(open(path))


# --------------------------------------------------------------------------
# chatbot/ server helpers + ingest chunkers
# --------------------------------------------------------------------------
def test_app_tool_call_parsing_and_prompt():
    import app
    assert app.CONTROL_CONTRACT_VERSION == CONTROL_CONTRACT_VERSION
    txt = ('try 2 channels.\n<TOOLCALL>{"name":"set_ui_value",'
           '"arguments":{"control_id":"segmentation.cellpose.number_of_channels",'
           '"value":2}}</TOOLCALL>')
    clean, calls = app.parse_tool_calls(txt)
    assert calls[0]["name"] == "set_ui_value" and "<TOOLCALL>" not in clean
    assert app.split_streamable(txt) == "try 2 channels."
    assert app.parse_tool_calls("x <TOOLCALL>{bad}</TOOLCALL>")[1] == []
    sp = app.build_system_prompt(
        {"current_step": "tracking", "step_schema": [], "parameters": {}, "metadata": {}, "queue": []},
        [{"title": "btrack", "text": "good for crowded cells"}], [])
    assert "BEHAV3D" in sp and "good for crowded cells" in sp
    assert "at most the actions needed for this one user turn" in sp


def test_backend_hides_bulk_metadata_tool_after_forms_exist():
    import app
    tools = [
        {"name": "set_ui_value"},
        {"name": "bulk_fill_metadata"},
        {"name": "fill_metadata_builder"},
    ]
    fresh = app.tools_for_context(tools, {"metadata": {"loaded": False}})
    assert {tool["name"] for tool in fresh} == {
        "set_ui_value", "bulk_fill_metadata", "fill_metadata_builder",
    }
    existing = app.tools_for_context(tools, {
        "metadata": {"loaded": True, "record_source": "metadata_builder_draft"},
        "metadata_builder": {"open": True},
        "ui_state": {"controls": [{"id": "metadata.samples.0.pixel_distance_xy"}]},
    })
    assert {tool["name"] for tool in existing} == {
        "set_ui_value", "fill_metadata_builder",
    }
    description = (
        "I have three movies with two immune cell types and collagen, four channels, "
        "1.15 um pixel size and 4 um z spacing. Help me set up the analysis."
    )
    assert app.should_force_bulk_metadata(
        {"metadata": {"loaded": False}}, description, fresh,
    )
    assert not app.should_force_bulk_metadata(
        {"metadata": {"loaded": True}}, description, existing,
    )
    assert not app.should_force_bulk_metadata(
        {"metadata": {"loaded": False}}, "What is metadata?", fresh,
    )
    forced_choice, forced_extra = app.model_tool_policy(True, True)
    assert forced_choice is None
    assert forced_extra == {"thinking": {"type": "disabled"}}
    assert app.model_tool_policy(False, True) == ("auto", None)
    assert app.model_tool_policy(False, False) == (None, None)
    sanitized = app.sanitize_bulk_metadata_arguments({"samples": [{
        "sample_name": "Movie_1", "dimension_order": "TCZYX",
        "pixel_distance_xy": 1.15, "pixel_distance_z": 4,
        "time_interval": 15, "time_unit": "s",
        "cell_types": {"T cells": {}},
    }]}, "I have three movies sampled every 15, but I do not know the unit")
    assert sanitized == {"samples": [{
        "pixel_distance_xy": 1.15, "pixel_distance_z": 4, "time_interval": 15,
    }]}
    recovered = app.recover_single_control_action(
        [],
        {"ui_state": {"controls": [{
            "id": "tracking.tcell.btrack.maximum_search_radius",
            "label": "Maximum search radius", "value": 100,
            "visible": True, "enabled": True,
        }]}},
        "Please adjust the maximum search radius.",
        "Would you like me to set the maximum search radius to **18 um**?",
    )
    assert recovered == [{
        "name": "set_ui_value",
        "arguments": {
            "control_id": "tracking.tcell.btrack.maximum_search_radius", "value": 18,
        },
    }]
    ambiguous = app.recover_single_control_action(
        [],
        {"ui_state": {"controls": [
            {"id": "one", "value": 1, "visible": True, "enabled": True},
            {"id": "two", "value": 2, "visible": True, "enabled": True},
        ]}},
        "Please adjust these.", "Set both to 3.",
    )
    assert ambiguous == []


def test_tool_call_parsing_tolerates_malformed_markers():
    import app
    # the exact malformed shape a small model emitted: no leading '<', no closing tag
    malformed = ('TOOLCALL>{"name": "set_ui_value", "arguments": '
                 '{"control_id": "data.metadata_csv", "value": "/p/m.csv"}} '
                 'TOOLCALL>{"name": "set_ui_value", "arguments": '
                 '{"control_id": "data.dimension_order_all", "value": "TZCYX"}}')
    clean, calls = app.parse_tool_calls(malformed)
    assert len(calls) == 2 and "TOOLCALL" not in clean
    assert calls[0]["arguments"]["control_id"] == "data.metadata_csv"
    # canonical form with nested args still parses and is stripped from display
    clean2, calls2 = app.parse_tool_calls(
        'Sure.\n<TOOLCALL>{"name":"set_ui_value","arguments":'
        '{"control_id":"tracking.tcell.method","value":"btrack"}}</TOOLCALL>')
    assert calls2[0]["arguments"]["value"] == "btrack" and clean2 == "Sure."
    # bare proposal form + list value
    _, calls3 = app.parse_tool_calls(
        'set_ui_value{"control_id":"features.tcell.feature_groups",'
        '"value":["movement","contact"]}')
    assert calls3[0]["arguments"]["value"] == ["movement", "contact"]
    clean4, calls4 = app.parse_tool_calls(
        'Done.\nfill_metadata_builder{"field":"n_samples","value":3}')
    assert clean4 == "Done."
    assert calls4 == [{"name": "fill_metadata_builder",
                       "arguments": {"field": "n_samples", "value": 3}}]


def test_to_openai_tools_injects_key_enum():
    import app
    from behav3d.napari._assistant_actions import TOOL_SCHEMA
    enum = ["data.metadata_csv", "tracking.tcell.method"]
    tools = app.to_openai_tools(TOOL_SCHEMA, key_enum=enum)
    assert all(t["type"] == "function" for t in tools)
    sp = next(t for t in tools if t["function"]["name"] == "set_ui_value")
    assert sp["function"]["parameters"]["properties"]["control_id"]["enum"] == enum
    # original schema is not mutated (deep-copied)
    assert "enum" not in TOOL_SCHEMA[0]["parameters"]["properties"]["control_id"]


def test_assemble_tool_calls_merges_streamed_fragments():
    import app
    frags = [
        {"index": 0, "name": "set_ui_value", "arguments": ""},
        {"index": 0, "name": None, "arguments": '{"control_id": "tracking.tcell'},
        {"index": 0, "name": None, "arguments": '.trackpy.search_range", "value": 40}'},
    ]
    calls = app.assemble_tool_calls(frags)
    assert calls == [{"name": "set_ui_value",
                      "arguments": {"control_id": "tracking.tcell.trackpy.search_range",
                                    "value": 40}}]
    # incomplete / unparseable arguments are dropped, not crashed on
    assert app.assemble_tool_calls([{"index": 0, "name": "set_ui_value",
                                     "arguments": "{bad"}]) == []


def test_doc_module_lists_in_sync():
    # app.py inlines _DOC_PY_MODULES (can't import ingest at top level in-container);
    # guard against it drifting from ingest's authoritative copy.
    import app, ingest
    if getattr(app, "modal", None) is not None:
        assert app._DOC_PY_MODULES == ingest._DOC_PY_MODULES


def test_ingest_chunkers():
    import ingest
    md = "# Title\nbody text\n## Section\nmore body"
    chunks = ingest.chunk_markdown(md, "README.md")
    assert len(chunks) >= 2 and all(c.source == "README.md" for c in chunks)
    cards = [{"key": "a.b", "step": "tracking", "category": None, "type": "int",
              "default": 1, "choices": None, "description": "desc"}]
    cc = ingest.cards_to_chunks(cards)
    assert "Parameter: a.b" in cc[0].text and cc[0].title == "a.b"


def test_build_actions_bulk_fill_metadata():
    # One bulk_fill_metadata call carries the whole form; build_actions keeps the
    # payload and a sample count, and rejects an empty samples list.
    acts = build_actions([
        {"name": "bulk_fill_metadata", "arguments": {
            "n_samples": 2, "n_immune": 1, "immune_names": ["tcell"],
            "samples": [
                {"sample_name": "s1", "raw_image_path": "/a.tif",
                 "cell_types": {"tcell": {"line": "Jurkat"}}},
                {"sample_name": "s2", "raw_image_path": "/b.tif"},
            ]}},
        {"name": "bulk_fill_metadata", "arguments": {"samples": []}},
    ], [], {})
    assert acts[0].kind == "bulk_fill_metadata" and acts[0].ok
    assert acts[0].data["sample_count"] == 2
    assert acts[0].data["immune_names"] == ["tcell"]
    assert "2 samples" in acts[0].preview
    assert not acts[1].ok                       # empty samples → rejected


def test_build_actions_select_segmentation_method():
    acts = build_actions([
        {"name": "select_segmentation_method", "arguments": {"value": "Cellpose"}},
        {"name": "select_segmentation_method", "arguments": {"value": ""}},
    ], [], {})
    assert acts[0].kind == "select_segmentation_method" and acts[0].ok
    assert acts[0].data["value"] == "Cellpose" and "Cellpose" in acts[0].preview
    assert not acts[1].ok                        # empty value → rejected


def test_new_tools_registered():
    import app
    from behav3d.napari._assistant_actions import TOOL_SCHEMA
    names = {t["name"] for t in TOOL_SCHEMA}
    assert {"set_ui_value", "bulk_fill_metadata", "show_track_length_distribution",
            "create_cell_type_group", "create_btrack_config_copy", "open_result",
            "recommend_edt", "summarize_track_counts"} <= names
    assert "set_parameter" not in names
    # the backend text-fallback parser must also recognise the new tool names
    assert "bulk_fill_metadata" in app._TOOL_NAMES
    assert "set_ui_value" in app._TOOL_NAMES
    assert "recommend_edt" in app._TOOL_NAMES
    assert "summarize_track_counts" in app._TOOL_NAMES


def test_build_actions_for_edt_and_track_count_queries():
    actions = build_actions([
        {"name": "recommend_edt", "arguments": {
            "cell_type": "tcell", "cell_diameter_um": 10,
        }},
        {"name": "summarize_track_counts", "arguments": {
            "cell_type": "tcell", "position_t": 200,
            "minimum_lengths": [200, 30, 30],
        }},
    ], [], {})
    assert actions[0].ok and actions[0].kind == "recommend_edt"
    assert actions[1].ok and actions[1].data["minimum_lengths"] == [30, 200]
    assert actions[1].data["position_t"] == 200


def test_prompt_data_prep_reconciles_state_and_lists_tools():
    import app
    sp = app.build_system_prompt(
        {"current_step": "data_preparation", "step_schema": [], "parameters": {},
         "metadata": {"loaded": False},
         "metadata_builder": {"n_samples": 22, "n_immune": 1, "open": True,
                              "immune_names": ["tcell"]},
         "queue": []},
        [], [])
    # The full builder state is represented as structured context.
    assert '"n_samples": 22' in sp
    assert '"loaded": false' in sp
    assert "Never ask for a value already present" in sp
    assert "n_samples" in sp  # present only as structured live context


def test_prompt_treats_metadata_builder_draft_as_current():
    import app
    sp = app.build_system_prompt(
        {"current_step": "segmentation", "metadata": {
            "loaded": True, "record_source": "metadata_builder_draft",
            "records": [{"sample_name": "Sample01", "dimension_order": "TCZYX"}],
            "validation": [], "save_required": True,
        }, "metadata_builder": {"open": True}, "queue": []},
        [], [])
    assert "supersede the last saved DataFrame" in sp
    assert "draft changes still need to be saved" in sp


def test_experiment_reference_discovers_notes_and_compacts_configuration(tmp_path=None):
    import tempfile
    from pathlib import Path

    root = Path(tmp_path or tempfile.mkdtemp())
    (root / "README_BEHAV3D_Exp085.md").write_text(
        "# Experiment\nWithin-well paired KO versus rescue comparison.",
        encoding="utf-8",
    )
    (root / "behav3d_parameters_clean.yml").write_text(
        "\n".join([
            "paths:",
            "  metadata_csv: /private/raw/metadata.csv",
            "pixel_classifier:",
            "  apoc_strategy: APOC Probability Map + Watershed",
            "apoc:",
            "  apoc_organoid1_prob_mask_threshold: 0.55",
            "  apoc_organoid1_prob_seed_threshold: 1.0",
            "  apoc_organoid1_feature_string: very-large-internal-value",
            "tracking:",
            "  organoid1:",
            "    method: propagation",
            "features:",
            "  organoid1:",
            "    features_choice: [intensity, contact, death]",
            "    dead_mask_percentage_threshold: 2.0",
            "track_filtering:",
            "  organoid1:",
            "    min_length_enabled: true",
            "    min_track_length: 100",
        ]),
        encoding="utf-8",
    )
    reference = _experiment_reference_context(str(root), {})
    assert reference["notes"][0]["source"] == "README_BEHAV3D_Exp085.md"
    settings = reference["saved_configurations"][0]["settings"]
    assert settings["tracking"]["organoid1"]["method"] == "propagation"
    assert settings["segmentation"]["apoc_by_cell_type"]["organoid1"] == {
        "prob_mask_threshold": 0.55,
        "prob_seed_threshold": 1.0,
    }
    assert settings["features"]["organoid1"]["dead_mask_percentage_threshold"] == 2.0
    serialized = str(reference)
    assert "/private/raw" not in serialized
    assert "very-large-internal-value" not in serialized
    assert "not proof that a module ran" in reference["configuration_caveat"]


def test_prompt_scopes_experiment_reference_and_requires_result_evidence():
    import app
    reference = {
        "notes": [{
            "source": "README_BEHAV3D_Exp010.md",
            "text": "Safety profiling with tumor and healthy organoids.",
        }],
        "saved_configurations": [{
            "source": "behav3d_parameters.yml",
            "settings": {"features": {"TEG": {
                "features_choice": ["invasiveness"],
            }}},
        }],
    }
    sp = app.build_system_prompt(
        {
            "current_step": "analysis",
            "metadata": {"loaded": True, "records": []},
            "experiment_reference": reference,
            "results": [],
        },
        [], [],
    )
    assert "for this dataset only" in sp
    assert "not proof that" in sp
    assert "Claim an output is available only when" in sp
    assert '"experiment_reference"' in sp
    assert "Safety profiling with tumor and healthy organoids" in sp


def test_prompt_encodes_pi_feedback_scenarios():
    import app
    sp = app.build_system_prompt(
        {"current_step": "tracking", "metadata": {
            "loaded": True,
            "records": [{"sample_name": "Movie1", "time_interval": 15,
                         "time_unit": "s"}],
        }, "ui_state": {"controls": []}},
        [], [])
    assert "call bulk_fill_metadata directly" in sp
    assert "Calling only open_builder" in sp
    assert "Never infer a dimension order" in sp
    assert "never rebuild the form from values mentioned earlier" in sp
    assert "explicit request to fill, set, update, fix, or adjust" in sp
    assert "Seed threshold is the main splitting lever" in sp
    assert "Mask threshold primarily defines the foreground contour" in sp
    assert "With plain Mask + EDT/Watershed, raise EDT" in sp
    assert "With Peak EDT/Watershed, lower EDT" in sp
    assert "do not infer motion from a label" in sp
    assert "btrack is the routine default" in sp
    assert "Reporter Propagation" in sp
    assert "60 um/min at 15 seconds per frame is 15 um/frame" in sp
    assert "Do not provide a numeric example speed" in sp
    assert "merely because its inherited thresholds contain defaults" in sp
    assert "Filtering must be run even when all filters are disabled" in sp
    assert "keep Start offset at 1" in sp
    assert "Average linkage is the default" in sp
    assert "Never call that combination contradictory" in sp
    assert "Do not call a chosen minimum reasonable" in sp


class _FakeSpin:
    def __init__(self, value, minimum=0, maximum=9999):
        self._value = value
        self._minimum = minimum
        self._maximum = maximum

    def value(self): return self._value
    def setValue(self, value): self._value = type(self._value)(value)
    def minimum(self): return self._minimum
    def maximum(self): return self._maximum
    def isEnabled(self): return True
    def isHidden(self): return False


class _FakeCheck:
    def __init__(self, checked=False): self._checked = checked
    def isChecked(self): return self._checked
    def setChecked(self, value): self._checked = bool(value)
    def isEnabled(self): return True
    def isHidden(self): return False


class _FakeLine:
    def __init__(self, value=""): self._value = value
    def text(self): return self._value
    def setText(self, value): self._value = str(value)
    def isEnabled(self): return True
    def isHidden(self): return False


class _FakeGroup:
    def __init__(self, title="Sample 1", checked=True):
        self._title, self._checked = title, checked
    def title(self): return self._title
    def isChecked(self): return self._checked


class _FakeCombo:
    def __init__(self, items, index=0): self.items, self.index = list(items), index
    def count(self): return len(self.items)
    def itemText(self, index): return self.items[index]
    def currentText(self): return self.items[self.index]
    def currentIndex(self): return self.index
    def setCurrentIndex(self, index): self.index = index
    def isEnabled(self): return True
    def isHidden(self): return False


class _FakeListItem:
    def __init__(self, text, selected=False):
        self._text, self._selected = text, selected
    def text(self): return self._text
    def isSelected(self): return self._selected
    def setSelected(self, selected): self._selected = bool(selected)


class _FakeList:
    def __init__(self, items, selected=()):
        selected = set(selected)
        self.items = [_FakeListItem(item, item in selected) for item in items]
    def count(self): return len(self.items)
    def item(self, index): return self.items[index]
    def selectedItems(self): return [item for item in self.items if item.isSelected()]
    def isEnabled(self): return True
    def isHidden(self): return False


class _FakeTabs:
    def __init__(self, panel): self.panel = panel
    def currentWidget(self): return self.panel


class _FakeWorkflowTabs:
    def __init__(self, index=0): self.index = index
    def currentIndex(self): return self.index
    def setCurrentIndex(self, index): self.index = index
    def tabText(self, index):
        return ["Data Preparation", "Visualization", "Segmentation"][index]


class _FakeTrackingPanel:
    def __init__(self, cell_type, radius):
        self.cell_type = cell_type
        self.combo_method = _FakeCombo(
            ["LAP", "TrackPy", "Propagation", "Reporter Propagation", "btrack"],
            4,
        )
        self.lap_merge_cost = _FakeSpin(0)
        self.lap_split_cost = _FakeSpin(0)
        self.tp_adaptive_stop = _FakeSpin(10.0)
        self.tp_adaptive_step = _FakeSpin(0.95)
        self.bt_max_search_radius = _FakeSpin(radius, 1, 9999)
        self.bt_use_visual_features = _FakeCheck(False)
        self.bt_use_optimize = _FakeCheck(False)
        self.bt_dist_thresh = _FakeSpin(40.0)
        self.bt_time_thresh = _FakeSpin(5)
        self.bt_hyp_checks = {"P_FP": _FakeCheck(True), "P_branch": _FakeCheck(False)}
        self._bt_unit_mgr = type("_Unit", (), {"physical": True})()
        self.persisted = 0

    def _persist(self): self.persisted += 1


def test_live_control_registry_targets_actual_cell_type_only():
    from types import SimpleNamespace
    tcell = _FakeTrackingPanel("tcell", 100)
    organoid = _FakeTrackingPanel("organoid1", 200)
    organoid.bt_use_optimize.setChecked(True)
    tracking = SimpleNamespace(
        panels={"tcell": tcell, "organoid1": organoid},
        cell_tabs=_FakeTabs(tcell),
    )
    main = SimpleNamespace(tracking_tab=tracking)
    controls = control_registry(main)
    ids = {item["id"] for item in controls}
    assert "tracking.tcell.btrack.maximum_search_radius" in ids
    assert "tracking.organoid1.btrack.maximum_search_radius" in ids
    assert "tracking.tcell.lap.merging_distance" in ids
    assert "tracking.tcell.trackpy.adaptive_step" in ids
    by_id = {item["id"]: item for item in controls}
    assert by_id["tracking.tcell.btrack.use_global_optimization"]["visible"] is True
    assert by_id["tracking.tcell.btrack.maximum_search_radius"]["unit"] == "um"
    assert "tracking.tcell.btrack.distance_threshold" not in by_id
    assert "tracking.tcell.btrack.hypotheses" not in by_id
    assert by_id["tracking.organoid1.btrack.distance_threshold"]["visible"] is True
    assert by_id["tracking.organoid1.btrack.hypotheses"]["visible"] is True
    assert active_cell_type(main, "tracking") == "tcell"
    assert apply_set_ui_value(main, "tracking.tcell.btrack.maximum_search_radius", 125)
    assert tcell.bt_max_search_radius.value() == 125
    assert organoid.bt_max_search_radius.value() == 200
    assert tcell.persisted == 1 and organoid.persisted == 0


def test_tracking_registry_separates_reporter_propagation_from_btrack():
    from types import SimpleNamespace
    panel = _FakeTrackingPanel("reporter", 100)
    panel.combo_method.setCurrentIndex(3)
    panel.rp_min_overlap_fraction = _FakeSpin(0.1, 0.0, 1.0)
    panel.rp_segment_size_min = _FakeSpin(100)
    main = SimpleNamespace(
        tracking_tab=SimpleNamespace(
            panels={"reporter": panel},
            cell_tabs=_FakeTabs(panel),
        )
    )
    controls = {item["id"]: item for item in control_registry(main)}
    assert controls[
        "tracking.reporter.reporter_propagation.minimum_overlap"
    ]["visible"] is True
    assert controls[
        "tracking.reporter.btrack.maximum_search_radius"
    ]["visible"] is False


def test_apoc_same_value_is_noop_and_method_specific():
    from types import SimpleNamespace
    method = _FakeCombo([
        "APOC (GPU)", "ConvPaint", "Pixel Classifier (Random Forest)",
        "Cellpose", "Import segmentation",
    ], 0)
    seg = SimpleNamespace(
        method_combo=method,
        apoc_page=SimpleNamespace(spin_examples=_FakeSpin(3)),
        convpaint_page=SimpleNamespace(spin_examples=_FakeSpin(3)),
        pixel_classifier_page=SimpleNamespace(spin_examples=_FakeSpin(3)),
        cellpose_page=SimpleNamespace(),
    )
    main = SimpleNamespace(segmentation_tab=seg)
    controls = control_registry(main)
    visible = {item["id"] for item in controls if item["visible"]}
    assert "segmentation.apoc.examples_per_sample" in visible
    assert "segmentation.random_forest.examples_per_sample" not in visible
    actions = build_actions([{
        "name": "set_ui_value",
        "arguments": {"control_id": "segmentation.apoc.examples_per_sample", "value": 3},
    }], [], {}, controls=controls)
    assert len(actions) == 1 and actions[0].data.get("no_op") is True
    assert not actions[0].ok  # no confirmation card should be rendered


def test_cellpose_sam_registry_exposes_live_core_controls():
    from types import SimpleNamespace
    persisted = []
    channels = [_FakeCheck(True), _FakeCheck(False)]
    sam = SimpleNamespace(
        cell_type_combo=_FakeCombo(["tcell", "organoid1"]),
        check_all_cell_types=_FakeCheck(False),
        check_process_all=_FakeCheck(True),
        combo_gpu_device=_FakeCombo(["Auto", "CUDA 0"]),
        btn_force_cpu=_FakeCheck(False),
        combo_model=_FakeCombo(["cpsam", "cpsam_v2"]),
        spin_diameter=_FakeSpin(0.0, 0.0, 1000.0),
        spin_flow_threshold=_FakeSpin(0.4, 0.0, 3.0),
        spin_cellprob=_FakeSpin(0.0, -6.0, 6.0),
        check_do_3d=_FakeCheck(True),
        spin_stitch=_FakeSpin(0.0, 0.0, 1.0),
        spin_batch_size=_FakeSpin(8, 1, 256),
        check_drop_2d=_FakeCheck(True),
        spin_size_min=_FakeSpin(0),
        spin_size_max=_FakeSpin(0),
        channel_panel=SimpleNamespace(_checkboxes=channels),
        _unit_mgr=SimpleNamespace(physical=True),
        _persist_params=lambda: persisted.append(True),
    )
    segmentation = SimpleNamespace(
        method_combo=_FakeCombo([
            "APOC (GPU)", "ConvPaint", "Pixel Classifier (Random Forest)",
            "Cellpose", "Cellpose-SAM (zero-shot)", "Import segmentation",
        ], 4),
        apoc_page=SimpleNamespace(spin_examples=_FakeSpin(3)),
        convpaint_page=SimpleNamespace(spin_examples=_FakeSpin(3)),
        pixel_classifier_page=SimpleNamespace(spin_examples=_FakeSpin(3)),
        cellpose_page=SimpleNamespace(),
        cellpose_sam_page=sam,
    )
    main = SimpleNamespace(segmentation_tab=segmentation)
    controls = {item["id"]: item for item in control_registry(main)}
    assert controls["segmentation.cellpose_sam.diameter"]["unit"] == "um"
    assert controls["segmentation.cellpose_sam.flow_threshold"]["visible"] is False
    assert controls["segmentation.cellpose_sam.stitch_threshold"]["visible"] is False
    channel_id = "segmentation.cellpose_sam.tcell.channels"
    assert controls[channel_id]["choices"] == ["Channel 0", "Channel 1"]
    assert active_cell_type(main, "segmentation") == "tcell"
    assert apply_set_ui_value(main, "segmentation.cellpose_sam.batch_size", 4)
    assert sam.spin_batch_size.value() == 4
    assert persisted


def test_metadata_records_are_complete_and_validated():
    import pandas as pd
    frame = pd.DataFrame([
        {"sample_name": "s1", "raw_image_path": "/missing/a.tif",
         "dimension_order": "TCZYX", "pixel_distance_xy": 0.5,
         "pixel_distance_z": 2.0, "time_interval": 5, "time_unit": "min",
         "im_tcell_line_condition": "Jurkat_control"},
        {"sample_name": "s2", "raw_image_path": "/missing/b.tif",
         "dimension_order": "ZCTXY", "pixel_distance_xy": 0.5,
         "pixel_distance_z": 2.0, "time_interval": 5, "time_unit": "min",
         "im_tcell_line_condition": "Jurkat_treated"},
    ])
    summary = summarize_metadata(frame)
    assert len(summary["records"]) == 2
    assert summary["records"][1]["im_tcell_line_condition"] == "Jurkat_treated"
    assert any(issue["field"] == "dimension_order" and issue["severity"] == "error"
               for issue in summary["validation"])


def test_metadata_builder_draft_validation_uses_live_form_value():
    from types import SimpleNamespace
    form = {
        "group": _FakeGroup(),
        "basic": {
            "sample_name": _FakeLine("Sample01"),
            "raw_image_path": _FakeLine("/missing/sample.tif"),
            "dimension_order": _FakeLine("TCZYX"),
            "pixel_distance_xy": _FakeSpin(0.5),
            "pixel_distance_z": _FakeSpin(2.0),
            "time_interval": _FakeSpin(5.0),
            "time_unit": _FakeCombo(["s", "m", "h"], 1),
        },
        "dead_channel": {}, "cell_types": {},
    }
    dp = SimpleNamespace(
        builder_grp=_FakeGroup(checked=True), _sample_forms=[form],
        n_samples_spin=_FakeSpin(1), n_organoid_spin=_FakeSpin(0),
        n_immune_spin=_FakeSpin(0), n_other_spin=_FakeSpin(0),
        include_dead_cb=_FakeCheck(False), _organoid_name_edits=[],
        _immune_name_edits=[], _other_name_edits=[],
    )
    state = _metadata_builder_state(dp)
    assert state["draft_records"][0]["dimension_order"] == "TCZYX"
    assert not any(issue["field"] == "dimension_order"
                   for issue in state["draft_validation"])


def test_context_keeps_builder_draft_authoritative_after_tab_switch():
    import pandas as pd
    from types import SimpleNamespace
    form = {
        "group": _FakeGroup(),
        "basic": {
            "sample_name": _FakeLine("Sample01"),
            "raw_image_path": _FakeLine("/missing/sample.tif"),
            "dimension_order": _FakeLine("TCZYX"),
            "pixel_distance_xy": _FakeSpin(0.5),
            "pixel_distance_z": _FakeSpin(2.0),
            "time_interval": _FakeSpin(5.0),
            "time_unit": _FakeCombo(["s", "m", "h"], 1),
        },
        "dead_channel": {}, "cell_types": {},
    }
    # The saved DataFrame still has the old blank value while the live draft has
    # already been corrected. This is the exact state from the reported loop.
    dp = SimpleNamespace(
        metadata=pd.DataFrame([{
            "sample_name": "Sample01", "raw_image_path": "/missing/sample.tif",
            "dimension_order": None, "pixel_distance_xy": 0.5,
            "pixel_distance_z": 2.0, "time_interval": 5.0, "time_unit": "m",
        }]),
        output_dir="", behav3d_parameters={},
        builder_grp=_FakeGroup(checked=True), _sample_forms=[form],
        n_samples_spin=_FakeSpin(1), n_organoid_spin=_FakeSpin(0),
        n_immune_spin=_FakeSpin(0), n_other_spin=_FakeSpin(0),
        include_dead_cb=_FakeCheck(False), _organoid_name_edits=[],
        _immune_name_edits=[], _other_name_edits=[],
    )
    context = build_context(SimpleNamespace(
        tabs=_FakeWorkflowTabs(2), data_prep_tab=dp,
    ))
    assert context["current_step"] == "segmentation"
    assert context["metadata"]["record_source"] == "metadata_builder_draft"
    assert context["metadata"]["records"][0]["dimension_order"] == "TCZYX"
    assert not any(issue["field"] == "dimension_order"
                   for issue in context["metadata"]["validation"])
    assert context["metadata_builder"]["open"] is True


def test_live_registry_covers_filtering_and_hmm_controls():
    from types import SimpleNamespace
    filter_panel = SimpleNamespace(
        en_exp_duration=_FakeCheck(True), spin_exp_duration=_FakeSpin(350),
        en_min_length=_FakeCheck(True), spin_min_length=_FakeSpin(30),
        en_max_length=_FakeCheck(True), spin_max_length=_FakeSpin(30),
        check_split_long_tracks=_FakeCheck(False),
        check_filter_min_size=_FakeCheck(False), spin_min_size_t1=_FakeSpin(1000),
        check_filter_dead_t0=_FakeCheck(False),
        combo_time_type=_FakeCombo(["frames", "hours"]),
    )
    state = SimpleNamespace(
        spin_window_size=_FakeSpin(5),
        spin_hmm_feature_smoothing_window=_FakeSpin(5),
        spin_hmm_n_states=_FakeSpin(4),
        combo_hmm_n_states_mode=_FakeCombo(["fixed", "auto"]),
        spin_hmm_k_min=_FakeSpin(2), spin_hmm_k_max=_FakeSpin(8),
        chk_hmm_sticky=_FakeCheck(False),
        _timepoint_checkboxes={}, _logscale_checkboxes={}, _bingrp_checkboxes={},
    )
    main = SimpleNamespace(
        filtering_tab=SimpleNamespace(panels={"tcell": filter_panel}),
        analysis_tab=SimpleNamespace(single_cell_tab=SimpleNamespace(
            cell_type_combo=_FakeCombo(["tcell"]), state_tab=state,
        )),
    )
    controls = {item["id"]: item for item in control_registry(main)}
    assert "filtering.tcell.minimum_initial_size.enabled" in controls
    assert "filtering.tcell.dead_at_first_timepoint.enabled" in controls
    assert "filtering.tcell.time_unit" in controls
    assert controls["filtering.tcell.maximum_length.enabled"]["label"].endswith(
        "Trim retained tracks to a common length"
    )
    assert controls["filtering.tcell.maximum_length.timepoints"]["label"].endswith(
        "Common output track length"
    )
    assert "filtering.tcell.maximum_length.split_long_tracks" in controls
    states = "analysis.state_classification.tcell.number_of_states"
    assert states in controls and controls[states]["value"] == 4
    assert controls[states]["visible"] is True


def test_active_killing_registry_targets_one_selected_type():
    from types import SimpleNamespace
    active = SimpleNamespace(
        immune_combo=_FakeCombo(["tcell"]),
        target_list=_FakeList(
            ["organoid1", "organoid2"],
            selected=["organoid1", "organoid2"],
        ),
        spin_obs_window=_FakeSpin(5),
        death_signal_combo=_FakeCombo([
            "percentage_dead_mask", "mean_dead_dye", "nr_dead_mask_pixels",
        ]),
        check_abs_threshold=_FakeCheck(False),
        spin_threshold_mult=_FakeSpin(1.5),
        spin_abs_threshold=_FakeSpin(25.0),
        spin_min_contact=_FakeSpin(1),
        spin_top_n=_FakeSpin(10),
    )
    main = SimpleNamespace(feature_extraction_tab=SimpleNamespace(
        panels={},
        active_killing_panel=active,
        _ak_toggle_btn=_FakeCheck(True),
    ))
    controls = {item["id"]: item for item in control_registry(main)}
    target_id = "features.active_killing.target_types"
    assert controls[target_id]["visible"] is True
    assert controls[target_id]["value"] == ["organoid1", "organoid2"]
    assert apply_set_ui_value(main, target_id, ["organoid2"])
    assert controls["features.active_killing.absolute_threshold"]["visible"] is False
    assert active_cell_type(main, "feature_extraction") == "tcell"
    assert [item.text() for item in active.target_list.selectedItems()] == ["organoid2"]


def test_analysis_registry_exposes_only_active_trajectory_controls():
    from types import SimpleNamespace
    persisted = []
    track = SimpleNamespace(
        spin_traj_size=_FakeSpin(100),
        spin_n_clusters=_FakeSpin(6),
        combo_linkage=_FakeCombo(["average", "complete", "single"]),
        combo_trim=_FakeCombo(["first", "last"], 1),
        chk_split_long_tracks=_FakeCheck(False),
        chk_parallel=_FakeCheck(True),
        chk_save_dist=_FakeCheck(False),
        chk_use_original=_FakeCheck(False),
        spin_seed=_FakeSpin(12345),
        _persist_track_cfg=lambda cell_type: persisted.append(cell_type),
    )
    state = SimpleNamespace(
        spin_window_size=_FakeSpin(5),
        combo_hmm_n_states_mode=_FakeCombo(["fixed", "auto"]),
        chk_hmm_sticky=_FakeCheck(False),
        _timepoint_checkboxes={},
        _logscale_checkboxes={},
        _bingrp_checkboxes={},
    )
    single = SimpleNamespace(
        cell_type_combo=_FakeCombo(["tcell"]),
        state_tab=state,
        track_tab=track,
        _stack=_FakeCombo(["overview", "settings"], 1),
        inner_tabs=_FakeCombo(["Behavioral State", "State Trajectory"], 1),
    )
    main = SimpleNamespace(
        analysis_tab=SimpleNamespace(single_cell_tab=single)
    )
    controls = {item["id"]: item for item in control_registry(main)}
    trajectory = "analysis.state_trajectory.tcell.linkage"
    hmm = "analysis.state_classification.tcell.window_size"
    assert controls[trajectory]["visible"] is True
    assert controls[hmm]["visible"] is False
    assert apply_set_ui_value(main, trajectory, "complete")
    assert track.combo_linkage.currentText() == "complete"
    assert persisted == ["tcell"]


def test_edt_recommendations_use_per_sample_xy_resolution():
    result = calculate_edt_recommendations([
        {"sample_name": "Well_A1", "pixel_distance_xy": 0.5,
         "distance_unit": "um"},
        {"sample_name": "Well_A2", "pixel_distance_xy": 1000,
         "distance_unit": "nm"},
    ])
    assert result["rows"][0]["object_diameter_px"] == 20
    assert result["rows"][0]["edt_candidates_px"] == [4, 5, 6]
    assert result["rows"][1]["object_diameter_px"] == 10
    assert result["rows"][1]["edt_candidates_px"] == [2, 2.5, 3]
    text = format_edt_recommendations(result, "tcell")
    assert "Well_A1" in text and "global starting value" in text
    assert "higher EDT values generally split touching objects more" in text
    assert "On Peak EDT/Watershed, the direction is reversed" in text
    assert "Lower EDT values split touching objects more aggressively" not in text

    organoid = calculate_edt_recommendations(
        [{"sample_name": "Org", "pixel_distance_xy": 0.5,
          "distance_unit": "um"}],
        organoid_cells_across=5,
    )
    assert organoid["object_diameter_um"] == 50
    assert organoid["rows"][0]["object_diameter_px"] == 100


def test_streamed_assistant_text_uses_researcher_facing_labels():
    text = researcher_facing_text(
        "pixel_distance_xy is 0.5; position_t is 10 for sample_name and TrackID."
    )
    assert text == "XY pixel size is 0.5; timepoint is 10 for sample name and track ID."


def test_segmentation_registry_exposes_exact_edt_control():
    from types import SimpleNamespace
    edt = _FakeSpin(2.5, 0.0, 50.0)
    random_forest = SimpleNamespace(
        spin_examples=_FakeSpin(3),
        param_widgets={"tcell": {"edt": edt}},
        _persist_params=lambda: None,
    )
    segmentation = SimpleNamespace(
        method_combo=_FakeCombo([
            "APOC (GPU)", "ConvPaint", "Pixel Classifier (Random Forest)",
            "Cellpose", "Import segmentation",
        ], 2),
        apoc_page=SimpleNamespace(spin_examples=_FakeSpin(3)),
        convpaint_page=SimpleNamespace(spin_examples=_FakeSpin(3)),
        pixel_classifier_page=random_forest,
        cellpose_page=SimpleNamespace(),
    )
    main = SimpleNamespace(segmentation_tab=segmentation)
    controls = {item["id"]: item for item in control_registry(main)}
    control_id = "segmentation.random_forest.tcell.edt_threshold"
    assert control_id in controls and controls[control_id]["visible"] is True
    assert apply_set_ui_value(main, control_id, 5.0)
    assert edt.value() == 5.0


def test_segmentation_registry_exposes_probability_thresholds_and_strategy():
    from types import SimpleNamespace
    persisted = []
    strategy = _FakeCombo([
        "APOC Probability Map + Watershed",
        "APOC Mask + EDT/Watershed Resegmentation",
    ])
    tab = SimpleNamespace(
        _per_tab_strategy_combo=strategy,
        prob_mask_threshold_spin=_FakeSpin(0.4, 0.0, 1.0),
        prob_seed_threshold_spin=_FakeSpin(0.7, 0.0, 1.0),
        edt_threshold_spin=None,
    )
    training = SimpleNamespace(tabs={"Tcell_cmtrm": tab}, strategy_combo=None)
    apoc = SimpleNamespace(
        spin_examples=_FakeSpin(3), _training_widget=training,
        _collect_apoc_tab_config=lambda: {"ok": True},
        _save_apoc_params_to_yaml=lambda **kwargs: persisted.append(kwargs),
    )
    segmentation = SimpleNamespace(
        method_combo=_FakeCombo([
            "APOC (GPU)", "ConvPaint", "Pixel Classifier (Random Forest)",
            "Cellpose", "Import segmentation",
        ], 0),
        apoc_page=apoc,
        convpaint_page=SimpleNamespace(spin_examples=_FakeSpin(3)),
        pixel_classifier_page=SimpleNamespace(spin_examples=_FakeSpin(3)),
        cellpose_page=SimpleNamespace(),
    )
    main = SimpleNamespace(segmentation_tab=segmentation)
    controls = {item["id"]: item for item in control_registry(main)}
    base = "segmentation.apoc.Tcell_cmtrm"
    assert controls[f"{base}.mask_threshold"]["strategy"] == strategy.currentText()
    assert controls[f"{base}.seed_threshold"]["value"] == 0.7
    assert f"{base}.edt_threshold" not in controls
    assert apply_set_ui_value(main, f"{base}.mask_threshold", 0.45)
    assert apply_set_ui_value(main, f"{base}.seed_threshold", 0.75)
    assert tab.prob_mask_threshold_spin.value() == 0.45
    assert tab.prob_seed_threshold_spin.value() == 0.75
    assert len(persisted) == 2


def test_active_segmentation_cell_type_uses_method_subtab():
    from types import SimpleNamespace
    training = SimpleNamespace(
        tab_widget=_FakeCombo(["tcell", "organoid1"], 1),
        _tab_cell_types=["tcell", "organoid1"],
    )
    main = SimpleNamespace(segmentation_tab=SimpleNamespace(
        method_combo=_FakeCombo(["APOC", "ConvPaint"], 0),
        apoc_page=SimpleNamespace(_training_widget=training),
    ))
    assert active_cell_type(main, "segmentation") == "organoid1"


def test_track_count_table_filters_full_track_length_at_requested_timepoint():
    import pandas as pd
    rows = []
    for timepoint in range(20):
        rows.append({"sample_name": "Well_A1", "TrackID": 1,
                     "position_t": timepoint})
    for timepoint in range(50):
        rows.append({"sample_name": "Well_A1", "TrackID": 2,
                     "position_t": timepoint})
    for timepoint in range(1, 201):
        rows.append({"sample_name": "Well_A1", "TrackID": 3,
                     "position_t": timepoint})
    rows.extend([
        {"sample_name": "Well_A2", "TrackID": 1, "position_t": 0},
        {"sample_name": "Well_A2", "TrackID": 1, "position_t": 0},
    ])
    frame = pd.DataFrame(rows)
    table, details = calculate_track_count_table(
        frame, minimum_lengths=[20, 50, 100, 200], position_t=0
    )
    a1 = table.loc[table["sample_name"] == "Well_A1"].iloc[0]
    assert a1["Cell count at timepoint 0 (not filtered)"] == 2
    assert a1[">=20 timepoints"] == 2
    assert a1[">=50 timepoints"] == 1
    assert a1[">=100 timepoints"] == 0
    assert details["position_present"] is True

    at_200, _ = calculate_track_count_table(
        frame, minimum_lengths=[30, 200], position_t=200
    )
    a1_200 = at_200.loc[at_200["sample_name"] == "Well_A1"].iloc[0]
    assert a1_200[">=200 timepoints"] == 1


def test_track_count_summary_is_saved_under_quality_control(tmp_path=None):
    import pandas as pd
    import tempfile
    from pathlib import Path
    root = Path(tmp_path or tempfile.mkdtemp())
    source = root / "analysis" / "tcell" / "track_features"
    source.mkdir(parents=True)
    pd.DataFrame([
        {"sample_name": "Well_A1", "TrackID": 1, "position_t": timepoint}
        for timepoint in range(30)
    ]).to_csv(source / "BEHAV3D_tcell_combined_track_features.csv", index=False)
    result = generate_track_count_summary(
        root, "tcell", minimum_lengths=[20, 50, 100, 200], position_t=0
    )
    assert result["output_path"].is_file()
    assert result["output_path"].parent.name == "quality_control"
    assert result["output_path"].name == "BEHAV3D_tcell_track_count_summary.csv"
    markdown = format_track_count_summary(result, "tcell")
    assert "| Sample |" in markdown
    assert "sample_name" not in markdown
    query = generate_track_count_summary(
        root, "tcell", minimum_lengths=[30], position_t=20
    )
    assert query["output_path"].is_file()
    assert query["output_path"] != result["output_path"]
    catalog_entry = next(
        item for item in scan_outputs(root)
        if item.path == result["output_path"]
    )
    assert catalog_entry.category == "filtering"
    assert "timepoint 0" in catalog_entry.description


def test_feedback_guidance_fixtures():
    import json
    from pathlib import Path
    from guidance import KNOWLEDGE_VERSION, select_guidance_cards

    fixture = Path(__file__).parent / "fixtures" / "assistant_feedback_transcripts.json"
    cases = json.loads(fixture.read_text(encoding="utf-8"))
    assert KNOWLEDGE_VERSION == "2026.07.23.2"
    for case in cases:
        cards = select_guidance_cards(
            {"current_step": case["step"]}, case["user"], case.get("intent"))
        text = " ".join(card["text"] for card in cards).lower()
        for required in case.get("required", []):
            assert required.lower() in text, (case["id"], required)
        for forbidden in case.get("forbidden", []):
            assert forbidden.lower() not in text, (case["id"], forbidden)
    method_cards = select_guidance_cards(
        {
            "current_step": "segmentation",
            "segmentation": {"method": "Cellpose-SAM (zero-shot)"},
        },
        "My objects are touching; how should I tune this method?",
    )
    method_ids = {card["id"] for card in method_cards}
    assert "cellpose_sam" in method_ids and "apoc" not in method_ids
    experiment_cards = select_guidance_cards(
        {
            "current_step": "analysis",
            "experiment_reference": {"notes": [{"text": "paired design"}]},
        },
        "What should I compare?",
    )
    assert "experiment_design" in {card["id"] for card in experiment_cards}


def test_assistant_has_no_hidden_model_continuation():
    from pathlib import Path
    source = (Path(__file__).parents[1] / "behav3d" / "napari" / "_assistant.py").read_text()
    assert "def _auto_continue" not in source
    assert "_dispatch_proactive" not in source
    assert "Checking your setup" not in source
    assert CONTROL_CONTRACT_VERSION == "2.6"


# --------------------------------------------------------------------------
if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    passed = 0
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
        passed += 1
    print(f"\n{passed}/{len(fns)} tests passed")
