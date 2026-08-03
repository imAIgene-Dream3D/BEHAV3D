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
    ProposedAction, get_by_dotted, set_by_dotted, validate_value, build_actions,
    apply_action, apply_set_ui_value, normalize_metadata_line_value,
)
from behav3d.napari._assistant_context import (
    summarize_metadata, _diff_from_defaults, validate_metadata_records,
    _metadata_builder_state, _experiment_reference_context,
    _image_dimensions_state, _current_log_state, _segmentation_state,
    _interface_capabilities, _feature_extraction_state, build_context,
)
from behav3d.napari._assistant_controls import (
    CONTROL_CONTRACT_VERSION, active_cell_type, control_registry,
)
from behav3d.napari._assistant_recommendations import (
    calculate_edt_recommendations, format_edt_recommendations,
)
from behav3d.napari._assistant import (
    researcher_facing_text, streaming_transcript_block, transcript_block_role,
)
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


def test_api_streaming_text_sanitizes_split_internal_names():
    import app

    visible, pending = app.split_researcher_stream_buffer(
        "Use percentage_"
    )
    assert visible == "Use "
    visible2, pending = app.split_researcher_stream_buffer(
        pending + "dead_mask and killing_"
    )
    visible3, pending = app.split_researcher_stream_buffer(
        pending + "efficiency now."
    )
    visible4, pending = app.split_researcher_stream_buffer(pending, final=True)
    assert visible2 + visible3 + visible4 == (
        "dead-mask percentage and killing efficiency now."
    )
    assert pending == ""
    assert app.researcher_facing_text(
        "Use `percentage_dead_mask` and `killing_efficiency`."
    ) == "Use dead-mask percentage and killing efficiency."
    assert app.researcher_facing_text(
        "Use time_interval with time_unit."
    ) == "Use time interval with time unit."
    assert app.researcher_facing_text("Use `APOC`.") == "Use APOC."


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
    assert forced_choice == "required"
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
    absent = app.sanitize_bulk_metadata_arguments({"samples": [{
        "cell_types": {
            "T cells": {"line": "GD2_CART"},
            "Macrophages": {"line": None},
        },
    }]}, "T cells are GD2_CART and macrophages were not added")
    assert absent["samples"][0]["cell_types"]["Macrophages"]["line"] == "not_added"
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


def test_metadata_pixel_size_action_reuses_user_supplied_values():
    import app

    controls = []
    for index in range(2):
        controls.extend([
            {
                "id": f"metadata.samples.{index}.pixel_distance_xy",
                "value": 0.5,
            },
            {
                "id": f"metadata.samples.{index}.pixel_distance_z",
                "value": 2.0,
            },
        ])
    result = app.metadata_pixel_size_action(
        {
            "current_step": "data_preparation",
            "ui_state": {"controls": controls},
        },
        [
            {
                "role": "user",
                "content": "XY pixel size is 1.15 um and Z spacing is 4 um.",
            },
            {"role": "assistant", "content": "I can update the metadata."},
            {"role": "user", "content": "Can you fill in the correct pixel size?"},
        ],
    )

    changed = {
        call["arguments"]["control_id"]: call["arguments"]["value"]
        for call in result["calls"]
    }
    assert changed == {
        "metadata.samples.0.pixel_distance_xy": 1.15,
        "metadata.samples.0.pixel_distance_z": 4.0,
        "metadata.samples.1.pixel_distance_xy": 1.15,
        "metadata.samples.1.pixel_distance_z": 4.0,
    }
    assert "time unit is unchanged" in result["text"]


def test_metadata_pixel_size_action_does_not_intercept_fresh_setup():
    import app

    result = app.metadata_pixel_size_action(
        {
            "current_step": "data_preparation",
            "metadata": {"loaded": False},
            "metadata_builder": {"open": False},
            "ui_state": {"controls": []},
        },
        [{
            "role": "user",
            "content": (
                "I want to set up an analysis with three movies. XY pixel size is "
                "1.15 um and Z spacing is 4 um. Help me set up the analysis."
            ),
        }],
    )

    assert result is None


def test_backend_exposes_metadata_persistence_only_when_available():
    import app

    tools = [
        {"name": "save_metadata"},
        {"name": "load_metadata"},
        {"name": "set_ui_value"},
    ]
    unavailable = app.tools_for_context(tools, {
        "metadata": {"loaded": False},
        "metadata_builder": {"actions": {
            "save_available": False, "load_available": False,
        }},
    })
    assert {item["name"] for item in unavailable} == {"set_ui_value"}

    available = app.tools_for_context(tools, {
        "metadata": {"loaded": False},
        "metadata_builder": {"actions": {
            "save_available": True, "load_available": True,
        }},
    })
    assert {item["name"] for item in available} == {
        "save_metadata", "load_metadata", "set_ui_value",
    }


def test_organoid_lines_require_processing_group_confirmation():
    import app

    context = {"metadata": {"loaded": False}, "metadata_builder": {}}
    messages = [{
        "role": "user",
        "content": (
            "I co-culture T cells and macrophages with different organoid lines. "
            "Can you set up the metadata?"
        ),
    }]
    response = app.organoid_processing_question(context, messages)
    assert "separate organoid types" in response
    assert "one organoid type" in response
    assert "segmented or tracked separately" in response

    messages.append({"role": "assistant", "content": response})
    messages.append({
        "role": "user",
        "content": "Treat them as one organoid type and record the line per sample.",
    })
    assert app.organoid_processing_question(context, messages) is None


def test_analysis_question_on_metadata_tab_is_not_hijacked_by_organoid_setup():
    import app

    context = {
        "current_step": "data_preparation",
        "metadata": {"loaded": False, "records": []},
        "metadata_builder": {"sample_forms_created": False},
    }
    messages = [{
        "role": "user",
        "content": (
            "I have an experiment with T cells, Macrophages and Organoids of "
            "different lines. What analysis would be possible for this data?"
        ),
    }]
    assert app.organoid_processing_question(context, messages) is None
    result = app.deterministic_turn_response(context, messages, [])
    text = result["text"]
    assert "No metadata is loaded" in text
    assert "Interaction Analysis" in text
    assert "Invasiveness Analysis" in text
    assert "Contact-Based Grouping" in text
    assert "Contact State-Shift Analysis" in text
    assert "Before I build the metadata" not in text
    assert "separate organoid types" not in text


def test_metadata_completion_uses_mandatory_well_and_line_fields():
    import app

    context = {
        "metadata": {
            "records": [{"sample_name": "Sample01"}],
            "validation": [
                {"severity": "error", "message": "Sample01: mandatory well is missing."},
                {
                    "severity": "error",
                    "message": "Sample01: mandatory T-cells line is missing.",
                },
            ],
        },
        "metadata_builder": {"sample_forms_created": True},
    }
    response = app.metadata_completion_summary(
        context, [{"role": "user", "content": "Is this all that is needed?"}]
    )
    assert "mandatory well" in response
    assert "mandatory T-cells line" in response
    assert "condition" in response and "optional" in response
    assert "1" in response


def test_metadata_identifier_inferences_require_confirmation():
    import app

    response = app.metadata_identifier_confirmation_question(
        {"metadata_builder": {"sample_forms_created": True}},
        [{"role": "user", "content": (
            "No well identifier. M21 and M23 are macrophage lines and "
            "T-cells are GD2_CART."
        )}],
    )
    assert "mandatory" in response
    assert "condition is optional" in response
    assert "not_added" in response
    assert "confirm" in response
    assert "filenames" in response
    assert "M21" not in response and "M23" not in response
    assert "GD2" not in response

    missing_lines = app.metadata_identifier_confirmation_question(
        {"metadata_builder": {"sample_forms_created": True}},
        [{"role": "user", "content": "I have not specified lines yet."}],
    )
    assert "line for every configured" in missing_lines
    assert "not_added" in missing_lines


def test_metadata_time_conversion_updates_every_sample_without_stale_action():
    import app

    controls = []
    for index in range(3):
        controls.extend([
            {
                "id": f"metadata.samples.{index}.time_interval",
                "value": 2.0, "visible": True, "enabled": True,
            },
            {
                "id": f"metadata.samples.{index}.time_unit",
                "value": "s", "visible": True, "enabled": True,
            },
        ])
    action = app.metadata_time_conversion_action(
        {"ui_state": {"controls": controls}},
        [{"role": "user", "content": "The time unit is s while I set 2 minutes."}],
    )
    assert "120 seconds" in action["text"]
    assert len(action["calls"]) == 3
    assert all(call["arguments"]["value"] == 120 for call in action["calls"])
    assert all(
        call["arguments"]["control_id"].endswith(".time_interval")
        for call in action["calls"]
    )


def test_metadata_save_request_uses_real_action_or_reports_blocker():
    import app

    available = app.metadata_persistence_action(
        {"metadata_builder": {"actions": {"save_available": True}}},
        [{"role": "user", "content": "Yeah save it please"}],
    )
    assert available["calls"] == [{"name": "save_metadata", "arguments": {}}]
    assert "confirm" in available["text"].lower()

    blocked = app.metadata_persistence_action(
        {"metadata_builder": {
            "actions": {"save_available": False},
            "draft_validation": [{
                "severity": "error",
                "message": "Sample01: mandatory well is missing.",
            }],
        }},
        [{"role": "user", "content": "Save the metadata"}],
    )
    assert blocked["calls"] == []
    assert "cannot save" in blocked["text"].lower()
    assert "mandatory well" in blocked["text"]


def test_analysis_choose_explains_and_named_view_opens_directly():
    import app

    summary = app.analysis_choice_summary(
        {
            "assistant_session": {"intent": "choose_analysis"},
            "metadata": {
                "loaded": True,
                "n_samples": 8,
                "cell_types": {
                    "organoid": ["Organoids"],
                    "immune": ["Macrophages", "T-cells"],
                    "other": [],
                },
                "records": [{
                    "or_Organoids_line_condition": "DO7",
                    "im_T-cells_line_condition": "GD2_CART",
                    "im_Macrophages_line_condition": "M21",
                    "dead_channel": 3,
                }],
            },
        },
        [{"role": "user", "content": "Can you help me pick what analysis is nice?"}],
    )
    for name in (
        "Death Dynamics", "Interaction Analysis", "Invasiveness Analysis",
        "Active Killing", "Behavioral State", "State Trajectory",
        "Contact-Based Grouping", "Contact State-Shift Analysis", "Backprojection",
    ):
        assert name in summary
    assert "8 samples" in summary
    assert "DO7" in summary and "GD2_CART" in summary and "M21" in summary
    assert "Questions suggested by your metadata" in summary
    assert "Categorical DTW" in summary

    action = app.analysis_navigation_action(
        {"current_step": "analysis", "analysis": {"view": "behavioral_state"}},
        [{"role": "user", "content": "Can you take me to Death Dynamics?"}],
    )
    assert action["calls"] == [{
        "name": "open_analysis_view",
        "arguments": {"view": "death_dynamics"},
    }]
    already_open = app.analysis_navigation_action(
        {"current_step": "analysis", "analysis": {"view": "death_dynamics"}},
        [{"role": "user", "content": "Open Death Dynamics"}],
    )
    assert already_open["calls"] == []
    assert "already" in already_open["text"].lower()
    ambiguous = app.recover_single_control_action(
        [],
        {"ui_state": {"controls": [
            {"id": "one", "value": 1, "visible": True, "enabled": True},
            {"id": "two", "value": 2, "visible": True, "enabled": True},
        ]}},
        "Please adjust these.", "Set both to 3.",
    )
    assert ambiguous == []


def test_generic_tracking_guide_asks_for_motion_before_method_context():
    import app

    context = {
        "current_step": "tracking",
        "active_cell_type": "Tcell_HIV",
    }
    question = app.tracking_motion_question(
        context, [{"role": "user", "content": "Guide tracking"}],
    )
    assert "Tcell_HIV" in question
    assert "consecutive frames" in question
    assert not any(name in question.lower() for name in (
        "btrack", "propagation", "trackpy", "laptrack",
    ))
    exact = app.tracking_motion_question(
        context, [{"role": "user", "content": "Which method?"}],
    )
    assert "consecutive frames" in exact
    assert "segmentation" not in exact.lower()
    assert app.tracking_motion_question(
        context,
        [{"role": "user", "content": (
            "The cells move about 12 microns per frame. Help me choose a tracking method."
        )}],
    ) is None
    assert app.tracking_motion_question(
        {"current_step": "filtering"},
        [{"role": "user", "content": "Guide tracking"}],
    ) is None


def test_segmentation_method_gate_requires_signal_cleanliness():
    import app

    context = {
        "current_step": "segmentation",
        "assistant_session": {"intent": "compare_segmentation_methods"},
    }
    question = app.segmentation_signal_question(
        context, [{"role": "user", "content": "Choose a method"}],
    )
    assert "clean" in question
    assert "bleed-through" in question
    assert "same channel" in question
    assert app.segmentation_signal_question(
        context,
        [{"role": "user", "content": (
            "The target is isolated in a clean channel. Which segmentation method?"
        )}],
    ) is None


def test_stalled_operation_gate_requests_log_without_diagnosing():
    import app

    context = {
        "current_step": "segmentation",
        "current_log": {
            "recent_lines": ["Loading training data..."],
            "has_explicit_error": False,
        },
    }
    question = app.missing_log_error_question(
        context,
        [{"role": "user", "content": (
            "I clicked Generate Training Data but nothing appeared. What is wrong?"
        )}],
    )
    assert "copy and paste" in question.lower()
    assert "exact message" in question.lower()
    assert not any(term in question.lower() for term in (
        "file path", "dimension order", "gpu memory",
    ))


def test_tracking_radius_is_calculated_from_speed_and_cadence():
    import app

    result = app.tracking_radius_action(
        {
            "current_step": "tracking",
            "active_cell_type": "Tcell_cmtrm",
            "metadata": {"records": [{
                "time_interval": 15, "time_unit": "s",
            }]},
            "ui_state": {"controls": [{
                "id": "tracking.Tcell_cmtrm.btrack.maximum_search_radius",
                "cell_type": "Tcell_cmtrm",
                "visible": True,
                "enabled": True,
            }]},
        },
        [{"role": "user", "content": (
            "The fastest cells move 60 um per minute. "
            "Please adjust the maximum search radius."
        )}],
    )
    assert result["calls"][0]["arguments"]["value"] == 18
    assert "15 µm per frame" in result["text"]
    assert "20% margin" in result["text"]


def test_active_killing_action_is_cadence_grounded_and_explicitly_proposed():
    import app

    controls = [
        {
            "id": control_id,
            "visible": True,
            "enabled": True,
        }
        for control_id in (
            "features.active_killing.target_types",
            "features.active_killing.observation_window",
            "features.active_killing.death_signal",
            "features.active_killing.use_absolute_threshold",
            "features.active_killing.absolute_threshold",
        )
    ]
    result = app.active_killing_action(
        {
            "current_step": "feature_extraction",
            "metadata": {"records": [{
                "time_interval": 2, "time_unit": "min",
            }]},
            "ui_state": {"controls": controls},
        },
        [{"role": "user", "content": (
            "Configure active killing for tcell against organoid1 only. "
            "I expect killing within 10 minutes. "
            "Use an absolute threshold of 30 dead pixels."
        )}],
    )
    values = {
        call["arguments"]["control_id"]: call["arguments"]["value"]
        for call in result["calls"]
    }
    assert values["features.active_killing.target_types"] == ["organoid1"]
    assert values["features.active_killing.observation_window"] == 5
    assert values["features.active_killing.death_signal"] == "Dead-mask pixel count"
    assert values["features.active_killing.use_absolute_threshold"] is True
    assert values["features.active_killing.absolute_threshold"] == 30
    assert "proposing" in result["text"]
    assert "require your confirmation" in result["text"]


def test_active_killing_action_refuses_incomplete_absolute_mode():
    import app

    result = app.active_killing_action(
        {
            "current_step": "feature_extraction",
            "metadata": {"records": [{
                "time_interval": 2, "time_unit": "min",
            }]},
            "ui_state": {"controls": []},
        },
        [{"role": "user", "content": (
            "Configure active killing against 27t and mdo within 10 minutes."
        )}],
    )
    assert result["calls"] == []
    assert "positive minimum dead-mask pixel increase" in result["text"]
    assert "remains 0" in result["text"]


def test_active_killing_acceptance_proposes_every_agreed_field():
    import app

    controls = []
    specs = {
        "target_types": {
            "value": ["27t", "mdo"], "choices": ["27t", "mdo"],
        },
        "observation_window": {"value": 5},
        "death_signal": {
            "value": "Dead-mask percentage",
            "choices": [
                "Dead-mask percentage", "Mean dead-dye intensity",
                "Dead-mask pixel count",
            ],
        },
        "use_absolute_threshold": {"value": False},
        "absolute_threshold": {"value": 0.0, "active": False},
        "minimum_contact_duration": {"value": 1},
    }
    for suffix, values in specs.items():
        controls.append({
            "id": f"features.active_killing.{suffix}",
            "visible": True,
            "enabled": True,
            **values,
        })
    result = app.active_killing_confirmation_action(
        {
            "current_step": "feature_extraction",
            "ui_state": {"controls": controls},
        },
        [
            {"role": "assistant", "content": (
                "Active Killing configuration for tcell against 27t and mdo: "
                "Death signal: Dead-mask pixel count. Absolute threshold: "
                "30 dead pixels. Observation window: 5 timepoints. "
                "Minimum contact duration: 1 frame."
            )},
            {"role": "user", "content": "Ok, these settings seem ok"},
        ],
    )
    values = {
        call["arguments"]["control_id"]: call["arguments"]["value"]
        for call in result["calls"]
    }
    assert values["features.active_killing.target_types"] == ["27t", "mdo"]
    assert values["features.active_killing.death_signal"] == "Dead-mask pixel count"
    assert values["features.active_killing.use_absolute_threshold"] is True
    assert values["features.active_killing.absolute_threshold"] == 30
    assert values["features.active_killing.observation_window"] == 5
    assert values["features.active_killing.minimum_contact_duration"] == 1
    assert "complete agreed Active Killing setup" in result["text"]
    assert "independently" in result["text"]
    assert "not ready until every action card" in result["text"]


def test_active_killing_readiness_rejects_zero_absolute_threshold():
    import app

    summary = app.active_killing_readiness_summary(
        {
            "current_step": "feature_extraction",
            "feature_extraction": {"active_killing": {
                "setup_ready": False,
                "setup_issues": [
                    "Absolute signal-increase threshold must be greater than 0."
                ],
            }},
        },
        [
            {"role": "assistant", "content": "Active Killing setup"},
            {"role": "user", "content": "Is it ready?"},
        ],
    )
    assert "not ready yet" in summary
    assert "greater than 0" in summary


def test_hmm_movement_guidance_lists_every_live_option_before_editing():
    import app

    prefix = "analysis.state_classification.tcell."
    context = {
        "current_step": "analysis",
        "analysis": {
            "view": "behavioral_state",
            "selected_cell_type": "tcell",
        },
        "ui_state": {"controls": [
            {
                "id": prefix + "timepoint_features",
                "visible": True,
                "enabled": True,
                "choices": [
                    "speed", "displacement", "cumulative_displacement",
                    "directional_persistence", "mean_dead_dye",
                ],
                "value": ["speed"],
            },
            {
                "id": prefix + "window_features",
                "visible": True,
                "enabled": True,
                "choices": [
                    "net_displacement", "straightness",
                    "mean_square_displacement",
                ],
                "value": ["net_displacement"],
            },
            {
                "id": prefix + "binary_feature_groups",
                "visible": True,
                "enabled": True,
                "choices": ["27t_contact", "mdo_contact"],
                "value": ["27t_contact", "mdo_contact"],
            },
        ]},
    }
    text = app.hmm_movement_feature_guidance(
        context,
        [{"role": "user", "content": "Only movement features of relevance"}],
    )
    for label in (
        "Speed", "Displacement", "Cumulative displacement",
        "Directional persistence", "Net displacement", "Straightness",
        "Mean square displacement",
    ):
        assert label in text
    assert "Mean dead dye" not in text
    assert "use all available movement features" in text
    assert "applied after the HMM" in text


def test_hmm_all_movement_choice_updates_both_feature_lists():
    import app

    prefix = "analysis.state_classification.tcell."
    context = {
        "current_step": "analysis",
        "analysis": {
            "view": "behavioral_state",
            "selected_cell_type": "tcell",
        },
        "ui_state": {"controls": [
            {
                "id": prefix + "timepoint_features",
                "visible": True,
                "enabled": True,
                "choices": [
                    "speed", "displacement", "directional_persistence",
                    "mean_dead_dye",
                ],
                "value": ["speed"],
            },
            {
                "id": prefix + "window_features",
                "visible": True,
                "enabled": True,
                "choices": ["net_displacement", "straightness"],
                "value": ["net_displacement"],
            },
            {
                "id": prefix + "binary_feature_groups",
                "visible": True,
                "enabled": True,
                "choices": ["27t_contact", "mdo_contact"],
                "value": ["27t_contact", "mdo_contact"],
            },
        ]},
    }
    result = app.hmm_movement_setup_action(
        context,
        [
            {"role": "assistant", "content": (
                "Choose names or say use all available movement features."
            )},
            {"role": "user", "content": "Use all available movement features"},
        ],
    )
    values = {
        call["arguments"]["control_id"]: call["arguments"]["value"]
        for call in result["calls"]
    }
    assert values[prefix + "timepoint_features"] == [
        "speed", "displacement", "directional_persistence",
    ]
    assert values[prefix + "window_features"] == [
        "net_displacement", "straightness",
    ]
    assert "mean_dead_dye" not in str(values)


def _tcell_hmm_context():
    prefix = "analysis.state_classification.T-cells."
    return {
        "current_step": "analysis",
        "active_cell_type": "T-cells",
        "analysis": {
            "view": "behavioral_state",
            "selected_cell_type": "T-cells",
        },
        "ui_state": {"controls": [
            {
                "id": prefix + "timepoint_features",
                "visible": True, "enabled": True,
                "choices": ["speed", "elongation"],
                "value": ["speed"],
            },
            {
                "id": prefix + "window_features",
                "visible": True, "enabled": True,
                "choices": ["net_displacement", "straightness"],
                "value": ["net_displacement"],
            },
            {
                "id": prefix + "binary_feature_groups",
                "visible": True, "enabled": True,
                "choices": [
                    "Organoid_contact", "Macrophages_contact", "dead",
                ],
                "value": [],
            },
            {
                "id": prefix + "n_states",
                "visible": True, "enabled": True, "value": 6,
            },
        ]},
    }


def test_hmm_guidance_uses_live_selected_cell_type():
    import app

    context = _tcell_hmm_context()
    setup = app.hmm_setup_guidance(
        context,
        [{"role": "user", "content": (
            "I want to do behavioral analysis, can you take me through the steps?"
        )}],
    )
    assert "T-cells" in setup
    assert "currently have" in setup
    assert "rename" in setup.lower() and "merge" in setup.lower()

    groups = app.hmm_binary_group_guidance(
        context,
        [{"role": "user", "content": (
            "Would it be worth adding macrophage contact?"
        )}],
    )
    assert "T-cells" in groups
    assert "a cell from T-cells is directly touching Macrophages".lower() in groups.lower()
    assert "not a different population" in groups


def test_hmm_binary_group_edit_targets_selected_cell_control():
    import app

    context = _tcell_hmm_context()
    result = app.hmm_binary_group_action(
        context,
        [{"role": "user", "content": "Add organoid contact and also add dead"}],
    )
    assert "T-cells" in result["text"]
    assert result["calls"] == [{
        "name": "set_ui_value",
        "arguments": {
            "control_id": (
                "analysis.state_classification.T-cells.binary_feature_groups"
            ),
            "value": ["Organoid_contact", "dead"],
        },
    }]


def test_hmm_states_can_be_merged_by_reusing_names():
    import app

    response = app.hmm_state_merge_guidance(
        _tcell_hmm_context(),
        [{"role": "user", "content": (
            "If I have 6 states, can I select which ones to keep?"
        )}],
    )
    assert "T-cells" in response
    assert "same name" in response
    assert "merge" in response
    assert "Rename Primary Dynamic State Clusters" in response
    assert "outside BEHAV3D" not in response


def test_chat_transcript_shows_waiting_block_and_tracks_role_colors():
    assert streaming_transcript_block(True, "") == (
        "**BEHAV3D Assistant**\n\n*Preparing a response...*"
    )
    assert streaming_transcript_block(False, None) is None
    role = transcript_block_role("You")
    assert role == "user"
    assert transcript_block_role("My question", role) == "user"
    assert transcript_block_role("BEHAV3D Assistant", role) == "assistant"


def test_metadata_absence_values_are_csv_safe():
    import app

    assert normalize_metadata_line_value(None) == "not_added"
    assert normalize_metadata_line_value("None") == "not_added"
    assert normalize_metadata_line_value("(not_added)") == "not_added"
    assert normalize_metadata_line_value("M21") == "M21"
    controls = [{
        "id": "metadata.samples.0.cell_types.Macrophages.line",
        "label": "Sample 1, Macrophages: line",
        "value": "",
        "visible": True,
        "enabled": True,
    }]
    action = build_actions([{
        "name": "set_ui_value",
        "arguments": {
            "control_id": "metadata.samples.0.cell_types.Macrophages.line",
            "value": None,
        },
    }], [], {}, controls=controls)[0]
    assert action.ok
    assert action.data["value"] == "not_added"
    result = app.metadata_absence_action(
        {
            "current_step": "data_preparation",
            "ui_state": {"controls": controls},
        },
        [{"role": "user", "content": (
            "Macrophages were not added in Sample 1; set that line."
        )}],
    )
    assert "not added" in result["text"]
    assert result["calls"][0]["arguments"]["value"] == "not_added"


def test_equal_track_filter_review_explains_valid_fixed_length_workflow():
    import app

    controls = [
        {"id": "filtering.tcell.minimum_length.enabled", "value": True},
        {"id": "filtering.tcell.minimum_length.timepoints", "value": 30},
        {"id": "filtering.tcell.maximum_length.enabled", "value": True},
        {"id": "filtering.tcell.maximum_length.timepoints", "value": 30},
    ]
    summary = app.equal_track_filter_summary(
        {
            "current_step": "filtering",
            "active_cell_type": "tcell",
            "ui_state": {"controls": controls},
        },
        [{"role": "user", "content": "Review filters"}],
    )
    assert "valid" in summary
    assert "removes tracks shorter" in summary
    assert "trims retained longer tracks" in summary
    assert "uniform window" in summary


def test_merged_probability_watershed_guidance_explains_seed_confidence():
    import app

    result = app.merged_probability_watershed_guidance(
        {
            "current_step": "segmentation",
            "active_cell_type": "tcell",
            "segmentation": {"apoc": {"cell_type_strategies": {
                "tcell": "APOC Probability Map + Watershed",
            }}},
            "ui_state": {"controls": [
                {"id": "segmentation.apoc.tcell.seed_threshold", "value": 0.7},
                {"id": "segmentation.apoc.tcell.mask_threshold", "value": 0.4},
            ]},
        },
        [{"role": "user", "content": "Touching cells are still not split."}],
    )
    assert "raise it in small increments" in result
    assert "higher-confidence seed cores" in result
    assert "less confidence" not in result


def test_metadata_channel_mapping_stays_out_of_metadata_builder():
    import app

    messages = [
        {"role": "user", "content": (
            "Two replicates have the immune channels swapped. How should I set that up?"
        )},
        {"role": "assistant", "content": "Tell me where you want to map them."},
        {"role": "user", "content": (
            "Walk me through exactly where I name them and set the channels."
        )},
    ]
    result = app.metadata_channel_mapping_guidance({}, messages)
    assert "does not map raw channel indices" in result
    assert "line" in result.lower()
    assert "convpaint" in result.lower()
    assert "shared multiclass model" in result.lower()
    assert "not independent per sample" in result.lower()
    assert "i should not" not in result.lower()


def test_apoc_channel_guidance_rejects_invalid_index_and_dead_negative_class():
    import app

    context = {
        "current_step": "segmentation",
        "metadata": {"image_dimensions": [{"channel_count": 4}]},
    }
    invalid = app.apoc_channel_selection_guidance(
        context,
        [{"role": "user", "content": (
            "13T is ch1, blue ch0, green ch4. Which channels do I pick for APOC?"
        )}],
    )
    assert "indexed 0-3" in invalid
    assert "conflicts with the data" in invalid
    assert "before recommending" in invalid

    corrected = app.apoc_channel_selection_guidance(
        context,
        [{"role": "user", "content": (
            "13T is ch1, blue ch0, green ch3, dead ch2. "
            "Which channels do I pick for APOC?"
        )}],
    )
    assert "13t -> channel 1" in corrected.lower()
    assert "blue -> channel 0" in corrected.lower()
    assert "green -> channel 3" in corrected.lower()
    assert "do not automatically add **channel 2**" in corrected.lower()
    assert "negative/background class" in corrected.lower()


def test_apoc_channel_guidance_waits_for_training_controls():
    import app

    result = app.apoc_channel_selection_guidance(
        {
            "current_step": "segmentation",
            "active_cell_type": "tcell",
            "segmentation": {"apoc": {
                "training_data_loaded": False,
                "cell_types": {"tcell": {"channel_controls_ready": False}},
            }},
        },
        [{"role": "user", "content": (
            "How do I choose the APOC image channel inputs?"
        )}],
    )
    assert "Generate Training Data" in result
    assert "does not mean APOC uses every channel" in result
    assert "not a reason to switch" in result


def test_apoc_feature_grid_guidance_names_real_filters_and_scope():
    import app

    result = app.apoc_feature_grid_guidance(
        {
            "current_step": "segmentation",
            "segmentation": {"method": "APOC (GPU)"},
            "ui_state": {"controls": [
                {
                    "id": "segmentation.apoc.tcell.feature_scales",
                    "value": "1, 2, 5",
                },
                {
                    "id": "segmentation.apoc.tcell.feature_filters",
                    "value": ["Gaussian blur at sigma 1 px"],
                },
            ]},
        },
        [{"role": "user", "content": (
            "Can I tune the actual APOC feature values and filter sigmas?"
        )}],
    )
    for phrase in (
        "feature scales in pixels", "Gaussian blur", "Difference of Gaussians",
        "Laplacian of Gaussian", "Sobel-of-Gaussian",
        "not a structure tensor", "Small structures", "1, 2, 5, 10, 25",
        "original intensity", "Feature Extraction",
    ):
        assert phrase.lower() in result.lower()
    assert "present in the current live apoc controls" in result.lower()
    for phrase in (
        "Show classifier statistics", "greener, more informative",
        "redder, less informative", "broader candidate set",
    ):
        assert phrase.lower() in result.lower()


def test_apoc_feature_grid_follow_up_recommends_organoid_preset():
    import app

    result = app.apoc_feature_grid_guidance(
        {
            "current_step": "segmentation",
            "segmentation": {"method": "APOC (GPU)"},
            "metadata": {"cell_types": {"organoid": ["MDO"]}},
            "ui_state": {"controls": [
                {
                    "id": "segmentation.apoc.MDO.feature_scales",
                    "cell_type": "MDO",
                    "value": "1, 2, 5",
                },
                {
                    "id": "segmentation.apoc.MDO.feature_filters",
                    "cell_type": "MDO",
                    "value": ["Gaussian blur at sigma 1 px"],
                },
            ]},
        },
        [
            {"role": "user", "content": "Explain how to tune the APOC features."},
            {"role": "assistant", "content": "The Tune Features panel controls filters."},
            {"role": "user", "content": "What do you recommend for MDO?"},
        ],
    )

    for phrase in (
        "Large structures", "1, 2, 5, 10, and 25 pixels", "Retrain",
        "probability-map preview",
    ):
        assert phrase.lower() in result.lower()


def test_apoc_feature_preset_action_targets_segmentation_organoids():
    import app

    controls = []
    for cell_type in ("27T", "MDO"):
        controls.extend([
            {
                "id": f"segmentation.apoc.{cell_type}.feature_preset",
                "cell_type": cell_type,
                "value": "Medium structures",
                "visible": True,
                "enabled": True,
            },
            {
                "id": f"segmentation.apoc.{cell_type}.show_feature_tuning",
                "cell_type": cell_type,
                "value": False,
                "visible": True,
                "enabled": True,
            },
        ])
    result = app.apoc_feature_preset_action(
        {
            "current_step": "segmentation",
            "segmentation": {"method": "APOC (GPU)"},
            "metadata": {"cell_types": {
                "organoid": ["27T", "MDO"], "immune": ["tcell"],
            }},
            "ui_state": {"controls": controls},
        },
        [{"role": "user", "content": (
            "Fill in the correct features for the MDO and the 27T for me."
        )}],
    )
    changed = {
        call["arguments"]["control_id"]: call["arguments"]["value"]
        for call in result["calls"]
    }
    assert changed["segmentation.apoc.27T.feature_preset"] == "Large structures"
    assert changed["segmentation.apoc.MDO.feature_preset"] == "Large structures"
    assert changed["segmentation.apoc.27T.show_feature_tuning"] is True
    assert changed["segmentation.apoc.MDO.show_feature_tuning"] is True
    assert "Segmentation > APOC > Tune Features" in result["text"]
    assert "Feature Extraction or instance post-processing" in result["text"]


def test_apoc_feature_preset_action_honors_explicit_small_preset():
    import app

    result = app.apoc_feature_preset_action(
        {
            "current_step": "segmentation",
            "active_cell_type": "tcell",
            "segmentation": {"method": "APOC (GPU)"},
            "metadata": {"cell_types": {}},
            "ui_state": {"controls": [{
                "id": "segmentation.apoc.tcell.feature_preset",
                "cell_type": "tcell",
                "value": "Large structures",
            }]},
        },
        [{
            "role": "user",
            "content": (
                "These are individual T cells. Please tune the APOC features using "
                "the normal preset for Small structures."
            ),
        }],
    )

    changed = {
        call["arguments"]["control_id"]: call["arguments"]["value"]
        for call in result["calls"]
    }
    assert changed == {
        "segmentation.apoc.tcell.feature_preset": "Small structures",
    }
    assert "1, 2, and 5 pixels" in result["text"]


def test_btrack_step2_action_enables_then_calibrates_gap_controls():
    import app

    base_controls = [
        {
            "id": "tracking.tcell.btrack.use_global_optimization",
            "cell_type": "tcell",
            "value": False,
            "visible": True,
            "enabled": True,
        },
        {
            "id": "tracking.tcell.btrack.distance_threshold",
            "cell_type": "tcell",
            "value": 60,
            "unit": "um",
            "visible": False,
            "enabled": False,
        },
        {
            "id": "tracking.tcell.btrack.time_threshold",
            "cell_type": "tcell",
            "value": 3,
            "unit": "frames",
            "visible": False,
            "enabled": False,
        },
    ]
    enable = app.btrack_step2_action(
        {
            "current_step": "tracking",
            "active_cell_type": "tcell",
            "ui_state": {"controls": base_controls},
        },
        [{"role": "user", "content": "I mean Step 2 of tracking."}],
    )
    assert enable["calls"] == [{
        "name": "set_ui_value",
        "arguments": {
            "control_id": "tracking.tcell.btrack.use_global_optimization",
            "value": True,
        },
    }]
    assert "Global Hypothesis Optimizer" in enable["text"]
    assert "not the organoid Propagation" in enable["text"]
    assert "saved but inactive" in enable["text"]

    enabled_controls = [dict(control) for control in base_controls]
    enabled_controls[0]["value"] = True
    for control in enabled_controls[1:]:
        control["visible"] = True
        control["enabled"] = True
    enabled_controls.append({
        "id": "tracking.tcell.btrack.hypotheses",
        "cell_type": "tcell",
        "value": ["P_FP", "P_init", "P_term", "P_link"],
        "visible": True,
        "enabled": True,
    })
    tune = app.btrack_step2_action(
        {
            "current_step": "tracking",
            "active_cell_type": "tcell",
            "ui_state": {"controls": enabled_controls},
        },
        [{"role": "user", "content": (
            "For Step 2, bridge gaps up to 4 frames and 40 um."
        )}],
    )
    changed = {
        call["arguments"]["control_id"]: call["arguments"]["value"]
        for call in tune["calls"]
    }
    assert changed["tracking.tcell.btrack.distance_threshold"] == 40
    assert changed["tracking.tcell.btrack.time_threshold"] == 4
    assert "P_FP, P_init, P_term, P_link" in tune["text"]


def test_segmentation_minimum_size_uses_half_estimated_cell_volume():
    import app

    result = app.segmentation_minimum_size_action(
        {
            "current_step": "segmentation",
            "active_cell_type": "tcell",
            "metadata": {"records": [{
                "pixel_distance_xy": 1.0,
                "pixel_distance_z": 2.0,
                "distance_unit": "um",
            }]},
            "ui_state": {"controls": [{
                "id": "segmentation.apoc.tcell.minimum_size",
                "cell_type": "tcell",
                "value": 10,
                "unit": "voxels",
                "visible": True,
                "enabled": True,
            }]},
        },
        [{"role": "user", "content": (
            "My T cells are about 10 um across. Set the minimum size."
        )}],
    )
    assert "full object volume of about **524 µm³**" in result["text"]
    assert "start Minimum size at 50%" in result["text"]
    assert "**131 voxels**" in result["text"]
    assert result["calls"] == [{
        "name": "set_ui_value",
        "arguments": {
            "control_id": "segmentation.apoc.tcell.minimum_size",
            "value": 131,
        },
    }]


def test_edt_direction_guidance_is_strategy_specific_and_mentions_fallback():
    import app

    base = {
        "current_step": "segmentation",
        "active_cell_type": "13T",
        "segmentation": {"apoc": {"cell_type_strategies": {
            "13T": "APOC Mask + EDT/Watershed Resegmentation",
        }}},
    }
    mask_edt = app.edt_direction_guidance(
        base,
        [{"role": "user", "content": (
            "Does higher EDT mean more splitting? Explain the direction."
        )}],
    )
    assert "higher edt generally means more splitting" in mask_edt.lower()
    assert "raise for more splitting, lower for less" in mask_edt.lower()
    assert "falls back to the original unsplit component" in mask_edt.lower()

    peak = app.edt_direction_guidance(
        {
            **base,
            "segmentation": {"apoc": {"cell_type_strategies": {
                "13T": "APOC Peak EDT/Watershed",
            }}},
        },
        [{"role": "user", "content": "Which EDT direction gives more splitting?"}],
    )
    assert "lower edt generally means more splitting" in peak.lower()
    assert "minimum peak-height filter" in peak.lower()


def test_contact_distance_guidance_uses_xy_pixel_not_diagonal():
    import app

    result = app.feature_threshold_guidance(
        {
            "current_step": "feature_extraction",
            "active_cell_type": "tcell",
            "metadata": {"records": [{"pixel_distance_xy": 1.01}]},
            "ui_state": {"controls": [{
                "id": "features.tcell.contact_distance",
                "cell_type": "tcell",
                "value": 1.01,
            }]},
        },
        [{"role": "user", "content": (
            "How can I set the contact and dead-mask percentage threshold correctly?"
        )}],
    )
    assert "0 µm means strict mask touching" in result
    assert "one xy pixel" in result.lower()
    assert "not a voxel diagonal" in result.lower()
    assert "one-pixel xy gap" in result.lower()
    assert "green is below" in result.lower()
    assert "universal numeric range" in result.lower()


def test_dead_threshold_guidance_prioritizes_the_live_viewer_preview():
    import app

    result = app.feature_threshold_guidance(
        {
            "current_step": "feature_extraction",
            "active_cell_type": "tcell",
            "results": [{
                "id": "analysis/tcell/BEHAV3D_dead_dye_distribution.pdf",
                "label": "Dead dye distribution",
                "viewable": True,
            }],
        },
        [{"role": "user", "content": (
            "How should I set the dead-mask percentage threshold for the first time?"
        )}],
    )
    text = result.lower()
    assert "preview dead threshold in viewer" in text
    assert "green is below the threshold" in text
    assert "red is above it" in text
    assert "hovering" in text
    assert "universal numeric range" in text
    assert "let me open" not in text
    assert "result pdf" in text


def test_result_opening_correction_stops_repeated_opening_claims():
    import app

    result = app.result_opening_correction(
        {"current_step": "feature_extraction"},
        [
            {"role": "user", "content": "How do I set the dead threshold?"},
            {"role": "assistant", "content": "Let me open the result PDF."},
            {"role": "user", "content": "I think you cannot open it."},
        ],
    )
    text = result.lower()
    assert "no result was opened" in text
    assert "listed as viewable does not mean it has opened" in text
    assert "preview dead threshold in viewer" in text
    assert "let me" not in text


def test_interface_capabilities_report_current_exposure_without_cross_module_mapping():
    capabilities = _interface_capabilities([
        {"id": "segmentation.apoc.tcell.feature_scales"},
        {"id": "segmentation.apoc.tcell.feature_filters"},
    ])
    metadata = capabilities["metadata_builder"]
    apoc = capabilities["segmentation"]["apoc"]
    assert "raw-channel-to-cell-type mapping" in metadata["does_not_own"]
    assert apoc["custom_feature_scales_currently_exposed"] is True
    assert apoc["custom_feature_filters_currently_exposed"] is True
    assert "Sobel-of-Gaussian" in apoc["feature_filters"]


def test_safety_profile_summary_preserves_reference_without_claiming_execution():
    import app

    summary = app.safety_profile_summary(
        {
            "current_step": "analysis",
            "metadata": {"records": [{
                "time_interval": 2, "time_unit": "min",
            }]},
            "experiment_reference": {"notes": [{"text": (
                "Exp010 is multi-organoid safety profiling with 27T and MDO. "
                "Active killing is a contact-associated rise to 1.5 times baseline "
                "within 5 frames after one frame of contact."
            )}]},
            "results": [],
        },
        [{"role": "user", "content": (
            "How should I frame the safety comparison, and what exactly does "
            "active killing mean here?"
        )}],
    )
    assert "TEG cells are the immune effector" in summary
    assert "tumor 27T and healthy MDO" in summary
    assert "5 frames (10 minutes)" in summary
    assert "not a completed analysis" in summary
    assert "independent analysis for each target plus a combined analysis" in summary
    assert "configured" not in summary.lower()


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
            "recommend_edt", "summarize_track_counts", "save_metadata",
            "load_metadata", "open_analysis_view"} <= names
    assert "set_parameter" not in names
    # the backend text-fallback parser must also recognise the new tool names
    assert "bulk_fill_metadata" in app._TOOL_NAMES
    assert "set_ui_value" in app._TOOL_NAMES
    assert "recommend_edt" in app._TOOL_NAMES
    assert "summarize_track_counts" in app._TOOL_NAMES
    assert "save_metadata" in app._TOOL_NAMES
    assert "load_metadata" in app._TOOL_NAMES
    assert "open_analysis_view" in app._TOOL_NAMES


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


def test_build_actions_for_metadata_persistence_and_analysis_view():
    actions = build_actions([
        {"name": "save_metadata", "arguments": {}},
        {"name": "load_metadata", "arguments": {}},
        {"name": "open_analysis_view", "arguments": {"view": "death_dynamics"}},
        {"name": "open_analysis_view", "arguments": {"view": "unknown"}},
    ], [], {})
    assert actions[0].ok and actions[0].kind == "save_metadata"
    assert actions[1].ok and actions[1].kind == "load_metadata"
    assert actions[2].ok and actions[2].data["view"] == "death_dynamics"
    assert not actions[3].ok


def test_metadata_persistence_actions_call_existing_ui_handlers():
    from types import SimpleNamespace

    invoked = []
    main = SimpleNamespace(data_prep_tab=SimpleNamespace(
        _on_save_metadata=lambda: invoked.append("save") or True,
        _on_load_metadata=lambda: invoked.append("load") or True,
    ))
    save = ProposedAction("save_metadata")
    load = ProposedAction("load_metadata")
    assert apply_action(main, save)
    assert apply_action(main, load)
    assert invoked == ["save", "load"]
    assert "saved and activated" in save.data["result_markdown"]
    assert "loading has started" in load.data["result_markdown"]


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


def test_experiment_reference_discovers_calcium_and_compacts_historical_settings(
    tmp_path=None,
):
    import tempfile
    from pathlib import Path

    root = Path(tmp_path or tempfile.mkdtemp())
    (root / "README_Calcium_example.md").write_text(
        "# Calcium experiment\nNear-static reporter cells.",
        encoding="utf-8",
    )
    (root / "behav3d_parameters.yml").write_text(
        "\n".join([
            "paths:",
            "  output_dir: /private/results",
            "pixel_classifier:",
            "  convpaint_strategy: ConvPaint Probability + Watershed",
            "  convpaint_fe_alias: vgg",
            "  tcell_edt_threshold: 2.5",
            "  tcell_segment_size_min: 12",
            "track_filtering:",
            "  tcell:",
            "    min_length_enabled: true",
            "    min_track_length: 15",
            "  active_killing:",
            "    observation_window: 7",
            "    killing_threshold_multiplier: 2",
            "state_classification:",
            "  tcell:",
            "    n_states: 3",
            "    feature_smoothing_window: 5",
            "    selected_features: [speed, displacement]",
            "track_classification:",
            "  tcell:",
            "    behavioral_trajectory_size: 15",
            "    n_clusters: 3",
            "    linkage: average",
            "cell_type_groups:",
            "  Tcells: [tcell, tcell_control]",
            "single_cell:",
            "  selected_cell_type: Tcells",
        ]),
        encoding="utf-8",
    )

    reference = _experiment_reference_context(str(root), {})
    assert reference["notes"][0]["source"] == "README_Calcium_example.md"
    settings = reference["saved_configurations"][0]["settings"]
    assert settings["segmentation"]["convpaint_fe_alias"] == "vgg"
    assert settings["segmentation"]["saved_instance_by_cell_type"]["tcell"] == {
        "edt_threshold": 2.5,
        "segment_size_min": 12,
    }
    assert settings["active_killing"] == {
        "observation_window": 7,
        "killing_threshold_multiplier": 2,
    }
    assert settings["state_classification"]["by_cell_type"]["tcell"] == {
        "n_states": 3,
        "selected_features": ["speed", "displacement"],
        "feature_smoothing_window": 5,
    }
    assert settings["track_classification"]["by_cell_type"]["tcell"] == {
        "behavioral_trajectory_size": 15,
        "n_clusters": 3,
        "linkage": "average",
    }
    assert settings["cell_type_groups"] == {
        "Tcells": ["tcell", "tcell_control"],
    }
    assert settings["single_cell"]["selected_cell_type"] == "Tcells"
    assert "live metadata" in reference["source_policy"]
    assert "/private/results" not in str(reference)


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


def test_historical_reference_profiles_are_selected_and_never_auto_applied():
    from guidance import select_guidance_cards
    import app

    cards = select_guidance_cards(
        {"current_step": "tracking"},
        "Do you have values from a previous experiment for btrack?",
    )
    profile = next(item for item in cards if item["id"] == "reference_examples")
    for phrase in (
        "never defaults", "IVM HIV", "CD4/CD8-13T",
        "never edit a live control from an example alone",
    ):
        assert phrase.lower() in profile["text"].lower()

    prompt = app.build_system_prompt(
        {
            "current_step": "tracking",
            "metadata": {"loaded": True, "records": []},
        },
        cards,
        [{"name": "set_ui_value"}],
    )
    assert "HISTORICAL REFERENCE PROFILES are examples" in prompt
    assert "Never issue a form action from a historical value alone" in prompt
    assert "README notes for study intent" in prompt


def test_historical_reference_guidance_preserves_provenance_and_calibration():
    import app

    btrack = app.historical_reference_guidance(
        {
            "metadata": {"records": [{
                "time_interval": 30,
                "time_unit": "s",
            }]},
        },
        [{
            "role": "user",
            "content": (
                "Do you have example values from a previous experiment for btrack "
                "on T cells? Can I copy them?"
            ),
        }],
    )
    for phrase in (
        "IVM HIV", "CD4/CD8-13T", "should not be copied directly",
        "one-frame displacement", "largest spatial gap",
        "largest missing-frame gap",
    ):
        assert phrase.lower() in btrack["text"].lower()
    assert btrack["calls"] == []

    calcium = app.historical_reference_guidance(
        {},
        [{
            "role": "user",
            "content": (
                "Was there a previous experiment with near-static calcium reporter "
                "cells? What method and example values did it use?"
            ),
        }],
    )
    for phrase in (
        "Cellpose-SAM", "imported", "Reporter Propagation", "100-voxel",
        "10% overlap", "example values, not defaults",
    ):
        assert phrase.lower() in calcium["text"].lower()
    assert calcium["calls"] == []

    microglia = app.historical_reference_guidance(
        {},
        [{
            "role": "user",
            "content": (
                "What settings and design were used in the previous microglia "
                "Exp91 experiment?"
            ),
        }],
    )
    for phrase in (
        "Exp91", "eight wells", "1.77 µm", "120-second", "None_None",
        "APOC Probability Map + Watershed", "maximum search radius 150",
        "absolute increase of 30", "Start offset **0**", "README",
        "n=1", "historical values, not defaults",
    ):
        assert phrase.lower() in microglia["text"].lower()
    assert microglia["calls"] == []


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
    assert "btrack Step 2 is the Global Hypothesis Optimizer" in sp
    assert "offer to enable Step 2 for a complete setup" in sp
    assert "do not present disabled defaults as recommendations" in sp
    assert "Filtering must be run even when all filters are disabled" in sp
    assert "keep Start offset at 1" in sp
    assert "Average linkage is the default" in sp
    assert "Never call that combination contradictory" in sp
    assert "Do not call a chosen minimum reasonable" in sp
    assert "Never use 'already selected' as evidence" in sp
    assert "Bleed-through is not recorded in metadata" in sp
    assert "Fluorophore names such as GFP/RFP are not needed" in sp
    assert "APOC has separate Image Channel Inputs" in sp
    assert "does not mean APOC uses all channels" in sp
    assert "'Tune Features' means the classifier feature preset" in sp
    assert "ask them to copy and paste the latest error lines" in sp
    assert "Intensity and Contact are required for every cell type" in sp
    assert "never suggest dropping it from a T-cell population" in sp
    assert "automatically produces an independent analysis" in sp
    assert "Never apply only the mode checkbox" in sp
    assert "setup_ready true" in sp
    assert "enumerate all movement choices" in sp
    assert "currently selected speed and net displacement" in sp
    assert "read analysis.selected_cell_type" in sp
    assert "assigning the same name to multiple primary clusters merges them" in sp
    assert "Death Dynamics, Interaction Analysis, Invasiveness Analysis" in sp
    assert "line to the literal CSV-safe value 'not_added'" in sp
    assert "Never narrate internal rules" in sp
    assert "Opening a result requires an open_result call" in sp
    assert "A result merely listed as viewable has not been opened" in sp


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
    def __init__(self, checked=False, enabled=True):
        self._checked = checked
        self._enabled = enabled
    def isChecked(self): return self._checked
    def setChecked(self, value): self._checked = bool(value)
    def isEnabled(self): return self._enabled
    def isHidden(self): return False


class _FakeNamedCheck(_FakeCheck):
    def __init__(self, text, checked=False):
        super().__init__(checked)
        self._text = text
    def text(self): return self._text


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
    def setCurrentText(self, text):
        if text in self.items:
            self.index = self.items.index(text)
    def isEnabled(self): return True
    def isHidden(self): return False


class _FakeLog:
    def __init__(self, text): self._text = text
    def toPlainText(self): return self._text


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
    assert by_id["tracking.tcell.btrack.distance_threshold"]["visible"] is False
    assert by_id["tracking.tcell.btrack.hypotheses"]["visible"] is False
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


def test_metadata_validation_matches_required_well_and_population_lines():
    records = [{
        "sample_name": "Sample01",
        "exp_nr": 91,
        "well": "",
        "raw_image_path": "/missing/sample.zarr",
        "pixel_distance_xy": 1.77,
        "pixel_distance_z": 1.77,
        "distance_unit": "um",
        "time_interval": 120,
        "time_unit": "s",
        "cell_types": {
            "organoid": {"line": "DO7", "condition": ""},
            "T-cells": {"line": "", "condition": ""},
            "Macrophages": {"line": "not_added", "condition": ""},
        },
    }]
    issues = validate_metadata_records(records)
    errors = [item["message"] for item in issues if item["severity"] == "error"]
    assert any("mandatory well" in message for message in errors)
    assert any("mandatory T-cells line" in message for message in errors)
    assert not any("condition" in message for message in errors)
    assert not any("Macrophages line" in message for message in errors)


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


def test_unchanged_loaded_metadata_copy_does_not_require_save():
    import pandas as pd
    from types import SimpleNamespace

    form = {
        "group": _FakeGroup(),
        "basic": {
            "sample_name": _FakeLine("Sample01"),
            "raw_image_path": _FakeLine("/missing/sample.zarr"),
            "dimension_order": _FakeLine("TCZYX"),
            "pixel_distance_xy": _FakeSpin(0.5),
            "pixel_distance_z": _FakeSpin(2.0),
            "time_interval": _FakeSpin(5.0),
            "time_unit": _FakeCombo(["s", "m", "h"], 1),
        },
        "dead_channel": {}, "cell_types": {},
    }
    dp = SimpleNamespace(
        metadata=pd.DataFrame([{
            "sample_name": "Sample01",
            "raw_image_path": "/missing/sample.zarr",
            "dimension_order": "TCZYX",
            "pixel_distance_xy": 0.5,
            "pixel_distance_z": 2.0,
            "time_interval": 5.0,
            "time_unit": "m",
        }]),
        output_dir="", behav3d_parameters={},
        builder_grp=_FakeGroup(checked=True), _sample_forms=[form],
        n_samples_spin=_FakeSpin(1), n_organoid_spin=_FakeSpin(0),
        n_immune_spin=_FakeSpin(0), n_other_spin=_FakeSpin(0),
        include_dead_cb=_FakeCheck(False), _organoid_name_edits=[],
        _immune_name_edits=[], _other_name_edits=[],
        _metadata_builder_dirty=False,
    )
    state = _metadata_builder_state(dp)
    assert state["record_source"] == "loaded_metadata_copy"
    assert state["save_required"] is False
    context = build_context(SimpleNamespace(
        tabs=_FakeWorkflowTabs(2), data_prep_tab=dp,
    ))
    assert context["metadata"].get("record_source") != "metadata_builder_draft"
    assert context["metadata_builder"]["save_required"] is False


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


def test_active_killing_registry_keeps_dependent_threshold_editable():
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
    signal_id = "features.active_killing.death_signal"
    assert controls[target_id]["visible"] is True
    assert controls[target_id]["value"] == ["organoid1", "organoid2"]
    assert controls[signal_id]["value"] == "Dead-mask percentage"
    assert controls[signal_id]["choices"] == [
        "Dead-mask percentage", "Mean dead-dye intensity", "Dead-mask pixel count",
    ]
    assert apply_set_ui_value(main, signal_id, "Dead-mask pixel count")
    assert active.death_signal_combo.currentText() == "nr_dead_mask_pixels"
    assert apply_set_ui_value(main, target_id, ["organoid2"])
    threshold = controls["features.active_killing.absolute_threshold"]
    assert threshold["visible"] is True
    assert threshold["enabled"] is True
    assert threshold["active"] is False
    assert apply_set_ui_value(
        main, "features.active_killing.absolute_threshold", 30.0
    )
    assert active.spin_abs_threshold.value() == 30.0
    assert active_cell_type(main, "feature_extraction") == "tcell"
    assert [item.text() for item in active.target_list.selectedItems()] == ["organoid2"]


def test_feature_context_reports_active_killing_setup_readiness():
    from types import SimpleNamespace

    active = SimpleNamespace(
        immune_combo=_FakeCombo(["tcell"]),
        target_list=_FakeList(["27t", "mdo"], selected=["27t", "mdo"]),
        spin_obs_window=_FakeSpin(5),
        death_signal_combo=_FakeCombo([
            "percentage_dead_mask", "mean_dead_dye", "nr_dead_mask_pixels",
        ], 2),
        check_abs_threshold=_FakeCheck(True),
        spin_threshold_mult=_FakeSpin(1.5),
        spin_abs_threshold=_FakeSpin(0.0),
        spin_min_contact=_FakeSpin(1),
    )
    main = SimpleNamespace(feature_extraction_tab=SimpleNamespace(
        active_killing_panel=active,
        _ak_toggle_btn=_FakeCheck(True),
    ))
    state = _feature_extraction_state(main)["active_killing"]
    assert state["setup_ready"] is False
    assert "greater than 0" in state["setup_issues"][0]
    active.spin_abs_threshold.setValue(30.0)
    state = _feature_extraction_state(main)["active_killing"]
    assert state["setup_ready"] is True
    assert state["setup_issues"] == []


def test_feature_registry_marks_and_preserves_mandatory_groups():
    from types import SimpleNamespace

    checks = {
        "movement": _FakeCheck(True, enabled=False),
        "intensity": _FakeCheck(True, enabled=False),
        "morphology": _FakeCheck(True),
        "contact": _FakeCheck(True, enabled=False),
        "invasiveness": _FakeCheck(True),
        "death": _FakeCheck(True, enabled=False),
    }
    panel = SimpleNamespace(
        feature_checks=checks,
        contact_threshold=_FakeSpin(0.0),
        spin_dead_threshold=_FakeSpin(5.0),
        spin_workers=_FakeSpin(1),
        _persist=lambda: None,
    )
    main = SimpleNamespace(feature_extraction_tab=SimpleNamespace(
        panels={"tcell": panel},
        active_killing_panel=None,
    ))
    control = next(
        item for item in control_registry(main)
        if item["id"] == "features.tcell.feature_groups"
    )
    assert control["required_choices"] == [
        "movement", "intensity", "contact", "death",
    ]
    assert control["editable_choices"] == ["morphology", "invasiveness"]

    action = build_actions(
        [{
            "name": "set_ui_value",
            "arguments": {
                "control_id": control["id"],
                "value": ["movement", "contact", "death"],
            },
        }],
        cards=[],
        params={},
        controls=[control],
    )[0]
    assert action.ok is False
    assert "must keep:" in action.message
    assert "intensity" in action.message


def test_feature_group_guidance_keeps_intensity_for_dead_dye():
    import app

    text = app.feature_group_requirement_guidance(
        {
            "current_step": "feature_extraction",
            "active_cell_type": "tcell",
            "ui_state": {"controls": [{
                "id": "features.tcell.feature_groups",
                "visible": True,
                "enabled": True,
                "choices": [
                    "movement", "intensity", "morphology", "contact",
                    "invasiveness", "death",
                ],
                "required_choices": [
                    "movement", "intensity", "contact", "death",
                ],
                "value": [
                    "movement", "intensity", "morphology", "contact",
                    "invasiveness", "death",
                ],
            }]},
        },
        [{"role": "user", "content": "Now adjust T cells"}],
    )
    assert "Intensity" in text
    assert "mean dead-dye intensity" in text
    assert "will not suggest removing" in text
    assert "Morphology, Invasiveness" in text


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
        "pixel_distance_xy is 0.5; position_t is 10 for sample_name and TrackID. "
        "Use summarize_track_counts with nr_dead_mask_pixels, percentage_dead_mask, "
        "mean_dead_dye, and killing_efficiency."
    )
    assert text == (
        "XY pixel size is 0.5; timepoint is 10 for sample name and track ID. "
        "Use track-count preview with dead-mask pixel count, dead-mask percentage, "
        "mean dead-dye intensity, and killing efficiency."
    )


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


def test_apoc_registry_exposes_channel_and_feature_controls():
    from types import SimpleNamespace

    persisted = []
    channels = [
        _FakeNamedCheck("Channel 0", True),
        _FakeNamedCheck("Channel 1", True),
        _FakeNamedCheck("Channel 2", True),
    ]
    filters = {
        ("gaussian_blur", "1"): _FakeCheck(True),
        ("difference_of_gaussian", "1"): _FakeCheck(True),
        ("sobel_of_gaussian_blur", "2"): _FakeCheck(False),
    }
    tab = SimpleNamespace(
        channel_checkboxes=channels,
        feature_combo=_FakeCombo([
            "small_preset", "medium_preset", "large_preset", "custom",
        ], 2),
        tune_group=_FakeCheck(False),
        sigma_input=_FakeLine("1, 2, 5"),
        _feat_sigma_checks=filters,
        consider_original_cb=_FakeCheck(True),
        max_depth_spin=_FakeSpin(5),
        num_ensembles_spin=_FakeSpin(100),
        _per_tab_strategy_combo=_FakeCombo([
            "APOC Probability Map + Watershed",
        ]),
        prob_mask_threshold_spin=_FakeSpin(0.5, 0.0, 1.0),
        prob_seed_threshold_spin=_FakeSpin(0.8, 0.0, 1.0),
        _update_preview=lambda: None,
        _on_update_grid=lambda: None,
    )
    training = SimpleNamespace(tabs={"27t": tab}, strategy_combo=None)
    apoc = SimpleNamespace(
        spin_examples=_FakeSpin(3),
        _training_widget=training,
        _collect_apoc_tab_config=lambda: {"apoc_27t_channels": []},
        _save_apoc_params_to_yaml=lambda **kwargs: persisted.append(kwargs),
    )
    main = SimpleNamespace(segmentation_tab=SimpleNamespace(
        method_combo=_FakeCombo([
            "APOC (GPU)", "ConvPaint", "Pixel Classifier (Random Forest)",
            "Cellpose", "Cellpose-SAM (zero-shot)", "Import segmentation",
        ]),
        apoc_page=apoc,
        convpaint_page=SimpleNamespace(spin_examples=_FakeSpin(3)),
        pixel_classifier_page=SimpleNamespace(spin_examples=_FakeSpin(3)),
        cellpose_page=SimpleNamespace(),
    ))

    controls = {item["id"]: item for item in control_registry(main)}
    base = "segmentation.apoc.27t"
    assert controls[f"{base}.input_channels"]["choices"] == [
        "Channel 0", "Channel 1", "Channel 2",
    ]
    assert controls[f"{base}.feature_preset"]["value"] == "Large structures"
    assert controls[f"{base}.feature_preset"]["choices"] == [
        "Small structures", "Medium structures",
        "Large structures", "Custom feature selection",
    ]
    assert "Gaussian blur at sigma 1 px" in controls[
        f"{base}.feature_filters"
    ]["choices"]

    assert apply_set_ui_value(main, f"{base}.input_channels", [
        "Channel 1", "Channel 2",
    ])
    assert apply_set_ui_value(main, f"{base}.feature_preset", "Small structures")
    assert apply_set_ui_value(main, f"{base}.feature_filters", [
        "Gaussian blur at sigma 1 px",
        "Sobel edge at sigma 2 px",
    ])
    assert [check.isChecked() for check in channels] == [False, True, True]
    assert tab.feature_combo.currentText() == "custom"
    assert tab.tune_group.isChecked() is True
    assert filters[("gaussian_blur", "1")].isChecked() is True
    assert filters[("difference_of_gaussian", "1")].isChecked() is False
    assert filters[("sobel_of_gaussian_blur", "2")].isChecked() is True
    assert len(persisted) == 3


def test_segmentation_context_reports_apoc_training_state(tmp_path=None):
    import tempfile
    from pathlib import Path
    from types import SimpleNamespace

    root = Path(tmp_path or tempfile.mkdtemp())
    classifier = root / "PixelClassifier_27t.cl"
    classifier.write_text("classifier", encoding="utf-8")
    tab = SimpleNamespace(
        channel_checkboxes=[
            _FakeNamedCheck("Channel 0", False),
            _FakeNamedCheck("Channel 1", True),
            _FakeNamedCheck("Channel 2", True),
        ],
        feature_combo=_FakeCombo([
            "small_preset", "medium_preset", "large_preset", "custom",
        ], 2),
        tune_group=_FakeCheck(False),
        sigma_input=_FakeLine("1, 2, 5, 10, 25"),
        _feat_sigma_checks={
            ("gaussian_blur", "1"): _FakeCheck(True),
        },
        consider_original_cb=_FakeCheck(True),
        _per_tab_strategy_combo=_FakeCombo([
            "APOC Probability Map + Watershed",
        ]),
        _get_clf_path=lambda: classifier,
    )
    training = SimpleNamespace(tabs={"27t": tab}, strategy_combo=None)
    main = SimpleNamespace(segmentation_tab=SimpleNamespace(
        method_combo=_FakeCombo([
            "APOC (GPU)", "ConvPaint", "Pixel Classifier (Random Forest)",
            "Cellpose", "Cellpose-SAM (zero-shot)", "Import segmentation",
        ]),
        apoc_page=SimpleNamespace(
            _training_widget=training, _is_session_active=True,
        ),
    ))
    state = _segmentation_state(main)
    apoc = state["apoc"]
    assert apoc["training_data_loaded"] is True
    assert apoc["cell_types"]["27t"]["available_input_channels"] == [
        "Channel 0", "Channel 1", "Channel 2",
    ]
    assert apoc["cell_types"]["27t"]["selected_input_channels"] == [
        "Channel 1", "Channel 2",
    ]
    assert apoc["cell_types"]["27t"]["trained_classifier_found"] is True


def test_image_dimensions_and_current_log_are_serialized():
    from types import SimpleNamespace

    class _Item:
        def __init__(self, value): self.value = value
        def text(self): return self.value

    class _Table:
        def rowCount(self): return 1
        def item(self, row, column):
            return _Item(["Sample01", "(20, 4, 8, 256, 256)"][column])
        def cellWidget(self, row, column):
            return _FakeCombo(["TCZYX"])

    dimensions = _image_dimensions_state(SimpleNamespace(dim_table=_Table()))
    assert dimensions == [{
        "sample_name": "Sample01",
        "shape": "(20, 4, 8, 256, 256)",
        "dimension_order": "TCZYX",
        "channel_count": 4,
        "timepoint_count": 20,
    }]
    log_state = _current_log_state(
        SimpleNamespace(segmentation_tab=SimpleNamespace(log=_FakeLog(
            "[13:24:25] Loading training data...\n"
            "[13:24:26] Running _load_training_images...\n"
            "[13:24:27] Error: image could not be opened"
        ))),
        "segmentation",
    )
    assert log_state["has_explicit_error"] is True
    assert log_state["errors"] == [
        "[13:24:27] Error: image could not be opened",
    ]


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
    assert KNOWLEDGE_VERSION == "2026.07.27.28"
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
    assert CONTROL_CONTRACT_VERSION == "3.2"


# --------------------------------------------------------------------------
if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    passed = 0
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
        passed += 1
    print(f"\n{passed}/{len(fns)} tests passed")
