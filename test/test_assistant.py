"""
Tests for the BEHAV3D napari co-pilot assistant (pure-logic layer).

Runs under pytest *or* standalone:  python test/test_assistant.py
Requires the `behav3d` conda env (numpy, the behav3d package). No GPU / Modal /
Qt needed — only the dependency-free logic is exercised here.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "chatbot"))

from behav3d.napari._assistant_schema import (
    flatten_config_to_cards, cards_for_step, dump_cards_json,
)
from behav3d.napari._assistant_actions import (
    get_by_dotted, set_by_dotted, validate_value, build_actions,
)
from behav3d.napari._assistant_context import summarize_metadata, _diff_from_defaults


# --------------------------------------------------------------------------
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
    txt = ('try ~12px.\n<TOOLCALL>{"name":"set_parameter",'
           '"arguments":{"key":"cellpose.number_of_channels","value":2}}</TOOLCALL>')
    clean, calls = app.parse_tool_calls(txt)
    assert calls[0]["name"] == "set_parameter" and "<TOOLCALL>" not in clean
    assert app.split_streamable(txt) == "try ~12px."
    assert app.parse_tool_calls("x <TOOLCALL>{bad}</TOOLCALL>")[1] == []
    sp = app.build_system_prompt(
        {"current_step": "tracking", "step_schema": [], "parameters": {}, "metadata": {}, "queue": []},
        [{"title": "btrack", "text": "good for crowded cells"}], [])
    assert "BEHAV3D" in sp and "good for crowded cells" in sp
    sp2 = app.build_system_prompt(
        {"current_step": "segmentation", "step_schema": [],
         "required_params_at_default": [
             {"key": "pixel_classifier.examples_per_sample", "default": 3, "description": "desc"},
         ],
         "assistant_session": {"confirmed_parameter_keys": ["pixel_classifier.examples_per_sample"]},
         "parameters": {}, "metadata": {}, "queue": []},
        [], [])
    assert "CONFIRMED VALUES" in sp2
    assert "PARAMS STILL AT DEFAULT" not in sp2


def test_tool_call_parsing_tolerates_malformed_markers():
    import app
    # the exact malformed shape a small model emitted: no leading '<', no closing tag
    malformed = ('TOOLCALL>{"name": "set_parameter", "arguments": '
                 '{"key": "paths.metadata_csv", "value": "/p/m.csv"}} '
                 'TOOLCALL>{"name": "set_parameter", "arguments": '
                 '{"key": "dim_order.default_apply_all", "value": "TZCYX"}}')
    clean, calls = app.parse_tool_calls(malformed)
    assert len(calls) == 2 and "TOOLCALL" not in clean
    assert calls[0]["arguments"]["key"] == "paths.metadata_csv"
    # canonical form with nested args still parses and is stripped from display
    clean2, calls2 = app.parse_tool_calls(
        'Sure.\n<TOOLCALL>{"name":"set_parameter","arguments":'
        '{"key":"tracking.immune.method","value":"btrack"}}</TOOLCALL>')
    assert calls2[0]["arguments"]["value"] == "btrack" and clean2 == "Sure."
    # bare proposal form + list value
    _, calls3 = app.parse_tool_calls('set_parameter{"key":"x","value":[100,100,10]}')
    assert calls3[0]["arguments"]["value"] == [100, 100, 10]
    clean4, calls4 = app.parse_tool_calls(
        'Done.\nfill_metadata_builder{"field":"n_samples","value":3}')
    assert clean4 == "Done."
    assert calls4 == [{"name": "fill_metadata_builder",
                       "arguments": {"field": "n_samples", "value": 3}}]


def test_to_openai_tools_injects_key_enum():
    import app
    from behav3d.napari._assistant_actions import TOOL_SCHEMA
    enum = ["paths.metadata_csv", "tracking.immune.method"]
    tools = app.to_openai_tools(TOOL_SCHEMA, key_enum=enum)
    assert all(t["type"] == "function" for t in tools)
    sp = next(t for t in tools if t["function"]["name"] == "set_parameter")
    assert sp["function"]["parameters"]["properties"]["key"]["enum"] == enum
    # original schema is not mutated (deep-copied)
    assert "enum" not in TOOL_SCHEMA[0]["parameters"]["properties"]["key"]


def test_assemble_tool_calls_merges_streamed_fragments():
    import app
    frags = [
        {"index": 0, "name": "set_parameter", "arguments": ""},
        {"index": 0, "name": None, "arguments": '{"key": "tracking.immune'},
        {"index": 0, "name": None, "arguments": '.trackpy.search_range_px", "value": 40}'},
    ]
    calls = app.assemble_tool_calls(frags)
    assert calls == [{"name": "set_parameter",
                      "arguments": {"key": "tracking.immune.trackpy.search_range_px",
                                    "value": 40}}]
    # incomplete / unparseable arguments are dropped, not crashed on
    assert app.assemble_tool_calls([{"index": 0, "name": "set_parameter",
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


# --------------------------------------------------------------------------
if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    passed = 0
    for fn in fns:
        fn()
        print(f"PASS {fn.__name__}")
        passed += 1
    print(f"\n{passed}/{len(fns)} tests passed")
