"""
BEHAV3D assistant — Modal service (thin DeepSeek proxy + RAG + tool calling).

Modal is a CPU-only middle layer: it authenticates the napari client, retrieves
BEHAV3D knowledge (RAG), then calls the **DeepSeek API** with native function
calling and relays the stream back. The DeepSeek key never leaves the server (it
lives in a Modal secret); clients only ever see the Modal URL + a shared token.

Deploy:   modal deploy chatbot/app.py
Dev:      modal serve chatbot/app.py        # hot-reloading local proxy
Ingest:   modal run chatbot/app.py::ingest   # (re)build the RAG index

Endpoints (FastAPI, behind a Bearer token):
  GET  /health                  -> {"ok": true}
  POST /chat                    -> text/event-stream of:
        {"type":"token","text":...}        streamed assistant text
        {"type":"tool_calls","calls":[...]}  proposed actions (client confirms)
        {"type":"done"} | {"type":"error","message":...}

Tool calls use DeepSeek's native function-calling; `set_parameter.key` is
constrained to the real BEHAV3D parameter keys (enum) so the model can't invent
them. `parse_tool_calls` is retained only as a fallback for a call accidentally
embedded in the content text.

The pure helpers (:func:`build_system_prompt`, :func:`to_openai_tools`,
:func:`assemble_tool_calls`, :func:`parse_tool_calls`) are unit-tested without
Modal/network.

NOTE: do **not** add ``from __future__ import annotations`` here. The FastAPI
``/chat`` route is defined inside the nested ``web()`` function and relies on the
real ``Request`` type annotation being resolvable; stringized annotations make
FastAPI mistake ``request`` for a query parameter and return HTTP 422.
"""
import json
import re

# ===========================================================================
# Pure logic (no Modal / network) — unit-testable
# ===========================================================================
_TOOLCALL_RE = re.compile(r"<TOOLCALL>\s*(\{.*?\})\s*</TOOLCALL>", re.DOTALL)

_TOOL_NAMES = ("set_parameter", "navigate_to_step", "add_queue_step", "fill_metadata_builder")


_CATEGORY_LABELS = {
    "immune": "immune cells",
    "organoid": "organoids",
    "other": "other cells",
    "all_organoids": "all organoids",
}

_LEAF_LABELS = {
    "metadata_csv": "metadata CSV",
    "output_dir": "output directory",
    "default_apply_all": "dimension order",
    "examples_per_sample": "examples per sample",
    "sample_specific_classifier": "sample-specific classifier",
    "workers": "worker count",
    "use_all_timepoints": "use all timepoints",
    "tp_start": "first timepoint",
    "tp_end": "last timepoint",
    "number_of_channels": "number of channels",
    "labels_mode": "channel-label setup",
    "channel_labels": "channel labels",
    "method": "method",
    "overwrite": "overwrite existing outputs",
    "search_range_px": "search range",
    "memory_frames": "memory frames",
    "track_cost_px": "frame-to-frame linking distance",
    "gap_close_cost_px": "gap-closing distance",
    "gap_close_max_frames": "maximum gap length",
    "max_search_radius": "maximum search radius",
    "features_choice": "feature groups",
    "contact_threshold": "contact distance",
    "n_workers": "worker count",
    "dead_perc_threshold": "dead-pixel threshold",
    "exp_duration": "experiment duration",
    "min_track_length": "minimum track length",
    "max_track_length": "maximum track length",
    "umap_min_dist": "UMAP minimum distance",
    "umap_n_neighbors": "UMAP neighbors",
    "nr_of_clusters": "number of clusters",
}


def label_for_key(key) -> str:
    """Small local copy of the client-side labeler; Modal only has JSON cards."""
    if not key:
        return "Parameter"
    parts = str(key).split(".")
    leaf = parts[-1]
    label = _LEAF_LABELS.get(leaf, leaf.replace("_", " "))
    category = next((p for p in parts if p in _CATEGORY_LABELS), None)
    method = next((p for p in parts if p in {"lap", "trackpy", "btrack"}), None)
    top = parts[0] if parts else ""
    prefixes = []
    if category:
        prefixes.append(_CATEGORY_LABELS[category])
    if method and leaf != "method":
        prefixes.append({"lap": "LAP", "trackpy": "TrackPy", "btrack": "btrack"}[method])
    elif top == "pixel_classifier":
        prefixes.append("pixel classifier")
    elif top == "cellpose":
        prefixes.append("Cellpose")
    elif top == "features":
        prefixes.append("feature extraction")
    elif top == "track_filtering":
        prefixes.append("track filtering")
    elif top == "analysis":
        prefixes.append("analysis")
    return ": ".join(prefixes + [label]) if prefixes else label.capitalize()


def build_system_prompt(context: dict, retrieved: list[dict], tools: list[dict]) -> str:
    """Assemble the system prompt from persona + context + retrieved docs."""
    parts = [
        "You are the BEHAV3D co-pilot, embedded in a napari plugin for analysing "
        "cell behaviour in 3D fluorescent imaging. Be concise and concrete. Ground "
        "every recommendation in the BEHAV3D KNOWLEDGE below and the user's CONTEXT. "
        "If you are unsure, ask one clarifying question about their data.\n\n"
        "FORMATTING (responses render as markdown in a narrow side panel): lead "
        "with a one-sentence answer; use short paragraphs separated by blank lines; "
        "use bullet/numbered lists for any set of items or steps; bold key terms "
        "and researcher-facing labels; use `code` only for literal values; avoid "
        "dense walls of text.\n\n"
        "USER-FACING LANGUAGE RULE: Internal names such as `n_samples`, "
        "`metadata.loaded`, `tracking.immune.btrack.max_search_radius`, and tool "
        "names are for tool calls only. In normal text, use visible labels such as "
        "\"Number of samples\", \"metadata loaded\", and \"immune cells: maximum "
        "search radius\". Do not mention Python variables, dotted keys, JSON, or "
        "tool-call syntax unless the user explicitly asks for implementation details.\n\n"
        "TOOLS — always call a tool instead of asking the user to click manually:\n"
        "- `fill_metadata_builder`: guide through the Metadata Builder step-by-step. "
        "  Use whenever the user asks to build metadata or be walked through setup.\n"
        "- `set_parameter`: propose a parameter change using a dotted key from "
        "  PARAMETERS FOR THIS STEP. NEVER use set_parameter for the segmentation "
        "  method selector — tell the user to click the Method dropdown instead.\n"
        "- `navigate_to_step`: switch the active tab. Use proactively when a "
        "  prerequisite is missing (e.g. metadata not loaded → navigate to data_preparation).\n"
        "- `add_queue_step`: add a processing step to the queue.\n\n"
        "TOOL-CALL RULE (CRITICAL): When the user gives you a value during a guided flow, "
        "you MUST emit the fill_metadata_builder or set_parameter tool call in THAT SAME "
        "response — never acknowledge an answer in text only. "
        "Pattern: user answers → you call the tool → you ask the next question. "
        "If you skip the tool call, the UI will not update.\n\n"
        "ONE-QUESTION RULE: Ask exactly ONE question per turn. After the tool call, "
        "ask the NEXT question. Never bundle multiple questions.\n\n"
        "PREREQ RULE: If the current step or the step the user is asking about has "
        "a blocker, address ONLY that blocker first. Call navigate_to_step to the "
        "required tab when the user is on the wrong tab and explain why.\n\n"
        "The app may fill empty fields immediately and asks before overwriting "
        "existing values. In text, say what you are filling or asking next.",
    ]

    step = context.get("current_step", "?")
    md = context.get("metadata", {})
    assistant_session = context.get("assistant_session", {}) or {}
    confirmed_keys = set(assistant_session.get("confirmed_parameter_keys", []) or [])
    active_ct = context.get("active_cell_type_tab")
    ct_line = f"\n- Active cell-type sub-tab: {active_ct}" if active_ct else ""
    parts.append(
        "## CONTEXT\n"
        f"- Current step/tab: {context.get('current_tab_label', step)}\n"
        f"- Output dir set: {context.get('output_dir_set')}\n"
        f"- Samples loaded: {md.get('n_samples', 0) if md.get('loaded') else 0}\n"
        f"- Cell types: {md.get('cell_types', {})}\n"
        f"- Non-default parameters: {context.get('parameters', {})}\n"
        f"- Queue: {[s.get('type') for s in context.get('queue', [])]}"
        + ct_line
    )

    # Step readiness — show blockers so the model can redirect proactively.
    readiness = context.get("step_readiness", {})
    current_ready = readiness.get(step)
    if current_ready is not None:
        status = "ready" if current_ready.get("ready", True) else "not ready"
        blockers = ", ".join(current_ready.get("blockers", [])) or "none"
        parts.append(
            "## CURRENT STEP READINESS\n"
            f"- {step}: {status}; blockers: {blockers}"
        )
    downstream = {
        k: v for k, v in readiness.items()
        if k != step and not v.get("ready", True)
    }
    if downstream:
        lines = [
            f"- {k}: {', '.join(v.get('blockers', ['unknown blocker']))}"
            for k, v in downstream.items()
        ]
        parts.append(
            "## DOWNSTREAM BLOCKERS\n"
            "Use these only when the user asks about that later step:\n" + "\n".join(lines)
        )

    mb = context.get("metadata_builder", {})
    if mb:
        org = mb.get("organoid_names", [])
        imm = mb.get("immune_names", [])
        oth = mb.get("other_names", [])
        sample_forms = mb.get("sample_forms", []) or []
        sample_lines = []
        for sample in sample_forms[:3]:
            sample_lines.append(
                f"- {sample.get('label', 'Sample')}: basic={sample.get('basic', {})}; "
                f"cell_types={list((sample.get('cell_types') or {}).keys())}"
            )
        samples_text = "\n".join(sample_lines) if sample_lines else "- none"
        parts.append(
            "## METADATA BUILDER STATE\n"
            f"- Builder open: {mb.get('open')}\n"
            f"- n_samples={mb.get('n_samples')}, n_organoids={mb.get('n_organoids')}, "
            f"n_immune={mb.get('n_immune')}, n_other={mb.get('n_other')}, "
            f"include_dead={mb.get('include_dead')}\n"
            f"- Cell types configured: {mb.get('cell_types_configured')} "
            f"(organoids={org}, immune={imm}, other={oth})\n"
            f"- Sample forms created: {mb.get('sample_forms_created')} "
            f"(count={mb.get('sample_form_count', 0)})\n"
            f"- Visible sample forms (first 3):\n{samples_text}\n\n"
            "Use `fill_metadata_builder` tool calls to fill the form. "
            "In user-facing text, use visible labels such as Number of samples, "
            "Organoid types, Immune types, Image path, Pixel XY, Pixel Z, and "
            "Time interval. "
            "REQUIRED ORDER: open_builder → n_samples / n_organoids / n_immune / n_other / "
            "include_dead → configure_cell_types → organoid_name / immune_name / other_name "
            "(index 0, 1, …; use immune_multicolor and immune_multicolor_channels when needed) "
            "→ create_sample_forms → per-sample fields (sample_name, exp_nr, well, "
            "raw_image_path, pixel_distance_xy, pixel_distance_z, time_interval, time_unit, "
            "dead_channel_number, dead_mask_path). "
            "For per-sample cell-type rows, use cell_line, cell_condition, "
            "cell_segments_image_path, cell_tracks_image_path, and cell_tracks_csv_path "
            "with the exact visible cell_type name. "
            "After filling Sample 1, call fill_down to copy shared values to all others. "
            "Never skip configure_cell_types or create_sample_forms — they are required actions."
        )

    seg = context.get("segmentation", {}) or {}
    if seg:
        parts.append(
            "## SEGMENTATION UI STATE\n"
            f"- Current method: {seg.get('method')}\n"
            f"- Available methods: {seg.get('available_methods', [])}\n"
            "APOC (GPU), ConvPaint, Pixel Classifier (Random Forest), Cellpose, "
            "and Import segmentation are distinct choices. If the user says APOC, "
            "guide the APOC page; do not call it Pixel Classifier."
        )

    # Params still at default — surface the most actionable missing configs.
    at_default = context.get("required_params_at_default", [])
    if confirmed_keys:
        at_default = [p for p in at_default if p.get("key") not in confirmed_keys]
        parts.append(
            "## CONFIRMED VALUES\n"
            "The user has already accepted these current values; do not ask about "
            "them again unless they ask to change them:\n"
            + "\n".join(f"- {label_for_key(k)} [internal key: `{k}`]" for k in sorted(confirmed_keys))
        )
    if at_default:
        lines = [
            f"- {label_for_key(p['key'])} [internal key: `{p['key']}`] "
            f"(default {p['default']!r}): {p.get('description', '')}"
            for p in at_default[:6]
        ]
        parts.append(
            "## PARAMS STILL AT DEFAULT\n"
            "These parameters have not been changed from their defaults — ask about "
            "the most important one if guiding the user:\n" + "\n".join(lines)
        )

    schema = context.get("step_schema", [])
    if schema:
        lines = [f"- {label_for_key(c['key'])} [internal key: `{c['key']}`] "
                 f"(default {c['default']!r}"
                 + (f", choices {c['choices']}" if c.get("choices") else "") + ")"
                 for c in schema[:40]]
        parts.append(
            "## PARAMETERS FOR THIS STEP\n"
            "Use the internal key only inside set_parameter calls. Use the label in text.\n"
            + "\n".join(lines)
        )

    # Per-step guided flow instructions.
    _step_guides = {
        "data_preparation": (
            "Check whether the output directory is set first. If it is not set, ask for "
            "the output directory path; when the user answers, call set_parameter with "
            "internal key paths.output_dir. "
            "If metadata is not loaded, guide through the builder in this exact sequence:\n"
            "1. Call fill_metadata_builder(open_builder)\n"
            "2. Ask 'How many samples?' → user answers N → call fill_metadata_builder(n_samples=N)\n"
            "3. Ask 'How many organoid types?' → call fill_metadata_builder(n_organoids=N)\n"
            "4. Ask 'How many immune cell types?' → call fill_metadata_builder(n_immune=N)\n"
            "5. Ask 'How many other cell types?' → call fill_metadata_builder(n_other=N)\n"
            "6. Call fill_metadata_builder(configure_cell_types) (no user input needed)\n"
            "7. For each organoid type: ask its name → call fill_metadata_builder(organoid_name=X, index=i)\n"
            "8. Repeat for immune_name and other_name\n"
            "9. Call fill_metadata_builder(create_sample_forms) (no user input needed)\n"
            "10. For each sample: ask image path → call fill_metadata_builder(raw_image_path=P, index=i); "
            "ask pixel size → fill pixel_distance_xy/z; ask time interval → fill time_interval/time_unit\n"
            "11. Ask for line/condition only when the metadata needs it; fill those with "
            "cell_line/cell_condition and the visible cell_type name\n"
            "12. Call fill_metadata_builder(fill_down) to copy shared values when appropriate\n"
            "CRITICAL: steps that say 'user answers N → call ...' require a tool call in that same turn. "
            "Never skip it."
        ),
        "segmentation": (
            "Verify metadata is loaded and output_dir is set before anything else. "
            "Ask about the global Method dropdown first. The visible choices are APOC (GPU), "
            "ConvPaint (DL pixel classifier), Pixel Classifier (Random Forest), Cellpose "
            "(Deep Learning), and Import segmentation. APOC is not the same thing as "
            "Pixel Classifier (Random Forest). Do NOT call set_parameter for the Method "
            "dropdown — tell the user the exact visible option to select. Once the user "
            "accepts the current/default value for a setting, move to the next setting; "
            "do not ask about it again. For Cellpose: ask number_of_channels, then "
            "labels_mode. For Pixel Classifier: ask examples_per_sample, workers, "
            "use_all_timepoints. For APOC or ConvPaint: guide the user through generating "
            "training data, labeling, training, and then batch segmentation."
        ),
        "tracking": (
            "Read cell types from context. For each immune type: default=btrack, state it and ask "
            "to confirm, then ask max_search_radius. For organoid: default=propagation. "
            "For other: default=lap. Propose set_parameter calls for method and key params, "
            "one at a time. After each type is configured, move to the next."
        ),
        "feature_extraction": (
            "Check metadata.columns for *_tracks_image_path — missing means tracking is not done. "
            "Use features.<category>.contact_threshold for the contact distance (in µm). "
            "Use features.<category>.n_workers for worker count. "
            "Check metadata.columns for a dead_channel column to detect viability staining."
        ),
    }
    guide = _step_guides.get(step)
    if guide:
        parts.append(f"## GUIDED FLOW — {step.upper()}\n{guide}")

    if retrieved:
        docs = "\n\n".join(f"[{d.get('title','doc')}] {d.get('text','')}" for d in retrieved)
        parts.append("## BEHAV3D KNOWLEDGE (retrieved)\n" + docs)

    return "\n\n".join(parts)


def _balanced_json_spans(text: str):
    """Yield (start, end, obj) for every balanced, string-aware {...} that parses
    as JSON. Independent of the TOOLCALL markers, so it survives the malformed
    output small models emit (missing `<`, missing closing tag, nested objects)."""
    i, n = 0, len(text)
    while i < n:
        if text[i] != "{":
            i += 1
            continue
        depth, in_str, esc, j = 0, False, False, i
        while j < n:
            c = text[j]
            if in_str:
                if esc:
                    esc = False
                elif c == "\\":
                    esc = True
                elif c == '"':
                    in_str = False
            elif c == '"':
                in_str = True
            elif c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
                if depth == 0:
                    try:
                        yield i, j + 1, json.loads(text[i:j + 1])
                    except Exception:
                        pass
                    break
            j += 1
        i = j + 1


def parse_tool_calls(text: str) -> tuple[str, list[dict]]:
    """Extract tool calls from model output, tolerating malformed markers.

    Accepts the canonical ``<TOOLCALL>{...}</TOOLCALL>`` as well as the common
    small-model variants: ``TOOLCALL>{...}`` (missing `<`/closing tag) and bare
    ``set_parameter{...}`` where the object holds the arguments. Returns
    (clean_text_for_display, calls)."""
    calls: list[dict] = []
    for start, _end, obj in _balanced_json_spans(text):
        if not isinstance(obj, dict):
            continue
        name = obj.get("name")
        if name in _TOOL_NAMES:
            calls.append({"name": name, "arguments": obj.get("arguments", {}) or {}})
        elif any(k in obj for k in ("key", "step", "step_type", "field")):
            # bare arguments object — infer the tool from the preceding token
            pre = text[max(0, start - 48):start]
            m = re.search(
                r"(set_parameter|navigate_to_step|add_queue_step|fill_metadata_builder)"
                r"[^A-Za-z0-9_]*$",
                pre,
            )
            if m:
                calls.append({"name": m.group(1), "arguments": obj})
    return split_streamable(text), calls


def split_streamable(text: str) -> str:
    """Return only the human-visible prefix — text before any tool-call syntax,
    whether well-formed (`<TOOLCALL>`), malformed (`TOOLCALL>`), or bare
    (`set_parameter{`/`set_parameter(`)."""
    idxs = [text.find(mk) for mk in ("<TOOLCALL", "TOOLCALL>") if text.find(mk) != -1]
    m = re.search(
        r"(?:set_parameter|navigate_to_step|add_queue_step|fill_metadata_builder)\s*[\({]",
        text,
    )
    if m:
        idxs.append(m.start())
    return text[:min(idxs)].rstrip() if idxs else text.rstrip()


_QUEUE_STEP_TYPES = ["segment", "train", "track", "feature_extract", "filter", "active_killing"]


def to_openai_tools(tool_schema: list[dict], key_enum=None) -> list[dict]:
    """Wrap our TOOL_SCHEMA (name/description/parameters) into OpenAI `tools`
    format. If `key_enum` is given, constrain `set_parameter.key` to those values
    so the model cannot invent parameter keys."""
    out = []
    for t in tool_schema or []:
        params = json.loads(json.dumps(t.get("parameters", {})))  # deep copy
        if t.get("name") == "set_parameter" and key_enum:
            props = params.setdefault("properties", {})
            key_prop = props.setdefault("key", {"type": "string"})
            key_prop["enum"] = list(key_enum)
        if t.get("name") == "add_queue_step":
            props = params.setdefault("properties", {})
            st_prop = props.setdefault("step_type", {"type": "string"})
            st_prop["enum"] = _QUEUE_STEP_TYPES
        out.append({
            "type": "function",
            "function": {
                "name": t.get("name"),
                "description": t.get("description", ""),
                "parameters": params,
            },
        })
    return out


def assemble_tool_calls(fragments: list[dict]) -> list[dict]:
    """Merge streamed OpenAI `delta.tool_calls` fragments into our
    `[{name, arguments}]` shape. Each fragment carries an `index`, an optional
    `name`, and an `arguments` string delta; we concatenate by index and JSON-parse
    the accumulated arguments."""
    by_index: dict = {}
    order: list = []
    for frag in fragments or []:
        idx = frag.get("index", 0)
        if idx not in by_index:
            by_index[idx] = {"name": None, "args": ""}
            order.append(idx)
        if frag.get("name"):
            by_index[idx]["name"] = frag["name"]
        if frag.get("arguments"):
            by_index[idx]["args"] += frag["arguments"]
    calls = []
    for idx in order:
        entry = by_index[idx]
        if not entry["name"]:
            continue
        try:
            args = json.loads(entry["args"]) if entry["args"].strip() else {}
        except json.JSONDecodeError:
            continue
        calls.append({"name": entry["name"], "arguments": args})
    return calls


# ===========================================================================
# Modal app (imported lazily so the pure helpers import without modal installed)
# ===========================================================================
try:
    import modal
except Exception:  # allows unit-testing the pure helpers without modal
    modal = None

if modal is not None:
    import os as _os

    VOLUME_NAME = "behav3d-assistant-index"
    INDEX_DIR = "/index"
    DEEPSEEK_BASE_URL = "https://api.deepseek.com"

    # Source modules whose docstrings seed the knowledge base. Inlined here (NOT
    # imported from ingest) because Modal imports app.py inside *every* container,
    # including images that don't carry ingest.py — a top-level `import ingest`
    # there raises ModuleNotFoundError. Keep this list in sync with
    # ingest._DOC_PY_MODULES.
    _DOC_PY_MODULES = [
        "behav3d/preprocessing/segmentation/cellpose_prediction.py",
        "behav3d/preprocessing/segmentation/apoc_segment.py",
        "behav3d/preprocessing/tracking/btrack_tracking.py",
        "behav3d/preprocessing/tracking/laptracking.py",
        "behav3d/features/timepoint_features.py",
        "behav3d/analysis/interaction_analysis.py",
    ]

    # Single CPU image for ingest + web. No GPU, no vLLM, no model weights — Modal
    # is just a proxy that retrieves docs and calls the DeepSeek API.
    service_image = (
        modal.Image.debian_slim(python_version="3.12")
        .pip_install(
            "fastapi[standard]", "sentence-transformers==3.2.1",
            "numpy", "huggingface_hub", "openai>=1.40",
        )
        .add_local_file("chatbot/embeddings.py", "/root/embeddings.py")
        .add_local_file("chatbot/ingest.py", "/root/ingest.py")
        .add_local_file("chatbot/schema_cards.json", "/root/schema_cards.json")
    )
    # Knowledge sources — paths are relative to CWD, so run modal from the repo
    # root. README is required; WIKI.md is optional (indexed if present).
    for _md in ("README.md", "WIKI.md"):
        if _os.path.exists(_md):
            service_image = service_image.add_local_file(_md, f"/root/repo/{_md}")
    for _rel in _DOC_PY_MODULES:
        if _os.path.exists(_rel):
            service_image = service_image.add_local_file(_rel, f"/root/repo/{_rel}")

    app = modal.App("behav3d-assistant")
    volume = modal.Volume.from_name(VOLUME_NAME, create_if_missing=True)
    # Shared bearer token the napari client must present.
    auth_secret = modal.Secret.from_name("behav3d-assistant-auth")  # {"BEHAV3D_ASSISTANT_TOKEN": "..."}
    # DeepSeek API key — stays server-side, never sent to the client.
    deepseek_secret = modal.Secret.from_name("deepseek-api-key")    # {"DEEPSEEK_API_KEY": "..."}

    # -- Index build -------------------------------------------------------
    @app.function(image=service_image, volumes={INDEX_DIR: volume}, timeout=900)
    def ingest():
        import ingest as ing  # /root/ingest.py
        n = ing.build_and_save(repo_root="/root/repo",
                               schema_cards_path="/root/schema_cards.json",
                               out_dir=INDEX_DIR)
        volume.commit()
        print(f"Built RAG index with {n} chunks.")

    # -- Web endpoint (CPU) ------------------------------------------------
    # scaledown_window keeps a warm container 5 min after the last request. For a
    # zero-cold-start endpoint add `min_containers=1` — cheap on CPU (~$15-45/mo),
    # vs ~$0 idle with a ~few-second cold start.
    @app.function(image=service_image, volumes={INDEX_DIR: volume},
                  secrets=[auth_secret, deepseek_secret], scaledown_window=300,
                  min_containers=1)
    @modal.asgi_app()
    def web():
        import os
        from fastapi import FastAPI, Request, HTTPException
        from fastapi.responses import StreamingResponse
        from openai import OpenAI
        from embeddings import Embedder, VectorIndex

        api = FastAPI(title="BEHAV3D Assistant")
        embedder = Embedder()
        try:
            index = VectorIndex.load(INDEX_DIR)
        except Exception:
            index = VectorIndex()
        # Valid parameter keys → constrain set_parameter so the model can't invent.
        try:
            with open("/root/schema_cards.json", encoding="utf-8") as f:
                _key_enum = [c["key"] for c in json.load(f)
                             if not c["key"].startswith("calculated_features.")]
        except Exception:
            _key_enum = []

        expected_token = os.environ.get("BEHAV3D_ASSISTANT_TOKEN", "")
        # deepseek-v4-flash is the explicit V4 Flash id (tool-calls + streaming).
        # The older `deepseek-chat` alias also maps to it but is being deprecated.
        # `deepseek-v4-pro` is the stronger/pricier option. Override via DEEPSEEK_MODEL.
        deepseek_model = os.environ.get("DEEPSEEK_MODEL", "deepseek-v4-flash")
        client = OpenAI(base_url=DEEPSEEK_BASE_URL,
                        api_key=os.environ.get("DEEPSEEK_API_KEY", ""))

        def _auth(request: Request):
            if not expected_token:
                return
            if request.headers.get("authorization", "") != f"Bearer {expected_token}":
                raise HTTPException(status_code=401, detail="bad token")

        @api.get("/health")
        def health():
            return {"ok": True, "chunks": len(index.chunks), "model": deepseek_model}

        @api.post("/chat")
        async def chat(request: Request):
            _auth(request)
            body = await request.json()
            messages = body.get("messages", [])
            context = body.get("context", {})
            tools = body.get("tools", [])

            user_msg = next((m["content"] for m in reversed(messages)
                             if m.get("role") == "user"), "")
            query = f"{context.get('current_step','')} {user_msg}"
            try:
                retrieved = index.search(embedder.encode([query])[0], k=6)
            except Exception:
                retrieved = []

            system = build_system_prompt(context, retrieved, tools)
            convo = [{"role": "system", "content": system}]
            convo += [m for m in messages if m.get("role") != "system"]
            oai_tools = to_openai_tools(tools, key_enum=_key_enum)

            def sse(obj):
                return "data: " + json.dumps(obj) + "\n\n"

            def event_stream():
                content = ""
                tool_frags = []
                try:
                    stream = client.chat.completions.create(
                        model=deepseek_model, messages=convo,
                        tools=oai_tools or None,
                        tool_choice="auto" if oai_tools else None,
                        temperature=0.3, stream=True,
                    )
                    for chunk in stream:
                        delta = chunk.choices[0].delta
                        if getattr(delta, "content", None):
                            content += delta.content
                            yield sse({"type": "token", "text": delta.content})
                        for tc in (getattr(delta, "tool_calls", None) or []):
                            fn = getattr(tc, "function", None)
                            tool_frags.append({
                                "index": getattr(tc, "index", 0),
                                "name": getattr(fn, "name", None) if fn else None,
                                "arguments": getattr(fn, "arguments", None) if fn else None,
                            })
                except Exception as e:
                    yield sse({"type": "error", "message": f"DeepSeek API error: {e}"})
                    yield sse({"type": "done"})
                    return

                calls = assemble_tool_calls(tool_frags)
                if not calls:                       # fallback: call embedded in text
                    _, calls = parse_tool_calls(content)
                if calls:
                    yield sse({"type": "tool_calls", "calls": calls})
                yield sse({"type": "done"})

            return StreamingResponse(event_stream(), media_type="text/event-stream")

        return api
