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
        "and parameter names; use `code` for values/keys; avoid dense walls of text.\n\n"
        "To change a setting, CALL the set_parameter / navigate_to_step / "
        "add_queue_step tools (the user confirms before anything applies, so never "
        "claim a value is already changed — say you are proposing it). Only use "
        "parameter keys from PARAMETERS FOR THIS STEP.",
    ]

    step = context.get("current_step", "?")
    md = context.get("metadata", {})
    parts.append(
        "## CONTEXT\n"
        f"- Current step/tab: {context.get('current_tab_label', step)}\n"
        f"- Output dir set: {context.get('output_dir_set')}\n"
        f"- Samples loaded: {md.get('n_samples', 0) if md.get('loaded') else 0}\n"
        f"- Cell types: {md.get('cell_types', {})}\n"
        f"- Non-default parameters: {context.get('parameters', {})}\n"
        f"- Queue: {[s.get('type') for s in context.get('queue', [])]}"
    )

    schema = context.get("step_schema", [])
    if schema:
        lines = [f"- {c['key']} (default {c['default']!r}"
                 + (f", choices {c['choices']}" if c.get("choices") else "") + ")"
                 for c in schema[:40]]
        parts.append("## PARAMETERS FOR THIS STEP\n" + "\n".join(lines))

    if retrieved:
        docs = "\n\n".join(f"[{d.get('title','doc')}] {d.get('text','')}" for d in retrieved)
        parts.append("## BEHAV3D KNOWLEDGE (retrieved)\n" + docs)

    return "\n\n".join(parts)


_TOOL_NAMES = ("set_parameter", "navigate_to_step", "add_queue_step")


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
        elif any(k in obj for k in ("key", "step", "step_type")):
            # bare arguments object — infer the tool from the preceding token
            pre = text[max(0, start - 48):start]
            m = re.search(r"(set_parameter|navigate_to_step|add_queue_step)[^A-Za-z0-9_]*$", pre)
            if m:
                calls.append({"name": m.group(1), "arguments": obj})
    return split_streamable(text), calls


def split_streamable(text: str) -> str:
    """Return only the human-visible prefix — text before any tool-call syntax,
    whether well-formed (`<TOOLCALL>`), malformed (`TOOLCALL>`), or bare
    (`set_parameter{`/`set_parameter(`)."""
    idxs = [text.find(mk) for mk in ("<TOOLCALL", "TOOLCALL>") if text.find(mk) != -1]
    m = re.search(r"(?:set_parameter|navigate_to_step|add_queue_step)\s*[\({]", text)
    if m:
        idxs.append(m.start())
    return text[:min(idxs)].rstrip() if idxs else text.rstrip()


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
                  secrets=[auth_secret, deepseek_secret], scaledown_window=300)
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
