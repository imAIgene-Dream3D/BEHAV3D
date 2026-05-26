"""
BEHAV3D assistant — Modal service (LLM + RAG + tool calling).

Deploy:   modal deploy chatbot/app.py
Dev:      modal serve chatbot/app.py        # hot-reloading local proxy
Ingest:   modal run chatbot/app.py::ingest   # (re)build the RAG index

Endpoints (FastAPI, behind a Bearer token):
  GET  /health                  -> {"ok": true}
  POST /chat                    -> text/event-stream of:
        {"type":"token","text":...}        streamed assistant text
        {"type":"tool_calls","calls":[...]}  proposed actions (client confirms)
        {"type":"done"}

The LLM is a self-hosted open instruct model on a GPU (vLLM). Tool calls use a
simple, model-agnostic ``<TOOLCALL>{json}</TOOLCALL>`` protocol that we parse out
of the generated text — robust across open models and trivial to stream.

The pure helpers (:func:`build_system_prompt`, :func:`parse_tool_calls`,
:func:`split_streamable`) are unit-tested without Modal/GPU.

NOTE: do **not** add ``from __future__ import annotations`` here. The FastAPI
``/chat`` route is defined inside the nested ``web()`` function and relies on the
real ``Request`` type annotation being resolvable; stringized annotations make
FastAPI mistake ``request`` for a query parameter and return HTTP 422.
"""
import json
import re

# ===========================================================================
# Pure logic (no Modal / GPU) — unit-testable
# ===========================================================================
_TOOLCALL_RE = re.compile(r"<TOOLCALL>\s*(\{.*?\})\s*</TOOLCALL>", re.DOTALL)

TOOL_PROTOCOL = (
    "When you want to PROPOSE a UI action, emit one line per action in the exact "
    "form <TOOLCALL>{\"name\": <tool>, \"arguments\": {...}}</TOOLCALL> at the end "
    "of your reply. The user must confirm before anything is applied, so never "
    "claim you have changed a value — say you are proposing it. Available tools: "
    "set_parameter{key,value}, navigate_to_step{step}, add_queue_step{step_type,params}."
)


def build_system_prompt(context: dict, retrieved: list[dict], tools: list[dict]) -> str:
    """Assemble the system prompt from persona + context + retrieved docs + tools."""
    parts = [
        "You are the BEHAV3D co-pilot, embedded in a napari plugin for analysing "
        "cell behaviour in 3D fluorescent imaging. Be concise and concrete. Ground "
        "every recommendation in the BEHAV3D KNOWLEDGE below and the user's CONTEXT. "
        "If you are unsure, ask one clarifying question about their data.\n\n"
        "FORMATTING (responses render as markdown in a narrow side panel): lead "
        "with a one-sentence answer; use short paragraphs separated by blank lines; "
        "use bullet/numbered lists for any set of items or steps; bold key terms "
        "and parameter names; use `code` for values/keys; avoid dense walls of text.",
        TOOL_PROTOCOL,
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


# ===========================================================================
# Modal app (imported lazily so the pure helpers import without modal installed)
# ===========================================================================
try:
    import modal
except Exception:  # allows unit-testing the pure helpers without modal
    modal = None

if modal is not None:
    MODEL_NAME = "Qwen/Qwen2.5-7B-Instruct"
    # Pin a commit for reproducibility (recommended). "main" tracks the latest.
    MODEL_REVISION = "main"
    VOLUME_NAME = "behav3d-assistant-index"
    INDEX_DIR = "/index"
    MODEL_CACHE_DIR = "/models"   # HF cache lives here, on a persistent Volume

    import os as _os

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

    # Light image for ingest + web (no vLLM → fast builds). It carries the local
    # python files and all knowledge sources.
    service_image = (
        modal.Image.debian_slim(python_version="3.12")
        .pip_install(
            "fastapi[standard]", "sentence-transformers==3.2.1",
            "numpy", "huggingface_hub",
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
    # Source modules whose docstrings seed the knowledge base.
    for _rel in _DOC_PY_MODULES:
        if _os.path.exists(_rel):
            service_image = service_image.add_local_file(_rel, f"/root/repo/{_rel}")

    # Heavy image only for the GPU inference engine. HF_HOME points at the mounted
    # weights Volume so vLLM reads cached weights instead of re-downloading from
    # HuggingFace on every cold start. hf_xet accelerates the one-time pull.
    # vLLM 0.6.3 pins torch 2.4 but leaves transformers unbounded; a too-new
    # transformers imports `torch.distributed.tensor.device_mesh` (torch >=2.5)
    # and crashes `import vllm`. Pin transformers to the 0.6.3-compatible 4.45.x
    # (still supports Qwen2.5). torch is left to vLLM's own pin.
    inference_image = (
        modal.Image.debian_slim(python_version="3.12")
        .pip_install("vllm==0.6.3", "transformers==4.45.2", "huggingface_hub", "hf_xet")
        .env({"HF_HOME": f"{MODEL_CACHE_DIR}/hf", "HF_XET_HIGH_PERFORMANCE": "1"})
    )

    app = modal.App("behav3d-assistant")
    volume = modal.Volume.from_name(VOLUME_NAME, create_if_missing=True)
    # Persistent cache for the (~15 GB) model weights — downloaded once, then read
    # from the Volume on every cold start.
    model_cache = modal.Volume.from_name("behav3d-model-cache", create_if_missing=True)
    # A shared secret the napari client must present as a Bearer token.
    auth_secret = modal.Secret.from_name("behav3d-assistant-auth")  # {"BEHAV3D_ASSISTANT_TOKEN": "..."}

    # -- One-time weight download (run: modal run chatbot/app.py::download_model) --
    @app.function(image=inference_image, volumes={MODEL_CACHE_DIR: model_cache}, timeout=1800)
    def download_model():
        from huggingface_hub import snapshot_download
        snapshot_download(MODEL_NAME, revision=MODEL_REVISION)
        model_cache.commit()
        print(f"Cached {MODEL_NAME}@{MODEL_REVISION} to the model-cache Volume.")

    # -- Index build -------------------------------------------------------
    @app.function(image=service_image, volumes={INDEX_DIR: volume}, timeout=900)
    def ingest():
        import ingest as ing  # /root/ingest.py
        n = ing.build_and_save(repo_root="/root/repo",
                               schema_cards_path="/root/schema_cards.json",
                               out_dir=INDEX_DIR)
        volume.commit()
        print(f"Built RAG index with {n} chunks.")

    # -- Inference engine --------------------------------------------------
    # Mounts the weights Volume (not the index). scaledown_window keeps a warm
    # container for 5 min after the last request. For zero cold starts during a
    # working session, add `min_containers=1` (you pay for an idle L4 while set).
    @app.cls(image=inference_image, gpu="L4", volumes={MODEL_CACHE_DIR: model_cache},
             scaledown_window=300, timeout=600)
    class Engine:
        @modal.enter()
        def load(self):
            # Weights are read from the mounted Volume (HF_HOME), so no network
            # download happens here — only the disk→GPU load.
            from vllm import LLM, SamplingParams
            self.LLM = LLM(model=MODEL_NAME, revision=MODEL_REVISION,
                           max_model_len=8192, gpu_memory_utilization=0.90)
            self.SamplingParams = SamplingParams

        @modal.method()
        def generate(self, messages: list[dict]) -> str:
            params = self.SamplingParams(temperature=0.3, top_p=0.9, max_tokens=900)
            # vLLM applies the model's chat template for [{"role","content"}]
            outputs = self.LLM.chat(messages, params)
            return outputs[0].outputs[0].text

    # -- Web endpoint ------------------------------------------------------
    @app.function(image=service_image, volumes={INDEX_DIR: volume}, secrets=[auth_secret],
                  scaledown_window=300)
    @modal.asgi_app()
    def web():
        import os
        from fastapi import FastAPI, Request, HTTPException
        from fastapi.responses import StreamingResponse
        from embeddings import Embedder, VectorIndex

        api = FastAPI(title="BEHAV3D Assistant")
        embedder = Embedder()
        try:
            index = VectorIndex.load(INDEX_DIR)
        except Exception:
            index = VectorIndex()
        expected_token = os.environ.get("BEHAV3D_ASSISTANT_TOKEN", "")

        def _auth(request: Request):
            if not expected_token:
                return
            auth = request.headers.get("authorization", "")
            if auth != f"Bearer {expected_token}":
                raise HTTPException(status_code=401, detail="bad token")

        @api.get("/health")
        def health():
            return {"ok": True, "chunks": len(index.chunks)}

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
            retrieved = []
            try:
                qv = embedder.encode([query])[0]
                retrieved = index.search(qv, k=6)
            except Exception:
                retrieved = []

            system = build_system_prompt(context, retrieved, tools)
            convo = [{"role": "system", "content": system}]
            convo += [m for m in messages if m.get("role") != "system"]

            def event_stream():
                try:
                    full = Engine().generate.remote(convo)
                except Exception as e:
                    # Surface engine failures (crash-loop, OOM, timeout) immediately
                    # instead of letting the client hang until its own timeout.
                    yield "data: " + json.dumps(
                        {"type": "error", "message": f"Inference engine error: {e}"}
                    ) + "\n\n"
                    yield "data: " + json.dumps({"type": "done"}) + "\n\n"
                    return
                visible = split_streamable(full)
                # chunked pseudo-streaming for a responsive feel
                for i in range(0, len(visible), 24):
                    piece = visible[i:i + 24]
                    yield "data: " + json.dumps({"type": "token", "text": piece}) + "\n\n"
                _, calls = parse_tool_calls(full)
                if calls:
                    yield "data: " + json.dumps({"type": "tool_calls", "calls": calls}) + "\n\n"
                yield "data: " + json.dumps({"type": "done"}) + "\n\n"

            return StreamingResponse(event_stream(), media_type="text/event-stream")

        return api
