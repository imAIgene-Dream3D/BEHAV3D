# BEHAV3D Assistant — Modal service

Self-hosted LLM + RAG backend for the BEHAV3D napari co-pilot. The napari client
(`behav3d/napari/_assistant*.py`) talks to this over HTTPS.

## What's here
| File | Role |
|------|------|
| `app.py` | Modal app: vLLM open model on a GPU, RAG retrieval, tool-calling, SSE `/chat` endpoint (Bearer-token auth). |
| `ingest.py` | Builds the RAG index from `README.md`, `WIKI.md` (if present), `schema_cards.json`, and module docstrings. |
| `embeddings.py` | sentence-transformers embedder + tiny numpy cosine index (persisted to a Modal Volume). |
| `schema_cards.json` | Decoupled snapshot of the 195 BEHAV3D parameter cards (so the cloud side never imports the heavy `behav3d` package). |

Model: `Qwen/Qwen2.5-7B-Instruct` on an `L4` GPU (edit `MODEL_NAME` / `gpu=` in `app.py` to scale up).

## One-time setup
```bash
pip install modal && modal token new          # authenticate Modal
# shared secret the napari client must present as its Bearer token:
modal secret create behav3d-assistant-auth BEHAV3D_ASSISTANT_TOKEN=$(openssl rand -hex 16)
```

## Build the index, cache weights, then deploy
```bash
modal run   chatbot/app.py::download_model   # one-time: caches ~15 GB of weights to a Volume (minutes)
modal run   chatbot/app.py::ingest           # build/refresh the RAG index (re-run after editing WIKI.md)
modal deploy chatbot/app.py                  # deploy /chat + /health; prints the public URL
```

### Cold starts
`download_model` makes every subsequent cold start read weights from the Volume
(`HF_HOME` on `/models`) instead of re-downloading from HuggingFace — cutting cold
start from minutes to ~30–60 s (container boot + disk→GPU load). The `Engine` stays
warm for `scaledown_window` (300 s) after the last request. For **zero** cold starts
during a working session, add `min_containers=1` to the `Engine`'s `@app.cls(...)`
(you pay for an idle L4 while it's set). Pin `MODEL_REVISION` to a commit hash for
reproducible weights.

## Point the napari client at it
Copy the deployed URL and the secret token into `napari/assistant_config.json`
(gitignored; see `napari/assistant_config.example.json`):
```json
{ "endpoint": "https://<your-app>.modal.run", "token": "<the-secret-token>", "timeout": 60 }
```
Env vars `BEHAV3D_ASSISTANT_ENDPOINT` / `BEHAV3D_ASSISTANT_TOKEN` override the file.

If no endpoint is configured or the service is unreachable, the napari dock
automatically falls back to **offline mode** (schema-grounded parameter answers,
no LLM).

## Local development
```bash
modal serve chatbot/app.py             # hot-reloading proxy; gives a temporary URL
curl -N -X POST <url>/chat -H "Authorization: Bearer <token>" \
  -H "Content-Type: application/json" \
  -d '{"messages":[{"role":"user","content":"what diameter for my T-cells?"}],"context":{"current_step":"segmentation"}}'
```

## Refreshing the parameter schema
`schema_cards.json` is generated from the live BEHAV3D config:
```bash
conda run -n behav3d python -c "from behav3d.napari._assistant_schema import dump_cards_json; dump_cards_json('chatbot/schema_cards.json')"
```
Re-run `modal run chatbot/app.py::ingest` afterwards.
