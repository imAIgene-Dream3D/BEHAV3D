# BEHAV3D Assistant — Modal service (DeepSeek proxy)

A **CPU-only** Modal service for the BEHAV3D napari co-pilot. It authenticates the
client, retrieves BEHAV3D knowledge (RAG), then calls the **DeepSeek API** with
native function-calling and relays the stream back. No GPU, no model weights. The
DeepSeek key stays server-side (Modal secret); clients only ever see the Modal URL
+ a shared bearer token.

## What's here
| File | Role |
|------|------|
| `app.py` | Modal app: RAG retrieval + DeepSeek call, SSE `/chat` + `/health` (Bearer auth). |
| `ingest.py` | Builds the RAG index from `README.md`, `WIKI.md` (if present), `schema_cards.json`, and module docstrings. |
| `embeddings.py` | sentence-transformers embedder + tiny numpy cosine index (persisted to a Modal Volume). |
| `schema_cards.json` | Decoupled snapshot of the 195 BEHAV3D parameter cards; also the `enum` that constrains `set_parameter.key`. |

Model: DeepSeek **`deepseek-v4-flash`** by default (V4 Flash; tool-calls +
streaming). Override with the `DEEPSEEK_MODEL` env var — e.g. `deepseek-v4-pro`
for the stronger/pricier model — by adding it to the `deepseek-api-key` secret.
(The old `deepseek-chat` alias maps to v4-flash too but is being deprecated.)

## One-time setup
```bash
pip install modal && python -m modal token new          # authenticate Modal
# shared bearer token the napari client must present:
python -m modal secret create behav3d-assistant-auth BEHAV3D_ASSISTANT_TOKEN=$(openssl rand -hex 16)
# your DeepSeek API key (stays server-side):
python -m modal secret create deepseek-api-key DEEPSEEK_API_KEY=sk-...   # optionally DEEPSEEK_MODEL=...
```

## Build the index, then deploy
```bash
python -m modal run    chatbot/app.py::ingest    # build/refresh the RAG index (re-run after editing WIKI.md)
python -m modal deploy chatbot/app.py            # fast — CPU image only; prints the public URL
```
No weight download and no GPU cold start: the endpoint is responsive immediately
(only the DeepSeek API latency). It stays warm for `scaledown_window` (300 s); for a
zero-cold-start endpoint add `min_containers=1` to `web`'s `@app.function(...)`
(cheap on CPU, ~$15–45/mo vs ~$0 idle).

## Point the napari client at it
Copy the deployed URL and the bearer token into `napari/assistant_config.json`
(gitignored; see `napari/assistant_config.example.json`):
```json
{ "endpoint": "https://<your-app>.modal.run", "token": "<the-bearer-token>", "timeout": 60 }
```
Env vars `BEHAV3D_ASSISTANT_ENDPOINT` / `BEHAV3D_ASSISTANT_TOKEN` override the file.
If no endpoint is configured or the service is unreachable, the napari dock falls
back to **offline mode** (schema-grounded parameter answers, no LLM).

## Local development
```bash
python -m modal serve chatbot/app.py             # hot-reloading proxy; gives a temporary URL
curl -N -X POST <url>/chat -H "Authorization: Bearer <token>" \
  -H "Content-Type: application/json" \
  -d '{"messages":[{"role":"user","content":"set the trackpy search range to 40 for immune"}],"context":{"current_step":"tracking"}}'
```

## Refreshing the parameter schema
`schema_cards.json` is generated from the live BEHAV3D config:
```bash
conda run -n behav3d python -c "from behav3d.napari._assistant_schema import dump_cards_json; dump_cards_json('chatbot/schema_cards.json')"
```
Re-run `python -m modal run chatbot/app.py::ingest` afterwards.

## Tearing down the old GPU deployment
The earlier Qwen/vLLM setup is gone from the code. Remove its leftover cloud
resources once (they cost nothing idle, but clean them up):
```bash
python -m modal app list                                  # find the old app if separate
python -m modal volume delete behav3d-model-cache         # the ~15 GB weights cache
```
Re-`python -m modal deploy chatbot/app.py` replaces the running `behav3d-assistant` app with
the CPU version (no GPU containers remain).
