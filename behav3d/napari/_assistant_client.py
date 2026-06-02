"""
BEHAV3D assistant — network client to the Modal-hosted LLM service.

Runs the request on a background QThread worker (QObject + moveToThread, matching
the pattern used by ``_segment_editor._EditWorker``) so streamed tokens never
block the napari GUI thread.

Wire protocol (client → server):  POST {endpoint}/chat
    {
      "messages": [{"role": "user"|"assistant", "content": "..."}],
      "context":  <build_context() dict>,
      "tools":    <TOOL_SCHEMA>,
    }
Server → client: an SSE / newline-delimited-JSON stream of events:
    {"type": "token",      "text": "..."}        # streamed assistant text
    {"type": "tool_calls", "calls": [ {name, arguments}, ... ]}
    {"type": "done"}
    {"type": "error",      "message": "..."}

If no endpoint is configured or the server is unreachable, the client falls back
to :class:`LocalResponder` (degraded mode): schema-grounded answers to
"what is X" / "explain X" with no LLM. This honours the "requires internet"
design while never hard-failing.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Optional

from qtpy.QtCore import QObject, Signal


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
def _config_path() -> Path:
    # napari/assistant_config.json lives next to behav3d_env.json (repo /napari)
    return Path(__file__).resolve().parents[2] / "napari" / "assistant_config.json"


def load_client_config() -> dict:
    """Resolve endpoint + token. Env vars override the JSON file."""
    cfg: dict = {}
    p = _config_path()
    if p.exists():
        try:
            cfg = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            cfg = {}
    endpoint = os.environ.get("BEHAV3D_ASSISTANT_ENDPOINT", cfg.get("endpoint", ""))
    token = os.environ.get("BEHAV3D_ASSISTANT_TOKEN", cfg.get("token", ""))
    return {"endpoint": endpoint.rstrip("/"), "token": token,
            "timeout": cfg.get("timeout", 60)}


# ---------------------------------------------------------------------------
# Degraded / offline responder
# ---------------------------------------------------------------------------
class LocalResponder:
    """Answers parameter questions from the local schema cards — no network."""

    def __init__(self):
        try:
            from behav3d.napari._assistant_schema import flatten_config_to_cards
            self.cards = flatten_config_to_cards()
        except Exception:
            self.cards = []

    def answer(self, message: str) -> str:
        q = (message or "").lower().strip()
        # find cards whose key/title appears in the question
        hits = [c for c in self.cards
                if c["title"].lower() in q or c["key"].lower() in q]
        if not hits:
            # token overlap fallback
            toks = {t for t in q.replace(".", " ").split() if len(t) > 2}
            scored = []
            for c in self.cards:
                hay = (c["key"] + " " + c["description"]).lower()
                scored.append((sum(t in hay for t in toks), c))
            scored.sort(key=lambda x: x[0], reverse=True)
            hits = [c for s, c in scored[:3] if s > 0]
        if not hits:
            return ("⚠️ Offline mode (no assistant server reachable). I can only "
                    "answer parameter questions locally. Try naming a parameter, "
                    "e.g. “what is search_range_px?”")
        lines = ["⚠️ *Offline mode* — answering from the local parameter schema:\n"]
        for c in hits[:3]:
            dv = c.get("default")
            ch = f" Choices: {c['choices']}." if c.get("choices") else ""
            lines.append(f"**{c['key']}** (default `{dv}`)\n{c['description']}{ch}\n")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Background worker
# ---------------------------------------------------------------------------
class ChatWorker(QObject):
    """Streams one assistant turn. Lives in a QThread (moveToThread)."""

    token = Signal(str)            # incremental assistant text
    tool_calls = Signal(object)    # list[dict] of proposed calls
    error = Signal(str)
    degraded = Signal(str)         # emitted (full text) when answered offline
    finished = Signal()            # always last; quits the thread

    def __init__(self, messages: list[dict], context: dict, tools: list[dict]):
        super().__init__()
        self._messages = messages
        self._context = context
        self._tools = tools

    def run(self):
        cfg = load_client_config()
        endpoint = cfg.get("endpoint")
        if not endpoint:
            self._fallback("No assistant endpoint configured.")
            return
        try:
            self._stream(cfg)
        except Exception as e:  # network / parse failure → degrade
            self._fallback(f"Could not reach assistant server ({e}).")
        finally:
            self.finished.emit()

    # -- networking ---------------------------------------------------------
    def _stream(self, cfg: dict):
        import requests  # available in the behav3d env
        headers = {"Content-Type": "application/json"}
        payload = {"messages": self._messages, "context": self._context, "tools": self._tools}
        with requests.post(f"{cfg['endpoint']}/chat", json=payload, headers=headers,
                           stream=True, timeout=cfg.get("timeout", 60)) as resp:
            resp.raise_for_status()
            for raw in resp.iter_lines(decode_unicode=True):
                if not raw:
                    continue
                line = raw[6:] if raw.startswith("data: ") else raw  # tolerate SSE prefix
                try:
                    evt = json.loads(line)
                except json.JSONDecodeError:
                    continue
                etype = evt.get("type")
                if etype == "token":
                    self.token.emit(evt.get("text", ""))
                elif etype == "tool_calls":
                    self.tool_calls.emit(evt.get("calls", []))
                elif etype == "error":
                    self.error.emit(evt.get("message", "Unknown server error."))
                    return
                elif etype == "done":
                    return

    # -- degraded mode ------------------------------------------------------
    def _fallback(self, reason: str):
        user_msg = ""
        for m in reversed(self._messages):
            if m.get("role") == "user":
                user_msg = m.get("content", "")
                break
        text = LocalResponder().answer(user_msg)
        self.degraded.emit(text)
