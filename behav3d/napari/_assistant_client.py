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
    {"type": "status",     "stage": "provider", "message": "..."}
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
import time
from pathlib import Path

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
            "timeout": cfg.get("timeout", 45)}


def _configured_read_timeout(config: dict) -> float:
    timeout = config.get("timeout", 45)
    if isinstance(timeout, (tuple, list)):
        timeout = timeout[-1]
    try:
        return max(1.0, float(timeout))
    except (TypeError, ValueError):
        return 45.0


def classify_request_failure(
    error: Exception,
    *,
    stage: str = "connecting",
    read_timeout: float = 45,
) -> dict:
    """Turn transport exceptions into component-specific, user-facing status."""
    error_name = type(error).__name__
    response = getattr(error, "response", None)
    status_code = getattr(response, "status_code", None)

    if error_name == "ConnectTimeout":
        return {
            "level": "offline",
            "stage": "modal",
            "component": "modal",
            "code": "modal_connect_timeout",
            "message": (
                "Modal did not respond within 10 seconds. The service may be "
                "starting or temporarily unavailable."
            ),
        }
    if error_name == "ReadTimeout":
        seconds = int(read_timeout) if float(read_timeout).is_integer() else read_timeout
        if stage in {"provider", "streaming"}:
            message = (
                f"DeepSeek stopped responding for {seconds} seconds. Modal is reachable, "
                "but the model response timed out."
            )
            component = "deepseek"
            code = "deepseek_timeout"
        else:
            message = (
                f"Modal was reached but did not send an update for {seconds} seconds."
            )
            component = "modal"
            code = "modal_read_timeout"
        return {
            "level": "error",
            "stage": stage,
            "component": component,
            "code": code,
            "message": message,
        }
    if error_name == "ConnectionError":
        return {
            "level": "offline",
            "stage": "modal",
            "component": "modal",
            "code": "modal_unreachable",
            "message": (
                "Could not reach Modal. Check the internet connection and assistant endpoint."
            ),
        }
    if status_code is not None:
        return {
            "level": "error",
            "stage": "modal",
            "component": "modal",
            "code": "modal_http_error",
            "http_status": int(status_code),
            "message": f"Modal returned HTTP {status_code} instead of an assistant response.",
        }
    if isinstance(error, (json.JSONDecodeError, ValueError)):
        return {
            "level": "error",
            "stage": stage,
            "component": "modal",
            "code": "invalid_response",
            "message": "The assistant service returned an invalid response.",
        }
    return {
        "level": "error",
        "stage": stage,
        "component": "unknown",
        "code": "request_failed",
        "message": f"The assistant request failed ({error_name}).",
    }


def request_failure_is_retryable(info: dict) -> bool:
    """Return True for transient failures that are safe to retry before output."""
    code = str((info or {}).get("code") or "")
    if code == "modal_http_error":
        return int((info or {}).get("http_status") or 0) in {429, 502, 503, 504}
    return code in {
        "modal_connect_timeout",
        "modal_read_timeout",
        "modal_unreachable",
        "deepseek_timeout",
        "provider_error",
    }


def diagnose_assistant_service(
    config: dict | None = None,
    *,
    probe_provider: bool = False,
    requester=None,
) -> dict:
    """Check Modal/RAG health and optionally make a minimal DeepSeek request.

    ``requester`` is injectable so the response mapping can be unit-tested without
    a live endpoint. The normal startup check is lightweight; the explicit UI
    action sets ``probe_provider=True`` to verify DeepSeek as well.
    """
    config = dict(config or load_client_config())
    endpoint = str(config.get("endpoint") or "").rstrip("/")
    if not endpoint:
        return {
            "level": "offline",
            "stage": "configuration",
            "component": "client",
            "code": "endpoint_missing",
            "modal": "not_configured",
            "deepseek": "not_checked",
            "message": "Assistant endpoint is not configured.",
        }

    if requester is None:
        import requests
        requester = requests.get

    started = time.monotonic()
    try:
        response = requester(
            f"{endpoint}/health",
            params={"probe_provider": "true"} if probe_provider else None,
            timeout=(5, min(_configured_read_timeout(config), 30)),
        )
        response.raise_for_status()
        payload = response.json()
        if not isinstance(payload, dict):
            raise ValueError("health response is not an object")
    except Exception as error:
        result = classify_request_failure(
            error, stage="modal", read_timeout=min(_configured_read_timeout(config), 30)
        )
        result.update({
            "modal": "offline",
            "deepseek": "not_checked",
            "latency_ms": round((time.monotonic() - started) * 1000),
            "endpoint": endpoint,
        })
        return result

    service = payload.get("service") or {}
    retrieval = payload.get("retrieval") or {}
    provider = payload.get("provider") or {}
    modal_status = service.get("status", "online" if payload.get("ok") else "error")
    retrieval_status = retrieval.get(
        "status", "ready" if int(payload.get("chunks") or 0) > 0 else "empty"
    )
    provider_status = provider.get("status", "not_checked")
    level = "online"
    if modal_status != "online" or provider_status == "error":
        level = "error"
    elif retrieval_status != "ready":
        level = "degraded"

    if provider_status == "online":
        message = "Modal, BEHAV3D guidance, and DeepSeek are online."
    elif provider_status == "error":
        message = "Modal is online, but the DeepSeek health check failed."
    elif retrieval_status != "ready":
        message = "Modal is online, but the BEHAV3D guidance index is unavailable."
    else:
        message = "Modal and BEHAV3D guidance are online. DeepSeek was not tested."

    return {
        "level": level,
        "stage": "health",
        "component": "service",
        "code": "health_check",
        "message": message,
        "modal": modal_status,
        "deepseek": provider_status,
        "retrieval": retrieval_status,
        "chunks": retrieval.get("chunks", payload.get("chunks")),
        "model": payload.get("model"),
        "knowledge_version": payload.get("knowledge_version"),
        "provider_latency_ms": provider.get("latency_ms"),
        "provider_error": provider.get("error_type"),
        "latency_ms": round((time.monotonic() - started) * 1000),
        "endpoint": endpoint,
    }


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
            return ("Local help is limited to parameter definitions. Try asking, "
                    "for example, “What does search range mean?”")
        lines = ["I can still answer from the local parameter reference:\n"]
        from behav3d.napari._assistant_actions import humanize_parameter_key
        for c in hits[:3]:
            dv = c.get("default")
            ch = f" Choices: {c['choices']}." if c.get("choices") else ""
            lines.append(
                f"**{humanize_parameter_key(c['key'])}** (default `{dv}`)\n"
                f"{c['description']}{ch}\n"
            )
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Background worker
# ---------------------------------------------------------------------------
class ChatWorker(QObject):
    """Streams one assistant turn. Lives in a QThread (moveToThread)."""

    token = Signal(str)            # incremental assistant text
    tool_calls = Signal(object)    # list[dict] of proposed calls
    status = Signal(object)        # stage/component/message diagnostics
    error = Signal(str)
    degraded = Signal(str)         # emitted (full text) when answered offline
    finished = Signal()            # always last; quits the thread

    def __init__(self, messages: list[dict], context: dict, tools: list[dict]):
        super().__init__()
        self._messages = messages
        self._context = context
        self._tools = tools
        self._last_stage = "connecting"
        self._stream_started = False
        self._received_output = False

    def run(self):
        cfg = load_client_config()
        endpoint = cfg.get("endpoint")
        try:
            if not endpoint:
                self._fallback({
                    "level": "offline",
                    "stage": "configuration",
                    "component": "client",
                    "code": "endpoint_missing",
                    "message": "Assistant endpoint is not configured.",
                })
                return
            for attempt in range(2):
                self._last_stage = "connecting"
                self._stream_started = False
                try:
                    failure = self._stream(cfg)
                except Exception as error:
                    failure = classify_request_failure(
                        error,
                        stage=self._last_stage,
                        read_timeout=_configured_read_timeout(cfg),
                    )

                if failure is None:
                    return
                can_retry = (
                    attempt == 0
                    and not self._received_output
                    and request_failure_is_retryable(failure)
                )
                if can_retry:
                    self._emit_status({
                        "level": "working",
                        "stage": "retrying",
                        "component": failure.get("component", "service"),
                        "code": "automatic_retry",
                        "message": (
                            f"{failure.get('message', 'The request was interrupted')} "
                            "Retrying automatically in 2 seconds..."
                        ),
                        "attempt": 2,
                        "request_id": failure.get("request_id"),
                    })
                    time.sleep(2)
                    continue

                if failure.get("server_error") or self._received_output:
                    self._emit_status(failure)
                    self.error.emit(failure.get("message", "The assistant request failed."))
                else:
                    self._fallback(failure)
                return
        finally:
            self.finished.emit()

    def _emit_status(self, info: dict):
        info = dict(info or {})
        if info.get("stage"):
            self._last_stage = str(info["stage"])
        self.status.emit(info)

    # -- networking ---------------------------------------------------------
    def _stream(self, cfg: dict):
        import requests  # available in the behav3d env
        headers = {"Content-Type": "application/json"}
        payload = {"messages": self._messages, "context": self._context, "tools": self._tools}
        # Split timeout: cap connect at 10s, read at the configured value, so a
        # stalled server doesn't leave the dock "thinking" (and buttons greyed) for
        # the whole window. A warm Modal container (min_containers=1) makes 10s ample.
        timeout = cfg.get("timeout", 45)
        ct = timeout if isinstance(timeout, (tuple, list)) else (10, timeout)
        self._emit_status({
            "level": "working",
            "stage": "connecting",
            "component": "modal",
            "message": "Connecting to Modal...",
        })
        with requests.post(f"{cfg['endpoint']}/chat", json=payload, headers=headers,
                           stream=True, timeout=ct) as resp:
            resp.raise_for_status()
            self._emit_status({
                "level": "working",
                "stage": "connected",
                "component": "modal",
                "message": "Modal connected. Preparing the response...",
            })
            for raw in resp.iter_lines(decode_unicode=True):
                if not raw:
                    continue
                line = raw[6:] if raw.startswith("data: ") else raw  # tolerate SSE prefix
                try:
                    evt = json.loads(line)
                except json.JSONDecodeError:
                    continue
                etype = evt.get("type")
                if etype == "status":
                    if evt.get("stage") == "streaming":
                        self._stream_started = True
                    self._emit_status(evt)
                elif etype == "token":
                    self._received_output = self._received_output or bool(evt.get("text"))
                    if not self._stream_started:
                        self._stream_started = True
                        local_guidance = self._last_stage == "local_guidance"
                        self._emit_status({
                            "level": "working",
                            "stage": "streaming",
                            "component": "modal" if local_guidance else "deepseek",
                            "message": (
                                "Receiving BEHAV3D guidance..."
                                if local_guidance else
                                "Receiving the assistant response..."
                            ),
                            "request_id": evt.get("request_id"),
                        })
                    self.token.emit(evt.get("text", ""))
                elif etype == "tool_calls":
                    self._received_output = self._received_output or bool(evt.get("calls"))
                    if not self._stream_started:
                        self._stream_started = True
                        self._emit_status({
                            "level": "working",
                            "stage": "streaming",
                            "component": "deepseek",
                            "message": "Receiving proposed form changes...",
                            "request_id": evt.get("request_id"),
                        })
                    self.tool_calls.emit(evt.get("calls", []))
                elif etype == "error":
                    return {
                        "level": "error",
                        "stage": evt.get("stage", self._last_stage),
                        "component": evt.get("component", "server"),
                        "code": evt.get("code", "server_error"),
                        "message": evt.get("message", "Unknown server error."),
                        "request_id": evt.get("request_id"),
                        "error_type": evt.get("error_type"),
                        "http_status": evt.get("http_status"),
                        "server_error": True,
                    }
                elif etype == "done":
                    return
            raise ConnectionError("Assistant stream ended before a done event")

    # -- degraded mode ------------------------------------------------------
    def _fallback(self, reason: dict | str):
        if isinstance(reason, str):
            reason = {
                "level": "offline",
                "stage": self._last_stage,
                "component": "unknown",
                "code": "request_failed",
                "message": reason,
            }
        self._emit_status(reason)
        user_msg = ""
        for m in reversed(self._messages):
            if m.get("role") == "user":
                user_msg = m.get("content", "")
                break
        local_answer = LocalResponder().answer(user_msg)
        text = (
            "**Assistant offline**\n\n"
            f"**Reason:** {reason.get('message', 'The assistant service is unavailable.')}\n\n"
            f"{local_answer}"
        )
        self.degraded.emit(text)


class HealthWorker(QObject):
    """Runs a service health check without blocking the napari GUI thread."""

    result = Signal(object)
    finished = Signal()

    def __init__(self, *, probe_provider: bool = False):
        super().__init__()
        self._probe_provider = bool(probe_provider)

    def run(self):
        try:
            self.result.emit(diagnose_assistant_service(
                probe_provider=self._probe_provider
            ))
        finally:
            self.finished.emit()
