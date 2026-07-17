"""Reusable widgets for the Analysis tab's Guided / Advanced mode.

This module is deliberately free of any analysis logic — it only renders the
explanation-first "Guided" on-ramp and the small segmented switch that flips a
sub-tab between Guided and its existing (Advanced) settings form.

Pieces
------
- :class:`ModeSwitch`      — segmented "Guided | Advanced" control (a signal).
- :class:`AnalysisExplainer` — one analysis block: real name + plain-language
  subtitle, a collapsed explainer behind a tinted "See what this does" cue, an
  honest "what you'll decide" list, an "Ask the assistant" button, and an
  always-visible "Start" button.
- :class:`GuidedPanel`     — the scrollable list of explainers for one sub-tab.

Copy lives in :mod:`behav3d.napari.analysis_guided_copy`; the two sub-tab
modules own the Start wiring and the QStackedWidget that holds Guided (page 0)
and the existing form (page 1).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Callable

from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QToolButton,
    QFrame, QScrollArea, QSizePolicy,
)
from qtpy.QtCore import Qt, Signal


# ── Mode persistence ────────────────────────────────────────────────────────
def _ui_config_path() -> Path:
    # Sits next to napari/assistant_config.json (see _assistant_client).
    return Path(__file__).resolve().parents[2] / "napari" / "analysis_ui_config.json"


def load_guided_mode(default: bool = True) -> bool:
    """Return the saved Guided/Advanced preference (True = Guided)."""
    p = _ui_config_path()
    if p.exists():
        try:
            cfg = json.loads(p.read_text(encoding="utf-8"))
            return bool(cfg.get("analysis_guided_mode", default))
        except Exception:
            return default
    return default


def save_guided_mode(guided: bool) -> None:
    """Persist the Guided/Advanced preference. Failures are non-fatal."""
    p = _ui_config_path()
    try:
        cfg = {}
        if p.exists():
            try:
                cfg = json.loads(p.read_text(encoding="utf-8"))
            except Exception:
                cfg = {}
        cfg["analysis_guided_mode"] = bool(guided)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(cfg, indent=2), encoding="utf-8")
    except Exception:
        pass


# ── Assistant hand-off ──────────────────────────────────────────────────────
def find_main_widget(w):
    """Walk the Qt parent chain to the BEHAV3DWidget (holds ``assistant``)."""
    node = w
    while node is not None:
        if hasattr(node, "assistant") and hasattr(node, "_assistant_dock"):
            return node
        node = node.parent() if hasattr(node, "parent") else None
    return None


def ask_assistant(origin: QWidget, prompt: str) -> None:
    """Reveal the assistant dock and send ``prompt`` on the user's behalf."""
    main = find_main_widget(origin)
    if main is None or getattr(main, "assistant", None) is None:
        return
    try:
        main.ask_assistant(prompt)
    except Exception:
        pass


# ── Tag colouring for the "what you'll decide" list ─────────────────────────
_TAG_COLORS = {
    "required": "#ff8a80",    # only you know this
    "suggested": "#90caf9",   # suggested set, adjustable
    "estimated": "#90caf9",   # can be estimated
    "default": "#888888",     # tested defaults, optional
}


# ── Mode switch ─────────────────────────────────────────────────────────────
class ModeSwitch(QWidget):
    """Segmented 'Guided | Advanced' control. Emits ``modeChanged(bool guided)``."""

    modeChanged = Signal(bool)

    def __init__(self, guided: bool = True, parent=None):
        super().__init__(parent)
        self._guided = guided

        lay = QHBoxLayout(self)
        lay.setContentsMargins(6, 4, 6, 4)
        lay.setSpacing(0)
        lay.addStretch()

        self._btn_guided = QPushButton("Guided")
        self._btn_adv = QPushButton("Advanced")
        for b in (self._btn_guided, self._btn_adv):
            b.setCheckable(True)
            b.setCursor(Qt.PointingHandCursor)
            b.setFixedHeight(26)
        self._btn_guided.clicked.connect(lambda: self.set_guided(True))
        self._btn_adv.clicked.connect(lambda: self.set_guided(False))
        lay.addWidget(self._btn_guided)
        lay.addWidget(self._btn_adv)

        self._restyle()

    def is_guided(self) -> bool:
        return self._guided

    def set_guided(self, guided: bool, *, emit: bool = True):
        guided = bool(guided)
        changed = guided != self._guided
        self._guided = guided
        self._restyle()
        if emit and changed:
            self.modeChanged.emit(guided)

    def _restyle(self):
        self._btn_guided.setChecked(self._guided)
        self._btn_adv.setChecked(not self._guided)
        on = (
            "QPushButton { background:#1e88e5; color:white; font-weight:bold; "
            "border:1px solid #1e88e5; padding:2px 16px; font-size:12px; }"
        )
        off = (
            "QPushButton { background:transparent; color:#9aa0a6; "
            "border:1px solid #555; padding:2px 16px; font-size:12px; }"
            "QPushButton:hover { color:#ddd; }"
        )
        # Round only the outer edges of the segmented pair.
        left_on = on.replace("padding:2px 16px;", "padding:2px 16px; border-top-left-radius:13px; border-bottom-left-radius:13px;")
        left_off = off.replace("padding:2px 16px;", "padding:2px 16px; border-top-left-radius:13px; border-bottom-left-radius:13px;")
        right_on = on.replace("padding:2px 16px;", "padding:2px 16px; border-top-right-radius:13px; border-bottom-right-radius:13px;")
        right_off = off.replace("padding:2px 16px;", "padding:2px 16px; border-top-right-radius:13px; border-bottom-right-radius:13px;")
        self._btn_guided.setStyleSheet(left_on if self._guided else left_off)
        self._btn_adv.setStyleSheet(right_on if not self._guided else right_off)


# ── One analysis explainer ──────────────────────────────────────────────────
class AnalysisExplainer(QWidget):
    """A single analysis block for the Guided page.

    ``spec`` is one of the dicts in :mod:`analysis_guided_copy`. ``on_start`` is
    called (no args) when the user presses Start.
    """

    def __init__(self, spec: dict, on_start: Callable[[], None], parent=None):
        super().__init__(parent)
        self._spec = spec
        self._on_start = on_start

        card = QFrame(self)
        card.setObjectName("guidedCard")
        card.setStyleSheet(
            "#guidedCard { background:#2a2a2e; border:1px solid #3a3a40; "
            "border-radius:10px; }"
        )
        wrap = QVBoxLayout(self)
        wrap.setContentsMargins(0, 0, 0, 6)
        wrap.addWidget(card)

        lay = QVBoxLayout(card)
        lay.setContentsMargins(12, 10, 12, 10)
        lay.setSpacing(8)

        # Header: colored marker + real name, plain-language subtitle underneath.
        title = QLabel(
            f"<span style='color:{spec['color']}; font-size:15px;'>●</span>"
            f"&nbsp;<span style='font-weight:bold; font-size:14px; color:#eaeaea;'>"
            f"{spec['title']}</span>"
        )
        title.setTextFormat(Qt.RichText)
        lay.addWidget(title)
        subtitle = QLabel(spec["subtitle"])
        subtitle.setWordWrap(True)
        subtitle.setStyleSheet("color:#9aa0a6; font-size:12px; margin-left:16px;")
        lay.addWidget(subtitle)

        # Tinted "See what this does" cue that toggles the explainer body.
        self._cue = QToolButton()
        self._cue.setCheckable(True)
        self._cue.setChecked(False)
        self._cue.setCursor(Qt.PointingHandCursor)
        self._cue.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self._cue.setArrowType(Qt.RightArrow)
        self._cue.setText("ⓘ  New to this? See what this analysis does")
        self._cue.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self._cue.setStyleSheet(
            "QToolButton { background:#1f3a5f; color:#bbdefb; border:1px solid "
            "#3f6ea5; border-radius:6px; padding:6px 10px; text-align:left; "
            "font-size:12px; font-weight:bold; }"
            "QToolButton:hover { background:#274a78; }"
        )
        self._cue.clicked.connect(self._toggle_body)
        lay.addWidget(self._cue)

        # Explainer body (hidden by default).
        self._body = QFrame()
        self._body.setVisible(False)
        body = QVBoxLayout(self._body)
        body.setContentsMargins(4, 4, 4, 2)
        body.setSpacing(8)

        what = QLabel(spec["what_does"])
        what.setWordWrap(True)
        what.setStyleSheet("color:#dcdcdc; font-size:12px;")
        body.addWidget(what)

        if spec.get("concept"):
            c = spec["concept"]
            callout = QLabel(
                f"<span style='color:#ffd54f;'>💡</span> "
                f"<b style='color:#eaeaea;'>{c['term']}</b> "
                f"<span style='color:#c5cae9;'>{c['text']}</span>"
            )
            callout.setTextFormat(Qt.RichText)
            callout.setWordWrap(True)
            callout.setStyleSheet(
                "background:#20304a; border-radius:6px; padding:8px 10px; "
                "font-size:12px;"
            )
            body.addWidget(callout)

        get = QLabel(f"<b>What you'll get:</b>&nbsp;{spec['what_get']}")
        get.setTextFormat(Qt.RichText)
        get.setWordWrap(True)
        get.setStyleSheet("color:#c8c8c8; font-size:12px;")
        body.addWidget(get)

        decide_hdr = QLabel("WHAT YOU'LL DECIDE")
        decide_hdr.setStyleSheet(
            "color:#8a8f98; font-size:10px; font-weight:bold; letter-spacing:1px;"
        )
        body.addWidget(decide_hdr)
        hint = QLabel(
            "These depend on your biology — no setting can pick them perfectly. "
            "The assistant can go deeper on any of them."
        )
        hint.setWordWrap(True)
        hint.setStyleSheet("color:#7f858d; font-size:11px;")
        body.addWidget(hint)

        for item in spec.get("decide", []):
            color = _TAG_COLORS.get(item.get("kind", "default"), "#888888")
            row = QLabel(
                f"<span style='color:#cfcfcf; font-size:12px;'>• {item['label']}</span> "
                f"<span style='color:{color}; font-size:11px;'>{item['tag']}</span>"
            )
            row.setTextFormat(Qt.RichText)
            row.setWordWrap(True)
            body.addWidget(row)

        ask = QPushButton("💬  Ask the assistant")
        ask.setCursor(Qt.PointingHandCursor)
        ask.setStyleSheet(
            "QPushButton { background:transparent; color:#bbdefb; border:1px "
            "solid #3f6ea5; border-radius:6px; padding:5px 12px; font-size:12px; }"
            "QPushButton:hover { background:#20304a; }"
        )
        ask.clicked.connect(self._on_ask)
        ask_row = QHBoxLayout()
        ask_row.addWidget(ask)
        ask_row.addStretch()
        body.addLayout(ask_row)

        lay.addWidget(self._body)

        # Always-visible primary action.
        start = QPushButton("Start — open the settings  ▸")
        start.setCursor(Qt.PointingHandCursor)
        start.setStyleSheet(
            "QPushButton { background:#1e88e5; color:white; font-weight:bold; "
            "border:none; border-radius:6px; padding:7px 16px; font-size:13px; }"
            "QPushButton:hover { background:#1976d2; }"
        )
        start.clicked.connect(lambda: self._on_start())
        start_row = QHBoxLayout()
        start_row.addWidget(start)
        start_row.addStretch()
        lay.addLayout(start_row)

    def _toggle_body(self):
        show = self._cue.isChecked()
        self._body.setVisible(show)
        self._cue.setArrowType(Qt.DownArrow if show else Qt.RightArrow)

    def _on_ask(self):
        ask_assistant(self, self._spec.get("seed", ""))


# ── Guided page for one sub-tab ─────────────────────────────────────────────
class GuidedPanel(QScrollArea):
    """Scrollable list of :class:`AnalysisExplainer` for one sub-tab.

    ``analyses`` is a list of copy specs; ``start_cb`` is called with the
    analysis ``id`` when its Start button is pressed.
    """

    def __init__(self, analyses: list[dict], start_cb: Callable[[str], None],
                 intro: str = "", parent=None):
        super().__init__(parent)
        self.setWidgetResizable(True)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setFrameShape(QFrame.NoFrame)

        content = QWidget()
        lay = QVBoxLayout(content)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.setSpacing(6)

        if intro:
            intro_lbl = QLabel(intro)
            intro_lbl.setWordWrap(True)
            intro_lbl.setStyleSheet(
                "color:#9aa0a6; font-size:11px; padding:2px 4px 6px 4px;"
            )
            lay.addWidget(intro_lbl)

        for spec in analyses:
            aid = spec["id"]
            lay.addWidget(
                AnalysisExplainer(spec, on_start=lambda a=aid: start_cb(a))
            )

        lay.addStretch()
        self.setWidget(content)
