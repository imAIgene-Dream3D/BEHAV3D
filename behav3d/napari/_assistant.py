"""
BEHAV3D assistant — the right-side co-pilot dock widget.

``AssistantDock`` is added as a *separate* napari dock (right area) and holds a
reference to the main ``BEHAV3DWidget`` so it can read live workflow context
(:func:`_assistant_context.build_context`) and apply proposed parameter changes
(:mod:`_assistant_actions`) after the user confirms.

Layout (top → bottom):
  • context bar          — 📍 step · N samples · output dir · M queued (live)
  • chat transcript      — QTextBrowser, markdown, streamed assistant tokens
  • proposed-action tray — confirm cards: "old -> new  [Fill it in] [Dismiss]"
  • quick actions        — Explain this screen · Tell me about my data
  • input row            — text box + Send

All network I/O happens on a QThread worker (``ChatWorker``); the GUI never blocks.
"""
from __future__ import annotations

from typing import Optional

from qtpy.QtCore import Qt, QThread, QTimer, Signal
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QLineEdit,
    QTextBrowser, QFrame, QSizePolicy,
)

from behav3d.napari._assistant_context import build_context, context_summary_line
from behav3d.napari._assistant_actions import (
    build_actions, apply_action, TOOL_SCHEMA, ProposedAction, humanize_parameter_key,
)
from behav3d.napari._assistant_client import ChatWorker, load_client_config


_UNSET = object()  # sentinel for "parameter not yet in config"
_METADATA_FIELD_DEFAULTS = {
    "n_samples": 1,
    "n_organoids": 0,
    "n_immune": 0,
    "n_other": 0,
    "include_dead": False,
    "exp_nr": 1,
    "pixel_distance_xy": 0.5,
    "pixel_distance_z": 2.0,
    "time_interval": 1.0,
    "time_unit": "s",
}


def _coerce_bool(value) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in ("1", "true", "yes", "y", "on")
    return bool(value)


def _metadata_field_has_value(dp, field: str, index: int, cell_type: str | None = None) -> bool:
    """Return True if the metadata builder field already has a non-default value set."""
    if dp is None:
        return False
    try:
        from qtpy.QtWidgets import QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox
        if field == "n_samples":
            return dp.n_samples_spin.value() != _METADATA_FIELD_DEFAULTS["n_samples"]
        if field == "n_organoids":
            return dp.n_organoid_spin.value() != _METADATA_FIELD_DEFAULTS["n_organoids"]
        if field == "n_immune":
            return dp.n_immune_spin.value() != _METADATA_FIELD_DEFAULTS["n_immune"]
        if field == "n_other":
            return dp.n_other_spin.value() != _METADATA_FIELD_DEFAULTS["n_other"]
        if field == "include_dead":
            return dp.include_dead_cb.isChecked() != _METADATA_FIELD_DEFAULTS["include_dead"]
        if field in ("immune_multicolor", "immune_multicolor_channels"):
            flags = getattr(dp, "_immune_multicolor_flags", [])
            counts = getattr(dp, "_immune_multicolor_counts", [])
            if field == "immune_multicolor":
                return index < len(flags) and flags[index].isChecked()
            return index < len(counts) and counts[index].value() != 2
        for attr, cat in (("_organoid_name_edits", "organoid_name"),
                          ("_immune_name_edits",   "immune_name"),
                          ("_other_name_edits",    "other_name")):
            if field == cat:
                edits = getattr(dp, attr, [])
                return index < len(edits) and bool(edits[index].text().strip())
        # Per-sample fields
        forms = getattr(dp, "_sample_forms", [])
        if index >= len(forms):
            return False
        if field in ("dead_channel_number", "dead_mask_path"):
            dead_key = "number" if field == "dead_channel_number" else "mask_path"
            w = forms[index].get("dead_channel", {}).get(dead_key)
        elif field in (
            "cell_line", "cell_condition", "cell_segments_image_path",
            "cell_tracks_image_path", "cell_tracks_csv_path",
        ):
            field_map = {
                "cell_line": "line",
                "cell_condition": "condition",
                "cell_segments_image_path": "segments_image_path",
                "cell_tracks_image_path": "tracks_image_path",
                "cell_tracks_csv_path": "tracks_csv_path",
            }
            cell_fields = forms[index].get("cell_types", {}).get(str(cell_type), {})
            w = cell_fields.get(field_map[field])
        else:
            w = forms[index]["basic"].get(field)
        if w is None:
            return False
        if isinstance(w, QLineEdit):
            return bool(w.text().strip())
        if isinstance(w, (QSpinBox, QDoubleSpinBox)):
            val = w.value()
            if field in _METADATA_FIELD_DEFAULTS:
                return val != _METADATA_FIELD_DEFAULTS[field]
            return val != 0 and val != w.minimum()
        if isinstance(w, QComboBox):
            default = _METADATA_FIELD_DEFAULTS.get(field)
            if default is not None:
                return w.currentText() != str(default)
            return w.currentIndex() > 0
    except Exception:
        pass
    return False


def _fill_value_matches_current(
    dp,
    field: str,
    index: int,
    proposed,
    cell_type: str | None = None,
) -> bool:
    """Return True if the proposed fill value is the same as what's already in the widget.
    Used to auto-apply no-op repeated calls without showing a confirm card."""
    if dp is None:
        return False
    try:
        from qtpy.QtWidgets import QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox
        if field == "n_samples":
            return dp.n_samples_spin.value() == int(proposed)
        if field == "n_organoids":
            return dp.n_organoid_spin.value() == int(proposed)
        if field == "n_immune":
            return dp.n_immune_spin.value() == int(proposed)
        if field == "n_other":
            return dp.n_other_spin.value() == int(proposed)
        if field == "include_dead":
            return dp.include_dead_cb.isChecked() == _coerce_bool(proposed)
        if field in ("immune_multicolor", "immune_multicolor_channels"):
            flags = getattr(dp, "_immune_multicolor_flags", [])
            counts = getattr(dp, "_immune_multicolor_counts", [])
            if field == "immune_multicolor":
                return index < len(flags) and flags[index].isChecked() == _coerce_bool(proposed)
            return index < len(counts) and counts[index].value() == int(proposed)
        for attr, cat in (("_organoid_name_edits", "organoid_name"),
                          ("_immune_name_edits",   "immune_name"),
                          ("_other_name_edits",    "other_name")):
            if field == cat:
                edits = getattr(dp, attr, [])
                return (index < len(edits) and
                        edits[index].text().strip() == str(proposed).strip())
        forms = getattr(dp, "_sample_forms", [])
        if index >= len(forms):
            return False
        if field in ("dead_channel_number", "dead_mask_path"):
            dead_key = "number" if field == "dead_channel_number" else "mask_path"
            w = forms[index].get("dead_channel", {}).get(dead_key)
        elif field in (
            "cell_line", "cell_condition", "cell_segments_image_path",
            "cell_tracks_image_path", "cell_tracks_csv_path",
        ):
            field_map = {
                "cell_line": "line",
                "cell_condition": "condition",
                "cell_segments_image_path": "segments_image_path",
                "cell_tracks_image_path": "tracks_image_path",
                "cell_tracks_csv_path": "tracks_csv_path",
            }
            cell_fields = forms[index].get("cell_types", {}).get(str(cell_type), {})
            w = cell_fields.get(field_map[field])
        else:
            w = forms[index]["basic"].get(field)
        if w is None:
            return False
        if isinstance(w, QLineEdit):
            return w.text().strip() == str(proposed).strip()
        if isinstance(w, (QSpinBox, QDoubleSpinBox)):
            return w.value() == float(proposed)
        if isinstance(w, QComboBox):
            return w.currentText() == str(proposed)
    except Exception:
        pass
    return False


_SYSTEM_PRIMER = (
    "Answer as the BEHAV3D Assistant. Use researcher-facing labels, read the live "
    "interface context before asking for values, and keep answers concise."
)


_RESEARCHER_LABELS = {
    "pixel_distance_xy": "XY pixel size",
    "pixel_distance_z": "Z pixel size",
    "distance_unit": "distance unit",
    "raw_image_path": "raw image path",
    "dimension_order": "dimension order",
    "sample_name": "sample name",
    "time_interval": "time interval",
    "time_unit": "time unit",
    "position_t": "timepoint",
    "TrackID": "track ID",
    "track_id": "track ID",
    "cell_type": "cell type",
    "organoid_cells_across": "organoid size in cell widths",
    "minimum_lengths": "minimum track lengths",
}


def researcher_facing_text(text: str) -> str:
    """Replace context/schema field names that should never reach the chat UI."""
    result = str(text or "")
    for technical, label in _RESEARCHER_LABELS.items():
        result = result.replace(technical, label)
    return result


class _ActionCard(QFrame):
    """A confirm card for a single proposed action."""

    confirmed = Signal(object)   # ProposedAction
    dismissed = Signal(object)

    def __init__(self, action: ProposedAction, parent=None):
        super().__init__(parent)
        self.action = action
        self.setStyleSheet(
            "QFrame { background:#20303a; border:1px solid #3a5a6a; border-radius:4px; }"
        )
        lay = QVBoxLayout(self)
        lay.setContentsMargins(8, 6, 8, 6)
        lay.setSpacing(4)

        title = QLabel(action.preview or action.kind)
        title.setWordWrap(True)
        title.setStyleSheet("color:#e0f0ff; font-size:11px;")
        lay.addWidget(title)

        if not action.ok:
            warn = QLabel(action.message)
            warn.setWordWrap(True)
            warn.setStyleSheet("color:#ffb3b3; font-size:10px;")
            lay.addWidget(warn)

        row = QHBoxLayout()
        row.addStretch(1)
        if action.ok:
            btn_apply = QPushButton("Apply change")
            btn_apply.setStyleSheet("font-size:10px;")
            btn_apply.clicked.connect(lambda: self.confirmed.emit(self.action))
            row.addWidget(btn_apply)
        btn_dismiss = QPushButton("Keep current")
        btn_dismiss.setStyleSheet("font-size:10px;")
        btn_dismiss.clicked.connect(lambda: self.dismissed.emit(self.action))
        row.addWidget(btn_dismiss)
        lay.addLayout(row)


class AssistantDock(QWidget):
    """Right-side co-pilot chat dock for BEHAV3D."""

    def __init__(self, main_widget, parent=None):
        super().__init__(parent)
        self.main_widget = main_widget
        self.setMinimumWidth(300)
        self._history: list[dict] = []        # [{role, content}] sent to the model
        self._threads: list[tuple] = []       # keep (thread, worker) refs alive
        self._md_log: list[str] = []          # finalised transcript blocks (markdown)
        self._streaming_text = None           # in-progress assistant text, or None
        self._guided_flow_active: bool = False  # True while a step-by-step guide is running
        # Deterministic metadata wizard state (no LLM): _md_current is the question
        # awaiting an answer; _md_queue holds the remaining questions in this phase.
        self._md_current: dict | None = None
        self._md_queue: list[dict] = []
        self._md_phase: str | None = None
        self._current_intent: str = "free_form"
        self._request_active = False
        # Coalesce streamed-token re-renders. Re-rendering the whole transcript on
        # every token saturates the GUI thread on long answers (the napari window
        # freezes until the stream ends). This single-shot timer batches tokens into
        # at most one render per interval; the expensive full-document style walk is
        # applied only on the final render of the turn.
        self._render_timer = QTimer(self)
        self._render_timer.setSingleShot(True)
        self._render_timer.setInterval(90)
        self._render_timer.timeout.connect(self._render_streaming)

        # Cohesive dark-theme styling for the whole dock (matches napari surfaces).
        self.setStyleSheet("""
            AssistantDock { background-color: #2a2d30; }
            QTextBrowser {
                background-color: #232629;
                color: #e6e6e6;
                border: 1px solid #3a3f44;
                border-radius: 4px;
                padding: 10px;
                font-family: "Helvetica Neue", "Segoe UI", Arial, sans-serif;
                font-size: 13px;
                selection-background-color: #3d6b8a;
            }
            QLineEdit {
                background-color: #2f3338; color: #e6e6e6;
                border: 1px solid #3a3f44; border-radius: 4px;
                padding: 7px 10px; font-size: 13px;
            }
            QLineEdit:focus { border: 1px solid #5285a6; }
            QPushButton {
                background-color: #353a3f; color: #e6e6e6;
                border: 1px solid #454b50; border-radius: 4px;
                padding: 7px 10px; font-size: 12px;
            }
            QPushButton:hover { background-color: #3f464c; }
            QPushButton:pressed { background-color: #2c3035; }
        """)

        root = QVBoxLayout(self)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(8)

        # --- context bar --------------------------------------------------
        self.context_bar = QLabel("BEHAV3D Assistant")
        self.context_bar.setStyleSheet(
            "background:#2f3338; color:#9fc6e0; padding:7px 10px;"
            "border:1px solid #3a3f44; border-radius:4px; font-size:11px; font-weight:600;"
        )
        self.context_bar.setWordWrap(True)
        root.addWidget(self.context_bar)

        self.status_label = QLabel("")
        self.status_label.setStyleSheet("color:#9aa4aa; font-size:11px; padding:0 2px;")
        self.status_label.setVisible(False)
        root.addWidget(self.status_label)

        # --- transcript ---------------------------------------------------
        self.transcript = QTextBrowser()
        self.transcript.setOpenExternalLinks(True)
        self.transcript.setFrameShape(QFrame.NoFrame)
        self.transcript.document().setDocumentMargin(2)
        root.addWidget(self.transcript, stretch=1)

        # --- proposed-action tray ----------------------------------------
        self.action_tray = QWidget()
        self.action_tray_layout = QVBoxLayout(self.action_tray)
        self.action_tray_layout.setContentsMargins(0, 0, 0, 0)
        self.action_tray_layout.setSpacing(4)
        root.addWidget(self.action_tray)

        # --- quick actions (rebuilt on tab switch) ------------------------
        self.quick_widget = QWidget()
        self._quick_layout = QVBoxLayout(self.quick_widget)
        self._quick_layout.setContentsMargins(0, 0, 0, 0)
        self._quick_layout.setSpacing(4)
        root.addWidget(self.quick_widget)

        # --- input row ----------------------------------------------------
        input_row = QHBoxLayout()
        self.input = QLineEdit()
        self.input.setPlaceholderText("Ask about parameters, methods, your data…")
        self.input.returnPressed.connect(self._on_send)
        input_row.addWidget(self.input, stretch=1)
        self.btn_send = QPushButton("Send")
        self.btn_send.clicked.connect(self._on_send)
        input_row.addWidget(self.btn_send)
        root.addLayout(input_row)

        self.refresh_context_bar()
        self._greet()
        try:
            initial_idx = self.main_widget.tabs.currentIndex()
            self.main_widget.tabs.currentChanged.connect(self._on_tab_changed)
        except Exception:
            initial_idx = 0
        self._set_quick_buttons(initial_idx)

    # ------------------------------------------------------------------
    # Context bar
    # ------------------------------------------------------------------
    def refresh_context_bar(self):
        try:
            ctx = build_context(self.main_widget)
            self.context_bar.setText(context_summary_line(ctx))
        except Exception:
            self.context_bar.setText("BEHAV3D Assistant")

    def _greet(self):
        cfg = load_client_config()
        mode = "" if cfg.get("endpoint") else "  *(offline mode — no server configured)*"
        self._append_md(
            f"**BEHAV3D Assistant**{mode}\n\nAsk me which method or parameter "
            "values suit your data. I can also fill the forms in for you "
            "(you confirm first). Try *“Explain this screen”* to start."
        )

    # ------------------------------------------------------------------
    # Transcript helpers
    # ------------------------------------------------------------------
    def _render(self, style: bool = True):
        """Re-render the whole transcript as markdown.

        QTextCursor.insertMarkdown does not exist in PyQt5, so we keep a log of
        markdown blocks and render the document with QTextEdit.setMarkdown (which
        does render bold/italic/lists/code), then pin the scrollbar to the bottom.

        ``style=False`` skips the full-document style walk — used for the throttled
        renders during token streaming so a long answer doesn't freeze the GUI.
        """
        # Any direct/full render supersedes a pending throttled one.
        self._render_timer.stop()
        blocks = list(self._md_log)
        if self._streaming_text:
            blocks.append(f"**BEHAV3D Assistant**\n\n{self._streaming_text}")
        try:
            # A thin rule between turns gives clear visual separation.
            self.transcript.setMarkdown("\n\n---\n\n".join(blocks))
            if style:
                self._style_blocks()
        except Exception:
            # Fallback for any Qt build lacking setMarkdown.
            self.transcript.setPlainText("\n\n".join(blocks))
        sb = self.transcript.verticalScrollBar()
        sb.setValue(sb.maximum())

    def _render_streaming(self):
        """Throttled render during streaming — skips the expensive style walk,
        which is applied once when the turn finishes (_on_finished → _render)."""
        self._render(style=False)

    def _schedule_render(self):
        """Request a render soon, coalescing bursts of streamed tokens into at
        most one render per timer interval."""
        if not self._render_timer.isActive():
            self._render_timer.start()

    def _style_blocks(self):
        """setMarkdown renders paragraphs/lists very tightly. Walk the document
        and add line spacing + paragraph spacing so responses are readable."""
        try:
            from qtpy.QtGui import QTextCursor, QTextBlockFormat
        except Exception:
            return
        doc = self.transcript.document()
        block = doc.begin()
        while block.isValid():
            cur = QTextCursor(block)
            bf = cur.blockFormat()
            bf.setLineHeight(140, QTextBlockFormat.ProportionalHeight)  # 1.4× line spacing
            bf.setTopMargin(2.0)
            bf.setBottomMargin(9.0)                                     # gap between paragraphs
            cur.setBlockFormat(bf)
            block = block.next()

    def _append_md(self, markdown: str):
        self._md_log.append(markdown)
        self._render()

    def _append_user(self, text: str):
        self._append_md(f"**You**\n\n{text}")

    # ------------------------------------------------------------------
    # Sending
    # ------------------------------------------------------------------
    def ask_about(self, prompt: str):
        """Send ``prompt`` to the assistant on the user's behalf.

        Used by the Analysis tab's Guided "Ask the assistant" buttons. If a
        request is already in flight, the text is dropped into the input box
        instead of being lost, so the user can send it when ready.
        """
        prompt = (prompt or "").strip()
        if not prompt:
            return
        if self._request_active:
            try:
                self.input.setText(prompt)
            except Exception:
                pass
            return
        self._guided_flow_active = False
        self._send_message(prompt, intent="free_form")

    def _on_send(self):
        text = self.input.text().strip()
        if not text:
            return
        self.input.clear()
        # Deterministic metadata wizard owns the input while it's running (the
        # questions and answers never touch the LLM, so it cannot stall or loop).
        if self._md_current is not None:
            self._md_handle_answer(text)
            return
        # Short replies such as "APOC", "3 is fine", or "next" remain part of a
        # button-started guide. Longer questions return to normal free-form mode.
        guided_reply = self._guided_flow_active and len(text) <= 100
        intent = self._current_intent if guided_reply else "free_form"
        if not guided_reply:
            self._guided_flow_active = False
        self._send_message(text, intent=intent)

    def _send_message(
        self,
        text: str,
        *,
        show_as_user: bool = True,
        display_text: str | None = None,
        reset_guidance: bool = True,
        intent: str = "free_form",
    ):
        if self._request_active:
            return
        self._current_intent = intent
        if reset_guidance:
            # Any LLM-bound send (free-form or another guide button) cancels the
            # deterministic metadata wizard.
            self._md_current = None
            self._md_queue = []
            self._md_phase = None
        if show_as_user:
            self._append_user(display_text or text)
        self._history.append({"role": "user", "content": text})
        # Keep only the latest 12 user/assistant exchanges. Session facts and live
        # values are carried separately in structured context.
        self._history = self._history[-24:]
        self._set_busy(True)
        self._streaming_text = ""        # opens a live "Assistant:" block
        self._render()

        ctx = {}
        try:
            ctx = build_context(self.main_widget)
            ctx["assistant_session"] = {
                "intent": self._current_intent,
            }
        except Exception:
            ctx = {}

        messages = [{"role": "system", "content": _SYSTEM_PRIMER}] + self._history
        worker = ChatWorker(messages=messages, context=ctx, tools=TOOL_SCHEMA)
        thread = QThread()
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.token.connect(self._on_token)
        worker.tool_calls.connect(self._on_tool_calls)
        worker.degraded.connect(self._on_degraded)
        worker.error.connect(self._on_error)
        worker.finished.connect(thread.quit)
        worker.finished.connect(self._on_finished)
        thread.finished.connect(lambda: self._cleanup_thread(thread, worker))
        self._threads.append((thread, worker))
        thread.start()

    def _set_busy(self, busy: bool):
        self._request_active = bool(busy)
        self.btn_send.setEnabled(not busy)
        self.status_label.setText("Thinking..." if busy else "")
        self.status_label.setVisible(bool(busy))
        # Leave the input box editable while a turn is in flight so the user can
        # type their next answer without waiting for the stream to finish.
        for i in range(self._quick_layout.count()):
            row = self._quick_layout.itemAt(i)
            sub = row.layout() if row else None
            if sub:
                for j in range(sub.count()):
                    w = sub.itemAt(j).widget() if sub.itemAt(j) else None
                    if w:
                        w.setEnabled(not busy)

    # ------------------------------------------------------------------
    # Worker signal handlers
    # ------------------------------------------------------------------
    def _on_token(self, chunk: str):
        if self._streaming_text is None:
            self._streaming_text = ""
        self._streaming_text += chunk
        # The model sees structured field names in live context. Enforce the
        # researcher-facing vocabulary even if it repeats one despite prompting.
        self._streaming_text = researcher_facing_text(self._streaming_text)
        self._schedule_render()

    def _on_degraded(self, full_text: str):
        self._streaming_text = researcher_facing_text(full_text)
        self._render()

    def _on_error(self, message: str):
        # close any open streaming block, then show the error
        self._finalize_streaming()
        self._append_md(message)
        # Defensive: ensure the dock never stays stuck "thinking" if signal
        # ordering ever changes (worker.finished still clears busy too).
        self._set_busy(False)

    def _on_tool_calls(self, calls: list):
        try:
            ctx = build_context(self.main_widget)
            cards = ctx.get("step_schema", [])
            # include all cards so cross-step suggestions still validate
            from behav3d.napari._assistant_schema import flatten_config_to_cards
            cards = flatten_config_to_cards()
            params = getattr(self.main_widget.data_prep_tab, "behav3d_parameters", {})
            actions = build_actions(calls, cards, params,
                                    controls=(ctx.get("ui_state", {}) or {}).get("controls", []))
        except Exception:
            actions = []
        if actions and all(self._is_noop_action(action) for action in actions):
            # A model may still narrate a future action despite the live value
            # already matching. Replace that stale promise with the verified fact.
            self._streaming_text = self._noop_response(actions)
        if actions and self._streaming_text is not None:
            # Tool calls arrive after streamed prose. Finalise that prose before
            # appending action outcomes so the transcript preserves causality.
            self._finalize_streaming()
            self._render()
        auto_previews = []
        # Suppress per-widget pulse glows when applying a batch — a flood of fills
        # would otherwise schedule hundreds of QTimers and thrash stylesheets,
        # which can freeze the GUI. One context-bar refresh happens at the end.
        from behav3d.napari._assistant_actions import set_pulse_suppressed
        batch = len(actions) > 1
        if batch:
            set_pulse_suppressed(True)
        try:
            for act in actions:
                if self._is_noop_action(act):
                    continue
                if self._should_auto_apply(act):
                    if act.ok:
                        # Classify BEFORE applying (widget still holds pre-apply state).
                        _silent_auto = self._is_silent_auto_apply(act)
                        try:
                            ok = apply_action(self.main_widget, act)
                        except Exception:
                            ok = False
                        if ok:
                            if act.data.get("result_markdown"):
                                self._append_action_result(act.data["result_markdown"])
                            elif act.kind == "bulk_fill_metadata":
                                n = (act.data.get("sample_count")
                                     or len(act.data.get("samples", []) or []))
                                self._append_md(
                                    f"Filled the Metadata Builder for {n} "
                                    f"sample{'s' if n != 1 else ''}. Review the values "
                                    "and tell me anything you'd like to change."
                                )
                            # Only count as a meaningful change if it wasn't a structural
                            # trigger or a same-value no-op re-call from the model.
                            elif not _silent_auto:
                                auto_previews.append(act.preview)
                            if not batch:
                                self.refresh_context_bar()
                else:
                    self._add_action_card(act)
        finally:
            if batch:
                set_pulse_suppressed(False)
                self.refresh_context_bar()

        if auto_previews:
            self._append_md("Filled in: " + " · ".join(auto_previews))
        # A no-op means the proposed value was already present. Treat that as
        # confirmation from the user's previous answer and move forward silently
        # instead of rendering a scary "Already set" card.
        # Deliberately no hidden continuation request. One user or button input maps
        # to exactly one model request, even when actions were applied or were no-ops.

    def _is_noop_action(self, action: ProposedAction) -> bool:
        return bool(action.data.get("no_op"))

    def _noop_response(self, actions: list[ProposedAction]) -> str:
        details = []
        metadata_draft = False
        for action in actions:
            label = action.data.get("label") or "That field"
            value = action.data.get("value")
            details.append(f"{label} is already set to **{value}**")
            metadata_draft = metadata_draft or str(
                action.data.get("control_id") or ""
            ).startswith("metadata.")
        message = ". ".join(details) + ". No change was needed."
        if metadata_draft:
            message += (
                " The value is in the Metadata Builder draft; save the metadata "
                "when ready."
            )
        return message

    def _should_auto_apply(self, action: ProposedAction) -> bool:
        """Auto-apply by default; only show a confirm card if a value is already set."""
        if not action.ok:
            return False

        if action.kind == "navigate_to_step":
            return True  # navigation never overwrites anything

        if action.kind == "add_queue_step":
            return False  # adding to the queue is a deliberate action

        if action.kind == "bulk_fill_metadata":
            return False  # one confirmation covers the complete multi-field change

        if action.kind == "select_segmentation_method":
            # Apply only when it changes the current selection.
            seg = getattr(self.main_widget, "segmentation_tab", None)
            combo = getattr(seg, "method_combo", None) if seg is not None else None
            if combo is None:
                return False
            cur = combo.currentText()
            val = str(action.data.get("value", "")).strip()
            return not (cur == val or cur.startswith(val) or (val and val in cur))

        if action.kind == "fill_metadata_builder":
            field = action.data.get("field")
            # Structural triggers never overwrite anything.
            if field in ("open_builder", "configure_cell_types",
                         "create_sample_forms", "fill_down"):
                return True
            dp = getattr(self.main_widget, "data_prep_tab", None)
            idx = int(action.data.get("index", 0))
            cell_type = action.data.get("cell_type")
            if not _metadata_field_has_value(dp, field, idx, cell_type):
                return True   # blank field — auto-apply
            # Field already has a value; auto-apply only if the proposed value is the
            # same (bot redundantly repeating itself). Different value → confirm card.
            return _fill_value_matches_current(
                dp, field, idx, action.data.get("value"), cell_type
            )

        if action.kind == "set_parameter":
            try:
                from behav3d.napari._assistant_actions import get_by_dotted
                from behav3d.napari._assistant_schema import flatten_config_to_cards
                key = action.data.get("key")
                params = getattr(
                    self.main_widget.data_prep_tab, "behav3d_parameters", {})
                current = get_by_dotted(params, key, _UNSET)
                if current is _UNSET:
                    return True  # never been set
                for card in flatten_config_to_cards():
                    if card["key"] == key:
                        return current == card["default"]
            except Exception:
                pass
            return True

        if action.kind == "set_ui_value":
            old = action.data.get("old_value")
            # Empty text fields can be filled immediately. Existing values require
            # the explicit Apply change card.
            return old is None or (isinstance(old, str) and not old.strip())

        if action.kind in (
            "show_track_length_distribution", "open_result", "recommend_edt",
            "summarize_track_counts",
        ):
            return True

        if action.kind in ("create_cell_type_group", "create_btrack_config_copy"):
            return False

        return True  # unknown action kind: auto-apply

    def _is_silent_auto_apply(self, action: ProposedAction) -> bool:
        """True when the action auto-applies but should NOT count as a new user answer.

        Covers two cases:
        - Structural triggers (open_builder, configure_cell_types, etc.) — no user data.
        - Same-value no-ops — model re-confirming a value that is already set.
          Detected BEFORE applying, while the widget still holds its current value.
        """
        if action.kind != "fill_metadata_builder":
            return False
        field = action.data.get("field")
        if field in ("open_builder", "configure_cell_types",
                     "create_sample_forms", "fill_down"):
            return True
        # Same-value no-op: proposed equals current, checked pre-apply.
        # This covers both "field at a meaningful value" and "field still at default"
        # (e.g. n_samples=1 → 1==1) — both are no-ops that should not trigger auto-continue.
        dp = getattr(self.main_widget, "data_prep_tab", None)
        idx = int(action.data.get("index", 0))
        return _fill_value_matches_current(
            dp, field, idx, action.data.get("value"), action.data.get("cell_type")
        )

    def _finalize_streaming(self):
        """Move the in-progress assistant text into the finalised transcript log."""
        if self._streaming_text:
            self._streaming_text = researcher_facing_text(self._streaming_text)
            self._md_log.append(f"**BEHAV3D Assistant**\n\n{self._streaming_text}")
            self._history.append({"role": "assistant", "content": self._streaming_text})
        self._streaming_text = None

    def _on_finished(self):
        had_text = bool(self._streaming_text)
        self._finalize_streaming()
        self._render()
        self._set_busy(False)
        if had_text:
            self._append_log()

    def _cleanup_thread(self, thread, worker):
        self._threads = [(t, w) for (t, w) in self._threads if t is not thread]
        worker.deleteLater()
        thread.deleteLater()

    # ------------------------------------------------------------------
    # Proposed-action tray
    # ------------------------------------------------------------------
    def _add_action_card(self, action: ProposedAction):
        card = _ActionCard(action)
        card.confirmed.connect(self._apply_action)
        card.dismissed.connect(lambda a, c=card: self._remove_card(c))
        self.action_tray_layout.addWidget(card)

    def _remove_card(self, card: _ActionCard):
        card.setParent(None)
        card.deleteLater()

    def _append_action_result(self, markdown: str):
        """Display a deterministic action result and retain it for follow-ups."""
        self._append_md(markdown)
        self._history.append({"role": "assistant", "content": markdown})
        self._history = self._history[-24:]

    def _apply_action(self, action: ProposedAction):
        ok = False
        try:
            ok = apply_action(self.main_widget, action)
        except Exception as e:
            self._append_md(f"\nCould not apply: {e}")
        if ok:
            if action.data.get("result_markdown"):
                self._append_action_result(action.data["result_markdown"])
            elif action.kind == "set_parameter" and not action.data.get("widget_updated"):
                # stored in config, but no matching field exists on the current screen
                label = action.data.get("label") or humanize_parameter_key(action.data.get("key"))
                self._append_md(
                    f"Saved to config: **{label}** = "
                    f"`{action.data.get('value')!r}`\n\n*This step has no form field "
                    "here yet — it'll be used when you reach that part of the pipeline.*"
                )
            else:
                self._append_md(f"Filled in: {action.preview}")
            self.refresh_context_bar()
        else:
            self._append_md(f"Could not apply: {action.preview}")
        # remove all cards for this action
        for i in reversed(range(self.action_tray_layout.count())):
            w = self.action_tray_layout.itemAt(i).widget()
            if isinstance(w, _ActionCard) and w.action is action:
                self._remove_card(w)

    def _md_apply_structural(self, field: str) -> bool:
        """Apply a structural metadata-builder step (configure_cell_types /
        create_sample_forms) directly, without an LLM round-trip."""
        try:
            act = ProposedAction("fill_metadata_builder", field=field, value=None, index=0)
            return bool(apply_action(self.main_widget, act))
        except Exception:
            return False

    # ------------------------------------------------------------------
    # Quick actions
    # ------------------------------------------------------------------
    def _send_intent(self, intent: str, label: str):
        """Send a short visible command plus a structured intent in context."""
        self._send_message(label, display_text=label, intent=intent)

    def _explain_screen(self):
        self._send_intent("explain_screen", "Explain this screen")

    def _start_interview(self):
        self._send_intent("review_data", "Review my data")

    def _start_metadata_guide(self):
        self._append_user("Build metadata")
        self._begin_metadata_wizard()

    # ------------------------------------------------------------------
    # Deterministic metadata wizard (no LLM — cannot stall or loop)
    # ------------------------------------------------------------------
    _MD_COUNT_STEPS = [
        {"kind": "count", "field": "n_samples",
         "question": "How many samples (fields of view) do you have?"},
        {"kind": "count", "field": "n_organoids",
         "question": "How many organoid types do you have?"},
        {"kind": "count", "field": "n_immune",
         "question": "How many immune cell types do you have?"},
        {"kind": "count", "field": "n_other",
         "question": "How many other (non-organoid, non-immune) cell types do you have?"},
    ]

    def _begin_metadata_wizard(self):
        """Drive the structural metadata setup (sample/cell-type counts and names)
        entirely on the frontend. Questions are posted directly and the typed answer
        is applied immediately — the LLM is never in the loop, so it can't stall,
        loop, or skip a field."""
        self._guided_flow_active = True
        try:
            nav = ProposedAction("navigate_to_step", step="data_preparation")
            apply_action(self.main_widget, nav)
        except Exception:
            pass
        if not self._md_apply_structural("open_builder"):
            self._guided_flow_active = False
            self._append_md("Metadata Builder was not opened.")
            return
        self._append_md("Opened the Metadata Builder")
        self.refresh_context_bar()
        dp = getattr(self.main_widget, "data_prep_tab", None)
        if (getattr(dp, "_edit_mode", False)
                and bool(getattr(dp, "_sample_forms", []) or [])):
            self._guided_flow_active = False
            self._current_intent = "review_data"
            self._append_md(
                "Your loaded metadata is ready to edit. I can review the current "
                "sample values or help change a specific field."
            )
            return
        self._md_phase = "counts"
        self._md_queue = [dict(s) for s in self._MD_COUNT_STEPS]
        self._md_current = None
        self._md_ask_next()

    def _md_ask_next(self):
        """Post the next question in the current phase, or advance phases."""
        if self._md_queue:
            self._md_current = self._md_queue.pop(0)
            self._append_md(f"**BEHAV3D Assistant**\n\n{self._md_current['question']}")
            return
        if self._md_phase == "counts":
            # Counts done → build the cell-type name fields, then ask each name.
            self._md_apply_structural("configure_cell_types")
            self.refresh_context_bar()
            self._md_phase = "names"
            self._md_queue = self._md_build_name_steps()
            self._md_ask_next()
            return
        if self._md_phase == "names":
            # Names done → build the per-sample forms and hand off.
            self._md_apply_structural("create_sample_forms")
            self.refresh_context_bar()
            self._md_finish()
            return
        self._md_finish()

    def _md_build_name_steps(self) -> list[dict]:
        dp = getattr(self.main_widget, "data_prep_tab", None)
        steps: list[dict] = []
        if dp is None:
            return steps
        spec = [
            ("n_organoid_spin", "organoid_name", "organoid", 'e.g. "tumor"'),
            ("n_immune_spin", "immune_name", "immune cell", 'e.g. "CD8"'),
            ("n_other_spin", "other_name", "other cell", ""),
        ]
        for spin_attr, field, label, hint in spec:
            spin = getattr(dp, spin_attr, None)
            n = spin.value() if spin is not None else 0
            for i in range(n):
                eg = f" ({hint})" if hint else ""
                steps.append({"kind": "name", "field": field, "index": i,
                              "question": f"What is the name of {label} type {i + 1}?{eg}"})
        return steps

    def _md_handle_answer(self, text: str):
        """Apply the typed answer to the current wizard question and advance."""
        item = self._md_current
        self._append_user(text)
        if self._md_help_requested(text):
            self._append_md(
                f"**BEHAV3D Assistant**\n\n{self._md_help_text(item)}\n\n"
                f"{item['question']}"
            )
            return
        if item["kind"] == "count":
            try:
                value = int("".join(ch for ch in text if (ch.isdigit() or ch == "-")))
            except ValueError:
                # Not a number — treat as a free-form question: exit the wizard and
                # hand the message to the assistant instead of forcing a number.
                self._md_current = None
                self._md_queue = []
                self._md_phase = None
                self._send_message(text, show_as_user=False)
                return
        else:  # name
            value = text.strip()
            if not value:
                self._append_md(f"**BEHAV3D Assistant**\n\n{item['question']}")
                return  # re-ask; keep _md_current
        self._md_current = None
        self._md_apply_value(item["field"], value, item.get("index", 0))
        self.refresh_context_bar()
        self._md_ask_next()

    @staticmethod
    def _md_help_requested(text: str) -> bool:
        normalized = (text or "").strip().lower()
        return any(phrase in normalized for phrase in (
            "help", "not sure", "what does", "what counts", "how do i know",
        ))

    @staticmethod
    def _md_help_text(item: dict) -> str:
        help_by_field = {
            "n_samples": (
                "A sample is one independently imaged field of view or acquisition "
                "that will become one metadata row. Count image files or fields of "
                "view, not the number of cell types inside them."
            ),
            "n_organoids": (
                "Count the distinct organoid populations you want to segment and "
                "track separately. Use 0 if your experiment has no organoids."
            ),
            "n_immune": (
                "Count the distinct immune-cell populations you want BEHAV3D to "
                "process separately, such as T cells and macrophages."
            ),
            "n_other": (
                "Count any additional cell populations that are neither organoids "
                "nor immune cells. Use 0 if there are none."
            ),
        }
        return help_by_field.get(
            item.get("field"),
            "Use the short label you want to see for this cell population throughout BEHAV3D.",
        )

    def _md_apply_value(self, field: str, value, index: int = 0):
        from behav3d.napari._assistant_actions import _metadata_builder_preview
        act = ProposedAction("fill_metadata_builder", field=field, value=value, index=index)
        try:
            ok = apply_action(self.main_widget, act)
        except Exception:
            ok = False
        if ok:
            try:
                preview = _metadata_builder_preview(field, value, index)
            except Exception:
                preview = f"{field} = {value}"
            self._append_md(f"Filled in: {preview}")

    def _md_finish(self):
        self._md_current = None
        self._md_queue = []
        self._md_phase = None
        self._guided_flow_active = False
        self._append_md(
            "All set — the cell-type fields and per-sample forms are ready below. "
            "For each sample, fill in the **image path** and acquisition settings; use "
            "**Fill All from Sample 1** to copy shared values (pixel size, time interval), "
            "then **Save Metadata CSV** and click **Load Metadata**. Ask me about any "
            "field you're unsure of."
        )

    # ------------------------------------------------------------------
    # Tab-change awareness and dynamic quick buttons
    # ------------------------------------------------------------------
    def _on_tab_changed(self, index: int):
        """Called whenever the user switches tabs in the main widget."""
        self.refresh_context_bar()
        self._set_quick_buttons(index)

    def _set_quick_buttons(self, tab_index: int):
        """Clear and repopulate quick-action buttons for the given tab."""
        lay = self._quick_layout
        while lay.count():
            item = lay.takeAt(0)
            sub = item.layout()
            if sub:
                while sub.count():
                    si = sub.takeAt(0)
                    w = si.widget()
                    if w:
                        w.setParent(None)
                        w.deleteLater()

        def _btn(label, handler, tip=""):
            b = QPushButton(label)
            b.setStyleSheet("font-size:11px; padding:5px 8px;")
            if tip:
                b.setToolTip(tip)
            b.clicked.connect(handler)
            return b

        def _row(*btns):
            r = QHBoxLayout()
            for b in btns:
                r.addWidget(b)
            lay.addLayout(r)

        if tab_index == 0:   # data_preparation
            _row(_btn("Walk through setup", self._start_setup_guide),
                 _btn("Check what's missing", self._check_data_prereqs))
            _row(_btn("Build metadata", self._start_metadata_guide))
        elif tab_index == 2:  # segmentation
            _row(_btn("Guide segmentation", self._start_segmentation_guide),
                 _btn("Choose a method", self._explain_seg_methods))
            _row(_btn("Estimate EDT", self._estimate_edt))
        elif tab_index == 3:  # tracking
            _row(_btn("Guide tracking", self._start_tracking_guide),
                 _btn("Which method?", self._explain_tracking_methods))
        elif tab_index == 4:  # feature_extraction
            _row(_btn("Guide setup", self._start_feature_guide),
                 _btn("Check prerequisites", self._check_feature_prereqs))
        elif tab_index == 5:  # filtering
            _row(_btn("Review filters", self._review_filters),
                 _btn("Track lengths", self._show_track_lengths))
            _row(_btn("Cell counts", self._show_track_counts))
        elif tab_index == 6:  # analysis
            _row(_btn("Choose analysis", self._choose_analysis),
                 _btn("Review results", self._review_results))
        else:
            _row(_btn("Explain this screen", self._explain_screen),
                 _btn("Tell me about my data", self._start_interview))

    # ------------------------------------------------------------------
    # Guided-flow handlers for per-tab quick buttons
    # ------------------------------------------------------------------
    def _start_setup_guide(self):
        self._guided_flow_active = True
        ctx = {}
        try:
            ctx = build_context(self.main_widget)
        except Exception:
            pass
        if ctx.get("output_dir_set") and not ctx.get("metadata", {}).get("loaded"):
            # Output dir is already set → run the deterministic metadata wizard
            # (same robust path as the "Build metadata" button).
            self._append_user("Walk through setup")
            self._begin_metadata_wizard()
        else:
            self._send_intent("guide_data_setup", "Walk through setup")

    def _check_data_prereqs(self):
        self._send_intent("check_data_setup", "Check what's missing")

    def _start_segmentation_guide(self):
        self._guided_flow_active = True
        self._send_intent("guide_segmentation", "Guide segmentation")

    def _explain_seg_methods(self):
        self._send_intent("compare_segmentation_methods", "Choose a method")

    def _estimate_edt(self):
        ctx = build_context(self.main_widget)
        cell_type = ctx.get("active_cell_type")
        if not cell_type:
            groups = (ctx.get("metadata", {}).get("cell_types", {}) or {})
            candidates = []
            for category in ("organoid", "immune", "other"):
                candidates.extend(groups.get(category, []) or [])
            candidates = list(dict.fromkeys(str(value) for value in candidates))
            if len(candidates) == 1:
                cell_type = candidates[0]
        if not cell_type:
            self._send_intent("recommend_edt", "Estimate EDT")
            return
        self._append_user("Estimate EDT")
        action = ProposedAction(
            "recommend_edt", cell_type=cell_type, cell_diameter_um=10.0,
            organoid_cells_across=None,
        )
        action.preview = f"Calculate EDT starting values for {cell_type}"
        if apply_action(self.main_widget, action) and action.data.get("result_markdown"):
            self._append_action_result(action.data["result_markdown"])
        else:
            self._append_md(f"I could not calculate EDT values for **{cell_type}**.")

    def _start_tracking_guide(self):
        self._guided_flow_active = True
        self._send_intent("guide_tracking", "Guide tracking")

    def _explain_tracking_methods(self):
        self._send_intent("compare_tracking_methods", "Which method?")

    def _start_feature_guide(self):
        self._guided_flow_active = True
        self._send_intent("guide_feature_extraction", "Guide setup")

    def _check_feature_prereqs(self):
        self._send_intent("check_feature_prerequisites", "Check prerequisites")

    def _review_filters(self):
        self._send_intent("review_filters", "Review filters")

    def _show_track_lengths(self):
        ctx = build_context(self.main_widget)
        cell_type = ctx.get("active_cell_type")
        if not cell_type:
            self._append_md("Select a cell type in Filtering first.")
            return
        action = ProposedAction("show_track_length_distribution", cell_type=cell_type)
        action.preview = f"Show track-length distributions for {cell_type}"
        if not apply_action(self.main_widget, action):
            self._append_md(f"I could not open the track-length preview for **{cell_type}**.")

    def _show_track_counts(self):
        ctx = build_context(self.main_widget)
        cell_type = ctx.get("active_cell_type")
        if not cell_type:
            self._append_md("Select a cell type in Filtering first.")
            return
        self._append_user("Cell counts by minimum track length")
        action = ProposedAction(
            "summarize_track_counts", cell_type=cell_type, position_t=0,
            minimum_lengths=[20, 50, 100, 200],
        )
        action.preview = f"Summarize track counts for {cell_type}"
        if apply_action(self.main_widget, action) and action.data.get("result_markdown"):
            self._append_action_result(action.data["result_markdown"])
        else:
            self._append_md(f"I could not calculate track counts for **{cell_type}**.")

    def _choose_analysis(self):
        self._send_intent("choose_analysis", "Choose analysis")

    def _review_results(self):
        self._send_intent("review_results", "Review results")

    # ------------------------------------------------------------------
    # Session replay log
    # ------------------------------------------------------------------
    def _append_log(self):
        import os
        dp = getattr(self.main_widget, "data_prep_tab", None)
        out_dir = getattr(dp, "output_dir", "") if dp else ""
        if not out_dir or not self._history:
            return
        try:
            path = os.path.join(out_dir, "assistant_log.md")
            last = self._history[-1]
            with open(path, "a", encoding="utf-8") as f:
                f.write(f"\n**{last['role']}:** {last['content']}\n")
        except Exception:
            pass
