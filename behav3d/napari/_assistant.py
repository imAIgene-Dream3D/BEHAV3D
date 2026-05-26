"""
BEHAV3D assistant — the right-side co-pilot dock widget.

``AssistantDock`` is added as a *separate* napari dock (right area) and holds a
reference to the main ``BEHAV3DWidget`` so it can read live workflow context
(:func:`_assistant_context.build_context`) and apply proposed parameter changes
(:mod:`_assistant_actions`) after the user confirms.

Layout (top → bottom):
  • context bar          — 📍 step · N samples · output dir · M queued (live)
  • chat transcript      — QTextBrowser, markdown, streamed assistant tokens
  • proposed-action tray — confirm cards: "old → new  [Fill it in ✨] [Dismiss]"
  • quick actions        — Explain this screen · Tell me about my data
  • input row            — text box + Send

All network I/O happens on a QThread worker (``ChatWorker``); the GUI never blocks.
"""
from __future__ import annotations

from typing import Optional

from qtpy.QtCore import Qt, QThread, Signal
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QLineEdit,
    QTextBrowser, QFrame, QSizePolicy,
)

from behav3d.napari._assistant_context import build_context, context_summary_line
from behav3d.napari._assistant_actions import (
    build_actions, apply_action, TOOL_SCHEMA, ProposedAction,
)
from behav3d.napari._assistant_client import ChatWorker, load_client_config


_SYSTEM_PRIMER = (
    "You are the BEHAV3D co-pilot, embedded in the napari plugin for analysing "
    "cell behaviour in 3D fluorescent imaging. Help the user choose methods and "
    "parameter values for their data. Ground answers in the provided step_schema "
    "and retrieved docs. When you recommend a concrete value, also emit a "
    "set_parameter tool call (the user confirms before it applies).\n\n"
    "FORMATTING — keep responses easy to read in a narrow side panel:\n"
    "- Keep it short: lead with a one-sentence answer, then detail only if useful.\n"
    "- Use short paragraphs (1–2 sentences) separated by a blank line.\n"
    "- Use bullet points or a numbered list for any set of items or steps.\n"
    "- Bold key terms and parameter names; use `code` for values/keys.\n"
    "- Avoid long dense paragraphs and walls of text."
)


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
            warn = QLabel("⚠️ " + action.message)
            warn.setWordWrap(True)
            warn.setStyleSheet("color:#ffb3b3; font-size:10px;")
            lay.addWidget(warn)

        row = QHBoxLayout()
        row.addStretch(1)
        if action.ok:
            btn_apply = QPushButton("Fill it in ✨")
            btn_apply.setStyleSheet("font-size:10px;")
            btn_apply.clicked.connect(lambda: self.confirmed.emit(self.action))
            row.addWidget(btn_apply)
        btn_dismiss = QPushButton("Dismiss")
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

        # Cohesive dark-theme styling for the whole dock (matches napari surfaces).
        self.setStyleSheet("""
            AssistantDock { background-color: #2a2d30; }
            QTextBrowser {
                background-color: #232629;
                color: #e6e6e6;
                border: 1px solid #3a3f44;
                border-radius: 8px;
                padding: 12px;
                font-family: -apple-system, "Segoe UI", "Helvetica Neue", Arial, sans-serif;
                font-size: 13px;
                selection-background-color: #3d6b8a;
            }
            QLineEdit {
                background-color: #2f3338; color: #e6e6e6;
                border: 1px solid #3a3f44; border-radius: 8px;
                padding: 7px 10px; font-size: 13px;
            }
            QLineEdit:focus { border: 1px solid #5285a6; }
            QPushButton {
                background-color: #353a3f; color: #e6e6e6;
                border: 1px solid #454b50; border-radius: 8px;
                padding: 7px 10px; font-size: 12px;
            }
            QPushButton:hover { background-color: #3f464c; }
            QPushButton:pressed { background-color: #2c3035; }
        """)

        root = QVBoxLayout(self)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(8)

        # --- context bar --------------------------------------------------
        self.context_bar = QLabel("📍 …")
        self.context_bar.setStyleSheet(
            "background:#2f3338; color:#9fc6e0; padding:7px 10px;"
            "border:1px solid #3a3f44; border-radius:8px; font-size:11px; font-weight:600;"
        )
        self.context_bar.setWordWrap(True)
        root.addWidget(self.context_bar)

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

        # --- quick actions -----------------------------------------------
        quick = QHBoxLayout()
        self.btn_explain = QPushButton("💡 Explain this screen")
        self.btn_explain.setStyleSheet("font-size:10px;")
        self.btn_explain.clicked.connect(self._explain_screen)
        quick.addWidget(self.btn_explain)
        self.btn_interview = QPushButton("🧬 Tell me about my data")
        self.btn_interview.setStyleSheet("font-size:10px;")
        self.btn_interview.clicked.connect(self._start_interview)
        quick.addWidget(self.btn_interview)
        root.addLayout(quick)

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

    # ------------------------------------------------------------------
    # Context bar
    # ------------------------------------------------------------------
    def refresh_context_bar(self):
        try:
            ctx = build_context(self.main_widget)
            self.context_bar.setText(context_summary_line(ctx))
        except Exception:
            self.context_bar.setText("📍 BEHAV3D Assistant")

    def _greet(self):
        cfg = load_client_config()
        mode = "" if cfg.get("endpoint") else "  *(offline mode — no server configured)*"
        self._append_md(
            f"**🤖 BEHAV3D Assistant**{mode}\n\nAsk me which method or parameter "
            "values suit your data. I can also fill the forms in for you "
            "(you confirm first). Try *“Explain this screen”* to start."
        )

    # ------------------------------------------------------------------
    # Transcript helpers
    # ------------------------------------------------------------------
    def _render(self):
        """Re-render the whole transcript as markdown.

        QTextCursor.insertMarkdown does not exist in PyQt5, so we keep a log of
        markdown blocks and render the document with QTextEdit.setMarkdown (which
        does render bold/italic/lists/code), then pin the scrollbar to the bottom.
        """
        blocks = list(self._md_log)
        if self._streaming_text is not None:
            blocks.append(f"**🤖 Assistant**\n\n{self._streaming_text or '…'}")
        try:
            # A thin rule between turns gives clear visual separation.
            self.transcript.setMarkdown("\n\n---\n\n".join(blocks))
            self._style_blocks()
        except Exception:
            # Fallback for any Qt build lacking setMarkdown.
            self.transcript.setPlainText("\n\n".join(blocks))
        sb = self.transcript.verticalScrollBar()
        sb.setValue(sb.maximum())

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
        self._append_md(f"**🧑 You**\n\n{text}")

    # ------------------------------------------------------------------
    # Sending
    # ------------------------------------------------------------------
    def _on_send(self):
        text = self.input.text().strip()
        if not text:
            return
        self.input.clear()
        self._send_message(text)

    def _send_message(self, text: str):
        self._append_user(text)
        self._history.append({"role": "user", "content": text})
        self._set_busy(True)
        self._streaming_text = ""        # opens a live "Assistant:" block
        self._render()

        ctx = {}
        try:
            ctx = build_context(self.main_widget)
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
        self.btn_send.setEnabled(not busy)
        self.input.setEnabled(not busy)

    # ------------------------------------------------------------------
    # Worker signal handlers
    # ------------------------------------------------------------------
    def _on_token(self, chunk: str):
        if self._streaming_text is None:
            self._streaming_text = ""
        self._streaming_text += chunk
        self._render()

    def _on_degraded(self, full_text: str):
        self._streaming_text = full_text
        self._render()

    def _on_error(self, message: str):
        # close any open streaming block, then show the error
        self._finalize_streaming()
        self._append_md(f"⚠️ {message}")

    def _on_tool_calls(self, calls: list):
        try:
            ctx = build_context(self.main_widget)
            cards = ctx.get("step_schema", [])
            # include all cards so cross-step suggestions still validate
            from behav3d.napari._assistant_schema import flatten_config_to_cards
            cards = flatten_config_to_cards()
            params = getattr(self.main_widget.data_prep_tab, "behav3d_parameters", {})
            actions = build_actions(calls, cards, params)
        except Exception:
            actions = []
        for act in actions:
            self._add_action_card(act)

    def _finalize_streaming(self):
        """Move the in-progress assistant text into the finalised transcript log."""
        if self._streaming_text:
            self._md_log.append(f"**🤖 Assistant**\n\n{self._streaming_text}")
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

    def _apply_action(self, action: ProposedAction):
        ok = False
        try:
            ok = apply_action(self.main_widget, action)
        except Exception as e:
            self._append_md(f"\n⚠️ Could not apply: {e}")
        if ok:
            if action.kind == "set_parameter" and not action.data.get("widget_updated"):
                # stored in config, but no matching field exists on the current screen
                self._append_md(
                    f"✅ Saved to config: `{action.data.get('key')}` = "
                    f"`{action.data.get('value')!r}`\n\n*This step has no form field "
                    "here yet — it'll be used when you reach that part of the pipeline.*"
                )
            else:
                self._append_md(f"✅ Filled in: {action.preview}")
            self.refresh_context_bar()
        else:
            self._append_md(f"⚠️ Couldn't apply: {action.preview}")
        # remove all cards for this action
        for i in reversed(range(self.action_tray_layout.count())):
            w = self.action_tray_layout.itemAt(i).widget()
            if isinstance(w, _ActionCard) and w.action is action:
                self._remove_card(w)

    # ------------------------------------------------------------------
    # Quick actions
    # ------------------------------------------------------------------
    def _explain_screen(self):
        ctx = {}
        try:
            ctx = build_context(self.main_widget)
        except Exception:
            pass
        step = ctx.get("current_tab_label") or ctx.get("current_step", "this step")
        self._send_message(
            f"Explain the {step} tab: what does it do, and what do its key "
            "parameters mean for someone configuring it for the first time?"
        )

    def _start_interview(self):
        self._send_message(
            "I'd like help configuring BEHAV3D for my data. Ask me a short series "
            "of questions about my experiment (imaging modality, cell types, frame "
            "interval, pixel size, expected cell diameter, motility) and then "
            "recommend sensible default parameters for the current step."
        )

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
