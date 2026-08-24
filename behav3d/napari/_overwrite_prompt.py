"""
BEHAV3D napari plugin – shared overwrite-prompt helper.

Provides a single ``prompt_overwrite`` function used by every tab so all
Overwrite / Skip / Cancel dialogs have the same body, default button and
return values.

Return values
-------------
The function returns one of the string tokens:

- ``"overwrite"`` — user wants to overwrite the existing data
- ``"skip"``      — user wants to skip the existing data (batch only)
- ``"cancel"``    — user dismissed the dialog
- a custom token  — when the caller passes ``extra_buttons`` and the user
  picked one of them (the token is the second element of the tuple).
"""
from __future__ import annotations

from typing import Iterable, Optional, Sequence, Tuple

from qtpy.QtWidgets import QMessageBox, QWidget


def prompt_overwrite(
    parent: Optional[QWidget],
    title: str,
    items: Sequence[str],
    *,
    body_prefix: str = "The following data already exists:",
    body_suffix: str = "What do you want to do?",
    overwrite_label: str = "Overwrite",
    skip_label: str = "Skip",
    allow_skip: bool = True,
    extra_buttons: Optional[Iterable[Tuple[str, str, QMessageBox.ButtonRole]]] = None,
) -> str:
    """Show a uniform overwrite-prompt dialog.

    Parameters
    ----------
    parent : QWidget or None
        Parent widget for modal anchoring.
    title : str
        Window title (e.g. "Overwrite Existing Features?").
    items : sequence of str
        Bullet-listed names of the existing artefacts that would be
        overwritten. Shown as ``"  • item"`` lines between
        ``body_prefix`` and ``body_suffix``.
    body_prefix, body_suffix : str
        Surrounding text around the bullet list.
    overwrite_label, skip_label : str
        Labels for the destructive / accept buttons. Use
        ``"Overwrite All"`` / ``"Skip Existing"`` for batch dialogs.
    allow_skip : bool
        When False, the Skip button is hidden (single-cell flows often
        treat Skip the same as Cancel; some callers prefer to remove
        the redundancy explicitly).
    extra_buttons : iterable of (label, token, role), optional
        Additional buttons inserted between Skip and Cancel. ``token``
        is what ``prompt_overwrite`` returns when that button is
        picked. ``role`` is a ``QMessageBox.ButtonRole`` controlling
        Qt's button ordering / default styling.

    Returns
    -------
    str
        One of ``"overwrite"``, ``"skip"``, ``"cancel"`` or one of the
        custom tokens supplied via ``extra_buttons``.
    """
    details = "\n".join(f"  • {it}" for it in items) if items else ""
    text_parts = [body_prefix]
    if details:
        text_parts.append("")
        text_parts.append(details)
    text_parts.append("")
    text_parts.append(body_suffix)

    box = QMessageBox(parent)
    box.setWindowTitle(title)
    box.setIcon(QMessageBox.Warning)
    box.setText("\n".join(text_parts))

    btn_overwrite = box.addButton(overwrite_label, QMessageBox.DestructiveRole)

    extra_mapped: list[tuple[object, str]] = []
    if extra_buttons:
        for label, token, role in extra_buttons:
            btn = box.addButton(label, role)
            extra_mapped.append((btn, token))

    btn_skip = box.addButton(skip_label, QMessageBox.AcceptRole) if allow_skip else None
    btn_cancel = box.addButton("Cancel", QMessageBox.RejectRole)
    box.setDefaultButton(btn_cancel)

    box.exec_()
    clicked = box.clickedButton()

    if clicked is btn_overwrite:
        return "overwrite"
    if btn_skip is not None and clicked is btn_skip:
        return "skip"
    for btn, token in extra_mapped:
        if clicked is btn:
            return token
    return "cancel"


def prompt_overwrite_single(
    parent: Optional[QWidget],
    title: str,
    items: Sequence[str],
    *,
    skip_label: str = "Skip",
    body_suffix: str = "What do you want to do?",
    extra_buttons: Optional[Iterable[Tuple[str, str, QMessageBox.ButtonRole]]] = None,
) -> str:
    """Convenience wrapper for single-cell-type runs.

    Uses ``"Overwrite"`` and ``"Skip"`` labels. Callers that can tell some of the
    listed outputs are incomplete pass ``skip_label="Skip & Resume"``, so the button
    says what it will actually do.
    """
    return prompt_overwrite(
        parent,
        title,
        items,
        overwrite_label="Overwrite",
        skip_label=skip_label,
        body_suffix=body_suffix,
        extra_buttons=extra_buttons,
    )


def prompt_overwrite_batch(
    parent: Optional[QWidget],
    title: str,
    items: Sequence[str],
    *,
    skip_label: str = "Skip Existing",
    body_prefix: str = "The following data already exists:",
    body_suffix: str = "What do you want to do?",
    extra_buttons: Optional[Iterable[Tuple[str, str, QMessageBox.ButtonRole]]] = None,
) -> str:
    """Convenience wrapper for batch (multi-cell-type) runs.

    Uses ``"Overwrite All"`` and ``"Skip Existing"`` labels. Callers that can tell
    some of the listed outputs are incomplete pass ``skip_label="Skip & Resume"``,
    so the button says what it will actually do.
    """
    return prompt_overwrite(
        parent,
        title,
        items,
        overwrite_label="Overwrite All",
        skip_label=skip_label,
        body_prefix=body_prefix,
        body_suffix=body_suffix,
        extra_buttons=extra_buttons,
    )
