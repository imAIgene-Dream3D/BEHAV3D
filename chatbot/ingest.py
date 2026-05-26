"""
BEHAV3D assistant — knowledge-base ingestion / RAG index builder.

Collects the BEHAV3D knowledge sources, chunks them, embeds them, and writes a
:class:`VectorIndex` to the Modal Volume that the chat endpoint reads at query
time.

Sources (per the agreed design):
  * README.md                         — workflow & method guidance
  * WIKI.md                           — curated guide (indexed automatically if present)
  * schema_cards.json                 — the 195 parameter cards (decoupled snapshot)
  * docstrings of key behav3d modules — extracted with ``ast`` (no heavy import)

The chunking functions are pure Python and unit-testable without Modal.
"""
from __future__ import annotations

import ast
import json
import os
from pathlib import Path

from embeddings import Chunk, Embedder, VectorIndex

_MAX_CHARS = 1100   # soft cap per chunk


# ---------------------------------------------------------------------------
# Pure chunkers
# ---------------------------------------------------------------------------
def chunk_markdown(text: str, source: str) -> list[Chunk]:
    """Split markdown into heading-anchored chunks, sub-splitting long sections."""
    chunks: list[Chunk] = []
    current_title = source
    buf: list[str] = []

    def flush():
        body = "\n".join(buf).strip()
        if not body:
            return
        # sub-split overly long sections on blank lines
        if len(body) <= _MAX_CHARS:
            chunks.append(Chunk(text=f"# {current_title}\n{body}", source=source, title=current_title))
            return
        piece, size = [], 0
        for para in body.split("\n\n"):
            if size + len(para) > _MAX_CHARS and piece:
                chunks.append(Chunk(text=f"# {current_title}\n" + "\n\n".join(piece),
                                    source=source, title=current_title))
                piece, size = [], 0
            piece.append(para)
            size += len(para)
        if piece:
            chunks.append(Chunk(text=f"# {current_title}\n" + "\n\n".join(piece),
                                source=source, title=current_title))

    for line in text.splitlines():
        if line.lstrip().startswith("#"):
            flush()
            buf = []
            current_title = line.lstrip("# ").strip() or source
        else:
            buf.append(line)
    flush()
    return chunks


def cards_to_chunks(cards: list[dict]) -> list[Chunk]:
    """One chunk per parameter card — the precise, authoritative parameter text."""
    out = []
    for c in cards:
        choices = f" Allowed values: {c['choices']}." if c.get("choices") else ""
        text = (
            f"Parameter: {c['key']}\n"
            f"Step: {c['step']}"
            + (f" (cell category: {c['category']})" if c.get("category") else "")
            + f"\nType: {c['type']}. Default: {c['default']!r}.{choices}\n"
            f"{c['description']}"
        )
        out.append(Chunk(text=text, source=f"schema:{c['step']}", title=c["key"]))
    return out


def python_docstring_chunks(py_path: Path) -> list[Chunk]:
    """Extract module/class/function docstrings from a .py file via ast."""
    try:
        tree = ast.parse(py_path.read_text(encoding="utf-8"))
    except Exception:
        return []
    out: list[Chunk] = []
    mod_doc = ast.get_docstring(tree)
    if mod_doc:
        out.append(Chunk(text=mod_doc, source=f"docstring:{py_path.name}", title=py_path.stem))
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            doc = ast.get_docstring(node)
            if doc and len(doc) > 60:   # skip trivial one-liners
                out.append(Chunk(text=f"{node.name}: {doc}",
                                 source=f"docstring:{py_path.name}", title=node.name))
    return out


# ---------------------------------------------------------------------------
# Collection
# ---------------------------------------------------------------------------
_DOC_PY_MODULES = [
    "behav3d/preprocessing/segmentation/cellpose_prediction.py",
    "behav3d/preprocessing/segmentation/apoc_segment.py",
    "behav3d/preprocessing/tracking/btrack_tracking.py",
    "behav3d/preprocessing/tracking/laptracking.py",
    "behav3d/features/timepoint_features.py",
    "behav3d/analysis/interaction_analysis.py",
]


def collect_chunks(repo_root: Path, schema_cards_path: Path) -> list[Chunk]:
    chunks: list[Chunk] = []

    for md_name in ("README.md", "WIKI.md"):
        p = repo_root / md_name
        if p.exists():
            chunks += chunk_markdown(p.read_text(encoding="utf-8"), md_name)

    if schema_cards_path.exists():
        cards = json.loads(schema_cards_path.read_text(encoding="utf-8"))
        chunks += cards_to_chunks(cards)

    for rel in _DOC_PY_MODULES:
        p = repo_root / rel
        if p.exists():
            chunks += python_docstring_chunks(p)

    return chunks


# ---------------------------------------------------------------------------
# Modal entrypoint (imported by app.py; safe to import without modal installed)
# ---------------------------------------------------------------------------
def build_and_save(repo_root: str, schema_cards_path: str, out_dir: str) -> int:
    chunks = collect_chunks(Path(repo_root), Path(schema_cards_path))
    index = VectorIndex.build(chunks, Embedder())
    index.save(out_dir)
    return len(chunks)


if __name__ == "__main__":
    # Local dry-run: build against the repo without embedding (counts only).
    here = Path(__file__).resolve().parent
    repo = here.parent
    cs = collect_chunks(repo, here / "schema_cards.json")
    print(f"Collected {len(cs)} chunks")
    from collections import Counter
    print("by source-kind:", dict(Counter(c.source.split(':')[0] for c in cs)))
