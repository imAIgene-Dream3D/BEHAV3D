"""
BEHAV3D assistant — embeddings + tiny vector index.

A deliberately dependency-light retrieval layer: sentence-transformers for
embeddings and an in-memory numpy matrix for cosine search. The index is small
(repo docs + ~195 parameter cards), so a flat numpy search is plenty fast and
avoids a FAISS/Chroma build dependency.

The index is persisted as two files in a Modal Volume:
    index.npy   float32 [N, D]   (L2-normalised embeddings)
    index.json  list[dict]       (chunk metadata: text, source, title)
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass

import numpy as np

EMBED_MODEL = "BAAI/bge-small-en-v1.5"   # 384-dim, fast, CPU-friendly


@dataclass
class Chunk:
    text: str
    source: str        # e.g. "README.md", "schema:tracking", "docstring:_tracking"
    title: str         # short human label / citation anchor


class Embedder:
    """Lazy-loaded sentence-transformers wrapper."""

    def __init__(self, model_name: str = EMBED_MODEL):
        self.model_name = model_name
        self._model = None

    @property
    def model(self):
        if self._model is None:
            from sentence_transformers import SentenceTransformer
            self._model = SentenceTransformer(self.model_name)
        return self._model

    def encode(self, texts: list[str]) -> np.ndarray:
        vecs = self.model.encode(texts, normalize_embeddings=True,
                                 convert_to_numpy=True, show_progress_bar=False)
        return vecs.astype("float32")


class VectorIndex:
    def __init__(self, matrix: np.ndarray | None = None, chunks: list[dict] | None = None):
        self.matrix = matrix
        self.chunks = chunks or []

    # -- persistence --------------------------------------------------------
    def save(self, directory: str) -> None:
        os.makedirs(directory, exist_ok=True)
        np.save(os.path.join(directory, "index.npy"), self.matrix)
        with open(os.path.join(directory, "index.json"), "w", encoding="utf-8") as f:
            json.dump(self.chunks, f)

    @classmethod
    def load(cls, directory: str) -> "VectorIndex":
        matrix = np.load(os.path.join(directory, "index.npy"))
        with open(os.path.join(directory, "index.json"), encoding="utf-8") as f:
            chunks = json.load(f)
        return cls(matrix=matrix, chunks=chunks)

    # -- build & search -----------------------------------------------------
    @classmethod
    def build(cls, chunks: list[Chunk], embedder: Embedder) -> "VectorIndex":
        texts = [c.text for c in chunks]
        matrix = embedder.encode(texts) if texts else np.zeros((0, 384), "float32")
        meta = [{"text": c.text, "source": c.source, "title": c.title} for c in chunks]
        return cls(matrix=matrix, chunks=meta)

    def search(self, query_vec: np.ndarray, k: int = 6) -> list[dict]:
        if self.matrix is None or len(self.chunks) == 0:
            return []
        # cosine similarity == dot product on normalised vectors
        scores = self.matrix @ query_vec.reshape(-1)
        top = np.argsort(-scores)[:k]
        out = []
        for i in top:
            c = dict(self.chunks[int(i)])
            c["score"] = float(scores[int(i)])
            out.append(c)
        return out
