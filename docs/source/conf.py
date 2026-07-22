"""Sphinx configuration for the BEHAV3D wiki."""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.abspath("../.."))

# ── Project information ────────────────────────────────────────────────────
project = "BEHAV3D EXPLORER"
author = "imAIgene-Dream3D"
copyright = "2026, imAIgene-Dream3D"
release = "0.1.0"

# ── General configuration ──────────────────────────────────────────────────
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "myst_parser",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinxcontrib.mermaid",
]

# Mock heavy / GUI / GPU dependencies so autodoc can import the package on a
# minimal CI runner that does not have cellpose, napari, torch, etc.
autodoc_mock_imports = [
    "napari",
    "qtpy",
    "magicgui",
    "cellpose",
    "torch",
    "torchvision",
    "torchaudio",
    "dask",
    "zarr",
    "btrack",
    "trackpy",
    "apoc",
    "pyclesperanto_prototype",
    "convpaint",
    "scanpy",
    "anndata",
    "umap",
    "igraph",
    "leidenalg",
    "dtaidistance",
    "sklearn",
    "skimage",
    "tifffile",
    "imageio",
    "imageio_ffmpeg",
    "joblib",
    "yaml",
    "pandas",
    "numpy",
    "scipy",
    "matplotlib",
    "seaborn",
    "plotly",
    "h5py",
    "czifile",
    "aicspylibczi",
    "readlif",
    "ilastik",
]

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "tasklist",
    "substitution",
    "attrs_inline",
    "attrs_block",
    "dollarmath",
]

# Auto-generate anchor slugs for headings (h1–h4) so in-page and cross-page
# "#heading-slug" links resolve. Without this, MyST does not create heading
# targets and those links warn as missing.
myst_heading_anchors = 4

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

master_doc = "index"
language = "en"

templates_path = ["_templates"]
exclude_patterns: list[str] = []

# Be lenient about cross-references while content is still landing.
nitpicky = False

# ── HTML output ────────────────────────────────────────────────────────────
html_theme = "sphinx_rtd_theme"
html_title = "BEHAV3D EXPLORER Wiki"
html_static_path = ["_static"]
html_css_files = ["custom.css"]

html_theme_options = {
    "navigation_depth": 4,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "titles_only": True,
}

# Public URL where the wiki is published once GitHub Pages is enabled.
# Used by Sphinx for sitemap / canonical link tags.
html_baseurl = "https://imAIgene-Dream3D.github.io/BEHAV3D/"

# ── Napoleon (Google-style docstrings) ─────────────────────────────────────
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False

# ── Intersphinx ────────────────────────────────────────────────────────────
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
}
