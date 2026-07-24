# BEHAV3D EXPLORER Wiki

**BEHAV3 EXPLORER** is a Python toolkit for the behavioural analysis of cells in 3D fluorescent time-lapse microscopy. It takes you all the way from raw confocal stacks to behavioural phenotypes: segmentation → tracking → feature extraction → filtering → analysis.

Although the pipeline can also be run from Python (see `notebooks/run_behav3d.ipynb` for advanced users who want to script steps or explore the code), **this wiki documents the napari GUI** — the graphical interface for day-to-day use in the lab.

## How you use BEHAV3D EXPLORER

BEHAV3D EXPLORER provides a full graphical interface through the napari viewer. To use it, launch the **napari GUI** (see [](getting_started) for setup) and work through the dock widget's tabs from left to right:


```{mermaid}
flowchart LR
    subgraph Pipeline ["BEHAV3D EXPLORER"]
        direction LR
        DP["📋 Data Preparation"] --> Seg["🦠 Segmentation"] --> Trk["📍 Tracking"] --> FE["🧪 Feature Extraction"] --> Filt["🧹 Filtering"] --> Ana["📊 Analysis"]
    end
    Vz["👁 Visualization"] -.inspect.-> Pipeline
    PQ["🛒 Processing Queue"] -.batches.-> Pipeline
```

The **Visualization** and the **Processing Queue** work across the whole pipeline: Visualization lets you visualize any sample and its outputs at any stage of the pipeline, and the Processing Queue batches segmentation / tracking / feature-extraction / filtering / analysis steps to run sequentially across all your samples. **Backprojection** — painting analysis results back onto the raw images — is built into the Single Cell workflows as their final step.

```{note}
The **Analysis** tab is implemented end-to-end: Death Dynamics, Interaction Analysis, and Single Cell (both behavioural-**state** and **track** classification), each Single Cell workflow ending in a **Backprojection** step.
```

## Where to next

```{grid} 2
:gutter: 3

::::{grid-item-card} 🚀 Installation
:link: installation
:link-type: doc

Auto installer for Windows / macOS / Linux, manual conda/mamba steps, GPU/CPU PyTorch options, troubleshooting.
::::

::::{grid-item-card} 🧭 Getting Started
:link: getting_started
:link-type: doc

Launch the plugin, take the 5-minute tour, learn the dock widget layout.
::::

::::{grid-item-card} 🧰 Plugin Essentials
:link: plugin_essentials/index
:link-type: doc

Visualization, the Processing Queue, and the canonical output folder layout — concepts that span every step.
::::

::::{grid-item-card} 🔬 Preprocessing
:link: preprocessing/index
:link-type: doc

Data preparation: building the metadata table, setting up the output folder, converting raw images to Zarr for fast access.
::::

::::{grid-item-card} ⚙️ Processing
:link: processing/index
:link-type: doc

Segmentation (6 methods, including zero-shot Cellpose-SAM) and Tracking (Propagation, Reporter Propagation, btrack, Import) — the heavy computation that turns images into trackable cells.
::::

::::{grid-item-card} 📈 Analysis
:link: analysis/index
:link-type: doc

Feature extraction, filtering, death dynamics, interaction analysis, behavioural-state classification, and backprojection.
::::

::::{grid-item-card} 🧠 Cellpose Training
:link: cellpose_training
:link-type: doc

How to train your own Cellpose model with the `train_behav3d_cellpose.ipynb` notebook when the pretrained models aren't enough.
::::

```

## Table of contents

```{toctree}
:caption: Introduction
:maxdepth: 2

installation
getting_started
```

```{toctree}
:caption: Plugin Essentials
:maxdepth: 2

plugin_essentials/index
```

```{toctree}
:caption: Pipeline
:maxdepth: 3

preprocessing/index
processing/index
analysis/index
```

```{toctree}
:caption: Tools
:maxdepth: 2

cellpose_training
recommended_settings
```

```{toctree}
:caption: Reference
:maxdepth: 2

api/index
```

## Indices and tables

- {ref}`genindex`
- {ref}`modindex`
- {ref}`search`
