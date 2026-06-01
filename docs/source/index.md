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

The **Visualization** and the **Processing Queue** work across the whole pipeline: Visualization lets you visualize any sample and its outputs at any stage of the pipeline, and the Processing Queue batches segmentation / tracking / feature-extraction / filtering / analysis steps to run sequentially across all your samples.


```{note}
The **Analysis** tab is currently a stub in the plugin ("🚧 Coming soon"). The corresponding pages of this wiki are skeletons that will be filled in as each analysis sub-tab is implemented.
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

Segmentation (5 methods) and Tracking (5 methods) — the heavy computation that turns images into trackable cells.
::::

::::{grid-item-card} 📈 Analysis
:link: analysis/index
:link-type: doc

Feature extraction, filtering, plus skeletons for the behavioural / death-morphology / interaction / backprojection pipelines.
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
