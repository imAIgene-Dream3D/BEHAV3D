# Processing

Processing turns the prepared raw data into a set of **labelled, tracked objects** that the analysis stage can compute on.

```{mermaid}
flowchart LR
    Seg["🦠 Segmentation<br/>5 methods"] --> Trk["📍 Tracking<br/>5 methods"]
```

## Sections

| Page | Tab in napari | What it covers |
|---|---|---|
| [Segmentation overview](segmentation/index) | 🦠 | When to pick APOC vs ConvPaint vs Pixel Classifier vs Cellpose vs Import existing |
| [Tracking overview](tracking/index) | 📍 | When to pick LAP vs TrackPy vs Propagation vs btrack vs Import existing |

```{toctree}
:hidden:
:maxdepth: 2

segmentation/index
tracking/index
```
