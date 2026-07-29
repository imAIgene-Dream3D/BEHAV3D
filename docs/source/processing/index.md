# Processing

Processing turns the prepared raw data into a set of **labelled, tracked objects** that the analysis stage can compute on.

```{mermaid}
flowchart LR
    Seg["🦠 Segmentation<br/>6 methods"] --> Trk["📍 Tracking<br/>LAP · TrackPy · Propagation · Reporter Propagation · btrack · Import"]
```

## Sections

| Page | Tab in napari | What it covers |
|---|---|---|
| [Segmentation overview](segmentation/index) | 🦠 | When to pick APOC vs ConvPaint vs Pixel Classifier vs Cellpose vs Import existing |
| [Tracking overview](tracking/index) | 📍 | LAP, TrackPy, Propagation, Reporter Propagation, btrack, and Import existing |

```{toctree}
:hidden:
:maxdepth: 2

segmentation/index
tracking/index
```
