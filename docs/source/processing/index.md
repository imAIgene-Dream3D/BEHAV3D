# Processing

Processing turns the prepared raw data into a set of **labelled, tracked objects** that the analysis stage can compute on.

```{mermaid}
flowchart LR
    Seg["🦠 Segmentation<br/>6 methods"] --> Trk["📍 Tracking<br/>btrack · Fragmentation Propagation · Bounded Propagation · Reporter Propagation · TrackPy · LapTrack · Import"]
```

## Sections

| Page | Tab in napari | What it covers |
|---|---|---|
| [Segmentation overview](segmentation/index) | 🦠 | When to pick APOC vs ConvPaint vs Pixel Classifier vs Cellpose vs Cellpose-SAM vs Import existing |
| [Tracking overview](tracking/index) | 📍 | btrack, Fragmentation Propagation, Bounded Propagation, Reporter Propagation, TrackPy, LapTrack, and Import existing |

```{toctree}
:hidden:
:maxdepth: 2

segmentation/index
tracking/index
```
