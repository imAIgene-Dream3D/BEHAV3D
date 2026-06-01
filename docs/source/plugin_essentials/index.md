# Plugin Essentials

Some parts of the BEHAV3D EXPLORER plugin do not belong to a single pipeline step, you use them at every stage. This section documents those cross-tab tools.

| Page | What it covers |
|---|---|
| [Visualization](visualization) | The Visualization tab, how to load any sample (raw, segments, tracks) into napari layers. Used after every step to inspect outputs. |
| [Processing Queue](processing_queue) | The 🛒 panel at the bottom of the dock widget. Batches steps from segmentation / tracking / feature extraction / filtering and runs them sequentially across all samples. |
| [Output Directory & File Layout](output_layout) | The canonical folder tree that every tab writes to. Where to find segments, tracks, features, analysis results on disk. |

```{toctree}
:hidden:
:maxdepth: 1

visualization
processing_queue
output_layout
```
