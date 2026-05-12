# btrack — Bayesian Cell Tracking in BEHAV3D

btrack uses a **Kalman filter** to predict object motion frame-by-frame, and optionally a **global hypothesis model** to resolve track initiations, terminations, and gap-bridging in a second pass.

---

## Recommended Workflow

### Step 1 — Kalman filter linking (always runs)

Run **without the optimizer first** and validate the result before enabling Step 2.

| Parameter | What it controls | When to change |
|---|---|---|
<<<<<<< HEAD
| **Config preset** | Bundled motion + hypothesis model | *Cell* for T cells/organoids; *Particle* for small fast objects; *Custom* for your own JSON |
=======
| **Config preset** | Bundled tracking preset | *Cell* for motion-only cell tracking; *Particle* for small fast objects; *Custom* for your own JSON |
| **Use visual features** | Whether to use raw-image-derived appearance features during linking | Enable this when you want btrack to use `area`, `major_axis_length`, `minor_axis_length`, and per-channel mean intensities from `raw_image_path` |
>>>>>>> dev
| **Max search radius (px)** | Maximum per-frame displacement allowed | Increase if cells move fast; decrease to avoid false links |
| **Update method** | EXACT (default) vs APPROXIMATE | Switch to APPROXIMATE only for >1000 objects/frame |
| **Step size (frames)** | Batch size for iterative linking | Leave at 100 unless memory is an issue |

Check the result in napari before proceeding:
- Tracks jumping between unrelated cells → reduce *max search radius*
- Too many short track fragments → increase *max search radius* or `max_lost` in the config JSON

### Step 2 — Global hypothesis optimizer (opt-in)

Resolves false positives, track starts/ends, and gap-bridging globally.
Only enable once Step 1 looks correct.

| Parameter | What it controls |
|---|---|
| **Hypotheses** | Which biological events to model (see table below) |
| **Distance threshold** | Max spatial distance to consider when reconnecting fragments |
| **Time threshold** | Max frame gap to attempt reconnection |

**Available hypotheses:**

| Hypothesis | Meaning | Default |
|---|---|---|
| `P_FP` | False-positive detection | Always on |
| `P_init` | Track initiation (cells entering field of view) | On |
| `P_term` | Track termination (cells leaving field of view) | On |
| `P_link` | Standard frame-to-frame linking | On |
| `P_branch` | Track branching / cell division | Off |
| `P_dead` | Cell death | Off |
| `P_merge` | Track merging / cell fusion | Off |

---

## Configuration Files

Bundled JSON configs are in this folder (`behav3d/preprocessing/tracking/models/`):

| File | Preset name | Best for |
|---|---|---|
| `cell_config.json` | `cell` | T cells, tumour cells, organoids |
| `particle_config.json` | `particle` | Small, fast particles |

Copy and edit either file to create a custom config. Select *Custom JSON* in the GUI and browse to your file.

<<<<<<< HEAD
=======
When **Use visual features** is enabled, BEHAV3D adds these skimage-derived features to each object before tracking and then runs btrack with `tracking_updates=["motion", "visual"]`. Any valid btrack JSON can be used with or without visual features; the JSON itself does not decide that. The output tracks CSV stays the standard BEHAV3D schema.

>>>>>>> dev
### Key JSON fields (inside `"TrackerConfig"`)

```json
"MotionModel": {
  "accuracy": 10,           // Kalman filter measurement noise
  "prob_not_assign": 0.01,  // Probability that an object is not assigned
  "max_lost": 5             // Max frames an object can be missing
},
"HypothesisModel": {
  "hypotheses": ["P_FP", "P_init", "P_term", "P_link"],
  "dist_thresh": 60,        // Spatial threshold for hypothesis generation (px)
  "time_thresh": 3,         // Temporal threshold for hypothesis generation (frames)
  "relax": true             // Relax volume constraints at image boundaries
}
```

---

## Output

Each run produces two files per sample:

- `<sample>_tracks.csv` — BEHAV3D standard track table  
  Columns: `TrackID`, `SegmentID`, `position_t`, `position_x/y/z`, `pixel_position_x/y/z`
- `<sample>_tracked.zarr` — Segmentation image with labels replaced by `TrackID`

---

## References

- btrack documentation: https://btrack.readthedocs.io  
- Ulicna et al. (2021) *Automated deep lineage tree analysis using a Bayesian single cell tracking approach.* Frontiers in Computer Science.
