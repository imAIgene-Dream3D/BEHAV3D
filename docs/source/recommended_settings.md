# Deriving settings for a new dataset

BEHAV3D settings should be derived from the current acquisition and previewed. Cell
names and values from earlier experiments are not presets.

## Measurements to collect

Before recommending values, collect:

1. XY pixel size and Z spacing.
2. Time between consecutive frames.
3. Approximate object diameter, or number of cell widths across a multicellular structure.
4. Fastest plausible displacement between two consecutive frames.
5. Whether consecutive masks overlap, touch, merge, divide, enter, or leave the field.
6. Whether each signal is constitutive, switch-on, fluctuating, or intermittent.

Whenever a control uses frames, also report the corresponding physical duration.

## Segmentation scale

For a single cell, 10 micrometres is a transparent initial diameter assumption only
when the researcher has not supplied a measurement. Convert it to pixels per sample:

```text
diameter_pixels = diameter_micrometres / XY_pixel_size
```

For EDT, preview 20%, 25%, and 30% of that pixel diameter. For a multicellular
structure, first ask how many approximately 10 micrometre cell widths span it. A
10-30-cell range can be used as a sanity check, not as an automatic answer.

For a 3D minimum-size preview, estimate the full spherical volume and begin at 50%
to tolerate dim or incomplete masks:

```text
estimated_voxels = 0.5 * (4/3 * pi * (diameter_micrometres / 2)^3)
                   / (XY_pixel_size^2 * Z_spacing)
```

Elongated, protrusive, and confluent objects do not fit a diameter model. Use a
permissive threshold, annotate borders/background, and inspect the segmentation.

## Tracking method

Use frame-to-frame behavior:

| Observed behavior | Starting method |
|---|---|
| Detection loses overlap with its previous position | btrack |
| Detection overlaps and may fragment | Fragmentation Propagation |
| Detection overlaps/touches and IDs must remain within connected regions | Bounded Propagation |
| Physical object is static but its segmented reporter flickers | Reporter Propagation |

For btrack:

```text
maximum_search_radius = fastest_measured_one_frame_displacement * (1 + margin)
```

Use a modest 10-25% margin. Calibrate the global optimizer from the largest spatial
gap and missing-frame duration that should reconnect. Enable branching for division;
retain initialization and termination for field-of-view entry and exit. Lower Step
size only after an out-of-memory error.

For a moving object with a flickering reporter, track a constitutive channel and use
the reporter as an intensity feature, or use permissive segmentation plus btrack.

## Feature extraction and analysis

Choose features from the research question:

- Movement for motility and directional behavior.
- Contact for strict touching or calibrated proximity.
- Morphology for shape/plasticity questions.
- Channel intensity for switch-on or fluctuating reporters.
- Invasiveness for surface engagement between selected objects.

Death Dynamics measures the population fraction with a switch-on signal over time.
Active Killing measures a contact-associated signal rise in a selected target. These
names do not prove the signal represents death; preserve the researcher's biological
definition.

Behavioral State can classify any single-cell population. Ask whether states should
represent movement, contact, morphology, channel intensity, or a combination. Derive
window lengths from event duration and acquisition cadence rather than using a bare
frame default.

## Filtering

Inspect the track-length distribution before choosing a minimum. A minimum removes
shorter tracks; a maximum trims retained tracks to a common window, so equal minimum
and maximum values are valid. Filtering must still run with all filters disabled
because it creates the downstream feature table and interpolates missing timepoints.

Use the track-count preview for count questions. Counts must be calculated from the
current combined track-features CSV at the requested timepoint and threshold; never
estimate them from prose.
