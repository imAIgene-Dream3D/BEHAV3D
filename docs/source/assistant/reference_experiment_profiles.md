# Generalized phenotype guidance

The assistant must not choose a method or value from a biological name. A name can
be an optional example after the recommendation has been derived from measurable
properties of the current images.

## Describe the data before choosing a method

Record these properties for each processing population:

| Property | Evidence to request |
|---|---|
| Object scale | Approximate diameter in micrometres or cell widths |
| Sampling | XY pixel size, Z spacing, and time between frames |
| One-frame motion | Measured displacement or whether the object still overlaps itself |
| Shape | Stable, protrusive, elongated, dividing, or fragmenting |
| Topology | Isolated, touching, temporarily merged, or confluent |
| Signal | Constitutive, switch-on, fluctuating, or intermittently invisible |
| Density | Sparse enough to separate, or crowded/ambiguous |

"Fast" and "slow" are not fixed cell properties. An object is fast for tracking
when its displacement during one acquisition interval is large enough that its next
detection no longer overlaps the previous one. The same object can therefore be
slow at a short interval and fast at a long interval.

## General object profiles

These profiles describe image behavior, not biological identities.

### Small, shape-stable single objects

Often approximately one cell diameter (around 10 micrometres as an initial sizing
reference). Use observed frame-to-frame overlap to choose tracking:

- retained overlap: Fragmentation Tracking or Bounded Propagation may work;
- lost overlap: use btrack and derive the search radius from measured displacement;
- division: enable the btrack branching hypothesis;
- entry or exit at the field boundary: retain initialization and termination hypotheses.

### Large, plastic or protrusive single objects

Diameter alone is a poor mask model when the object extends processes or changes
shape. Use boundary annotations and a permissive size threshold. If touching
segmentations join while motion remains small, Bounded Propagation prevents one
track ID from spreading across disconnected regions. If detections lose overlap,
use btrack.

### Multicellular structures

Ask how many roughly 10 micrometre cell widths span the structure. Ten to thirty
cells across (roughly 100-300 micrometres) is a sanity-check range, not a default.
Use the actual answer and each sample's XY/Z sampling to derive segmentation scales,
minimum size, and EDT starting points.

For near-static or slowly drifting structures that overlap between frames, use
Fragmentation Tracking. Track multiple distinguishable structure populations in one
run when they coexist and need a shared overlap-based track space; preserve their
origin labels. Explants and much larger structures require their own measured scale.

### Static objects with an intermittent reporter

Reporter Propagation is valid only when both conditions hold:

1. the physical object is effectively static; and
2. the segmented reporter flickers or disappears.

If the object moves, track a constitutive channel and extract the reporter as an
intensity feature. If no constitutive channel exists, use permissive segmentation
that tolerates false positives and a moving-object tracker. Reporter Propagation
must not freeze a moving object's shape.

### Elongated objects or confluent sheets

Diameter-based minimum-size estimates are unreliable. Use a permissive seed or size
threshold, annotate borders and background explicitly, and inspect a preview. If
connected components repeatedly join while objects still overlap themselves, compare
Bounded Propagation with btrack.

## Tracking decision tree

1. Does the object remain spatially overlapping between consecutive frames?
   - No: use btrack.
   - Yes: continue.
2. Is the signal intermittent while the physical object is static?
   - Yes: use Reporter Propagation.
   - No: continue.
3. Can touching/merged masks cause one track ID to span disconnected regions?
   - Yes: use Bounded Propagation.
   - No: use Fragmentation Tracking.
4. When two methods are plausible, explain the topology difference and preview both.

For btrack, derive the maximum search radius from the fastest plausible measured
one-frame displacement plus a modest margin. Derive optimizer distance and time
thresholds from the largest spatial and missing-frame gaps that should be bridged.
Never import a numeric value from a named population or another experiment.

## Signal-to-analysis mapping

- A signal that rises and remains positive can drive Death Dynamics and, when tied
  to contact, Active Killing or Interaction Analysis. Call it death only when the
  reporter has been biologically validated as a death signal.
- A fluctuating signal can be a Behavioral State or State Trajectory feature.
- Contact, movement, morphology, and channel intensity can all define single-cell
  states. Ask what the researcher wants the state profile to represent.
- Convert every window or threshold expressed in frames into physical duration using
  metadata before explaining or recommending it.

## Metadata taxonomy

- A processing population is an object or signal distinguishable in the image and
  requiring its own segmentation and track IDs.
- Line records biological identity or source; Condition records treatment or state.
- Multicolor means one biological population was intentionally split across colors
  to reduce density, then segmented/tracked separately and recombined. Every color
  must share population, line, and condition.
- Different populations or biological identities are not Multicolor, even if color
  distinguishes them.
- With already acquired data, recognize the acquisition design that exists; do not
  recommend adding colors that were never acquired.
