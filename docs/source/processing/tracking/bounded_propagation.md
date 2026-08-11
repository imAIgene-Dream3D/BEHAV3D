# Bounded Propagation

Bounded Propagation is an overlap-based tracker for objects that remain spatially
overlapping between consecutive frames but may touch or form joined segmentation
regions. It uses watershed propagation like Fragmentation Tracking, with one extra
topology rule: a track ID cannot span more than one disconnected region.

## When to use it

Use Bounded Propagation when all of the following are true:

- objects retain substantial frame-to-frame overlap;
- masks can touch, merge, or change shape;
- an ID spreading across disconnected regions would be a meaningful tracking error.

If detections lose overlap because displacement is large relative to the frame
interval, use btrack instead. If masks remain separate and fragmentation is the main
failure mode, Fragmentation Tracking may be simpler. Preview both when topology makes
either interpretation plausible.

## How it works

Before watershed propagation, the current frame is split into connected regions.
Each existing track claims the region that overlaps its previous footprint most. A
track's seed is discarded outside that region. More than one track may claim the same
region, allowing watershed to split it; an unclaimed region becomes a new track.

## Parameters

| Control | What it does |
|---|---|
| **Minimum overlap fraction** | Minimum fractional overlap required before an existing track can claim a current connected region. Calibrate from a representative preview. |
| **Minimum segment size** | Removes small candidate regions before propagation. Derive it from measured object size and XY/Z sampling, then verify that real partial objects are retained. |

## Validation

Inspect at least one representative movie in Visualization. Verify that IDs remain
stable through touching events, that one ID does not cover disconnected objects, and
that small real objects are not removed by the size threshold.

## See also

- [Tracking overview](index)
- [Fragmentation Tracking](fragmentation_tracking)
- [btrack](btrack)
