# 🧱 Bounded Propagation

Bounded Propagation is an overlap-based tracker for static or slow-moving objects that remain spatially
overlapping between consecutive frames but may touch or form joined segmentation
regions. It uses watershed propagation like [Fragmentation Tracking](fragmentation_tracking), with one extra
topology rule: a track ID cannot span more than one disconnected region.

## When to use it

Use Bounded Propagation when all of the following are true:

- objects retain substantial frame-to-frame overlap (static/slow-moving);
- masks can touch, merge, or change shape but should keep separate identities;
- connected objects that separate during the experiment should be seen as separate tracked objects. 

If detections lose overlap because displacement is large relative to the frame
interval, use [btrack](btrack) instead. If spatially separate object should belong to the same track (dying objects that are falling apart) [Fragmentation Tracking](fragmentation_tracking) should be used. Preview both when topology makes either interpretation plausible.

## How it works

Before watershed propagation, the current frame's mask is split into **connected regions**. Each existing track keeps its seed only in the **one region its previous footprint overlaps most**; its markers everywhere else are discarded. Watershed then floods each region from the seeds it kept:

- a region only one track kept → flooded entirely to that track;
- a region several tracks kept → split between them by watershed (genuine merges/splits are preserved);
- a region no track kept → left unlabelled, then given a **single new track ID per disconnected leftover patch** (so one blob does not splinter into several new tracks).

Gap-closing dilation is deliberately **off** here, so a track ID can never leak across a background gap into a disconnected blob — that is what keeps each ID inside one connected region.

## Parameters

| Control | What it does |
|---|---|
| **Minimum overlap fraction** | Before a track may claim a connected region, its previous-frame footprint must cover at least this fraction **of that region's area**. `0` (default) = any shared pixel qualifies; raise it to stop a track claiming a region it barely grazes. |
| **Minimum segment size** | Minimum segment size in voxels, applied **after** watershed. It cleans up the small fragments watershed splitting itself can create (a shared region divided between tracks, or a leftover sliver spun into a new track), so it is normally set **much lower** than a segmentation-stage size filter (default `20`). Flat / 2-D segments (one voxel thick on any axis) are always removed regardless of this value. |

## Validation

Inspect at least one representative movie in Visualization. Verify that IDs remain
stable through touching events, that one ID does not cover disconnected objects, and
that small real objects are not removed by the size threshold.

## See also

- [Tracking overview](index)
- [Fragmentation Tracking](fragmentation_tracking)
- [btrack](btrack)
