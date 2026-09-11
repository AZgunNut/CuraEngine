# Looped Lightning experimental variants

This branch carries a runtime-selectable Lightning experiment set. Variant 0 is stock Cura Lightning. Variants 1-22 are intended to be selectable from Cura without restarting; change the selector and re-slice.

## Variants 1-10: one-change tests

1. Allow targets on the same rendered polyline, excluding source-adjacent segments.
2. Use actual tree leaves as closure sources instead of post-clipping polyline front points.
3. Search targets on actual Lightning tree branches via visitBranches().
4. Search/connect against pre-outline-clipping Lightning geometry.
5. Score targets using distance plus branch-direction alignment.
6. Replace straight closure with a simple quadratic arc.
7. Increase maximum closure reach from 12 mm to 20 mm.
8. Make closure reach proportional to source branch length.
9. Fall back to the second-best target when the nearest target has a poor return angle.
10. Reduce minimum closure distance from 2 extrusion widths to 0.75 extrusion width.

## Variants 11-20: two-change tests

11. Actual tree leaves + actual tree branch targets.
12. Same-polyline targets + angular target scoring.
13. Actual tree branch targets + angular target scoring.
14. Actual tree leaves + same-polyline targets.
15. Pre-clip geometry + actual tree leaves.
16. Same-polyline targets + quadratic arc closure.
17. Actual tree branch targets + quadratic arc closure.
18. Angular scoring + quadratic arc closure.
19. Dynamic reach + angular scoring.
20. Second-target fallback + quadratic arc closure.

## Variant 21: smoothed continuous Lightning path

After normal Lightning is fully generated for each layer, use that finished Lightning geometry as a guide rather than tracing every kink exactly. Build a locally averaged/smoothed guide, offset about four extrusion widths to a consistent side, and derive a replacement continuous route that preserves the overall Lightning form as closely as practical. The generated continuous route replaces the original Lightning output for this variant.

When two neighboring smoothed routes do not have enough room for their requested offsets, do not cross, merge, or bunch them. Through the cramped region place the two printable paths side-by-side around the median space between their Lightning guides at approximately one extrusion-width center-to-center spacing, then smoothly transition back to the requested offset when room opens.

Priority: continuity -> two printable side-by-side lines -> smoothness -> requested offset.

## Variant 22: smoothed sine continuous Lightning path

Same post-processing concept and collision fallback as variant 21, but modulate the one-sided offset from 0 to 4 extrusion widths and back to 0 with a 10 mm full sine wavelength along the smoothed guide. Phase is locked so every layer begins from the same side and same phase; do not randomly flip or alternate it between layers. Because the underlying Lightning evolves as it propagates down the model, the phase-locked sine should naturally walk/spiral through 3D as the form changes.

## Implementation constraints

- Keep Cura's internal parent/child Lightning topology acyclic; experimental cross-links are render-time geometry only.
- Variant 0 must remain a stock control.
- Runtime selection must not require restarting Cura; changing the selector and re-slicing is sufficient.
- Newly added closure geometry is not reused as a target unless a variant explicitly says so.
- Clip final experimental geometry to the normal infill outline.
