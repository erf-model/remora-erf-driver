# HANDOFF issue 10: spherical O2A SST remap

## What changed

- `Source/ERFRemoraMultiBlockContainer.cpp:2243-2244` documents the deliberate absence of O2A donor fill: uncovered ERF cells keep ERF's own `wrflowinp` SST.
- `Source/ERFRemoraMultiBlockContainer.cpp:2245-2352` adds the spherical branch of `BuildOceanToAtmosWeights()`.
- `Source/ERFRemoraMultiBlockContainer.cpp:2265-2269` builds a `SourceSearchGrid` over the staged REMORA rho-cell region using staged psi corner lon/lat.
- `Source/ERFRemoraMultiBlockContainer.cpp:2273-2276` forms each ERF destination cell polygon from `m_erf_node_coords`.
- `Source/ERFRemoraMultiBlockContainer.cpp:2295-2298` forms each candidate REMORA source cell polygon from the staged psi corners.
- `Source/ERFRemoraMultiBlockContainer.cpp:2299-2304` clips ERF destination polygons against wet REMORA candidate polygons, accumulates raw coverage, and records source overlaps.
- `Source/ERFRemoraMultiBlockContainer.cpp:2344-2348` stores `intersect_area / dst_area` weights and REMORA source indices in the same layout as the Cartesian O2A branch.
- `Source/ERFRemoraMultiBlockContainer.cpp:2353-2480` keeps the existing Cartesian O2A path under the non-spherical `else`; it still uses the replicated rho-mask lookup.
- `Source/ERFRemoraMultiBlockContainer.H:239-240` updates the `m_remora_staging_plan` member comment now that the spherical weight loop consumes it.

## Invariants

1. `covered_frac` is accumulated from raw destination-area fractions before normalization.

   Relevant line: `Source/ERFRemoraMultiBlockContainer.cpp:2303`

   ```cpp
   covered_frac += intersect_area / dst_area;
   ```

   This happens before weights are written at `Source/ERFRemoraMultiBlockContainer.cpp:2344-2348` and before `NormalizeWeightsAndCheckConservation()` runs at `Source/ERFRemoraMultiBlockContainer.cpp:2482`.

2. The coverage flag tests both a nonempty overlap set and the coverage threshold.

   Relevant lines: `Source/ERFRemoraMultiBlockContainer.cpp:2329-2330`

   ```cpp
   valid_arr(i, j, 0) =
       (!overlaps.empty() && covered_frac >= m_remap_min_coverage_frac) ? 1 : 0;
   ```

   There is also an earlier `if (overlaps.empty()) { continue; }` at `Source/ERFRemoraMultiBlockContainer.cpp:2308`, so cells with no wet REMORA overlaps remain at the zero-filled coverage value initialized at `Source/ERFRemoraMultiBlockContainer.cpp:2236`.

3. Dry REMORA rho cells are excluded before clipping, using the staged mask.

   Relevant line: `Source/ERFRemoraMultiBlockContainer.cpp:2293`

   ```cpp
   if (remora_wet(src_i, src_j, 0) <= Real(0.5)) { continue; }
   ```

   This `continue` occurs before the REMORA polygon is formed and before `ClipPolygon()` at `Source/ERFRemoraMultiBlockContainer.cpp:2295-2300`, so `covered_frac` is a wet-area fraction, not a geometric area fraction.

## Staged mask and corners

- ERF destination corners come from `m_erf_node_coords` at `Source/ERFRemoraMultiBlockContainer.cpp:2260`.
- Staged REMORA psi corners come from `m_remora_staging_plan->node_mf` at `Source/ERFRemoraMultiBlockContainer.cpp:2261`.
- The staged REMORA rho wet mask comes from `m_remora_staging_plan->wet_mf` at `Source/ERFRemoraMultiBlockContainer.cpp:2262`.
- The search source box is the staged wet-mask box for the current ERF box: `Source/ERFRemoraMultiBlockContainer.cpp:2263`.
- The spherical branch does not call `BuildReplicatedMaskLookup`; that remains only in the Cartesian `else` at `Source/ERFRemoraMultiBlockContainer.cpp:2366-2368`.

## Uncertainty

- Everything runtime is UNVERIFIED.
- Build behavior is UNVERIFIED; I did not build.
- Test behavior is UNVERIFIED; I did not run tests.
- There is NO test in either repository that can exercise this path. It needs OCCAMS-class spherical input.
- The correctness of using one `MFIter` index across `m_weight_o2a_mf`, `m_erf_node_coords`, and the staged REMORA MultiFabs rests on the same DistributionMapping/box-index reasoning already used by the A2O spherical staging consumer, not on runtime evidence.
- Polygon clipping, search-bin candidate selection, wet-mask exclusion, and coverage-threshold behavior in the new spherical branch rest on code inspection and mirroring of the existing A2O spherical and O2A Cartesian patterns, not on an executed spherical case.
- Antimeridian/polar longitude wrapping behavior remains UNVERIFIED; this follows the existing planar lon/lat clipping assumption in the A2O spherical path.
