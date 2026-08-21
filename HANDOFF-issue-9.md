# Issue 9 Handoff

## Changed

- `Source/ERFRemoraMultiBlockContainer.cpp:1673`: `BuildAtmosToOceanWeights` now derives `src_face_dir` from each A2O weight family, so `TAUX_FACE` maps to ERF x-faces and `TAUY_FACE` maps to ERF y-faces.
- `Source/ERFRemoraMultiBlockContainer.cpp:1701`: the spherical path now keeps the staged ERF node box and uses `FaceCellQuadFromNodes` for TauX/TauY source control volumes, while scalar/rho families keep `ExtractCellQuadFromCoords`.
- `Source/ERFRemoraMultiBlockContainer.cpp:1725`: spherical TauX/TauY source-face land checks use the two bracketing ERF cell-mask values.
- `Source/ERFRemoraMultiBlockContainer.cpp:1793`: spherical TauX/TauY source candidates are converted from cell candidates to face candidates by scanning the candidate-cell hull and extending the high index in the staggered direction.
- `Source/ERFRemoraMultiBlockContainer.cpp:1848`: the Cartesian path now allows the extra source-face index in the staggered direction and builds TauX/TauY source bounds with `FaceCellBoundsFromGeom`.

## Why

ERF packs `tau13` and `tau23` on x- and y-face layouts, but `BuildAtmosToOceanWeights` was storing ERF cell indices for the TauX/TauY A2O families. The stencil then dereferenced face-staggered source arrays with cell-centered indices, shifting wind stress by half an ERF cell and making the outermost source-face column unreachable. The change stores face indices for the face-staggered stress families.

## Maintainer Decisions

- Face land-mask convention in the spherical path: I chose land if either bracketing ERF cell has `erf_land >= 0.5`, matching REMORA's `msku`/`mskv` convention that a face is wet only when both adjacent cells are wet. At a domain edge, the bracketing lookup clamps to the single available ERF cell.
- Face candidates from `SourceSearchGrid`: I kept `SourceSearchGrid` cell-based, then for TauX/TauY scanned the min/max hull of returned cell candidates and extended the high index by one in the active face direction. Exact polygon clipping against the face control volumes still decides which faces enter the stencil.

## Deliberately Not Done

- Did not change the remap method.
- Did not add `driver.a2o_stress_method`.
- Did not implement bilinear wind-stress remapping.
- Did not touch anything under `Submodules/`.
- Did not touch `BuildOceanToAtmosWeights`.
- Did not change `Source/ERFRemoraMultiBlockContainer.H`.

## Unverified

Everything runtime is unverified. Per instruction, I did not build, did not run CMake, did not run MPI, and did not run tests. I only inspected the source diff.
