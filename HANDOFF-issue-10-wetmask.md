# Handoff: Issue 10 Wet Mask

## Changed

- `Source/ERFRemoraMultiBlockContainer.cpp:954` adds a flat-index helper for domain-sized lookups.
- `Source/ERFRemoraMultiBlockContainer.cpp:964` adds `BuildReplicatedMaskLookup`, which gathers REMORA mask values into a rank-local lookup used during O2A weight construction.
- `Source/ERFRemoraMultiBlockContainer.cpp:2047` fetches only REMORA's rho-centered `mskr` mask for the SST lane.
- `Source/ERFRemoraMultiBlockContainer.cpp:2052` builds a replicated rho-mask lookup when the mask is present and cell-centered.
- `Source/ERFRemoraMultiBlockContainer.cpp:2097` keeps `covered_frac` as the pre-normalization destination fraction, now restricted to wet REMORA rho-cell area.
- `Source/ERFRemoraMultiBlockContainer.cpp:2105` skips dry REMORA rho cells before clipping, coverage accumulation, and stencil insertion.
- `Source/ERFRemoraMultiBlockContainer.cpp:2146` preserves the nonempty-overlap plus `m_remap_min_coverage_frac` coverage guard.
- `Source/ERFRemoraMultiBlockContainer.H:276` updates `m_valid_o2a_mf` documentation to describe wet-source coverage.
- `Source/StandaloneSolverRuns.H:51` updates the `ERFApplyOceanSurfaceState` coverage contract.

## Normalization

`covered_frac` is accumulated at `Source/ERFRemoraMultiBlockContainer.cpp:2097` through `Source/ERFRemoraMultiBlockContainer.cpp:2120`, before stencil truncation writes and before `NormalizeWeightsAndCheckConservation` at `Source/ERFRemoraMultiBlockContainer.cpp:2171`. Normalization can still rescale kept weights to sum to one, but the coverage flag has already been computed from the unnormalized wet area fraction.

## Wetness Handling

Dry cells are excluded outright using `mskr > 0.5`. REMORA's mask contract is binary, with `1 = wet` and `0 = land`, so this is equivalent to weighting by source wetness for the current mask while avoiding dry cells in the SST stencil entirely. I did not introduce fractional wetness handling because the available mask does not represent partial wet fractions.

## Deliberately Not Done

- Did not implement the spherical O2A branch; the existing early warning and return remain unchanged.
- Did not add a bilinear O2A option.
- Did not touch `BuildAtmosToOceanWeights` or `BuildRemapWeights`.
- Did not touch `Submodules/`, `Exec/`, or `CMakeLists.txt`.

Everything runtime is UNVERIFIED. Per instruction, I did not build, run CMake, run MPI, or run tests.
