# Handoff: Issue 10 O2A Staging

## Changed

- `Source/ERFRemoraMultiBlockContainer.H:47` forward-declares `REMORAStagingPlan`.
- `Source/ERFRemoraMultiBlockContainer.H:170` declares `BuildREMORAStagingPlan()`.
- `Source/ERFRemoraMultiBlockContainer.H:239` adds `m_remora_staging_plan`.
- `Source/ERFRemoraMultiBlockContainer.cpp:351` defines the separate O2A staging structure.
- `Source/ERFRemoraMultiBlockContainer.cpp:1304` calls `BuildREMORAStagingPlan()` during spherical two-way setup.
- `Source/ERFRemoraMultiBlockContainer.cpp:1688` adds `BuildREMORAStagingPlan()`.
- `Source/ERFRemoraMultiBlockContainer.cpp:1703` fetches REMORA `lon_psi`/`lat_psi`.
- `Source/ERFRemoraMultiBlockContainer.cpp:1708` fetches REMORA rho mask `mskr`.
- `Source/ERFRemoraMultiBlockContainer.cpp:1736` builds REMORA source box lon/lat bounds directly from psi nodes.
- `Source/ERFRemoraMultiBlockContainer.cpp:1759` builds ERF destination box bounds from the existing ERF corner mesh.
- `Source/ERFRemoraMultiBlockContainer.cpp:1784` derives per-ERF-box REMORA rho staging regions.
- `Source/ERFRemoraMultiBlockContainer.cpp:1820` allocates staged REMORA node and wet-mask MultiFabs on the ERF DistributionMapping.
- `Source/ERFRemoraMultiBlockContainer.cpp:1836` creates fresh `MultiBlockCommMetaData` for REMORA node and rho copies.
- `Source/ERFRemoraMultiBlockContainer.cpp:1843` stages `lon_psi`, `lat_psi`, and `mskr`.
- `Source/ERFRemoraMultiBlockContainer.cpp:2244` documents cartesian replicated mask lookup versus spherical staged mask use.

## Difference From BuildERFStagingPlan

`BuildREMORAStagingPlan()` reverses the source and destination: REMORA rho/psi data are staged onto the ERF/O2A SST layout, while `BuildERFStagingPlan()` stages ERF data onto the REMORA layout.

The new plan does not reconstruct corners. REMORA `lon_psi` and `lat_psi` are already nodal cell-corner fields, so they are copied directly into `node_mf` components 0 and 1. It also stages only `mskr`, because O2A currently has exactly one SST lane at rho points. There is no A2O-style family loop and no REMORA face-family growth.

The per-box search remains intentionally coarse like `BuildERFStagingPlan()`: it reduces geographic bounds for all source and destination boxes, then unions whole intersecting source boxes into the staged region.

## A2O Verification

I verified A2O was not edited by inspecting `git diff -U0 -- Source/ERFRemoraMultiBlockContainer.cpp Source/ERFRemoraMultiBlockContainer.H`. The diff contains no modified lines inside `BuildERFStagingPlan()` or inside the `BuildAtmosToOceanWeights()` body; the only initialization change is the additional guarded `BuildREMORAStagingPlan()` call after the existing `BuildERFStagingPlan()`.

I also checked `rg` results: the A2O spherical consumer still reads only `m_staging_plan` in `BuildAtmosToOceanWeights()` at `Source/ERFRemoraMultiBlockContainer.cpp:1941`, and the new O2A staging state is a separate `m_remora_staging_plan` member.

## MultiBlockCommMetaData Uncertainty

The argument order and identity `MultiBlockIndexMapping` mirror `BuildERFStagingPlan()`: destination staging BoxArray/ERF DistributionMapping first, REMORA source domain and source BoxArray/DistributionMapping second. I am confident this is the intended mirror.

The runtime behavior is not exercised here. I added assertions that REMORA psi boxes are nodal conversions of the rho mask boxes, that `lon_psi` and `lat_psi` share box indexing and processor mapping, and that ERF corners are nodal conversions of the O2A weight boxes. If any of those assumptions are false in real solver layouts, this staging plan needs adjustment before the O2A weight loop consumes it.

Using the same node metadata for `lon_psi` and `lat_psi` should be correct under those assertions, but it remains MPI-runtime unverified.

## Runtime Status

Everything runtime is UNVERIFIED. I did not build, run CMake, run MPI, or run tests.
