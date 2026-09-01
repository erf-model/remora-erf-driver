# Installation and Reproducer Guide

## Quick Start

```bash
git clone --recursive git@github.com:erf-model/remora-erf-driver
cd remora-erf-driver
cmake -S . -B build
cmake --build build -j 8 --target erf-remora-driver
cd Exec/Upwelling
../../build/erf-remora-driver inputs_oneway max_step=10
```

Run coupled decks from their own directory: they import their solver decks with
relative `FILE =` paths.

## Prerequisites

### RRTMGP data files

ERF cases that activate radiation require the RRTMGP coefficient data files at
runtime. Set `RRTMGP_DIR` before running any ERF-based case that turns radiation on:

```bash
export RRTMGP_DIR=/path/to/rrtmgp-data-files
```

> **Build note:** `erf.rrtmgp_file_path` is only meaningful when ERF is built
> with RRTMGP support. In this coupled driver repository **MPI and RRTMGP are
> enabled by default** in the top-level CMake configuration. When you build ERF
> using cmake these are the additional flags needed for radiation
> `-DERF_ENABLE_RRTMGP=ON -DERF_ENABLE_NETCDF=ON -DERF_ENABLE_HDF5=ON` to be
> set in that build.

### `fcompare` executable

The superbuild can build `fcompare` when `ERF_REMORA_ENABLE_FCOMPARE=ON`
(this is the default when `ERF_REMORA_ENABLE_TESTS=ON`). After configure, the
top-level build populates `FCOMPARE_EXE` to the superbuild-owned
`amrex_fcompare` path.

If you want to use a different `fcompare` binary for manual workflows, provide
it via the environment:

```bash
export FCOMPARE_EXE=/path/to/fcompare
```

If `FCOMPARE_EXE` is not set for the helper script, it falls back to
`command -v fcompare` and fails if no `fcompare` is on `PATH`.

### Known limitations

- `driver.atm2ocn_mode = flux` now has a validated startup bootstrap:
  `ocean_step = 0` uses state lanes and later steps use imported flux lanes.
  For the Upwelling debug deck, use nonzero `erf.most.Cd`, `erf.most.Ch`,
  and `erf.most.Cq` so ERF produces live surface stress and heat/moisture
  fluxes.
- REMORA does not currently support vertical refinement in the z-direction in
  coupled mode; the z refinement ratio must be 1.

---

## Build

### 1. MPI-only build (default)

From the repository root:

```bash
cmake -S . -B build
cmake --build build -j 8 --target erf-remora-driver
```

If configured with `-DERF_REMORA_ENABLE_TESTS=ON`, the following CTest labels
are registered (`ERF_REMORA_ENABLE_TESTS` is a driver-repo CMake flag, not an
upstream ERF or REMORA flag):

```bash
ctest --test-dir build -L compile
ctest --test-dir build -L reproducer
```

### 1a. Local CTest workflow

This mirrors the current top-level test registration for local iteration.
It uses the lightweight no-radiation build that the wrapper scripts expect
for fast restart/MPI checks.

Configure:

```bash
cmake -S . -B build_ci \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DERF_REMORA_ENABLE_TESTS=ON \
  -DREMORA_ENABLE_PNETCDF=OFF \
  -DERF_ENABLE_NETCDF=OFF \
  -DERF_ENABLE_RRTMGP=OFF
```

Build:

```bash
cmake --build build_ci --parallel
```

Run restart tests:

```bash
ctest --test-dir build_ci -L restart --output-on-failure --no-tests=error
```

Run only the 2-rank MPI sweep (`np2` vs `np1` plus `np2` smoke):

```bash
ctest --test-dir build_ci -L parallel -R 'np2' --output-on-failure --no-tests=error
```

Run the full MPI sweep (`np2`, `np4`, `np6`, `np8` vs `np1`):

```bash
ctest --test-dir build_ci -L parallel --output-on-failure --no-tests=error
```

### 2. MPI + OpenMP build

Use a separate build directory to preserve the MPI-only build.

`ERF_REMORA_ENABLE_OPENMP` is a **driver-repo CMake flag**. When set to `ON`
the top-level `CMakeLists.txt` automatically force-enables
`Kokkos_ENABLE_OPENMP=ON` and `Kokkos_ENABLE_SERIAL=ON` alongside the
per-solver OpenMP flags, so you do not need to pass those Kokkos flags
manually:

```bash
cmake -S . -B build_omp \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_Fortran_COMPILER=mpif90 \
  -DERF_REMORA_ENABLE_OPENMP=ON
cmake --build build_omp -j 8 --target erf-remora-driver
```

### 3. Perlmutter (NERSC)

This is the build used for the OCCAMS cases, which need NetCDF and PNetCDF for
the WRF input and the REMORA grid/forcing files.

Environment:

```bash
source Submodules/ERF/Build/machines/perlmutter_erf.profile
module load cray-parallel-netcdf
```

The profile loads `gcc-native`, `cmake`, `cudatoolkit`, `cray-hdf5-parallel`,
`cray-netcdf-hdf5parallel` and `cray-libsci`. It has `cray-parallel-netcdf`
commented out, which is why that one is loaded separately — without it, REMORA's
PNetCDF reads of the `pnetcdf_compat` files will not link.

Configure:

```bash
cmake -S . -B build_ncutil \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=cc \
  -DCMAKE_CXX_COMPILER=CC
```

That is deliberately short. Everything else the coupled OCCAMS build needs is
already a default or is forced by this repository's own `CMakeLists.txt`:

- `ERF_ENABLE_MULTIBLOCK`, `ERF_ENABLE_TINY_PROFILE`,
  `REMORA_ENABLE_TINY_PROFILE` are forced ON, and `ERF_ENABLE_HDF5` forced OFF.
- `ERF_ENABLE_NETCDF`, `ERF_ENABLE_PNETCDF`, `REMORA_ENABLE_PNETCDF`,
  `ERF_ENABLE_RRTMGP` and `ERF_REMORA_ENABLE_MPI` all default ON here.
- `ERF_ENABLE_CRAY_AUTO_FIXES` defaults ON and detects Cray itself.

Build:

```bash
cmake --build build_ncutil -j12
```

Build once, and let a failure stop you. `cmake --build` leaves the previous
executable in place when it fails, so a broken rebuild followed by a run
silently exercises the old binary — which is the worst case when you are looking
for diagnostics that only the new code emits.

Notes on this configuration:

- `ERF_ENABLE_RRTMGP` is ON. The OCCAMS decks enable no radiation, so it is
  compiled but unused; turning it off would speed the build but changes nothing
  at runtime and forces a full rebuild.
- `ERF_ENABLE_NOAHMP` is OFF. The coupled-SST path works regardless, but the
  Noah-MP/REMORA-SST coexistence that the coupled-SST design exists to permit is
  not exercised by this build.

Run:

```bash
salloc -N 1 --constraint cpu -A m4304 -t 01:00:00 --qos interactive
cd Exec/OCCAMS
./run_occams.sh
```

After a pull, the cheap check first — one ocean step and two ERF steps, minutes
rather than a queue slot:

```bash
./run_occams.sh smoke
```

That reaches ERF WRF-init, REMORA NetCDF init, spherical remap weight
construction and one full exchange, with FPE trapping left on. It is not long
enough to expose dynamical instability, so it guards startup and coupling, not
the solution. Run `oneway` for that.

See `Exec/OCCAMS/run_occams.sh` for per-case runs and rank overrides, and
`Exec/README.md` for how to read the coupling diagnostics it prints.

---

## Reproducer Decks

The `driver.synthetic_only` control and the `Exec/CouplingReproducers` decks it
drove have been removed, as has `OceanSurf` — ERF now aborts on
`erf.land_surface_model = OceanSurf` and the coupled SST path is selected with
`erf.use_coupled_sst`, which the driver sets for you.

The decks that exist:

| Deck | Direction | Notes |
|---|---|---|
| `Exec/Upwelling/inputs_oneway` | `ERFtoREMORA` | Idealized; needs no external data. Run this first on a new build. |
| `Exec/Upwelling/inputs_twoway` | `two_way` | REMORA SST into ERF's surface layer. |
| `Exec/OCCAMS/inputs_OCCAMS_smoke` | `ERFtoREMORA` | Fast startup and exchange guard; minutes, not a queue slot. |
| `Exec/OCCAMS/inputs_OCCAMS_oneway` | `ERFtoREMORA` | Real OCCAMS domain. |
| `Exec/OCCAMS/inputs_OCCAMS_state` | `two_way` | State-passing contract. |
| `Exec/OCCAMS/inputs_OCCAMS_flux` | `two_way` | Flux-passing contract. |

Each coupled deck pulls its two solver decks in with `FILE =`; run the coupled
deck, not `inputs_ERF` or `inputs_REMORA` directly. The OCCAMS decks reference
absolute NERSC paths and are not runnable elsewhere without editing them.

`Exec/README.md` documents each deck and how to read the driver's coupling
diagnostics.

---

## Verification Runs

Most of what this section used to describe by hand is now in the CTest suite
(see "Local CTest workflow" above). Prefer the suite: it covers more cases, and
the manual procedures it replaced drifted out of date faster than they were run.

**Run the suite at `-j1`.** The harnesses `file(REMOVE_RECURSE)` their working
directories, so concurrent tests delete each other's output and fail in ways
that look like solver bugs.

### 1. Rank comparison — automated

```bash
ctest --test-dir build_ci -L parallel --output-on-failure --no-tests=error
```

64 tests: `np2`, `np4`, `np6` and `np8`, each as a smoke run and as an
`fcompare` against the serial baseline. This replaces the manual
`/tmp/upwelling_mpi` procedure that used to live here.

One test at a time:

```bash
ctest --test-dir build_ci \
  -R '^parallel_upwelling_twoway_state_matched_nopbl_np4_compare_to_np1$' \
  --output-on-failure
```

### 2. Restart validation — automated

```bash
ctest --test-dir build_ci -L restart --output-on-failure --no-tests=error
```

8 tests across direction, contract, and cadence. Each runs a cold start and a
restart and compares the two.

### 3. Coverage contract — automated

```bash
ctest --test-dir build_ci -R reproducer_coupled_sst_coverage --output-on-failure
```

Drives `inputs_twoway` four times and asserts on covered-cell counts rather than
temperature. The only test that exercises the coverage threshold.

### 4. Flux coupling — manual

```bash
mpirun -n 1 ./build/erf-remora-driver Exec/Upwelling/inputs_twoway \
  driver.atm2ocn_mode=flux \
  driver.plot_int_exchange=1
```

Expected behavior:

- step 0 bootstraps REMORA from ERF state lanes,
- step 1 exports live flux lanes including stress, sensible heat, latent heat,
  longwave radiation, and evaporation,
- the driver writes exchange diagnostics when `driver.plot_int_exchange=1`.

`Exec/Upwelling/inputs_ERF` already sets `erf.most.Cd = 0.0014`,
`erf.most.Ch = 0.0012` and `erf.most.Cq = 0.0012`, so no coefficient overrides
are needed. Earlier versions of this guide said those were zero and told you to
override them; that is no longer true, and passing the same values changes
nothing.

### 5. What the suite does not cover

The Upwelling decks are matched 1:1, all sea, with uniform SST and uniform
surface stress. They are structurally blind to the spherical remap, wet-mask
exclusion, non-conformal coverage, and donor fill. A green suite means "no
regression against the idealized case", not that those paths work.

Those first execute on OCCAMS:

```bash
cd Exec/OCCAMS
./run_occams.sh smoke      # startup and one exchange
./run_occams.sh            # oneway, state, flux
```

---

## Unified Inputs Example

AMReX `ParmParse` allows one inputs file to supply parameters to the driver and
both solver libraries simultaneously.

Key dispatch:

- `driver.*` — consumed by `erf-remora-driver`
- `erf.*` and shared mesh keys (`amr.*`, `geometry.*`) — consumed by ERF
- `remora.*` and shared mesh keys (`amr.*`, `geometry.*`) — consumed by REMORA

Example inputs file:

```text
# Import the two solver decks (relative to this file's directory)
FILE = inputs_REMORA
FILE = inputs_ERF

# Driver control
driver.coupling_direction = two_way
driver.atm2ocn_mode       = state

# Coupling cadence is derived, not declared: the driver computes
# atm_substeps_per_ocean_step = remora.fixed_dt / erf.fixed_dt and aborts
# unless that ratio rounds to an integer within tolerance.
erf.fixed_dt    = 3.0
remora.fixed_dt = 300.0

max_step = 864

# Domain topology lives on the solvers, not the driver
erf.n_cell    = 50 20 48
remora.n_cell = 50 20 40

# Optional driver controls
driver.remap_frame       = auto     # auto | spherical | cartesian
driver.remap_v           = 1        # remap geometry and coverage detail
driver.v                 = 0        # driver lifecycle detail
driver.plot_int_exchange = -1       # exchange-field plotfiles

# What a wet ocean cell with no atmosphere source cell receives.
# Spherical frame only; inert in the cartesian frame.
driver.remap_donor_fill  = withhold # withhold | coawst
driver.remap_donor_max_diagonals = 3.0   # coawst acceptance radius, coawst only
```

`driver.remap_donor_fill` decides what happens to a wet REMORA cell that no ERF
cell overlaps — because the two land masks disagree along the coast, or because
the ocean grid reaches past the atmosphere.

- `withhold` (the default) writes no stencil. The cell is reported uncovered and
  `ApplyConservativeRemap` assigns that lane's fallback. Read the lane list before
  assuming this is harmless: **only Patm (101325 Pa) and Tair (288.15 K) get
  standard-atmosphere values.** Uwind/Vwind get 0 m/s; Cloud, SWrad, LWrad and all
  eight flux lanes get a literal 0.0; Qair and Rain are not remapped at all and
  keep whatever `PrefillAtmosToOceanFallbacks` wrote. Withholding is what CMEPS,
  YAC, CDO and CESM do with an unmapped point, and it is the shipped default for
  that reason.
- `coawst` instead grafts a full-weight stencil pointing at the nearest non-land
  ERF cell, reproducing COAWST's `check_weights`. The resulting
  atmosphere-to-ocean operator is **not conservative**, and selecting the mode
  prints an unconditional warning saying so.

`driver.remap_donor_max_diagonals` bounds the `coawst` donor search to that many
of the donor cell's own control-volume diagonals. It is inert under `withhold`,
and the driver prints a note rather than letting a set value look effective. The
bound has no precedent in any production coupler; `Exec/README.md` explains why it
exists anyway and what the value does and does not control.

Retired keys, for anyone reading an old deck: `driver.synthetic_only`,
`driver.dt_atm`, `driver.atm_substeps_per_ocean_step`, `driver.ocean_steps`,
and the `driver.atm.*` / `driver.ocn.*` topology keys no longer exist. Domain
size is read from `erf.n_cell` and `remora.n_cell`. `ParmParse` ignores unknown
keys silently, so a stale deck will run and quietly disregard them.
`driver.coupled_one_step` is still parsed but deprecated, and the driver warns
that it is ignored.

Run with:

```bash
./build/erf-remora-driver /path/to/inputs
```
