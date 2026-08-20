# Installation and Reproducer Guide

## Quick Start

```bash
git clone --recursive git@github.com:erf-model/remora-erf-driver
cd remora-erf-driver
cmake -S . -B build
cmake --build build -j 8 --target erf-remora-driver
./build/erf-remora-driver Exec/CouplingReproducers/inputs_ERFtoREMORA_coupleToERF max_step=10
```

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
  -DREMORA_ENABLE_NETCDF=OFF \
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

---

## Ordered Reproducer Set

### Reproducer A — Coupling infrastructure control (`driver.synthetic_only = 1`)

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_coupling_infrastructure_check
```

Scientific intent: verify coupling cadence and subcycling bookkeeping; verify
cache slab ownership and update path without solver-physics confounders.

### Reproducer B — One-step standalone solver initialization (`driver.synthetic_only = 0`)

```bash
./build/erf-remora-driver \
  Exec/CouplingReproducers/inputs_one_step_solver_initialization_check \
  max_step=1
```

Scientific intent: verify both solvers initialize and advance prognostic state
in the shared executable; separate initialization issues from coupling-path
issues.

### Reproducer C — `ERFtoREMORA` state-exchange baseline

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_ERFtoREMORA_coupleToERF
```

Scientific intent: verify atmosphere-state transfer and ocean application at
each coupled step; establish baseline coupled behavior before enabling SST
feedback.

### Reproducer D — `two_way` SST feedback into ERF (`OceanSurf` active)

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_two_way_coupleToERF \
  max_step=600
```

Scientific intent: verify SST from REMORA is injected into ERF `OceanSurf`
`t_surf`; verify ERF surface-layer fluxes respond to ocean state (`Cd`, `Ch`,
`Cq` nonzero); produce measurable thermodynamic response in `temp` and `theta`
over long windows.

---

## Verification Runs

These runs assume you are in the repository root and that `RRTMGP_DIR` and
`FCOMPARE_EXE` are set (see Prerequisites above).

### 1. Upwelling smoke test

Runs one coupled ocean step using the Perlin upwelling case:

```bash
mpirun -n 1 ./build/erf-remora-driver Exec/Upwelling/inputs_perlin_ERF_REMORA \
  max_step=1 \
  erf.rrtmgp_file_path="${RRTMGP_DIR}"
```

### 2. Upwelling MPI rank comparison

The equal-spacing configuration uses `Exec/Upwelling/inputs_perlin_ERF_REMORA`:

- `driver.atm.n_cell = 50 20 48`
- `driver.ocn.n_cell = 50 20 40`

Create output directories:

```bash
mkdir -p /tmp/upwelling_mpi/n1 /tmp/upwelling_mpi/n2 \
         /tmp/upwelling_mpi/n4 /tmp/upwelling_mpi/n6 /tmp/upwelling_mpi/n8
```

Run at each MPI size. Shown for `n=1`; repeat with `-n 2`, `-n 4`, `-n 6`,
`-n 8`, changing only the output directory:

```bash
mpirun -n 1 ./build/erf-remora-driver Exec/Upwelling/inputs_perlin_ERF_REMORA \
  max_step=10 erf.rrtmgp_file_path="${RRTMGP_DIR}" \
  erf.plot_file_1=/tmp/upwelling_mpi/n1/erfplt \
  remora.plot_file=/tmp/upwelling_mpi/n1/remoraplt \
  driver.plot_file_exchange_a2o=/tmp/upwelling_mpi/n1/exchange_a2o \
  driver.plot_file_exchange_o2a=/tmp/upwelling_mpi/n1/exchange_o2a
```

Compare each run against the `n=1` baseline:

```bash
"${FCOMPARE_EXE}" -a /tmp/upwelling_mpi/n1/remoraplt00010 /tmp/upwelling_mpi/n2/remoraplt00010
"${FCOMPARE_EXE}" -a /tmp/upwelling_mpi/n1/remoraplt00010 /tmp/upwelling_mpi/n4/remoraplt00010
"${FCOMPARE_EXE}" -a /tmp/upwelling_mpi/n1/remoraplt00010 /tmp/upwelling_mpi/n6/remoraplt00010
"${FCOMPARE_EXE}" -a /tmp/upwelling_mpi/n1/remoraplt00010 /tmp/upwelling_mpi/n8/remoraplt00010
```

Expected result for the equal-spacing configuration: `PLOTFILE AGREE`.

For the complete rank-matrix automation see
`scripts/run_rank_matrix_and_fcompare.sh`.

### 3. Upwelling OpenMP comparison

```bash
export OMP_NUM_THREADS=4
mpirun -n 1 ./build_omp/erf-remora-driver Exec/Upwelling/inputs_perlin_ERF_REMORA \
  max_step=10 erf.rrtmgp_file_path="${RRTMGP_DIR}"
```

To compare the MPI-only and OpenMP results:

```bash
"${FCOMPARE_EXE}" -a /tmp/upwelling_mpi/n1/remoraplt00010 \
                     /tmp/upwelling_omp/n1/remoraplt00010
```

### 4. Restart validation

The restart test checks whether a cold-start run and a restarted run produce
the same solution at a later step.

Cold-start run:

```bash
mpirun -n 1 ./build/erf-remora-driver \
  Exec/Upwelling/inputs_debug_cold_start_substep \
  erf.rrtmgp_file_path="${RRTMGP_DIR}"
```

Restarted run (uses `erf.restart` to resume from the checkpoint written by the
cold-start run):

```bash
mpirun -n 1 ./build/erf-remora-driver \
  Exec/Upwelling/inputs_debug_restart_substep \
  erf.rrtmgp_file_path="${RRTMGP_DIR}"
```

Compare the final plotfiles:

```bash
cd Exec/Upwelling
"${FCOMPARE_EXE}" -a cold_remoraplt00002  remoraplt00002
"${FCOMPARE_EXE}" -a cold_erfplt00020     restart_erfplt00020
```

Expected behavior:

- REMORA plotfiles agree exactly (bitwise).
- ERF plotfiles agree to floating-point roundoff; very small differences in
  fields such as velocity or pressure are normal and expected.

### 5. Flux-coupling debug run

This short run exercises `driver.atm2ocn_mode=flux` on the Upwelling debug
deck. The coefficient overrides below are important for this setup because the
base `inputs_ERF` file sets `erf.most.Cd = 0`, `erf.most.Ch = 0`, and
`erf.most.Cq = 0`, which suppresses ERF surface stress and heat/moisture
fluxes.

```bash
mpirun -n 1 ./build/erf-remora-driver \
  Exec/Upwelling/inputs_debug_cold_start_substep \
  driver.atm2ocn_mode=flux \
  driver.plot_int_exchange=1 \
  erf.rrtmgp_file_path="${RRTMGP_DIR}" \
  erf.most.Cd=0.0014 \
  erf.most.Ch=0.0012 \
  erf.most.Cq=0.0012
```

Expected behavior:

- step 0 bootstraps REMORA from ERF state lanes
- step 1 exports live flux lanes including stress, sensible heat, latent heat,
  longwave radiation, and evaporation
- the driver writes exchange diagnostics when `driver.plot_int_exchange=1`

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
# Driver control
driver.synthetic_only              = 0
driver.coupling_direction          = two_way
driver.atm2ocn_mode                = state
driver.dt_atm                      = 2.0
driver.atm_substeps_per_ocean_step = 5
driver.ocean_steps                 = 1

# Driver domain topology
driver.atm.prob_lo = 0.0 0.0  0.0
driver.atm.prob_hi = 1.0 1.0  1.0
driver.atm.n_cell  = 32 32 32
driver.ocn.prob_lo = 0.0 0.0 -1.0
driver.ocn.prob_hi = 1.0 1.0  0.0
driver.ocn.n_cell  = 32 32 10
amr.max_level      = 0
amr.max_grid_size  = 32

# ERF controls
erf.fixed_dt    = 0.5
erf.les_type    = None
erf.use_gravity = false

# REMORA controls
remora.ic_type         = analytic
remora.grid_scale_type = constant
remora.sum_interval    = -1
```

Run with:

```bash
./build/erf-remora-driver /path/to/inputs
```
