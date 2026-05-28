# Installation and Reproducer Guide

## Quick Start
```bash
git clone --recursive git@github.com:erf-model/remora-erf-driver
cd remora-erf-driver
cmake -S . -B build && cmake --build build -j4 --target erf-remora-driver
./build/erf-remora-driver Exec/CouplingReproducers/inputs_ERFtoREMORA_coupleToERF max_step=10
```

## Coupled ABL Reproducer
```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_two_way_ABL
```

## Build



From repository root:

```bash
cmake -S . -B build
cmake --build build -j4 --target erf-remora-driver
```

If you update files in the ERF or REMORA source trees used by this build,
rebuild the same target:

```bash
cmake --build build -j4 --target erf-remora-driver
```

## Ordered Reproducer Set

### Reproducer A: Coupling infrastructure-only control (`driver.synthetic_only = 1`)

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_coupling_infrastructure_check
```

Scientific intent:

- verify stable atmosphere/ocean coupling cadence and subcycling bookkeeping,
- verify cache slab ownership/update path without solver-physics confounders,
- establish a clean control case for subsequent physical coupling runs.

### Reproducer B: One-step standalone solver initialization (`driver.synthetic_only = 0`)

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_one_step_solver_initialization_check \
  max_step=1
```

Scientific intent:

- verify both solvers can initialize and advance prognostic state in the shared executable,
- separate initialization/stability issues from coupling-path issues.

### Reproducer C: `ERFtoREMORA` state-exchange baseline

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_ERFtoREMORA_coupleToERF
```

Scientific intent:

- verify atmosphere-state transfer and ocean application at each coupled step,
- establish baseline coupled behavior before enabling SST feedback into ERF.

### Reproducer D: `two_way` SST feedback into ERF (`OceanSurf` active)

Use the dedicated OceanSurf deck plus longer runtime:

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_two_way_coupleToERF \
  max_step=600
```

Scientific intent:

- verify SST from REMORA is injected into persistent ERF OceanSurf `t_surf` (K),
- verify ERF surface-layer fluxes respond to ocean state (`Cd`, `Ch`, `Cq` nonzero),
- produce measurable thermodynamic response in ERF (`temp`, `theta`) over long windows.

### Quick validation diagnostics

```bash
~/codes/ERF/tools/fcompare_serial -d temp plt00000 plt00600
~/codes/ERF/tools/fcompare_serial -d theta plt00000 plt00600
```

Interpretation:

- machine-zero drift in `temp/theta` suggests negligible SST-to-ERF feedback,
- O(1e-1 K) scale over 600 s (case-dependent) indicates feedback is dynamically active.

### Dependency notes

Coupled decks require:

- REMORA source: `Exec/CoupleToERF/inputs`
- ERF source: `.Exec_dev/RemoraCoupling/inputs_shear_ERFtoREMORA` or
  `.Exec_dev/RemoraCoupling/inputs_shear_two_way`
- ERF source: `.Exec_dev/RemoraCoupling/input_sounding`

`driver.atm2ocn_mode=state` is the validated mode. `driver.atm2ocn_mode=flux`
is still disabled in this driver.

## Unified Inputs Example

AMReX `ParmParse` allows one inputs file to provide parameters for the driver
 and both solver libraries.

How keys are consumed:

- `driver.*` keys are read by `erf-remora-driver`.
- `erf.*` (and shared mesh keys like `amr.*`, `geometry.*`) are read by ERF.
- `remora.*` (and shared mesh keys like `amr.*`, `geometry.*`) are read by REMORA.

Example mixed inputs file:

```text
# Driver control
driver.synthetic_only = 0
driver.coupling_direction = two_way
driver.atm2ocn_mode = state
driver.dt_atm = 2.0
driver.atm_substeps_per_ocean_step = 5
driver.ocean_steps = 1

# Driver topology controls
driver.atm.prob_lo = 0.0 0.0 0.0
driver.atm.prob_hi = 1.0 1.0 1.0
driver.atm.n_cell = 32 32 32
driver.ocn.prob_lo = 0.0 0.0 -1.0
driver.ocn.prob_hi = 1.0 1.0 0.0
driver.ocn.n_cell = 32 32 10
amr.max_level = 0
amr.max_grid_size = 32

# ERF controls
erf.fixed_dt = 0.5
erf.les_type = None
erf.use_gravity = false

# REMORA controls
remora.ic_type = analytic
remora.grid_scale_type = constant
remora.sum_interval = -1
```

Run with:

```bash
./build/erf-remora-driver /path/to/inputs
```
