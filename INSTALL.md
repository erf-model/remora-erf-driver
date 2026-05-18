# Installation and Usage Guide

## Build the superbuild

From repository root:

```bash
cmake -S . -B build
cmake --build build -j4
```

## How To Run (Reviewer Reproducers)

### Reproducer A: Coupling infrastructure test (`driver.synthetic_only = 1`)

Run:

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_coupling_infrastructure_check
```

Optional (exercise flux-layer selection with the same cache allocation):

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_coupling_infrastructure_check \
  driver.atm2ocn_mode=flux
```

This reproducer verifies:

- shared AMReX initialization/finalization in the combined executable,
- multi-rate coupling clock progression with atmosphere subcycling,
- native AMReX interface cache ownership and controlled read/write paths,
- coupling infrastructure stability without executing full ERF/REMORA fluid physics.

Expected behavior:

- console diagnostics for initial and updated cache values,
- ocean-step completion messages with atmosphere tick progression,
- clean process exit.

### Reproducer B: One-step standalone solver initialization tests (`driver.synthetic_only = 0`)

Run:

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_one_step_solver_initialization_check \
  max_step=1 erf.v=0 remora.v=0
```

This reproducer verifies:

- ERF library path can read inputs, initialize, and advance one real timestep,
- REMORA library path can read inputs, initialize, and advance one real timestep,
- both solver libraries can run in one executable without crashing.

Expected outputs:

- `erfplt00000`, `erfplt00001`

Optional quick check:

```bash
sed -n '1,20p' erfplt00001/Header
```

### Reproducer C: Coupled ERF shear deck -> REMORA coupled apply/evolve loop

Run:

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_one_step_coupled_coupleToERF \
  max_step=10 \
  driver.dt_atm=1 \
  driver.atm_substeps_per_ocean_step=1 \
  remora.fixed_dt=1 \
  erf.plot_int_1=1 \
  remora.plot_int=1
```

This reproducer verifies:

- ERF and REMORA both advance through the driver-owned coupled loop,
- atmospheric state packing and ocean application occur at each coupled step,
- REMORA advances after receiving cache slabs from ERF.

Dependency note:

- The inputs file uses:
  - `FILE = external-src/ERF/.Exec_dev/RemoraCoupling/inputs_shear_oneway`
  - `erf.input_sounding_file = external-src/ERF/.Exec_dev/RemoraCoupling/input_sounding`
- Those files must exist in your ERF checkout.

### Flag status: expected now vs reserved

Expected for current reproducible checks:

- Coupling infrastructure check: `driver.synthetic_only = 1`.
- One-step standalone smoke: `driver.synthetic_only = 0`, `max_step=1`,
  `erf.prob_name = "UserDefined"`, `remora.ic_type = analytic`,
  `remora.fixed_dt` and `remora.fixed_ndtfast_ratio` set to finite values.
- The standalone one-step smoke deck keeps `remora.plot_int = -1` for stability
  with the current placeholder BlankProblem path.

Current coupled runtime behavior:

- `driver.atm2ocn_mode = state` is the validated coupling mode.
- `driver.atm2ocn_mode = flux` is disabled in the driver.
- The coupled loop is driver-owned: ERF advances one atmospheric step, REMORA
  applies the packed state, then REMORA advances one ocean step.
- The accepted smoke-test cadence is 1:1 with `driver.dt_atm = 1`,
  `driver.atm_substeps_per_ocean_step = 1`, and `remora.fixed_dt = 1`.

### Problem selection note for standalone smoke tests

When `driver.synthetic_only = 0`, each solver still chooses its own problem setup:

- ERF uses `erf.prob_name`.
- REMORA uses `remora.ic_type` plus its own init controls.

Current smoke deck notes:

- `Exec/CouplingReproducers/inputs_one_step_solver_initialization_check` keeps
  the standalone solver initialization smoke test separate from coupled runs.
- REMORA requires valid timestep controls in that path; the deck sets
  `remora.fixed_dt` and `remora.fixed_ndtfast_ratio` to avoid runaway `DT = 1e+20`.
- The coupled reproducer now uses `Exec/CouplingReproducers/inputs_one_step_coupled_coupleToERF`
  with `max_step=10`, `driver.dt_atm=1`, `driver.atm_substeps_per_ocean_step=1`,
  and `remora.fixed_dt=1`.

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
