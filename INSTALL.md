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

### Reproducer C: One-way coupled ERF shear deck -> REMORA one-step apply/evolve

Run:

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_one_step_coupled_erf_shear_oneway
```

This reproducer verifies:

- ERF one-step path can use the ERF `RemoraCoupling/inputs_shear_oneway` deck,
- ERF `PackAtmosphericStates` -> REMORA `ApplyAtmosphericStates` bridge runs in
  one process,
- REMORA can advance one step after applying cache slabs.

Dependency note:

- The inputs file uses:
  - `FILE = external-src/ERF/.Exec_dev/RemoraCoupling/inputs_shear_oneway`
  - `erf.input_sounding_file = external-src/ERF/.Exec_dev/RemoraCoupling/input_sounding`
- Those files must exist in your ERF checkout.

### Flag status: expected now vs prototype

Expected for current reproducible checks:

- Coupling infrastructure check: `driver.synthetic_only = 1`.
- One-step standalone smoke: `driver.synthetic_only = 0`, `max_step=1`,
  `erf.prob_name = "UserDefined"`, `remora.ic_type = analytic`,
  `remora.fixed_dt` and `remora.fixed_ndtfast_ratio` set to finite values.
- One-step smoke keeps `remora.plot_int = -1` for stability with the current
  placeholder BlankProblem path.

Prototype / staged (not full physics-coupled yet):

- `driver.atm2ocn_mode = flux` is currently for cache-layer routing checks;
  it does not yet represent a validated full flux-coupled physics path.
- In `driver.synthetic_only = 0` mode, default behavior remains standalone
  ERF then standalone REMORA smoke paths.
- A minimal one-way bridge path is available with `driver.coupled_one_step=1`:
  run one ERF step, pack cache slabs, apply them in REMORA, then run one
  REMORA step.
- Full multi-rate two-way runtime coupling in this mode is still not enabled yet.

### Problem selection note for standalone smoke tests

When `driver.synthetic_only = 0`, each solver still chooses its own problem setup:

- ERF uses `erf.prob_name`.
- REMORA uses `remora.ic_type` plus its own init controls.

Current smoke deck notes:

- `Exec/CouplingReproducers/inputs_one_step_solver_initialization_check` sets
  `driver.synthetic_only = 0` and sets
  `erf.prob_name = "UserDefined"` to use ERF's supported neutral placeholder
  problem path for this one-step initialization smoke test.
- REMORA requires valid timestep controls in this path; the deck sets
  `remora.fixed_dt` and `remora.fixed_ndtfast_ratio` to avoid runaway `DT = 1e+20`.
- This smoke input disables REMORA plotfile output (`remora.plot_int = -1`) so the
  one-step initialization test stays stable with the current placeholder BlankProblem
  initialization path.

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
