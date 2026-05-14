# Exec Reproducer Index

This directory contains inputs used to reproduce the current coupling and
standalone smoke milestones in the driver.

## Inputs

### `CouplingReproducers/inputs_coupling_infrastructure_check`

Purpose:

- coupling infrastructure reproducer through the driver
  (`driver.synthetic_only = 1`).

Run:

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_coupling_infrastructure_check
```

Verifies:

- coupling clock/subcycling progression,
- interface cache allocation and controlled slab updates,
- shared AMReX process lifecycle in the combined executable.

---

### `CouplingReproducers/inputs_one_step_solver_initialization_check`

Purpose:

- one-step uncoupled standalone smoke path
  (`driver.synthetic_only = 0`) where ERF runs then REMORA runs.

Run:

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_one_step_solver_initialization_check \
  max_step=1 erf.v=0 remora.v=0
```

Verifies:

- ERF and REMORA libraries can initialize and advance real solver state in the
  same executable,
- one-step non-crashing solver initialization/evolution path in one executable.

Expected outputs:

- `erfplt00000`, `erfplt00001`

Note:

- This smoke input intentionally sets `remora.plot_int = -1` for stability with
  current placeholder BlankProblem initialization.

---

### `CouplingReproducers/inputs_one_step_coupled_erf_shear_oneway`

Purpose:

- one-way coupled one-step path (`driver.coupled_one_step = 1`) where ERF runs
  one step with the ERF `RemoraCoupling/inputs_shear_oneway` setup, packs cache slabs,
  then REMORA applies the physical state and runs one step.

Run:

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_one_step_coupled_erf_shear_oneway
```

Dependency note:

- This input uses:
  - `FILE = external-src/ERF/.Exec_dev/RemoraCoupling/inputs_shear_oneway`
  - `erf.input_sounding_file = external-src/ERF/.Exec_dev/RemoraCoupling/input_sounding`
- It requires those files to exist in your ERF checkout (expected once that
  branch lands on ERF `development`).

---

### `CouplingReproducers/inputs_one_step_coupled_coupleToERF`

Purpose:

- one-way coupled one-step path (`driver.coupled_one_step = 1`) using the REMORA
  `CoupleToERF` problem setup.

Run:

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_one_step_coupled_coupleToERF
```

Dependency note:

- This input uses:
  - `FILE = external-src/REMORA/Exec/CoupleToERF/inputs`
  - `FILE = external-src/ERF/.Exec_dev/RemoraCoupling/inputs_shear_oneway`
  - `erf.input_sounding_file = external-src/ERF/.Exec_dev/RemoraCoupling/input_sounding`
- It requires those files to exist in your respective solver checkouts.

