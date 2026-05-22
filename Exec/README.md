# Exec Reproducer Index

This directory contains ordered reproducer inputs for coupling controls and
short physical experiments in the ERF-REMORA driver.

## Inputs

### `CouplingReproducers/inputs_coupling_infrastructure_check`

Purpose:

- coupling-infrastructure control experiment through the driver
  (`driver.synthetic_only = 1`).

Run:

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_coupling_infrastructure_check
```

What to verify:

- coupling clock/subcycling progression,
- interface cache allocation and controlled slab updates,
- shared AMReX process lifecycle in the combined executable.

---

### `CouplingReproducers/inputs_one_step_solver_initialization_check`

Purpose:

- one-step uncoupled initialization path
  (`driver.synthetic_only = 0`) where ERF runs then REMORA runs.

Run:

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_one_step_solver_initialization_check \
  max_step=1 erf.v=0 remora.v=0
```

What to verify:

- ERF and REMORA libraries can initialize and advance real solver state in the
  same executable,
- one-step standalone evolution path remains stable.

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
  then REMORA applies the transferred state and runs one step.

Run:

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_one_step_coupled_erf_shear_oneway
```

Dependency note:

- This input uses:
  - ERF source: `.Exec_dev/RemoraCoupling/inputs_shear_oneway`
  - ERF source: `.Exec_dev/RemoraCoupling/input_sounding`
- It requires those files to exist in your ERF source checkout (expected once that
  branch lands on ERF `development`).

---

### `CouplingReproducers/inputs_one_step_coupled_coupleToERF`

Purpose:

- one-way coupled one-step path (`driver.coupled_one_step = 1`) using the REMORA
  `CoupleToERF` setup as the baseline coupled experiment.

Run:

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_one_step_coupled_coupleToERF
```

Dependency note:

- This input uses:
  - REMORA source: `Exec/CoupleToERF/inputs`
  - ERF source: `.Exec_dev/RemoraCoupling/inputs_shear_oneway`
  - ERF source: `.Exec_dev/RemoraCoupling/input_sounding`
- It requires those files to exist in your respective solver source checkouts.

---

### `CouplingReproducers/inputs_coupled_coupleToERF_oceansurf`

Purpose:

- coupled ERF/REMORA run with OceanSurf enabled so SST from the ocean cache is
  copied into ERF LSM persistent `t_surf` and used by the ERF surface layer.

Run:

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_coupled_coupleToERF_oceansurf \
  max_step=600
```

What to verify:

- ocean SST cache is applied into OceanSurf state each ocean step,
- ERF surface-layer transfer is SST-sensitive (`erf.most.Cd/Ch/Cq` are set),
- ERF `temp/theta` show measurable evolution over long windows.

Dependency note:

- This input uses the same `FILE` dependencies as
  `inputs_one_step_coupled_coupleToERF`.
