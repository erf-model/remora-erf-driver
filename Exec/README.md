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
  max_step=1
```

What to verify:

- ERF and REMORA libraries can initialize and advance real solver state in the
  same executable,
- one-step standalone evolution path remains stable.

Expected outputs:

- `erfplt00000`, `erfplt00001`, `remoraplt00000`, `remoraplt00001`

---

### `CouplingReproducers/inputs_ERFtoREMORA_coupleToERF`

Purpose:

- `ERFtoREMORA` baseline (`driver.coupling_direction = ERFtoREMORA`) using the
  REMORA `CoupleToERF` setup. ERF advances one atmospheric step, packs cache
  slabs, then REMORA applies the transferred state and advances one ocean step
  without feeding REMORA SST back into active ERF surface physics.

Run:

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_ERFtoREMORA_coupleToERF
```

Dependency note:

- This input uses:
  - REMORA source: `Exec/CoupleToERF/inputs`
  - ERF source: `.Exec_dev/RemoraCoupling/inputs_shear_ERFtoREMORA`
  - ERF source: `.Exec_dev/RemoraCoupling/input_sounding`
- It requires those files to exist in your respective solver source checkouts.

---

### `CouplingReproducers/inputs_two_way_coupleToERF`

Purpose:

- `two_way` ERF/REMORA run with OceanSurf enabled so SST from the ocean cache
  is copied into ERF LSM persistent `t_surf` and used by the ERF surface layer.

Run:

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_two_way_coupleToERF \
  max_step=600
```

What to verify:

- ocean SST cache is applied into OceanSurf state each ocean step,
- ERF surface-layer transfer is SST-sensitive (`erf.most.Cd/Ch/Cq` are set in
  the ERF RemoraCoupling base deck),
- ERF `temp/theta` show measurable evolution over long windows.

Dependency note:

- This input uses the same `FILE` dependencies as
  `inputs_ERFtoREMORA_coupleToERF`, except that it selects the
  `inputs_shear_two_way` ERF base deck.

---

### `CouplingReproducers/inputs_two_way_ABL`

Purpose:

- two-way ABL reproducer (`driver.coupling_direction = two_way`) where REMORA
  SST feeds the ERF OceanSurf surface temperature path.

Run:

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_two_way_ABL
```

---

### `CouplingReproducers/inputs_two_way_couette`

Purpose:

- two-way Couette reproducer (`driver.coupling_direction = two_way`) for short
  SST-feedback validation against the idealized Couette atmosphere/ocean setup.

Run:

```bash
./build/erf-remora-driver Exec/CouplingReproducers/inputs_two_way_couette
```
