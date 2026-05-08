# ERF-REMORA Driver (Stage-0 Scaffold)

This repository hosts a top-level AMReX superbuild and a driver executable that
couples ERF and REMORA as libraries. The current sprint target is infrastructure
validation: 2D `MultiFab` cache ownership, read/write paths, and multi-rate clock
stability.

## Coupling Roadmap

### Stage-0 / Stage-1: Legacy COAWST State-Passing

The current interface contract follows the legacy COAWST state exchange defined in:

- Warner et al. (2010), *Development of a Coupled Ocean-Atmosphere-Wave-Sediment
  Transport (COAWST) Modeling System*, Fig. 5 (Block B semantics).

Current atmosphere -> ocean state list:

- `Uwind`, `Vwind` (10-m winds)
- `Patm`, `RH`, `Tair`, `cloud`, `rain`
- `SWrad`, `LWrad`

Current ocean -> atmosphere state list:

- `SST`

At this stage we use state-passing (legacy behavior): REMORA consumes atmospheric
state variables and computes surface fluxes internally via its bulk-physics path.

### Stage-2: Spatial Remapping (Uncommitted Scheme)

While Stage-0/1 relies on level-0, matched horizontal grids, Stage-2 will
introduce mismatched grids to support phenomena like coastal upwelling.

**Crucially, this architecture does not currently claim or target a specific
external interfacing or remapping scheme.** The design remains agnostic so we
can evaluate native AMReX multi-block capabilities (`ParallelCopy`,
`MultiBlockCommMetaData`) before deciding whether an external topological engine
is required for conservative area-weighted flux accounting.

### Stage-3 / Stage-4: Modern Flux-Passing + Wave Mediation

Planned transition:

- Add a runtime coupling-mode toggle (e.g., `ParmParse`) to switch between:
  - state-passing (legacy COAWST-like), and
  - direct flux-passing (ERF-computed stresses/heat/moisture fluxes).
- Extend exchanged fields for wave-mediated coupling (e.g., Stokes-drift-aware
  and wave-roughness-aware terms) when Stage-4 physics is activated.

Primary scientific references for this transition:

- Li et al. (2025), stochastic Ekman-Stokes coupling framework.
- Lattanzi et al. (ERF model papers, 2024/2025).
- Klion et al. (REMORA model paper, 2025).

## Architectural Rule

ERF and REMORA do not call each other directly. The driver owns:

- the coupling clock,
- the interface caches (`Vector<MultiFab*>`),
- and orchestration of pack/apply callbacks.

This mirrors the conductor pattern used in prior AMReX multiblock coupling work
and keeps solver-core changes minimal.

## Integration Mitigation Paths

When ERF and REMORA are linked into one executable, there are two practical
integration hazards:

- Header/translation-unit collisions from globally named symbols in public
  solver headers (for example utility functions, structs, and enums with the
  same names in both codebases).
- Link-time symbol collisions from duplicate global definitions.

Current mitigation path in this driver:

- Do not include both solver top-level headers in the same translation unit.
- Keep ERF-only and REMORA-only includes in separate `.cpp` files.
- Keep cross-solver exchange at a driver-owned API boundary (`Pack*`/`Apply*`).
- Keep duplicate global symbols file-local (`static`) where needed in solver
  repos for single-binary link safety.
- In library-only mode, isolate solver-specific `amrex_probinit` entry points.
  ERF and REMORA both declare this symbol as weak; in one executable, letting
  them share one global can route one solver to the wrong problem initializer.
- This repository now enforces explicit symbol prefixing in superbuild mode:
  - `ERFProblem`, `ERFProblemBase`, `ERFSolverChoice`, `erf_probinit`
  - `REMORAProblem`, `REMORAProblemBase`, `REMORASolverChoice`, `remora_probinit`

Policy details and PR-splitting guidance are documented in:

- `Docs/symbol-isolation.md`

Deferred/long-term mitigation path:

- Upstream namespace hygiene and broader symbol-scope cleanup in ERF/REMORA.
  This is intentionally deferred because it is high-impact and not required for
  Stage-0/1 infrastructure validation.

## Stage-0 Run Modes

- `driver.synthetic_only = 1` runs the driver-owned synthetic cache exchange
  only (no ERF/REMORA physics timestepping).
- `driver.synthetic_only = 0` runs an uncoupled smoke path:
  `RunERFStandalone()` followed by `RunREMORAStandalone()`.

Recommended commands:

```bash
# Synthetic cache/clock debug path only
./build/erf-remora-driver Exec/Stage0/inputs_driver_stage0

# Non-synthetic one-step smoke path (ERF then REMORA)
./build/erf-remora-driver Exec/Stage0/inputs_synthetic_stage0 \
  driver.synthetic_only=0 max_step=1 erf.v=0 remora.v=0
```

The second command should write:

- `erfplt00000`, `erfplt00001`
- `remoraplt00000`, `remoraplt00001`

## Problem Selection In Non-Synthetic Mode

In `driver.synthetic_only = 0`, each solver still chooses its own problem setup.

- ERF uses `erf.prob_name`.
- REMORA uses `remora.ic_type` plus its own init controls.

Notes for the current Stage-0 deck:

- `Exec/Stage0/inputs_synthetic_stage0` uses `erf.prob_name = "ERF-REMORA stage0"`.
  ERF prints a "problem name ... not known" warning and skips named perturbations;
  this is acceptable for the current one-step infrastructure smoke test.
- If you want a recognized ERF problem template, set `erf.prob_name` to a known
  ERF case (for example `"ABL"` or `"Ekman Spiral"`).
- REMORA requires a valid timestep in this path; `remora.fixed_dt` and
  `remora.fixed_ndtfast_ratio` are set in the Stage-0 deck to avoid a runaway
  `DT = 1e+20` condition.
- REMORA plotfiles require explicit plot variables; if `remora.plot_vars_3d`
  is empty, plotfile headers can have `nComp = 0` and AMRVis will fail to load.

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
driver.dt_atm = 2.0
driver.atm_substeps_per_ocean_step = 5
driver.ocean_steps = 1

# Shared AMReX mesh controls
amr.max_level = 0
amr.n_cell = 32 32 8
amr.max_grid_size = 32
geometry.prob_lo = 0.0 0.0 0.0
geometry.prob_hi = 1.0 1.0 1.0

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
