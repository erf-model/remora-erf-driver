# ERF-REMORA Driver (Initial Architecture Scaffold)

This repository hosts a top-level AMReX superbuild and a driver executable that
couples ERF and REMORA as libraries. The current focus is infrastructure
validation: 2D `MultiFab` cache ownership, read/write paths, multi-rate clock
stability, and physical state exchange.

## Coupling Roadmap

### Initial Architecture: Legacy COAWST State-Passing

The current interface contract follows the legacy COAWST state exchange defined in:

- Warner et al. (2010), *Development of a Coupled Ocean-Atmosphere-Wave-Sediment
  Transport (COAWST) Modeling System*, Fig. 5 (Block B semantics).

Current atmosphere -> ocean state list:

- `Uwind`, `Vwind` (10-m winds)
- `Patm`, `RH`, `Tair`, `cloud`, `rain`
- `SWrad`, `LWrad`

Current ocean -> atmosphere state list:

- `SST`

Note: REMORA does not currently support vertical refinement; the refinement ratio in the z-direction must be set to 1.

In the initial architecture we use state-passing (legacy behavior): REMORA consumes atmospheric
state variables and computes surface fluxes internally via its bulk-physics path.

Implementation note:

- The atmosphere -> ocean cache is now allocated as a superset layer vector
  (`Vector<MultiFab*>`) containing legacy state layers plus reserved direct-flux
  layers (`tau_x`, `tau_y`, `SHflux`, `LHflux`, `evap`).
- Coupling modes can swap active layer sets by index without changing the
  driver-owned cache allocation/ownership model.

### Future Spatial Remapping (Uncommitted Scheme)

While the current baseline uses level-0, matched horizontal grids, future work will
introduce mismatched grids to support phenomena like coastal upwelling.

**Crucially, this architecture does not currently claim or target a specific
external interfacing or remapping scheme.** The design remains agnostic so we
can evaluate native AMReX multi-block capabilities (`ParallelCopy`,
`MultiBlockCommMetaData`) before deciding whether an external topological engine
is required for conservative area-weighted flux accounting.

### Future Flux-Passing + Wave Mediation

Planned transition:

- Add a runtime coupling-mode toggle (e.g., `ParmParse`) to switch between:
  - state-passing (legacy COAWST-like), and
  - direct flux-passing (ERF-computed stresses/heat/moisture fluxes).
- Extend exchanged fields for wave-mediated coupling (e.g., Stokes-drift-aware
  and wave-roughness-aware terms) when extended coupling physics is activated.

Primary scientific references for this transition:

- Warner et al. (2010), COAWST baseline and stated need for consistent flux coupling.
- Zambon et al. (2014), COAWST hurricane application with direct stress/heat exchange.
- COAWST implementation anchors:
  - `Master/mct_roms_wrf.h` (`ATM2OCN_FLUXES` and `BULK_FLUXES` pathways)
  - `ROMS/Nonlinear/atm2ocn_flux.F` (direct imported flux path)
  - `ROMS/Nonlinear/bulk_flux.F` (state-driven bulk flux path)
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
  initial architecture and one-step testing validation.

## How To Run (Reviewer Reproducers)

For detailed build and run instructions, including reproducible smoke tests, see [INSTALL.md](INSTALL.md).

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
