# ERF-REMORA Driver

This repository builds a coupled driver that advances ERF (atmosphere) and
REMORA (ocean) in one executable. The driver owns the coupling clock, the
interface caches, and the pack/apply calls used to exchange state between the
two solvers.

For build and run instructions see [INSTALL.md](INSTALL.md).

## Architectural Rule

ERF and REMORA do not call each other directly. The driver owns:

- the coupling clock,
- the interface caches (`Vector<MultiFab*>`),
- and the orchestration of pack/apply callbacks.

This mirrors the conductor pattern used in AMReX multiblock coupling work and
keeps solver-core changes minimal.

## Coupling Contract

### Validated mode

`driver.atm2ocn_mode = state` is the only validated coupling path. In this
mode REMORA consumes atmospheric state variables and computes surface fluxes
internally via its bulk-physics path.

### Known limitations

- `driver.atm2ocn_mode = flux` is present in the driver but is **not** a
  documented or validated workflow.
- REMORA does not currently support vertical refinement in the z-direction in
  coupled mode; the z refinement ratio must be 1.
- Spatial remapping between non-matching horizontal grids is not yet
  implemented; the current baseline uses level-0, matched horizontal grids.

### Exchanged fields (COAWST state-passing baseline)

Atmosphere → ocean:
`Uwind`, `Vwind`, `Patm`, `Qair`, `Tair`, `cloud`, `rain`, `SWrad`, `LWrad`

Ocean → atmosphere:
`SST`

The atmosphere -> ocean cache is allocated as a superset layer vector that also
reserves direct-flux layers (`tau_x`, `tau_y`, `SHflux`, `LHflux`, `evap`)
for a future flux-passing mode without changing the driver-owned allocation
model.

## Test Problems

### CouplingReproducers

`Exec/CouplingReproducers` contains small coupled test cases for checking
driver behavior, state exchange, and two-way SST feedback:

- `inputs_coupling_infrastructure_check`
- `inputs_one_step_solver_initialization_check`
- `inputs_ERFtoREMORA_coupleToERF`
- `inputs_two_way_coupleToERF`

### Upwelling

`Exec/Upwelling` contains the coastal Perlin/upwelling problem used validation testing as well as for MPI
rank, OpenMP, and restart verification:

- `inputs_perlin_ERF_REMORA`
- `inputs_debug_cold_start_substep`
- `inputs_debug_restart_substep`
- `inputs_debug_cold_start`
- `inputs_debug_restart`

## Coupling Roadmap

### Current implementation

The current interface contract follows the legacy COAWST state exchange
(Warner et al. 2010, Fig. 5 Block B semantics). Validated behavior includes
driver-owned coupled stepping, state-passing exchange at each coupled step,
and 1:1 cadence smoke tests with `driver.dt_atm = 1` and
`remora.fixed_dt = 1`.

### Planned extensions

- Runtime toggle between state-passing and direct flux-passing modes.
- Wave-mediated coupling terms (Stokes-drift, wave-roughness-aware fluxes).
- Conservative area-weighted spatial remapping for non-matching grids.

Primary references: 
- Warner et al. (2010), COAWST baseline and stated need for consistent flux coupling.
- Zambon et al. (2014), COAWST hurricane application with direct stress/heat exchange.
- COAWST implementation anchors:
  - `Master/mct_roms_wrf.h` (`ATM2OCN_FLUXES` and `BULK_FLUXES` pathways)
  - `ROMS/Nonlinear/atm2ocn_flux.F` (direct imported flux path)
  - `ROMS/Nonlinear/bulk_flux.F` (state-driven bulk flux path)
- Lattanzi et al. (ERF model papers, 2024/2025).
- Klion et al. (REMORA model paper, 2025).

## Integration Mitigations

When ERF and REMORA are linked into one executable, two hazard classes arise:

1. **Header/translation-unit collisions** — globally named symbols present in
   both solver public headers (utility functions, structs, enums with identical
   names).
2. **Link-time symbol collisions** — duplicate global definitions.

Current mitigations:

- ERF-only and REMORA-only includes are kept in separate translation units.
- Cross-solver exchange is confined to a driver-owned API boundary
  (`Pack*`/`Apply*`).
- Duplicate global symbols are made file-local (`static`) in solver repos
  where needed.
- `amrex_probinit` entry points are explicitly prefixed:
  `erf_probinit` / `remora_probinit`.
- Problem and solver-choice types are handled appropriately
