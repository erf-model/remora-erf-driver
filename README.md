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

### Forcing lane ownership

Each atmosphere-to-ocean lane has one source, resolved at startup and printed as
a provenance table. ERF supplies a lane by default; a lane the REMORA deck
configures explicitly, or one ERF cannot fill, is left to REMORA. Override with
`driver.atm2ocn_lane_source.<lane> = erf | remora | auto`; see `Exec/README.md`.

### Known limitations

- `driver.atm2ocn_mode = flux` has a validated startup bootstrap (see
  `INSTALL.md`), but is less exercised than the state path.
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

### Upwelling

`Exec/Upwelling` contains the coastal upwelling problem of Perlin et al. (2007),
<https://doi.org/10.1175/JPO3104.1>, used for validation as well as for MPI
rank, OpenMP, and restart verification. It needs no external data:

- `inputs_oneway` — `ERFtoREMORA`
- `inputs_twoway` — `two_way`, with REMORA SST feeding ERF's surface layer

### OCCAMS

`Exec/OCCAMS` contains the real coupled domain. These decks reference absolute
NERSC paths for their WRF and REMORA input files and are not runnable elsewhere
without editing them:

- `inputs_OCCAMS_smoke` — fast startup and exchange guard
- `inputs_OCCAMS_oneway` — `ERFtoREMORA`
- `inputs_OCCAMS_state` — `two_way`, state-passing contract
- `inputs_OCCAMS_flux` — `two_way`, flux-passing contract

`Exec/OCCAMS/run_occams.sh` runs them and summarizes the coupling diagnostics.

See `Exec/README.md` for what each deck covers and how to read driver output,
and `INSTALL.md` for build and verification workflows.

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
