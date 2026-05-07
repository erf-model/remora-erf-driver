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
