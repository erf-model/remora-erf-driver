#ifndef ERF_REMORA_SOLVER_HOOKS_H_
#define ERF_REMORA_SOLVER_HOOKS_H_

#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

class ERF;
class REMORA;

namespace erf_remora::hooks
{
// Minimal-solver-intrusion boundary for ERF/REMORA coupling.
//
// Architectural mandate:
// - ERF and REMORA do not call each other directly.
// - Both solvers only pack/unpack AMReX slabs through this thin API.
// - The driver (`MultiBlockContainer`) owns orchestration, timing, and cache memory.
//
// Precedent:
// - `erf-amrwind-driver` conductor pattern (AMReX Multiblock indirection).
// - AMR-Wind/CFD orchestration style in `kynema-sgf/src/CFDSim.H`.
// - ERF side interface extension via methods declared in `erf-model/ERF/Source/ERF.H`.
//
// This keeps core fluid kernels isolated from coupling logic and reduces merge risk
// as Stage-1..4 capabilities are added.
ERF* CreateERF();
void DestroyERF(ERF* erf);

REMORA* CreateREMORA();
void DestroyREMORA(REMORA* remora);

void ERFPackAtmosphericStates(ERF& erf, amrex::Vector<amrex::MultiFab*>& states, amrex::Real time);
void ERFApplyOceanSurfaceState(ERF& erf, const amrex::Vector<amrex::MultiFab*>& state, amrex::Real time);

void REMORAPackSurfaceState(REMORA& remora, amrex::Vector<amrex::MultiFab*>& state, amrex::Real time);
void REMORAApplyAtmosphericStates(REMORA& remora, const amrex::Vector<amrex::MultiFab*>& states, amrex::Real time);
} // namespace erf_remora::hooks

#endif
