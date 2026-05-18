#ifndef ERF_REMORA_STANDALONE_SOLVER_RUNS_H_
#define ERF_REMORA_STANDALONE_SOLVER_RUNS_H_

#include "ERFRemoraMultiBlockContainer.H"

#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

namespace erf_remora
{
// Keep solver includes in separate translation units so ERF.H and REMORA.H
// never appear together in one TU.
void RunERFStandalone();
void RunREMORAStandalone();

// Persistent coupled-mode shims. Opaque handles keep ERF.H and REMORA.H out of
// driver container translation units.
void* CreateERFSolver(const ERFRemoraMultiBlockContainer& container);
void DestroyERFSolver(void* erf_handle) noexcept;
void InitERFSolver(void* erf_handle);
void ERFApplyOceanSurfaceState(void* erf_handle,
                               const amrex::Vector<amrex::MultiFab*>& state,
                               amrex::Real time);
void ERFPackAtmosphericStates(void* erf_handle,
                              amrex::Vector<amrex::MultiFab*>& states,
                              amrex::Real time);

void* CreateREMORASolver(const ERFRemoraMultiBlockContainer& container);
void DestroyREMORASolver(void* remora_handle) noexcept;
void InitREMORASolver(void* remora_handle);
void REMORAPackSurfaceState(void* remora_handle,
                            amrex::Vector<amrex::MultiFab*>& state,
                            amrex::Real time);
void REMORAApplyAtmosphericStates(void* remora_handle,
                                  const amrex::Vector<amrex::MultiFab*>& states,
                                  amrex::Real time);

// Coupled one-step path: container provides driver-owned topology so each solver
// is constructed via its explicit-topology constructor instead of parsing amr.*/geometry.*.
void RunERFOneStepAndPack(const ERFRemoraMultiBlockContainer& container,
                           amrex::Vector<amrex::MultiFab*>& states, amrex::Real time);
void RunREMORAApplyAndOneStep(const ERFRemoraMultiBlockContainer& container,
                               const amrex::Vector<amrex::MultiFab*>& states, amrex::Real time);
} // namespace erf_remora

#endif
