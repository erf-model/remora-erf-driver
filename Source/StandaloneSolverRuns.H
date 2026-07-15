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

const amrex::BoxArray&
ERF_Grids(void* erf_handle) noexcept;

const amrex::DistributionMapping&
ERF_DMap(void* erf_handle) noexcept;

// Persistent coupled-mode shims. Opaque handles keep ERF.H and REMORA.H out of
// driver container translation units.
void* CreateERFSolver(const ERFRemoraMultiBlockContainer& container);
void DestroyERFSolver(void* erf_handle) noexcept;
void ConfigureERFDriverAtmosToOceanCoupling(void* erf_handle,
                                            bool use_coupling_driver,
                                            bool use_two_way_coupling,
                                            bool use_state_contract);
void SetERFDriverAtmosToOceanStateContract(void* erf_handle,
                                           bool use_state_contract);
void InitERFSolver(void* erf_handle);
void ERFGetOceanToAtmosSurfaceLayout(void* erf_handle,
                                     amrex::BoxArray& ba,
                                     amrex::DistributionMapping& dm);
void ERFApplyOceanSurfaceState(void* erf_handle,
                               const amrex::Vector<amrex::MultiFab*>& state,
                               amrex::Real time);
void ERFPackAtmosphericStates(void* erf_handle,
                              amrex::Vector<amrex::MultiFab*>& states,
                              amrex::Real time);
amrex::Real ERFEvolveOneStep(void* erf_handle, amrex::Real time, amrex::Real dt_request);
void ERFWriteAtFinalTime(void* erf_handle);

void* CreateREMORASolver(const ERFRemoraMultiBlockContainer& container);
void DestroyREMORASolver(void* remora_handle) noexcept;
void ConfigureREMORADriverAtmosToOceanCoupling(void* remora_handle,
                                               bool use_coupling_driver,
                                               bool use_two_way_coupling,
                                               bool use_state_contract);
void SetREMORADriverAtmosToOceanStateContract(void* remora_handle,
                                              bool use_state_contract);
void InitREMORASolver(void* remora_handle);
void REMORAGetAtmosToOceanRhoLayout(void* remora_handle,
                                    amrex::BoxArray& ba,
                                    amrex::DistributionMapping& dm);
void REMORAGetAtmosToOceanUFaceLayout(void* remora_handle,
                                      amrex::BoxArray& ba,
                                      amrex::DistributionMapping& dm);
void REMORAGetAtmosToOceanVFaceLayout(void* remora_handle,
                                      amrex::BoxArray& ba,
                                      amrex::DistributionMapping& dm);
void REMORAPackSurfaceState(void* remora_handle,
                            amrex::Vector<amrex::MultiFab*>& state,
                            amrex::Real time);
void REMORAApplyAtmosphericStates(void* remora_handle,
                                  const amrex::Vector<amrex::MultiFab*>& states,
                                  amrex::Real time);
void REMORAApplyAtmosphericFluxes(void* remora_handle,
                                  const amrex::Vector<amrex::MultiFab*>& states,
                                  amrex::Real time);
amrex::Real REMORAEvolveOneStep(void* remora_handle, amrex::Real time, amrex::Real dt_request);
void REMORAWriteAtFinalTime(void* remora_handle);

// Coupled one-step path: container provides driver-owned topology so each solver
// is constructed via its explicit-topology constructor instead of parsing amr.*/geometry.*.
void RunERFOneStepAndPack(const ERFRemoraMultiBlockContainer& container,
                           amrex::Vector<amrex::MultiFab*>& states, amrex::Real time);
void RunREMORAApplyAndOneStep(const ERFRemoraMultiBlockContainer& container,
                               const amrex::Vector<amrex::MultiFab*>& states, amrex::Real time);
} // namespace erf_remora

#endif
