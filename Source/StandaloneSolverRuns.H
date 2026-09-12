#ifndef ERF_REMORA_STANDALONE_SOLVER_RUNS_H_
#define ERF_REMORA_STANDALONE_SOLVER_RUNS_H_

#include "ERFRemoraMultiBlockContainer.H"

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

namespace erf_remora
{
// Keep solver includes in separate translation units so ERF.H and REMORA.H
// never appear together in one TU.

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
void ERFGetOceanToAtmosCornerCoordinates(void* erf_handle,
                                         const amrex::MultiFab** x_corner,
                                         const amrex::MultiFab** y_corner);
// Land mask (1 = land, 0 = water) at ERF level-0, cell-centered on the rho grid.
void ERFGetLandMask(void* erf_handle, const amrex::iMultiFab** is_land);
// Which atmos->ocean lanes ERF's pack actually writes, as 0/1 in active-view
// slot order (9 state lanes, or 8 flux). A lane reported 0 would reach REMORA
// holding only the driver's prefill, which REMORA treats as authoritative.
// Defined next to ERF's pack so the two cannot drift apart.
void ERFGetFillableAtmosToOceanLanes(void* erf_handle,
                                     bool state_contract,
                                     amrex::Vector<int>& can_fill);
// erf_coverage is the per-cell ocean->atmos coverage flag on the layout returned by
// ERFGetOceanToAtmosSurfaceLayout: nonzero where the remap found enough wet
// REMORA rho-cell source area and at least one ocean donor for that ERF cell.
// Where it is zero, state[0] holds the remap's zero fill rather than an ocean
// value, and ERF keeps its own lower-boundary SST. Passing nullptr therefore
// means "nothing is covered", not "everything is".
void ERFApplyOceanSurfaceState(void* erf_handle,
                               const amrex::Vector<amrex::MultiFab*>& state,
                               amrex::Real time,
                               const amrex::iMultiFab* erf_coverage = nullptr);
// weight_mf_by_family/index_mf_by_family are keyed by WeightFamily slot order
// {SCALAR_CENTERED, TAUX_FACE, TAUY_FACE, WIND_FACE_TO_CENTER}.
void ERFPackAtmosphericStates(void* erf_handle,
                              amrex::Vector<amrex::MultiFab*>& states,
                              amrex::Real time,
                              const amrex::Vector<const amrex::MultiFab*>& weight_mf_by_family,
                              const amrex::Vector<const amrex::iMultiFab*>& index_mf_by_family,
                              int max_stencil_size,
                              const amrex::MultiFab* dst_mskr = nullptr,
                              const amrex::MultiFab* dst_msku = nullptr,
                              const amrex::MultiFab* dst_mskv = nullptr);
amrex::Real ERFEvolveOneStep(void* erf_handle, amrex::Real time, amrex::Real dt_request);
void ERFWriteAtFinalTime(void* erf_handle);

void* CreateREMORASolver(const ERFRemoraMultiBlockContainer& container);
void DestroyREMORASolver(void* remora_handle) noexcept;
void ConfigureREMORADriverAtmosToOceanCoupling(void* remora_handle,
                                               bool use_coupling_driver,
                                               bool use_two_way_coupling,
                                               bool use_state_contract);
void SetREMORAStartTime(void* remora_handle, amrex::Real time);
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
void REMORAGetAtmosToOceanPsiCoordinates(void* remora_handle,
                                         const amrex::MultiFab** x_psi,
                                         const amrex::MultiFab** y_psi);
// Spherical (degrees east/north) psi-point corner coordinates. Null when the
// REMORA grid file carried no lon_psi/lat_psi, which is the signal that the
// spherical remap frame is unavailable for this case.
void REMORAGetAtmosToOceanPsiLonLat(void* remora_handle,
                                    const amrex::MultiFab** lon_psi,
                                    const amrex::MultiFab** lat_psi);
// Wet/dry masks (1 = wet, 0 = land) at rho, u-face, and v-face points.
void REMORAGetMasks(void* remora_handle,
                    const amrex::MultiFab** mskr,
                    const amrex::MultiFab** msku,
                    const amrex::MultiFab** mskv);
void REMORAPackSurfaceState(void* remora_handle,
                            amrex::Vector<amrex::MultiFab*>& state,
                            amrex::Real time,
                            const amrex::MultiFab* weight_o2a_mf,
                            const amrex::iMultiFab* index_o2a_mf,
                            int max_stencil_size,
                            const amrex::iMultiFab* dst_land_mask = nullptr);
// Which atmos->ocean state lanes this run's REMORA deck configured explicitly,
// as 0/1 in state active-view slot order. Not recoverable from ParmParse: REMORA
// reads those keys with queryAdd, which materializes the defaults, so it records
// the answer during parsing instead.
void REMORAGetDeckConfiguredAtmosLanes(void* remora_handle,
                                       amrex::Vector<int>& deck_configured);
// Have REMORA consume the driver-supplied DOWNWARD longwave rather than
// computing it. Without this its default lwrad_type of "computed" ignores the
// lane and uses Berliand. Must precede REMORA InitData, which opens NetCDF
// forcing handles off the same selector this changes.
void SetREMORALongwaveFromDriver(void* remora_handle);
void REMORAApplyAtmosphericStates(void* remora_handle,
                                  const amrex::Vector<amrex::MultiFab*>& states,
                                  amrex::Real time);
void REMORAApplyAtmosphericFluxes(void* remora_handle,
                                  const amrex::Vector<amrex::MultiFab*>& states,
                                  amrex::Real time);
amrex::Real REMORAEvolveOneStep(void* remora_handle, amrex::Real time, amrex::Real dt_request);
void REMORAWriteAtFinalTime(void* remora_handle);
} // namespace erf_remora

#endif
