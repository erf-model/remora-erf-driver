#include "StandaloneSolverRuns.H"

#include <ERF.H>

namespace erf_remora
{
void*
CreateERFSolver(const ERFRemoraMultiBlockContainer& container)
{
    return static_cast<void*>(new ERF(container.GetAtmRealBox(),
                                      container.GetMaxLevel(),
                                      container.GetAtmNCell(),
                                      container.GetCoord(),
                                      container.GetRefRatio(),
                                      container.GetIsPeriodic(),
                                      "erf"));
}

const amrex::BoxArray&
ERF_Grids(void* erf_handle) noexcept
{
    ERF* my_erf = (ERF*)(erf_handle);
    return my_erf->boxArray(0);
}

const amrex::DistributionMapping&
ERF_DMap(void* erf_handle) noexcept
{
    ERF* my_erf = (ERF*)(erf_handle);
    return my_erf->DistributionMap(0);
}


void
DestroyERFSolver(void* erf_handle) noexcept
{
    delete static_cast<ERF*>(erf_handle);
}

void
ConfigureERFDriverAtmosToOceanCoupling(void* erf_handle,
                                       bool use_coupling_driver,
                                       bool use_two_way_coupling,
                                       bool use_state_contract)
{
    static_cast<ERF*>(erf_handle)->ConfigureDriverAtmosToOceanCoupling(
        use_coupling_driver, use_two_way_coupling, use_state_contract);
}

void
SetERFDriverAtmosToOceanStateContract(void* erf_handle,
                                      bool use_state_contract)
{
    static_cast<ERF*>(erf_handle)->SetDriverAtmosToOceanStateContract(use_state_contract);
}

void
InitERFSolver(void* erf_handle)
{
    static_cast<ERF*>(erf_handle)->InitData();
}

void
ERFGetOceanToAtmosSurfaceLayout(void* erf_handle,
                                amrex::BoxArray& ba,
                                amrex::DistributionMapping& dm)
{
    static_cast<ERF*>(erf_handle)->GetOceanToAtmosSurfaceLayout(ba, dm);
}

void
ERFGetOceanToAtmosCellCenterLonLat(void* erf_handle,
                                   const amrex::MultiFab** lon_cc,
                                   const amrex::MultiFab** lat_cc)
{
    const amrex::MultiFab* lon = nullptr;
    const amrex::MultiFab* lat = nullptr;
    static_cast<ERF*>(erf_handle)->GetOceanToAtmosCellCenterLonLat(lon, lat);
    if (lon_cc) { *lon_cc = lon; }
    if (lat_cc) { *lat_cc = lat; }
}

void
ERFGetLandMask(void* erf_handle, const amrex::iMultiFab** is_land)
{
    const amrex::iMultiFab* lmask = nullptr;
    static_cast<ERF*>(erf_handle)->GetLandMask(lmask);
    if (is_land) { *is_land = lmask; }
}

void
ERFGetFillableAtmosToOceanLanes(void* erf_handle,
                                bool state_contract,
                                amrex::Vector<int>& can_fill)
{
    if (erf_handle == nullptr) {
        // No atmosphere fills nothing. Claiming otherwise would have the driver
        // hand REMORA a view of prefill.
        can_fill.clear();
        return;
    }
    static_cast<ERF*>(erf_handle)->GetFillableAtmosToOceanLanes(state_contract, can_fill);
}

void
ERFApplyOceanSurfaceState(void* erf_handle,
                          const amrex::Vector<amrex::MultiFab*>& state,
                          amrex::Real time,
                          const amrex::iMultiFab* erf_coverage)
{
    static_cast<ERF*>(erf_handle)->ApplyOceanSurfaceState(state, time, erf_coverage);
}

void
ERFPackAtmosphericStates(void* erf_handle,
                         amrex::Vector<amrex::MultiFab*>& states,
                         amrex::Real time,
                         const amrex::Vector<const amrex::MultiFab*>& weight_mf_by_family,
                         const amrex::Vector<const amrex::iMultiFab*>& index_mf_by_family,
                         int max_stencil_size,
                         const amrex::MultiFab* dst_mskr,
                         const amrex::MultiFab* dst_msku,
                         const amrex::MultiFab* dst_mskv)
{
    static_cast<ERF*>(erf_handle)->PackAtmosphericStates(
        states, time, weight_mf_by_family, index_mf_by_family, max_stencil_size,
        dst_mskr, dst_msku, dst_mskv);
}

amrex::Real
ERFEvolveOneStep(void* erf_handle, amrex::Real time, amrex::Real dt_request)
{
    return static_cast<ERF*>(erf_handle)->EvolveOneStep(time, dt_request);
}

void
ERFWriteAtFinalTime(void* erf_handle)
{
    return static_cast<ERF*>(erf_handle)->WriteAtFinalTime();
}
} // namespace erf_remora
