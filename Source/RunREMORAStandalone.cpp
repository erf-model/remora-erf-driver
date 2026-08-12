#include "StandaloneSolverRuns.H"

#include <REMORA.H>

namespace erf_remora
{
void
RunREMORAStandalone()
{
    REMORA remora;
    remora.InitData();
    remora.Evolve();
}

void*
CreateREMORASolver(const ERFRemoraMultiBlockContainer& container)
{
    return static_cast<void*>(new REMORA(
        container.GetOcnRealBox(),
        container.GetMaxLevel(),
        container.GetOcnNCell(),
        container.GetCoord(),
        container.GetRefRatio(),
        container.GetIsPeriodic(),
        "remora"));
}

void
DestroyREMORASolver(void* remora_handle) noexcept
{
    delete static_cast<REMORA*>(remora_handle);
}

void
ConfigureREMORADriverAtmosToOceanCoupling(void* remora_handle,
                                          bool use_coupling_driver,
                                          bool use_two_way_coupling,
                                          bool use_state_contract)
{
    const auto mode = use_state_contract
        ? REMORA::DriverAtmosForcingMode::State
        : REMORA::DriverAtmosForcingMode::Flux;
    static_cast<REMORA*>(remora_handle)->ConfigureDriverAtmosToOceanCoupling(
        use_coupling_driver, use_two_way_coupling, mode);
}

void
SetREMORAStartTime(void* remora_handle, amrex::Real time)
{
    static_cast<REMORA*>(remora_handle)->SetStartTime(time);
}

void
SetREMORADriverAtmosToOceanStateContract(void* remora_handle,
                                         bool use_state_contract)
{
    const auto mode = use_state_contract
        ? REMORA::DriverAtmosForcingMode::State
        : REMORA::DriverAtmosForcingMode::Flux;
    static_cast<REMORA*>(remora_handle)->SetDriverAtmosToOceanForcingMode(mode);
}

void
InitREMORASolver(void* remora_handle)
{
    static_cast<REMORA*>(remora_handle)->InitData();
}

void
REMORAGetAtmosToOceanRhoLayout(void* remora_handle,
                               amrex::BoxArray& ba,
                               amrex::DistributionMapping& dm)
{
    static_cast<REMORA*>(remora_handle)->GetAtmosToOceanRhoLayout(ba, dm);
}

void
REMORAGetAtmosToOceanUFaceLayout(void* remora_handle,
                                 amrex::BoxArray& ba,
                                 amrex::DistributionMapping& dm)
{
    static_cast<REMORA*>(remora_handle)->GetAtmosToOceanUFaceLayout(ba, dm);
}

void
REMORAGetAtmosToOceanVFaceLayout(void* remora_handle,
                                 amrex::BoxArray& ba,
                                 amrex::DistributionMapping& dm)
{
    static_cast<REMORA*>(remora_handle)->GetAtmosToOceanVFaceLayout(ba, dm);
}

void
REMORAGetAtmosToOceanPsiCoordinates(void* remora_handle,
                                    const amrex::MultiFab** x_psi,
                                    const amrex::MultiFab** y_psi)
{
    const amrex::MultiFab* xp = nullptr;
    const amrex::MultiFab* yp = nullptr;
    static_cast<REMORA*>(remora_handle)->GetAtmosToOceanPsiCoordinates(xp, yp);
    if (x_psi) { *x_psi = xp; }
    if (y_psi) { *y_psi = yp; }
}

void
REMORAGetMasks(void* remora_handle,
              const amrex::MultiFab** mskr,
              const amrex::MultiFab** msku,
              const amrex::MultiFab** mskv)
{
    const amrex::MultiFab* r = nullptr;
    const amrex::MultiFab* u = nullptr;
    const amrex::MultiFab* v = nullptr;
    static_cast<REMORA*>(remora_handle)->GetLandSeaMasks(r, u, v);
    if (mskr) { *mskr = r; }
    if (msku) { *msku = u; }
    if (mskv) { *mskv = v; }
}

void
REMORAPackSurfaceState(void* remora_handle,
                       amrex::Vector<amrex::MultiFab*>& state,
                       amrex::Real time,
                       const amrex::MultiFab* weight_o2a_mf,
                       const amrex::iMultiFab* index_o2a_mf,
                       int max_stencil_size,
                       const amrex::iMultiFab* dst_land_mask)
{
    static_cast<REMORA*>(remora_handle)->PackSurfaceState(
        state, time, weight_o2a_mf, index_o2a_mf, max_stencil_size, dst_land_mask);
}

void
REMORAApplyAtmosphericStates(void* remora_handle,
                             const amrex::Vector<amrex::MultiFab*>& states,
                             amrex::Real time)
{
    static_cast<REMORA*>(remora_handle)->ApplyAtmosphericStates(states, time);
}

void
REMORAApplyAtmosphericFluxes(void* remora_handle,
                             const amrex::Vector<amrex::MultiFab*>& states,
                             amrex::Real time)
{
    static_cast<REMORA*>(remora_handle)->ApplyAtmosphericFluxes(states, time);
}

amrex::Real
REMORAEvolveOneStep(void* remora_handle, amrex::Real time, amrex::Real dt_request)
{
    return static_cast<REMORA*>(remora_handle)->EvolveOneStep(time, dt_request);
}

void
REMORAWriteAtFinalTime(void* remora_handle)
{
    return static_cast<REMORA*>(remora_handle)->WriteAtFinalTime();
}

void
RunREMORAApplyAndOneStep(const ERFRemoraMultiBlockContainer& container,
                          const amrex::Vector<amrex::MultiFab*>& states, amrex::Real time)
{
    static_cast<void>(container);
    REMORA remora;
    remora.InitData();
    remora.ApplyAtmosphericStates(states, time);
    remora.Evolve();
}
} // namespace erf_remora
