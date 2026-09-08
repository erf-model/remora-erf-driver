#include "StandaloneSolverRuns.H"

#include <REMORA.H>

#include <array>

namespace erf_remora
{
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
REMORAGetAtmosToOceanPsiLonLat(void* remora_handle,
                               const amrex::MultiFab** lon_psi,
                               const amrex::MultiFab** lat_psi)
{
    const amrex::MultiFab* lonp = nullptr;
    const amrex::MultiFab* latp = nullptr;
    static_cast<REMORA*>(remora_handle)->GetAtmosToOceanPsiLonLat(lonp, latp);
    if (lon_psi) { *lon_psi = lonp; }
    if (lat_psi) { *lat_psi = latp; }
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
REMORAGetDeckConfiguredAtmosLanes(void* remora_handle,
                                  amrex::Vector<int>& deck_configured)
{
    deck_configured.clear();
    if (remora_handle == nullptr) { return; }
    std::array<bool, AtmosState::NumTypes> configured{};
    static_cast<REMORA*>(remora_handle)->GetDeckConfiguredAtmosStateLanes(configured);
    deck_configured.resize(AtmosState::NumTypes, 0);
    for (int lane = 0; lane < AtmosState::NumTypes; ++lane) {
        deck_configured[lane] = configured[lane] ? 1 : 0;
    }
}

void
SetREMORALongwaveFromDriver(void* remora_handle)
{
    if (remora_handle == nullptr) { return; }
    static_cast<REMORA*>(remora_handle)->SetLongwaveFromDriver();
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
} // namespace erf_remora
