#include "StandaloneSolverRuns.H"

#include <ERF.H>

namespace erf_remora
{
void
RunERFStandalone()
{
    ERF erf;
    erf.InitData();
    erf.Evolve();
}

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
ERFApplyOceanSurfaceState(void* erf_handle,
                          const amrex::Vector<amrex::MultiFab*>& state,
                          amrex::Real time)
{
    static_cast<ERF*>(erf_handle)->ApplyOceanSurfaceState(state, time);
}

void
ERFPackAtmosphericStates(void* erf_handle,
                         amrex::Vector<amrex::MultiFab*>& states,
                         amrex::Real time)
{
    static_cast<ERF*>(erf_handle)->PackAtmosphericStates(states, time);
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

void
RunERFOneStepAndPack(const ERFRemoraMultiBlockContainer& container,
                      amrex::Vector<amrex::MultiFab*>& states, amrex::Real time)
{
    ERF erf(container.GetAtmRealBox(),
            container.GetMaxLevel(),
            container.GetAtmNCell(),
            container.GetCoord(),
            container.GetRefRatio(),
            container.GetIsPeriodic(),
            "erf");
    erf.InitData();
    erf.Evolve();
    erf.PackAtmosphericStates(states, time);
}
} // namespace erf_remora
