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

void
DestroyERFSolver(void* erf_handle) noexcept
{
    delete static_cast<ERF*>(erf_handle);
}

void
InitERFSolver(void* erf_handle)
{
    static_cast<ERF*>(erf_handle)->InitData();
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
