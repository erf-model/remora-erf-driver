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
InitREMORASolver(void* remora_handle)
{
    static_cast<REMORA*>(remora_handle)->InitData();
}

void
REMORAPackSurfaceState(void* remora_handle,
                       amrex::Vector<amrex::MultiFab*>& state,
                       amrex::Real time)
{
    static_cast<REMORA*>(remora_handle)->PackSurfaceState(state, time);
}

void
REMORAApplyAtmosphericStates(void* remora_handle,
                             const amrex::Vector<amrex::MultiFab*>& states,
                             amrex::Real time)
{
    static_cast<REMORA*>(remora_handle)->ApplyAtmosphericStates(states, time);
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
