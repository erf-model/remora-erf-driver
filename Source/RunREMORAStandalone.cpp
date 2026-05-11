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

void
RunREMORAApplyAndOneStep(const amrex::Vector<amrex::MultiFab*>& states, amrex::Real time)
{
    REMORA remora;
    remora.InitData();
    remora.ApplyAtmosphericStates(states, time);
    remora.Evolve();
}
} // namespace erf_remora
