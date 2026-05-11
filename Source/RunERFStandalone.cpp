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

void
RunERFOneStepAndPack(amrex::Vector<amrex::MultiFab*>& states, amrex::Real time)
{
    ERF erf;
    erf.InitData();
    erf.Evolve();
    erf.PackAtmosphericStates(states, time);
}
} // namespace erf_remora
