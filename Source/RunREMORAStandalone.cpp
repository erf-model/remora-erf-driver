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
RunREMORAApplyAndOneStep(const ERFRemoraMultiBlockContainer& container,
                          const amrex::Vector<amrex::MultiFab*>& states, amrex::Real time)
{
    REMORA remora(container.GetOcnRealBox(),
                  container.GetMaxLevel(),
                  container.GetOcnNCell(),
                  container.GetCoord(),
                  container.GetRefRatio(),
                  container.GetIsPeriodic(),
                  "remora");
    remora.InitData();
    remora.ApplyAtmosphericStates(states, time);
    remora.Evolve();
}
} // namespace erf_remora
