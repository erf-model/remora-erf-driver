#include "ERFRemoraSolverHooks.H"

#include <ERF.H>

namespace erf_remora::hooks
{
ERF*
CreateERF()
{
    return new ERF();
}

void
DestroyERF(ERF* erf)
{
    delete erf;
}

void
ERFPackAtmosphericStates(ERF& erf, amrex::Vector<amrex::MultiFab*>& states, amrex::Real time)
{
    erf.PackAtmosphericStates(states, time);
}

void
ERFApplyOceanSurfaceState(ERF& erf, const amrex::Vector<amrex::MultiFab*>& state, amrex::Real time)
{
    erf.ApplyOceanSurfaceState(state, time);
}
} // namespace erf_remora::hooks
