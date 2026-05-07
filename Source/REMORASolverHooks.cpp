#include "ERFRemoraSolverHooks.H"

#include <REMORA.H>

namespace erf_remora::hooks
{
REMORA*
CreateREMORA()
{
    return new REMORA();
}

void
DestroyREMORA(REMORA* remora)
{
    delete remora;
}

void
REMORAPackSurfaceState(REMORA& remora, amrex::Vector<amrex::MultiFab*>& state, amrex::Real time)
{
    remora.PackSurfaceState(state, time);
}

void
REMORAApplyAtmosphericStates(REMORA& remora, const amrex::Vector<amrex::MultiFab*>& states, amrex::Real time)
{
    remora.ApplyAtmosphericStates(states, time);
}
} // namespace erf_remora::hooks
