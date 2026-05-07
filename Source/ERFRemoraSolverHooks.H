#ifndef ERF_REMORA_SOLVER_HOOKS_H_
#define ERF_REMORA_SOLVER_HOOKS_H_

#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

class ERF;
class REMORA;

namespace erf_remora::hooks
{
ERF* CreateERF();
void DestroyERF(ERF* erf);

REMORA* CreateREMORA();
void DestroyREMORA(REMORA* remora);

void ERFPackAtmosphericStates(ERF& erf, amrex::Vector<amrex::MultiFab*>& states, amrex::Real time);
void ERFApplyOceanSurfaceState(ERF& erf, const amrex::Vector<amrex::MultiFab*>& state, amrex::Real time);

void REMORAPackSurfaceState(REMORA& remora, amrex::Vector<amrex::MultiFab*>& state, amrex::Real time);
void REMORAApplyAtmosphericStates(REMORA& remora, const amrex::Vector<amrex::MultiFab*>& states, amrex::Real time);
} // namespace erf_remora::hooks

#endif
