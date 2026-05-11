#ifndef ERF_REMORA_STANDALONE_SOLVER_RUNS_H_
#define ERF_REMORA_STANDALONE_SOLVER_RUNS_H_

#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

namespace erf_remora
{
// Keep solver includes in separate translation units so ERF.H and REMORA.H
// never appear together in one TU.
void RunERFStandalone();
void RunREMORAStandalone();
void RunERFOneStepAndPack(amrex::Vector<amrex::MultiFab*>& states, amrex::Real time);
void RunREMORAApplyAndOneStep(const amrex::Vector<amrex::MultiFab*>& states, amrex::Real time);
} // namespace erf_remora

#endif
