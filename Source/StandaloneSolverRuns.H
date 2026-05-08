#ifndef ERF_REMORA_STANDALONE_SOLVER_RUNS_H_
#define ERF_REMORA_STANDALONE_SOLVER_RUNS_H_

namespace erf_remora
{
// Keep solver includes in separate translation units so ERF.H and REMORA.H
// never appear together in one TU.
void RunERFStandalone();
void RunREMORAStandalone();
} // namespace erf_remora

#endif
