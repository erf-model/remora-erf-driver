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
} // namespace erf_remora
