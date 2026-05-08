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
} // namespace erf_remora
