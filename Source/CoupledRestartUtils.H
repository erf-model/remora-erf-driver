#ifndef ERF_REMORA_COUPLED_RESTART_UTILS_H_
#define ERF_REMORA_COUPLED_RESTART_UTILS_H_

#include <AMReX_REAL.H>

#include <string>

namespace erf_remora::restart
{
enum class SolverKind {
    ERF,
    REMORA
};

struct CheckpointHeaderMetadata
{
    std::string checkpoint_dir;
    SolverKind solver = SolverKind::ERF;
    int finest_level = 0;
    int step = 0;
    amrex::Real dt = amrex::Real(0.0);
    amrex::Real time = amrex::Real(0.0);
};

CheckpointHeaderMetadata ReadCheckpointHeaderMetadata(const std::string& checkpoint_dir,
                                                      SolverKind solver);

const char* SolverKindName(SolverKind solver) noexcept;
} // namespace erf_remora::restart

#endif
