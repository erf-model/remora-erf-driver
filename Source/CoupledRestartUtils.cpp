#include "CoupledRestartUtils.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>

#include <sstream>
#include <string>

namespace
{
void
GotoNextLine (std::istream& is)
{
    constexpr std::streamsize ignore_max {100000};
    is.ignore(ignore_max, '\n');
}

int
ComponentHeaderLines (erf_remora::restart::SolverKind solver) noexcept
{
    switch (solver) {
    case erf_remora::restart::SolverKind::ERF:
        return 4;
    case erf_remora::restart::SolverKind::REMORA:
        return 8;
    }

    return 0;
}
}

namespace erf_remora::restart
{
const char*
SolverKindName (SolverKind solver) noexcept
{
    switch (solver) {
    case SolverKind::ERF:
        return "ERF";
    case SolverKind::REMORA:
        return "REMORA";
    }

    return "Unknown";
}

CheckpointHeaderMetadata
ReadCheckpointHeaderMetadata (const std::string& checkpoint_dir, SolverKind solver)
{
    const std::string header_file = checkpoint_dir + "/Header";
    if (!amrex::FileExists(header_file)) {
        amrex::Abort(std::string("Missing restart header for ") + SolverKindName(solver) +
                     ": " + header_file);
    }

    amrex::Vector<char> file_chars;
    amrex::ParallelDescriptor::ReadAndBcastFile(header_file, file_chars);
    if (file_chars.empty()) {
        amrex::Abort(std::string("Empty restart header for ") + SolverKindName(solver) +
                     ": " + header_file);
    }

    std::string contents(file_chars.dataPtr(), file_chars.size());
    std::istringstream is(contents, std::istringstream::in);

    CheckpointHeaderMetadata meta;
    meta.checkpoint_dir = checkpoint_dir;
    meta.solver = solver;

    std::string line;
    std::getline(is, line); // title
    if (line.empty()) {
        amrex::Abort(std::string("Malformed restart header title for ") + SolverKindName(solver) +
                     ": " + header_file);
    }

    is >> meta.finest_level;
    if (!is.good()) {
        amrex::Abort(std::string("Malformed finest_level in restart header for ") + SolverKindName(solver) +
                     ": " + header_file);
    }
    GotoNextLine(is);

    for (int n = 0; n < ComponentHeaderLines(solver); ++n) {
        int chk_ncomp = 0;
        is >> chk_ncomp;
        if (!is.good()) {
            amrex::Abort(std::string("Malformed component-count section in restart header for ") +
                         SolverKindName(solver) + ": " + header_file);
        }
        GotoNextLine(is);
    }

    std::getline(is, line);
    {
        std::istringstream lis(line);
        if (!(lis >> meta.step)) {
            amrex::Abort(std::string("Missing istep[0] in restart header for ") +
                         SolverKindName(solver) + ": " + header_file);
        }
    }

    std::getline(is, line);
    {
        std::istringstream lis(line);
        if (!(lis >> meta.dt)) {
            amrex::Abort(std::string("Missing dt[0] in restart header for ") +
                         SolverKindName(solver) + ": " + header_file);
        }
    }

    std::getline(is, line);
    {
        std::istringstream lis(line);
        if (!(lis >> meta.time)) {
            amrex::Abort(std::string("Missing t_new[0] in restart header for ") +
                         SolverKindName(solver) + ": " + header_file);
        }
    }

    return meta;
}
} // namespace erf_remora::restart
