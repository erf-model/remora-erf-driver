#ifndef ERF_REMORA_COUPLING_CONTRACT_H_
#define ERF_REMORA_COUPLING_CONTRACT_H_

#include <AMReX_Vector.H>

#include <string>

namespace erf_remora::coupling
{
enum AtmosToOceanIndex : int {
    Uwind = 0,
    Vwind,
    Patm,
    RH,
    Tair,
    Cloud,
    Rain,
    SWrad,
    LWrad,
    NumAtmosToOceanStates
};

enum OceanToAtmosIndex : int {
    SST = 0,
    NumOceanToAtmosStates
};

inline const amrex::Vector<std::string>& AtmosToOceanNames ()
{
    static const amrex::Vector<std::string> names {
        "Uwind", "Vwind", "Patm", "RH", "Tair", "cloud", "rain", "SWrad", "LWrad"
    };
    return names;
}

inline const amrex::Vector<std::string>& OceanToAtmosNames ()
{
    static const amrex::Vector<std::string> names {"SST"};
    return names;
}
} // namespace erf_remora::coupling

#endif
