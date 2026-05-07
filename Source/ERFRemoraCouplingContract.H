#ifndef ERF_REMORA_COUPLING_CONTRACT_H_
#define ERF_REMORA_COUPLING_CONTRACT_H_

#include <AMReX_Vector.H>

#include <string>

namespace erf_remora::coupling
{
// Legacy contract reference:
//   Warner et al. (2010), "Development of a Coupled Ocean-Atmosphere-Wave-
//   Sediment Transport (COAWST) Modeling System", Fig. 5 / Block B.
//
// Stage-0 uses the same state-passing variable list as the WRF->ROMS path:
//  - 10-m winds: Uwind, Vwind [m s^-1]
//  - Atmospheric pressure: Patm [Pa]
//  - Relative humidity: RH [0-1 fraction]
//  - Air temperature: Tair [K]
//  - Cloud fraction: Cloud [0-1 fraction]
//  - Precipitation rate: Rain [kg m^-2 s^-1]  (project SI convention)
//    TODO(stage-1): verify/convert if upstream forcing uses legacy mm day^-1.
//  - Shortwave radiation: SWrad [W m^-2]
//  - Longwave radiation: LWrad [W m^-2]
//
// Ocean->atmosphere state:
//  - Sea-surface temperature: SST [K]
//
// Coupling mode note:
// - This enum is for legacy state-passing, where REMORA computes bulk fluxes
//   internally (COARE-style path).
// - TODO(stage-3/4): add a runtime toggle (e.g. ParmParse `driver.atm2ocn_mode`
//   or compile-time ATM2OCN_FLUXES) for direct flux-passing (tau_x, tau_y, heat/
//   moisture fluxes) to avoid cross-component flux inconsistencies.
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
