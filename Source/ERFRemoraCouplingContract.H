#ifndef ERF_REMORA_COUPLING_CONTRACT_H_
#define ERF_REMORA_COUPLING_CONTRACT_H_

#include <AMReX_Vector.H>

#include <string>

namespace erf_remora::coupling
{
// Legacy contract reference:
//   Warner et al. (2010), "Development of a Coupled Ocean-Atmosphere-Wave-
//   Sediment Transport (COAWST) Modeling System", Fig. 5 / Block B.
//   Zambon et al. (2014), Hurricane Ivan COAWST coupled configuration
//   (direct stress/heat flux transfer context for ATM2OCN_FLUXES-era coupling).
//
// COAWST implementation anchors for state-vs-flux split:
//   - Master/mct_roms_wrf.h (ATM2OCN_FLUXES vs BULK_FLUXES pathways)
//   - ROMS/Nonlinear/bulk_flux.F (state-driven flux computation path)
//   - ROMS/Nonlinear/atm2ocn_flux.F (direct imported flux application path)
//
// Initial architecture uses the same state-passing variable list as the WRF->ROMS path:
//  - 10-m winds: Uwind, Vwind [m s^-1]
//  - Atmospheric pressure: Patm [Pa]
//  - Relative humidity: RH [0-1 fraction]
//  - Air temperature: Tair [K]
//  - Cloud fraction: Cloud [0-1 fraction]
//  - Precipitation rate: Rain [kg m^-2 s^-1]  (project SI convention)
//    TODO(initial-step-testing): verify/convert if upstream forcing uses legacy mm day^-1.
//  - Shortwave radiation: SWrad [W m^-2]
//  - Longwave radiation: LWrad [W m^-2]
//
// Ocean->atmosphere state:
//  - Sea-surface temperature: SST [K]
//
// Coupling mode note:
// - This enum is for legacy state-passing, where REMORA computes bulk fluxes
//   internally (COARE-style path).
// - TODO(future-flux-passing): add a runtime toggle (e.g. ParmParse `driver.atm2ocn_mode`
//   or compile-time ATM2OCN_FLUXES) for direct flux-passing (tau_x, tau_y, heat/
//   moisture fluxes) to reduce cross-component flux inconsistencies seen when
//   each component computes interface fluxes independently.
enum class AtmosToOceanContract : int {
    StatePassing = 0,
    FluxPassing
};

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
    NumAtmosToOceanStateLayers,

    // Extra layers reserved for direct atmosphere->ocean flux exchange.
    // These are appended so legacy state indices remain unchanged.
    TauX = NumAtmosToOceanStateLayers,
    TauY,
    SHflux,
    LHflux,
    Evap,
    NumAtmosToOceanLayers,

    // Backward-compatible alias used by initial state-passing paths.
    NumAtmosToOceanStates = NumAtmosToOceanStateLayers
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

inline const amrex::Vector<std::string>& AtmosToOceanAllNames ()
{
    static const amrex::Vector<std::string> names {
        "Uwind", "Vwind", "Patm", "RH", "Tair", "cloud", "rain", "SWrad", "LWrad",
        "tau_x", "tau_y", "SHflux", "LHflux", "evap"
    };
    return names;
}

inline const amrex::Vector<int>& AtmosToOceanActiveIndices (AtmosToOceanContract contract)
{
    static const amrex::Vector<int> state_indices {
        Uwind, Vwind, Patm, RH, Tair, Cloud, Rain, SWrad, LWrad
    };
    static const amrex::Vector<int> flux_indices {
        TauX, TauY, SHflux, LHflux, SWrad, LWrad, Rain, Evap
    };
    return (contract == AtmosToOceanContract::FluxPassing) ? flux_indices : state_indices;
}

inline const amrex::Vector<std::string>& AtmosToOceanNames (AtmosToOceanContract contract)
{
    static const amrex::Vector<std::string> state_names {
        "Uwind", "Vwind", "Patm", "RH", "Tair", "cloud", "rain", "SWrad", "LWrad"
    };
    static const amrex::Vector<std::string> flux_names {
        "tau_x", "tau_y", "SHflux", "LHflux", "SWrad", "LWrad", "rain", "evap"
    };
    return (contract == AtmosToOceanContract::FluxPassing) ? flux_names : state_names;
}

inline int NumAtmosToOceanActiveLayers (AtmosToOceanContract contract)
{
    return static_cast<int>(AtmosToOceanActiveIndices(contract).size());
}

inline const amrex::Vector<std::string>& OceanToAtmosNames ()
{
    static const amrex::Vector<std::string> names {"SST"};
    return names;
}
} // namespace erf_remora::coupling

#endif
