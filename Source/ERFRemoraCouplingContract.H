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
// Condensed variable-match table
// (COAWST paper/code names -> this driver index -> ERF/REMORA symbols):
//
// Atmosphere -> Ocean cache
// (this repo: Vector<MultiFab*> m_atmos_to_ocean[_active])
//
//   idx  layer   COAWST names   ERF source/derive path    REMORA symbols (sink)                  units
//   ---  ------  -------------  ------------------------  -------------------------------------  ------------------------
//    0   Uwind   Uwind, U10     ERF Vars::xvel/yvel,      vec_uwind; Uwind_data_from_file       m s^-1
//                                eastward at cell centres
//    1   Vwind   Vwind, V10     ERF Vars::xvel/yvel,      vec_vwind; Vwind_data_from_file       m s^-1
//                                northward at cell centres
//    2   Patm    Pair, MSLP     derived "pressure" field  vec_Pair; solverChoice.Pair;          Pa (driver SI), mb in
//                                (base via p0_comp)       Pair_data_from_file                    some REMORA paths
//    3   RH      RH, Hair,      derived rel_humidity      solverChoice.Hair (RH) and/or         0-1 (driver)
//                RELH           from moisture state        vec_qair (specific humidity path)
//    4   Tair    Tair, T2       derived near-sfc temp     vec_Tair; solverChoice.Tair;          K (driver), C in
//                                (RhoTheta/moisture)      Tair_data_from_file                    REMORA internals/IO
//    5   Cloud   cloud, CLDFRA  ERF::cloud_fraction       vec_cloud; solverChoice.cloud;         0-1
//                                                          cloud_data_from_file
//    6   Rain    rain, RAIN     d(surface precip accum)/dt vec_rain; solverChoice.rain;         kg m^-2 s^-1
//                                from the microphysics     rain_data_from_file
//    7   SWrad   swrad, GSW     ERF rad/lsm fields        vec_srflx                              W m^-2
//    8   LWrad   lwrad, GLW     ERF rad/lsm fields        vec_longwave_down (state path)        W m^-2
//                                                          and/or vec_lrflx (flux field)
//    9   TauX    sustr, tau_x   eastward surface stress   vec_sustr (REMORA rotates and          N m^-2
//                                at cell centres           averages to u-faces)
//   10   TauY    svstr, tau_y   northward surface stress  vec_svstr (REMORA rotates and          N m^-2
//                                at cell centres           averages to v-faces)
//   11   SHflux  shflx, HFX     sensible-heat candidate   vec_shflx                              W m^-2
//   12   LHflux  lhflx, LH      latent-heat candidate     vec_lhflx                              W m^-2
//   13   Evap    evap, EVAP     evap export candidate     vec_evap; vec_EminusP                 kg m^-2 s^-1 (or m s^-1)
//
// Ocean -> Atmosphere cache
// (this repo: Vector<MultiFab*> m_ocean_to_atmos)
//
//   idx  layer   COAWST name    REMORA source symbols     ERF symbols (sink)                    units
//   ---  ------  -------------  ------------------------  -------------------------------------  ------------------------
//    0   SST     SST            surface tracer at         sst_lev and/or                         K in this driver
//                                Temp_comp (top cell);     lsm.Get_Data_Ptr(lev,0)
//                                no dedicated vec_sst
//
// Location and frame:
// - Every atmosphere->ocean lane crosses at cell centres (REMORA rho points).
// - Vectors (Uwind/Vwind, TauX/TauY) cross as eastward/northward components.
//   Each model turns them between its own grid axes and east/north: ERF from
//   its lat/lon, REMORA with its grid `angle`. The driver remaps every lane
//   with the same cell -> rho weights and rotates nothing.
//
// Unit notes:
// - RH may arrive as percent in some legacy forcing paths; normalize to 0-1 in the driver contract.
// - REMORA temperature state is often represented internally in Celsius.
// - Non-synthetic coupling must apply explicit conversions at pack/apply boundaries:
//   SST C->K for ERF-facing slabs; Tair/RH/qair conversions per active path.
//   Expected sites: REMORA::ApplyAtmosphericStates and ERF::ApplyOceanSurfaceState.
// - Flux mode requires explicit sign checks (upward-positive atmospheric vs
//   downward-positive ocean conventions) before writing SHflux/LHflux/rain/evap slabs.
//
// Why this table is asymmetric (9+ atmosphere->ocean vs 1 ocean->atmos):
// - This follows legacy COAWST Block-B state passing (Warner et al., 2010),
//   where atmosphere exports many near-surface states and ocean computes
//   stresses/heat internally (BULK_FLUXES-style workflow).
// - In that baseline contract, ocean feedback to atmosphere is minimal and
//   primarily SST; richer feedbacks are added in later coupling modes.
// - The modern COAWST direct-flux path (ATM2OCN_FLUXES family) shifts more
//   exchange to explicit flux layers and can justify expanding OCN->ATM too.
//
// Two independent extension axes (to avoid design confusion):
// - Axis 1: ATM->OCN contract mode (`driver.atm2ocn_mode = state|flux`).
//   This controls whether REMORA is driven by atmospheric states or direct
//   atmosphere-computed flux layers.
// - Axis 2: additional coupled physics (for example wave mediation).
//   Wave coupling can add variables in either direction (ATM<->WAV, OCN<->WAV,
//   and indirect ATM<->OCN impacts) and is not equivalent to choosing
//   `state` vs `flux` for ATM->OCN.
//
// Initial architecture uses the same state-passing variable list as the WRF->ROMS path:
//  - 10-m winds: Uwind, Vwind [m s^-1]
//  - Atmospheric pressure: Patm [Pa]
//  - Relative humidity: RH [0-1 fraction]
//  - Air temperature: Tair [K]
//  - Cloud fraction: Cloud [0-1 fraction]
//  - Precipitation rate: Rain [kg m^-2 s^-1], all phases as liquid water, >= 0.
//    ERF keeps only cumulative surface accumulations [kg m^-2 == mm], so the lane
//    is their change over the exchange interval divided by its length, the interval
//    mean (see IsAccumulationRateLane). 1 kg m^-2 s^-1 = 3600 mm h^-1.
//  - Shortwave radiation: SWrad [W m^-2]
//  - Longwave radiation: LWrad [W m^-2]
//
// Ocean->atmosphere state:
//  - Sea-surface temperature: SST [K]
//
// Coupling mode note:
// - Runtime toggle: ParmParse `driver.atm2ocn_mode = state|flux`.
// - `state`: legacy state-passing where REMORA computes bulk fluxes internally
//   (COARE-style path).
//   Active ATM->OCN layers: Uwind, Vwind, Patm, RH, Tair, cloud, rain, SWrad, LWrad.
// - `flux`: direct stress/heat/moisture exchange layers selected from the same
//   native AMReX Vector<MultiFab*> allocation.
//   Active ATM->OCN layers: tau_x, tau_y, SHflux, LHflux, SWrad, LWrad, rain, evap.
enum class AtmosToOceanContract : int {
    StatePassing = 0,
    FluxPassing
};

enum AtmosToOceanIndex : int {
    Uwind = 0,
    Vwind,
    Patm,
    Qair,
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
    RH = Qair,
    NumAtmosToOceanStates = NumAtmosToOceanStateLayers
};

enum OceanToAtmosIndex : int {
    SST = 0,
    NumOceanToAtmosStates
};

inline const amrex::Vector<std::string>& AtmosToOceanNames ()
{
    static const amrex::Vector<std::string> names {
        "Uwind", "Vwind", "Patm", "Qair", "Tair", "cloud", "rain", "SWrad", "LWrad"
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
        "Uwind", "Vwind", "Patm", "Qair", "Tair", "cloud", "rain", "SWrad", "LWrad"
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

inline const char* AtmosToOceanContractName (AtmosToOceanContract contract)
{
    return (contract == AtmosToOceanContract::FluxPassing) ? "flux" : "state";
}

// Where one atmos->ocean lane's data comes from, resolved once per lane.
//
// ERF by default, since coupling it to an ocean is the point. A lane the deck
// configured explicitly stays with the deck; a lane ERF cannot fill stays with
// REMORA whether the deck asked or not, because REMORA reads any non-null lane
// as authoritative, so handing over prefill is worse than handing over nothing.
enum class LaneSource : int {
    Erf = 0,        // ERF packs it; handed to REMORA
    RemoraDeck,     // deck named it; withheld, deck's own forcing stands
    RemoraDefault   // neither claims it; withheld, REMORA's default stands
};

inline const char* LaneSourceName (LaneSource source)
{
    switch (source) {
    case LaneSource::Erf:           return "ERF";
    case LaneSource::RemoraDeck:    return "REMORA";
    case LaneSource::RemoraDefault: return "REMORA";
    }
    return "?";
}

// One row per lane, replacing the separately-maintained name/index vectors.
// The provenance report and the driver.atm2ocn_lane_source.<name> override are
// keyed on it.
struct AtmosToOceanLaneSpec {
    int         index;              // AtmosToOceanIndex
    const char* name;               // as printed, and as spelled in overrides
    const char* remora_type_key;    // "remora.<key>", or nullptr if REMORA has none
    const char* units;
};

inline const amrex::Vector<AtmosToOceanLaneSpec>& AtmosToOceanLaneSpecs ()
{
    // Ordered by AtmosToOceanIndex, so specs[i].index == i.
    static const amrex::Vector<AtmosToOceanLaneSpec> specs {
        {Uwind,  "Uwind",  "uwind_type", "m s^-1"},
        {Vwind,  "Vwind",  "vwind_type", "m s^-1"},
        {Patm,   "Patm",   "pair_type",  "Pa"},
        // REMORA reads this as specific humidity or an RH fraction depending on
        // magnitude (REMORA_bulk_flux.cpp). Unresolved; one reason ERF's
        // humidity pack is off.
        {Qair,   "Qair",   "qair_type",  "kg kg^-1 or 0-1"},
        {Tair,   "Tair",   "tair_type",  "K"},
        {Cloud,  "cloud",  "cloud_type", "0-1"},
        {Rain,   "rain",   "rain_type",  "kg m^-2 s^-1"},
        {SWrad,  "SWrad",  "swrad_type", "W m^-2"},
        {LWrad,  "LWrad",  "lwrad_type", "W m^-2"},
        {TauX,   "tau_x",  nullptr,      "N m^-2"},
        {TauY,   "tau_y",  nullptr,      "N m^-2"},
        {SHflux, "SHflux", nullptr,      "W m^-2"},
        {LHflux, "LHflux", nullptr,      "W m^-2"},
        {Evap,   "evap",   nullptr,      "kg m^-2 s^-1"}
    };
    return specs;
}

// A lane ERF fills with the mean rate since its previous pack, differenced from a
// cumulative accumulation, rather than with an instantaneous value. Averaging one
// over an interval weights each sample by the ERF steps it spans, not by the
// trapezoid, which makes the result the interval's accumulation over its length.
inline bool IsAccumulationRateLane (int global_idx)
{
    return global_idx == Rain;
}

inline const AtmosToOceanLaneSpec& AtmosToOceanLane (int global_idx)
{
    return AtmosToOceanLaneSpecs()[global_idx];
}

inline const amrex::Vector<std::string>& OceanToAtmosNames ()
{
    static const amrex::Vector<std::string> names {"SST"};
    return names;
}
} // namespace erf_remora::coupling

#endif
