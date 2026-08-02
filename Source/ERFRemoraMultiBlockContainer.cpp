#include "ERFRemoraMultiBlockContainer.H"
#include "StandaloneSolverRuns.H"

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_String.H>

#include <cmath>
#include <string>

namespace erf_remora
{
namespace
{
void
DriverPrintHeader (int driver_verbose,
                   const char* title,
                   const char* detail = nullptr)
{
    if (driver_verbose < 2 || !amrex::ParallelDescriptor::IOProcessor()) { return; }
    amrex::Print() << "========================================================\n"
                   << "[driver] " << title << "\n";
    if (detail != nullptr && detail[0] != '\0') {
        amrex::Print() << detail << "\n";
    }
    amrex::Print() << "========================================================\n";
}

void
DriverPrintLine (int driver_verbose, const std::string& line)
{
    if (driver_verbose < 2 || !amrex::ParallelDescriptor::IOProcessor()) { return; }
    amrex::Print() << "[driver] " << line << "\n";
}

void
ValidateFluxActiveViewOrAbort (const amrex::Vector<amrex::MultiFab*>& active_view,
                               const amrex::Vector<std::string>& active_names)
{
    static const amrex::Vector<std::string> expected_names {
        "tau_x", "tau_y", "SHflux", "LHflux", "SWrad", "LWrad", "rain", "evap"
    };

    if (active_view.size() != expected_names.size()) {
        amrex::Abort("Flux contract active-view size mismatch: expected 8 active lanes.");
    }
    if (active_names.size() != expected_names.size()) {
        amrex::Abort("Flux contract active-name count mismatch: expected 8 names.");
    }

    for (int idx = 0; idx < static_cast<int>(expected_names.size()); ++idx) {
        if (active_view[idx] == nullptr) {
            amrex::Abort("Flux contract active-view contains a null lane pointer.");
        }
        if (active_names[idx] != expected_names[idx]) {
            amrex::Abort("Flux contract active-name order mismatch.");
        }
        const amrex::IntVect ixtype = active_view[idx]->boxArray().ixType().toIntVect();
        if (idx == 0 && ixtype != amrex::IntVect(1,0,0)) {
            amrex::Abort("Flux contract tau_x lane must be x-face centered.");
        }
        if (idx == 1 && ixtype != amrex::IntVect(0,1,0)) {
            amrex::Abort("Flux contract tau_y lane must be y-face centered.");
        }
        if (idx >= 2 && ixtype != amrex::IntVect::TheCellVector()) {
            amrex::Abort("Flux contract scalar flux lanes must remain cell-centered.");
        }
    }
}

std::string
JoinLaneNames (const amrex::Vector<std::string>& names)
{
    std::string joined;
    for (int idx = 0; idx < static_cast<int>(names.size()); ++idx) {
        if (idx != 0) {
            joined += ", ";
        }
        joined += names[idx];
    }
    return joined;
}

bool
PrintCouplingLayouts ()
{
    static bool initialized = false;
    static bool enabled = false;
    if (!initialized) {
        int flag = 0;
        amrex::ParmParse("driver").query("print_coupling_layouts", flag);
        enabled = (flag != 0);
        initialized = true;
    }
    return enabled;
}

void
PrintDistributionMap (const amrex::DistributionMapping& dm)
{
    const auto& pmap = dm.ProcessorMap();
    amrex::Print() << " pmap=[";
    for (int i = 0; i < static_cast<int>(pmap.size()); ++i) {
        amrex::Print() << (i == 0 ? "" : ",") << pmap[i];
    }
    amrex::Print() << "]";
}

void
PrintBoxArrayLayout (const char* label,
                     const amrex::BoxArray& ba,
                     const amrex::DistributionMapping& dm)
{
    if (!PrintCouplingLayouts()) { return; }
    amrex::Print() << "[coupling-layout] " << label
                   << " boxes=" << ba.size()
                   << " domain=" << ba.minimalBox();
    PrintDistributionMap(dm);
    amrex::Print() << "\n" << ba << "\n";
}

void
PrintMultiFabLayout (const char* label, const amrex::MultiFab& mf)
{
    if (!PrintCouplingLayouts()) { return; }
    amrex::Print() << "[coupling-layout] " << label
                   << " ncomp=" << mf.nComp()
                   << " ngrow=" << mf.nGrowVect();
    PrintDistributionMap(mf.DistributionMap());
    for (int comp = 0; comp < mf.nComp(); ++comp) {
        amrex::Print() << " c" << comp
                       << "_min=" << mf.min(comp)
                       << " c" << comp
                       << "_max=" << mf.max(comp);
    }
    amrex::Print() << "\n";
    PrintBoxArrayLayout(label, mf.boxArray(), mf.DistributionMap());
}

void
ValidateCouplingLayout (const amrex::BoxArray& fine_ba,
                        const amrex::DistributionMapping& fine_dm,
                        const amrex::BoxArray& coarse_ba,
                        const amrex::DistributionMapping& coarse_dm,
                        const amrex::IntVect& ratio,
                        const char* label,
                        bool domains_match = true)
{
    const amrex::IntVect fine_len = fine_ba.minimalBox().length();
    const amrex::IntVect coarse_len = coarse_ba.minimalBox().length();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "[coupling-layout] " << label
                       << " fine_boxes=" << fine_ba.size()
                       << " coarse_boxes=" << coarse_ba.size()
                       << " ratio=" << ratio
                       << " validation=checking\n";
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        fine_ba.size() == fine_dm.size(),
        (std::string(label) + ": fine BA/DM size mismatch.").c_str());
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        coarse_ba.size() == coarse_dm.size(),
        (std::string(label) + ": coarse BA/DM size mismatch.").c_str());
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        fine_ba.ixType() == coarse_ba.ixType(),
        (std::string(label) + ": fine/coarse BoxArray index types differ.").c_str());
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        fine_ba.minimalBox().smallEnd(0) == coarse_ba.minimalBox().smallEnd(0) &&
        fine_ba.minimalBox().smallEnd(1) == coarse_ba.minimalBox().smallEnd(1),
        (std::string(label) + ": fine/coarse slab origins are not aligned.").c_str());
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        fine_len[2] == 1 && coarse_len[2] == 1,
        "ValidateCouplingLayout: coupling slabs must remain one-cell thick in z.");

    if (domains_match) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            fine_len[0] == coarse_len[0] * ratio[0] &&
            fine_len[1] == coarse_len[1] * ratio[1],
            "ValidateCouplingLayout: for matching domains, fine/coarse slab lengths must match the requested ratio.");
    }

    amrex::BoxArray coarsened_fine = fine_ba;
    for (int i = 0; i < coarsened_fine.size(); ++i) {
        const amrex::Box original = coarsened_fine[i];
        amrex::Box coarsened = original;
        coarsened.coarsen(ratio);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            coarsened.refine(ratio) == original,
            (std::string(label) +
             ": fine box not evenly coarsenable by ratio at index " +
             std::to_string(i) +
             ". Direct average_down across this coupling layout is unsafe.").c_str());
    }
    coarsened_fine.coarsen(ratio);

    if (coarsened_fine == coarse_ba) {
        amrex::Print() << "[coupling-layout] " << label
                       << " ratio=" << ratio
                       << " direct_average_down_safe=yes\n";
    } else {
        amrex::Print() << "[coupling-layout] " << label
                       << " ratio=" << ratio
                       << " direct_average_down_safe=no"
                       << " warning=coarsen(fine_ba,ratio)!=coarse_ba"
                       << " action=pack_helper_will_average_down_to_tmp_then_parallelcopy\n";
    }
}

amrex::BoxArray
AtmosToOceanBoxArray (const amrex::BoxArray& ba_ocn, int idx)
{
    if (idx == coupling::TauX || idx == coupling::Uwind) {
        return amrex::convert(ba_ocn, amrex::IntVect(1,0,0));
    }
    if (idx == coupling::TauY || idx == coupling::Vwind) {
        return amrex::convert(ba_ocn, amrex::IntVect(0,1,0));
    }
    return ba_ocn;
}

amrex::IntVect
AtmosToOceanNGrow (int idx)
{
    // Flux ingestion in REMORA reads grown tiles around the valid region, so the
    // direct-flux lanes need ghost cells that match REMORA's current 2D stencil width.
    static constexpr amrex::IntVect flux_ngrow(2, 2, 0);
    if (idx >= coupling::TauX) {
        return flux_ngrow;
    }
    return amrex::IntVect(0, 0, 0);
}

void
CopyFieldToCellCenteredPlot (amrex::MultiFab& plot_mf, int comp, const amrex::MultiFab& field)
{
    using namespace amrex;

    const IntVect ixtype = field.boxArray().ixType().toIntVect();
    if (ixtype == IntVect::TheCellVector()) {
        MultiFab::Copy(plot_mf, field, 0, comp, 1, 0);
        return;
    }

    for (MFIter mfi(plot_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box bx = mfi.validbox();
        bx.makeSlab(2,0);
        auto const& src = field.const_array(mfi);
        auto dst = plot_mf.array(mfi, comp);

        if (ixtype == IntVect(1,0,0)) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                dst(i,j,k) = amrex::Real(0.5) * (src(i,j,k) + src(i+1,j,k));
            });
        } else if (ixtype == IntVect(0,1,0)) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                dst(i,j,k) = amrex::Real(0.5) * (src(i,j,k) + src(i,j+1,k));
            });
        } else {
            amrex::Abort("Unsupported non-cell-centered field in exchange plotfile writer.");
        }
    }
}
}

ERFRemoraMultiBlockContainer::~ERFRemoraMultiBlockContainer()
{
    if (m_erf != nullptr) {
        DestroyERFSolver(m_erf);
        m_erf = nullptr;
    }
    if (m_remora != nullptr) {
        DestroyREMORASolver(m_remora);
        m_remora = nullptr;
    }
}

ERFRemoraMultiBlockContainer::ERFRemoraMultiBlockContainer (
    const amrex::RealBox& atm_rb,
    const amrex::Vector<int>& n_cell_atm,
    const amrex::RealBox& ocn_rb,
    const amrex::Vector<int>& n_cell_ocn,
    int max_level,
    int coord,
    const amrex::Vector<amrex::IntVect>& ref_ratio,
    const amrex::Array<int,AMREX_SPACEDIM>& is_periodic,
    const amrex::IntVect& max_grid_size,
    amrex::Real erf_dt,
    int atm_substeps_per_ocean_step,
    coupling::AtmosToOceanContract atm2ocn_contract,
    bool use_two_way_coupling)
    : m_atm_rb(atm_rb), m_atm_n_cell(n_cell_atm),
      m_ocn_rb(ocn_rb), m_ocn_n_cell(n_cell_ocn),
      m_max_level(max_level), m_coord(coord),
      m_ref_ratio(ref_ratio), m_is_periodic(is_periodic),
      m_erf_dt(erf_dt),
      m_atm_substeps_per_ocean_step(atm_substeps_per_ocean_step),
      m_atmos_to_ocean_contract(atm2ocn_contract),
      m_use_two_way_coupling(use_two_way_coupling)
{
    // Validate containment: atmosphere must fully cover ocean
    const bool atmosphere_contains_ocean =
        (atm_rb.lo(0) <= ocn_rb.lo(0)) && (atm_rb.lo(1) <= ocn_rb.lo(1)) &&
        (atm_rb.hi(0) >= ocn_rb.hi(0)) && (atm_rb.hi(1) >= ocn_rb.hi(1));

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        atmosphere_contains_ocean,
        "Atmosphere domain must fully contain ocean domain.\n");

    // Determine if domains match horizontally
    m_domains_match =
        (atm_rb.lo(0) == ocn_rb.lo(0)) && (atm_rb.lo(1) == ocn_rb.lo(1)) &&
        (atm_rb.hi(0) == ocn_rb.hi(0)) && (atm_rb.hi(1) == ocn_rb.hi(1));

    // Compute horizontal resolution ratio
    amrex::IntVect horiz_ratio;
    if (m_domains_match) {
        // Matching domains: cell counts must have integer refinement ratio
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            n_cell_atm[0] % n_cell_ocn[0] == 0 &&
            n_cell_atm[1] % n_cell_ocn[1] == 0,
            "Matching domains: atmosphere cell count must be integer multiple of ocean cell count.");
        horiz_ratio = amrex::IntVect(n_cell_atm[0] / n_cell_ocn[0],
                                     n_cell_atm[1] / n_cell_ocn[1], 1);

        amrex::Print() << "Coupling configuration: MATCHING domains\n"
                       << "  Horizontal resolution ratio: " << horiz_ratio[0] << " x " << horiz_ratio[1] << "\n";
    } else {
        // Nonmatching domains: compute ratio from grid spacing
        // TODO: validate integer grid spacing ratio for conservative remapping
        horiz_ratio = amrex::IntVect(1, 1, 1);  // Placeholder for now

        amrex::Print() << "Coupling configuration: NONMATCHING domains\n"
                       << "  Coupling interface = ocean top surface\n";
    }

    // Coupling interface is the entire ocean top surface
    m_interface_rb = ocn_rb;
    m_interface_ncell = amrex::IntVect(n_cell_ocn[0], n_cell_ocn[1], n_cell_ocn[2]);

    // Define atmosphere domain Box
    amrex::Box domain(amrex::IntVect(0, 0, 0),
                      amrex::IntVect(n_cell_atm[0]-1, n_cell_atm[1]-1, 0));

    // Interface slab construction
    if (m_domains_match) {
        // Matching domains: use atmosphere extent, coarsen to ocean resolution
        m_ba.define(domain);
        m_ba.maxSize(max_grid_size);
        m_dm.define(m_ba);
        PrintBoxArrayLayout("Driver constructor ocean->atmos slab", m_ba, m_dm);

        // Derive ocean slab from atmosphere slab via coarsening
        m_ba_ocn = m_ba;
        m_ba_ocn.coarsen(horiz_ratio);
        m_dm_ocn.define(m_ba_ocn);
        PrintBoxArrayLayout("Driver constructor atmos->ocean slab", m_ba_ocn, m_dm_ocn);
        ValidateCouplingLayout(m_ba, m_dm, m_ba_ocn, m_dm_ocn,
                               horiz_ratio,
                               "ERF->REMORA horizontal coupling",
                               m_domains_match);
    } else {
        // Nonmatching domains: construct separate slabs for atmosphere and ocean extents
        // m_ba = atmosphere extent (for ocean->atmosphere cache)
        m_ba.define(domain);
        m_ba.maxSize(max_grid_size);
        m_dm.define(m_ba);
        PrintBoxArrayLayout("Driver constructor ocean->atmos slab (atm extent)", m_ba, m_dm);

        // m_ba_ocn = ocean extent (for atmosphere->ocean cache)
        amrex::Box ocn_domain(amrex::IntVect(0, 0, 0),
                              amrex::IntVect(n_cell_ocn[0]-1, n_cell_ocn[1]-1, 0));
        m_ba_ocn.define(ocn_domain);
        m_ba_ocn.maxSize(max_grid_size);
        m_dm_ocn.define(m_ba_ocn);
        PrintBoxArrayLayout("Driver constructor atmos->ocean slab (ocn extent)", m_ba_ocn, m_dm_ocn);

        // Note: No ValidateCouplingLayout for nonmatching - grids have different extents
    }

    amrex::Real dz_atm = (atm_rb.hi(2) - atm_rb.lo(2)) / n_cell_atm[2];
    amrex::RealBox slab_rb({atm_rb.lo(0), atm_rb.lo(1), 0.0},
                            {atm_rb.hi(0), atm_rb.hi(1), dz_atm});
    m_geom.define(domain, &slab_rb, amrex::CoordSys::cartesian, is_periodic.data());

    amrex::Real dz_ocn = (ocn_rb.hi(2) - ocn_rb.lo(2)) / n_cell_ocn[2];
    amrex::RealBox ocn_slab_rb({ocn_rb.lo(0), ocn_rb.lo(1), 0.0},
                                {ocn_rb.hi(0), ocn_rb.hi(1), dz_ocn});
    m_geom_ocn.define(m_ba_ocn.minimalBox(), &ocn_slab_rb, amrex::CoordSys::cartesian, is_periodic.data());

    m_atmos_to_ocean_state_active_indices =
        coupling::AtmosToOceanActiveIndices(coupling::AtmosToOceanContract::StatePassing);
    m_atmos_to_ocean_state_active_names =
        coupling::AtmosToOceanNames(coupling::AtmosToOceanContract::StatePassing);
    m_atmos_to_ocean_flux_active_indices =
        coupling::AtmosToOceanActiveIndices(coupling::AtmosToOceanContract::FluxPassing);
    m_atmos_to_ocean_flux_active_names =
        coupling::AtmosToOceanNames(coupling::AtmosToOceanContract::FluxPassing);

    amrex::ParmParse pp_driver("driver");
    pp_driver.query("v", m_driver_verbose);
    pp_driver.query("plot_int_exchange", m_plot_int_exchange);
    pp_driver.query("plot_file_exchange_o2a", m_plot_file_exchange_o2a);
    pp_driver.query("plot_file_exchange_a2o", m_plot_file_exchange_a2o);

    DriverPrintHeader(
        m_driver_verbose,
        "Construct ERFRemoraMultiBlockContainer",
        "Driver topology parsed; base slab geometry and contract metadata are being recorded.");
    DriverPrintLine(
        m_driver_verbose,
        std::string("configured_atm2ocn_mode=") +
            coupling::AtmosToOceanContractName(m_atmos_to_ocean_contract) +
            " two_way=" + (m_use_two_way_coupling ? "1" : "0"));
}

void
ERFRemoraMultiBlockContainer::InitializeSolvers ()
{
    DriverPrintHeader(
        m_driver_verbose,
        "Initialize Solvers",
        "Phase 1: construct ERF and REMORA objects. Phase 2: inject driver atm2ocn metadata before InitData. Phase 3: call InitData on both solvers.");

    if (m_erf != nullptr) {
        DriverPrintLine(m_driver_verbose, "destroying previous ERF solver handle");
        DestroyERFSolver(m_erf);
        m_erf = nullptr;
    }
    if (m_remora != nullptr) {
        DriverPrintLine(m_driver_verbose, "destroying previous REMORA solver handle");
        DestroyREMORASolver(m_remora);
        m_remora = nullptr;
    }

    DriverPrintLine(m_driver_verbose, "constructing ERF solver from driver-owned atmosphere topology");
    m_erf = CreateERFSolver(*this);
    DriverPrintLine(m_driver_verbose, "ERF constructor complete");
    DriverPrintLine(m_driver_verbose, "constructing REMORA solver from driver-owned ocean topology");
    m_remora = CreateREMORASolver(*this);
    DriverPrintLine(m_driver_verbose, "REMORA constructor complete");

    const bool bootstrap_uses_state_contract =
        (ActiveAtmosToOceanContract(/*ocean_step=*/0) == coupling::AtmosToOceanContract::StatePassing);
    DriverPrintLine(
        m_driver_verbose,
        std::string("injecting pre-InitData atm2ocn metadata: bootstrap_contract=") +
            (bootstrap_uses_state_contract ? "state" : "flux") +
            " configured_mode=" + coupling::AtmosToOceanContractName(m_atmos_to_ocean_contract) +
            " two_way=" + (m_use_two_way_coupling ? "1" : "0"));
    ConfigureERFDriverAtmosToOceanCoupling(
        m_erf, true, m_use_two_way_coupling, bootstrap_uses_state_contract);
    ConfigureREMORADriverAtmosToOceanCoupling(
        m_remora, true, m_use_two_way_coupling, bootstrap_uses_state_contract);

    DriverPrintLine(m_driver_verbose, "calling ERF InitData (enters ERF InitData_pre -> ERF InitData_post)");
    InitERFSolver(m_erf);
    DriverPrintLine(m_driver_verbose, "ERF InitData complete");
    DriverPrintLine(m_driver_verbose, "calling REMORA InitData");
    InitREMORASolver(m_remora);
    DriverPrintLine(m_driver_verbose, "REMORA InitData complete");
}

void
ERFRemoraMultiBlockContainer::PrefillAtmosToOceanFallbacks ()
{
    // Set default values in the atmos→ocean slab for fields that ERF may not
    // be able to provide (no moisture model, no radiation). Called by the driver
    // before PackAtmosphericStates so that PackAtmosphericStates can simply leave
    // unavailable slots untouched rather than guessing REMORA-specific constants.
    //
    // Under the current REMORA bulk-flux interpretation, values below 2 in the
    // humidity lane are treated as relative humidity [0-1]. Seed the fallback
    // from remora.air_humidity, then let ERF overwrite it with diagnosed RH
    // when atmospheric moisture is active.
    amrex::ParmParse pp_rem("remora");
    amrex::Real humidity_val = 0.0, cloud_val = 0.0, rain_val = 0.0;
    pp_rem.query("air_humidity", humidity_val);
    pp_rem.query("cloud", cloud_val);
    pp_rem.query("rain",  rain_val);

    for (int global_idx = 0; global_idx < static_cast<int>(m_atmos_to_ocean.size()); ++global_idx) {
        auto* mf = m_atmos_to_ocean[global_idx];
        if (mf == nullptr) {
            continue;
        }

        switch (global_idx) {
        case coupling::Qair:
            mf->setVal(humidity_val); // RH fallback [0-1] for REMORA bulk fluxes
            break;
        case coupling::Cloud:
            mf->setVal(cloud_val); // cloud fraction [0-1]
            break;
        case coupling::Rain:
            mf->setVal(rain_val);  // rain rate [kg m-2 s-1]
            break;
        case coupling::TauX:
            mf->setVal(0.0);       // surface stress [Pa]: reset each macro-step
            break;
        case coupling::TauY:
            mf->setVal(0.0);       // surface stress [Pa]: reset each macro-step
            break;
        case coupling::SHflux:
            mf->setVal(0.0);       // sensible heat flux [W m-2]
            break;
        case coupling::LHflux:
            mf->setVal(0.0);       // latent heat flux [W m-2]
            break;
        case coupling::SWrad:
            mf->setVal(0.0);       // shortwave [W m-2]: no radiation fallback
            break;
        case coupling::LWrad:
            mf->setVal(0.0);       // longwave  [W m-2]: no radiation fallback
            break;
        case coupling::Evap:
            mf->setVal(0.0);       // evaporation rate [kg m-2 s-1]
            break;
        default:
            break;
        }
    }
}

void
ERFRemoraMultiBlockContainer::SetAtmosphereTick (long long atm_tick)
{
    if (atm_tick < 0) {
        amrex::Abort("Cannot restore a negative atmosphere tick in the coupled driver.");
    }
    m_atm_tick = atm_tick;
}

void
ERFRemoraMultiBlockContainer::InitializeBlocks ()
{
    DriverPrintHeader(
        m_driver_verbose,
        "Initialize Coupling Caches",
        "Phase 4: query finalized post-InitData solver layouts. Phase 5: allocate driver-owned caches from solver layouts. Phase 6: seed bootstrap state and select active atm2ocn view.");
    m_atm_tick = 0;
    DriverPrintLine(m_driver_verbose, "resetting driver atmosphere tick to 0");
    AllocateCaches();
    DriverPrintLine(m_driver_verbose, "selecting bootstrap active atm2ocn cache view for ocean_step=0");
    SelectAtmosToOceanActiveView(/*ocean_step=*/0);
    if (m_erf != nullptr && m_remora != nullptr) {
        if (m_use_two_way_coupling && !m_ocean_to_atmos.empty()) {
            DriverPrintLine(m_driver_verbose, "packing REMORA surface state into ocean->atmos bootstrap cache");
            REMORAPackSurfaceState(m_remora, m_ocean_to_atmos, amrex::Real(0.0));
            FillOceanToAtmosGhosts();
        }
    } else {
        DriverPrintLine(m_driver_verbose, "synthetic-only path: filling synthetic coupling cache values");
        FillSyntheticData();
    }

    if (m_atmos_to_ocean_contract == coupling::AtmosToOceanContract::FluxPassing) {
        DriverPrintLine(m_driver_verbose, "validating steady-state flux active view shape and centering");
        ValidateFluxActiveViewOrAbort(m_atmos_to_ocean_flux_active, m_atmos_to_ocean_flux_active_names);
    }

    amrex::Print() << "Atmos->ocean contract mode: "
                       << coupling::AtmosToOceanContractName(m_atmos_to_ocean_contract)
                       << " (" << coupling::NumAtmosToOceanActiveLayers(m_atmos_to_ocean_contract)
                       << " active layers)\n";
    if (m_atmos_to_ocean_contract == coupling::AtmosToOceanContract::FluxPassing) {
        amrex::Print() << "Atmos->ocean bootstrap lanes (ocean_step=0): "
                       << JoinLaneNames(m_atmos_to_ocean_state_active_names) << "\n"
                       << "Atmos->ocean steady-state lanes (ocean_step>=1): "
                       << JoinLaneNames(m_atmos_to_ocean_flux_active_names) << "\n";
    } else {
        amrex::Print() << "Atmos->ocean active lanes: "
                       << JoinLaneNames(m_atmos_to_ocean_active_names) << "\n";
    }

    // TODO(initial-step-testing): implement a soft-start / ramp-in protocol for coupled forcing.
    // COAWST-style coupled systems can exhibit an initialization shock at t=0 when
    // atmosphere/ocean timescales are strongly separated. We should gradually ramp
    // exchanged tendencies across the first few ocean macro-steps.
}

void
ERFRemoraMultiBlockContainer::AllocateCaches ()
{
    if (m_erf != nullptr && m_remora != nullptr) {
        DriverPrintLine(m_driver_verbose, "allocating caches from solver-reported post-InitData layouts");
        AllocateCachesFromSolverLayouts();
    } else {
        DriverPrintLine(m_driver_verbose, "allocating synthetic fallback caches from driver-derived layouts");
        AllocateSyntheticCaches();
    }
}

void
ERFRemoraMultiBlockContainer::AllocateSyntheticCaches ()
{
    // Allocate a superset layer vector once (legacy state + optional flux layers).
    // Runtime coupling modes can then swap active layer views by index without
    // changing memory ownership or reallocating cache slabs.
    m_atmos_to_ocean_storage.resize(coupling::NumAtmosToOceanLayers);
    for (int idx = 0; idx < static_cast<int>(m_atmos_to_ocean_storage.size()); ++idx) {
        auto& mf_ptr = m_atmos_to_ocean_storage[idx];
        mf_ptr = std::make_unique<amrex::MultiFab>(
            AtmosToOceanBoxArray(m_ba_ocn, idx), m_dm_ocn, 1, AtmosToOceanNGrow(idx));
        mf_ptr->setVal(0.0);
        const std::string label = "Driver AllocateCaches atmos->ocean[" + std::to_string(idx) + "]";
        PrintMultiFabLayout(label.c_str(), *mf_ptr);
    }
    m_atmos_to_ocean = amrex::GetVecOfPtrs(m_atmos_to_ocean_storage);
    m_atmos_to_ocean_state_active.resize(m_atmos_to_ocean_state_active_indices.size());
    for (int local_idx = 0; local_idx < static_cast<int>(m_atmos_to_ocean_state_active_indices.size()); ++local_idx) {
        const int global_idx = m_atmos_to_ocean_state_active_indices[local_idx];
        m_atmos_to_ocean_state_active[local_idx] = m_atmos_to_ocean[global_idx];
    }
    m_atmos_to_ocean_flux_active.resize(m_atmos_to_ocean_flux_active_indices.size());
    for (int local_idx = 0; local_idx < static_cast<int>(m_atmos_to_ocean_flux_active_indices.size()); ++local_idx) {
        const int global_idx = m_atmos_to_ocean_flux_active_indices[local_idx];
        m_atmos_to_ocean_flux_active[local_idx] = m_atmos_to_ocean[global_idx];
    }

    if (m_use_two_way_coupling) {
        m_ocean_to_atmos_storage.resize(coupling::NumOceanToAtmosStates);
        const amrex::IntVect ocean_to_atmos_ng(1, 1, 0);
        for (int idx = 0; idx < static_cast<int>(m_ocean_to_atmos_storage.size()); ++idx) {
            auto& mf_ptr = m_ocean_to_atmos_storage[idx];
            // For nonmatching domains: use ocean extent for REMORA to pack into
            // For matching domains: m_ba == m_ba_ocn, so this works for both cases
            const amrex::BoxArray& ba_cache = m_domains_match ? m_ba : m_ba_ocn;
            const amrex::DistributionMapping& dm_cache = m_domains_match ? m_dm : m_dm_ocn;
            mf_ptr = std::make_unique<amrex::MultiFab>(ba_cache, dm_cache, 1, ocean_to_atmos_ng);
            mf_ptr->setVal(0.0);
            const std::string label = "Driver AllocateCaches ocean->atmos[" + std::to_string(idx) + "]";
            PrintMultiFabLayout(label.c_str(), *mf_ptr);
        }
    } else {
        m_ocean_to_atmos_storage.clear();
    }
    m_ocean_to_atmos = amrex::GetVecOfPtrs(m_ocean_to_atmos_storage);
}

void
ERFRemoraMultiBlockContainer::AllocateCachesFromSolverLayouts ()
{
    DriverPrintLine(m_driver_verbose, "querying REMORA rho/u/v atm2ocn interface layouts");
    amrex::BoxArray remora_rho_ba;
    amrex::DistributionMapping remora_rho_dm;
    amrex::BoxArray remora_u_ba;
    amrex::DistributionMapping remora_u_dm;
    amrex::BoxArray remora_v_ba;
    amrex::DistributionMapping remora_v_dm;

    REMORAGetAtmosToOceanRhoLayout(m_remora, remora_rho_ba, remora_rho_dm);
    REMORAGetAtmosToOceanUFaceLayout(m_remora, remora_u_ba, remora_u_dm);
    REMORAGetAtmosToOceanVFaceLayout(m_remora, remora_v_ba, remora_v_dm);

    if (m_use_two_way_coupling) {
        DriverPrintLine(m_driver_verbose, "querying ERF ocean-surface layout for ocean->atmos cache allocation");
        ERFGetOceanToAtmosSurfaceLayout(m_erf, m_ba, m_dm);
    }

    m_ba_ocn = remora_rho_ba;
    m_dm_ocn = remora_rho_dm;

    const amrex::Real dz_atm = (m_atm_rb.hi(2) - m_atm_rb.lo(2)) / m_atm_n_cell[2];
    const amrex::RealBox atm_slab_rb({m_atm_rb.lo(0), m_atm_rb.lo(1), 0.0},
                                     {m_atm_rb.hi(0), m_atm_rb.hi(1), dz_atm});
    m_geom.define(m_ba.minimalBox(), &atm_slab_rb, amrex::CoordSys::cartesian, m_is_periodic.data());

    const amrex::Real dz_ocn = (m_ocn_rb.hi(2) - m_ocn_rb.lo(2)) / m_ocn_n_cell[2];
    const amrex::RealBox ocn_slab_rb({m_ocn_rb.lo(0), m_ocn_rb.lo(1), 0.0},
                                     {m_ocn_rb.hi(0), m_ocn_rb.hi(1), dz_ocn});
    m_geom_ocn.define(m_ba_ocn.minimalBox(), &ocn_slab_rb, amrex::CoordSys::cartesian, m_is_periodic.data());

    DriverPrintLine(m_driver_verbose, "validating solver-reported coupling layouts against integer-ratio driver assumptions");
    ValidateCouplingLayoutsFromSolvers();

    m_atmos_to_ocean_storage.resize(coupling::NumAtmosToOceanLayers);
    for (int idx = 0; idx < static_cast<int>(m_atmos_to_ocean_storage.size()); ++idx) {
        auto& mf_ptr = m_atmos_to_ocean_storage[idx];
        const amrex::BoxArray* lane_ba = &remora_rho_ba;
        const amrex::DistributionMapping* lane_dm = &remora_rho_dm;
        if (idx == coupling::TauX) {
            lane_ba = &remora_u_ba;
            lane_dm = &remora_u_dm;
        } else if (idx == coupling::TauY) {
            lane_ba = &remora_v_ba;
            lane_dm = &remora_v_dm;
        }
        mf_ptr = std::make_unique<amrex::MultiFab>(*lane_ba, *lane_dm, 1, AtmosToOceanNGrow(idx));
        mf_ptr->setVal(0.0);
        const std::string label = "Driver AllocateCachesFromSolverLayouts atmos->ocean[" + std::to_string(idx) + "]";
        PrintMultiFabLayout(label.c_str(), *mf_ptr);
    }
    m_atmos_to_ocean = amrex::GetVecOfPtrs(m_atmos_to_ocean_storage);

    m_atmos_to_ocean_state_active.resize(m_atmos_to_ocean_state_active_indices.size());
    for (int local_idx = 0; local_idx < static_cast<int>(m_atmos_to_ocean_state_active_indices.size()); ++local_idx) {
        const int global_idx = m_atmos_to_ocean_state_active_indices[local_idx];
        m_atmos_to_ocean_state_active[local_idx] = m_atmos_to_ocean[global_idx];
    }
    m_atmos_to_ocean_flux_active.resize(m_atmos_to_ocean_flux_active_indices.size());
    for (int local_idx = 0; local_idx < static_cast<int>(m_atmos_to_ocean_flux_active_indices.size()); ++local_idx) {
        const int global_idx = m_atmos_to_ocean_flux_active_indices[local_idx];
        m_atmos_to_ocean_flux_active[local_idx] = m_atmos_to_ocean[global_idx];
    }

    if (m_use_two_way_coupling) {
        m_ocean_to_atmos_storage.resize(coupling::NumOceanToAtmosStates);
        const amrex::IntVect ocean_to_atmos_ng(1, 1, 0);
        for (int idx = 0; idx < static_cast<int>(m_ocean_to_atmos_storage.size()); ++idx) {
            auto& mf_ptr = m_ocean_to_atmos_storage[idx];
            // For nonmatching domains: use ocean extent (m_ba_ocn) for REMORA to pack into
            // For matching domains: m_ba == m_ba_ocn, so this works for both cases
            const amrex::BoxArray& ba_cache = m_domains_match ? m_ba : m_ba_ocn;
            const amrex::DistributionMapping& dm_cache = m_domains_match ? m_dm : m_dm_ocn;
            mf_ptr = std::make_unique<amrex::MultiFab>(ba_cache, dm_cache, 1, ocean_to_atmos_ng);
            mf_ptr->setVal(0.0);
            const std::string label = "Driver AllocateCachesFromSolverLayouts ocean->atmos[" + std::to_string(idx) + "]";
            PrintMultiFabLayout(label.c_str(), *mf_ptr);
        }
    } else {
        m_ocean_to_atmos_storage.clear();
    }
    m_ocean_to_atmos = amrex::GetVecOfPtrs(m_ocean_to_atmos_storage);
}

void
ERFRemoraMultiBlockContainer::ValidateCouplingLayoutsFromSolvers () const
{
    const amrex::IntVect horiz_ratio(m_atm_n_cell[0] / m_ocn_n_cell[0],
                                     m_atm_n_cell[1] / m_ocn_n_cell[1], 1);

    amrex::BoxArray remora_rho_ba;
    amrex::DistributionMapping remora_rho_dm;
    amrex::BoxArray remora_u_ba;
    amrex::DistributionMapping remora_u_dm;
    amrex::BoxArray remora_v_ba;
    amrex::DistributionMapping remora_v_dm;
    REMORAGetAtmosToOceanRhoLayout(m_remora, remora_rho_ba, remora_rho_dm);
    REMORAGetAtmosToOceanUFaceLayout(m_remora, remora_u_ba, remora_u_dm);
    REMORAGetAtmosToOceanVFaceLayout(m_remora, remora_v_ba, remora_v_dm);

    if (m_use_two_way_coupling) {
        amrex::BoxArray erf_surface_ba;
        amrex::DistributionMapping erf_surface_dm;
        ERFGetOceanToAtmosSurfaceLayout(m_erf, erf_surface_ba, erf_surface_dm);
        ValidateCouplingLayout(erf_surface_ba, erf_surface_dm,
                               remora_rho_ba, remora_rho_dm,
                               horiz_ratio,
                               "ERF surface <-> REMORA rho coupling",
                               m_domains_match);

        const amrex::BoxArray erf_xface_ba = amrex::convert(erf_surface_ba, amrex::IntVect(1,0,0));
        const amrex::BoxArray erf_yface_ba = amrex::convert(erf_surface_ba, amrex::IntVect(0,1,0));
        ValidateCouplingLayout(erf_xface_ba, erf_surface_dm,
                               remora_u_ba, remora_u_dm,
                               horiz_ratio,
                               "ERF x-face <-> REMORA u-face coupling",
                               m_domains_match);
        ValidateCouplingLayout(erf_yface_ba, erf_surface_dm,
                               remora_v_ba, remora_v_dm,
                               horiz_ratio,
                               "ERF y-face <-> REMORA v-face coupling",
                               m_domains_match);
    } else {
        ValidateCouplingLayout(m_ba, m_dm,
                               remora_rho_ba, remora_rho_dm,
                               horiz_ratio,
                               "Driver atmosphere slab <-> REMORA rho coupling",
                               m_domains_match);
    }
}

void
ERFRemoraMultiBlockContainer::FillSyntheticData ()
{
    for (int idx = 0; idx < static_cast<int>(m_atmos_to_ocean.size()); ++idx) {
        m_atmos_to_ocean[idx]->setVal(10.0 + static_cast<amrex::Real>(idx));
    }

    if (m_use_two_way_coupling && !m_ocean_to_atmos.empty()) {
        m_ocean_to_atmos[coupling::SST]->setVal(300.0);
        FillOceanToAtmosGhosts();
    }
}

void
ERFRemoraMultiBlockContainer::FillOceanToAtmosGhosts ()
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_geom.isPeriodic(0) == static_cast<bool>(m_is_periodic[0]) &&
        m_geom.isPeriodic(1) == static_cast<bool>(m_is_periodic[1]) &&
        m_geom.isPeriodic(2) == static_cast<bool>(m_is_periodic[2]),
        "Driver surface-slab geometry lost the requested periodic flags.");

    for (auto* mf : m_ocean_to_atmos) {
        if (mf != nullptr) {
            mf->FillBoundary(m_geom.periodicity());
        }
    }
}

void
ERFRemoraMultiBlockContainer::FillAtmosToOceanGhosts ()
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_geom_ocn.isPeriodic(0) == static_cast<bool>(m_is_periodic[0]) &&
        m_geom_ocn.isPeriodic(1) == static_cast<bool>(m_is_periodic[1]) &&
        m_geom_ocn.isPeriodic(2) == static_cast<bool>(m_is_periodic[2]),
        "Driver ocean-slab geometry lost the requested periodic flags.");

    for (auto* mf : m_atmos_to_ocean) {
        if (mf != nullptr && mf->nGrowVect().max() > 0) {
            mf->FillBoundary(m_geom_ocn.periodicity());
        }
    }
}

void
ERFRemoraMultiBlockContainer::InterpOceanToAtmos (
    const amrex::Vector<amrex::MultiFab*>& ocean_state,
    amrex::Vector<amrex::MultiFab*>& atm_state)
{
    BL_PROFILE("Driver::InterpOceanToAtmos");

    // 1. Initialize atmosphere state with fallback SST everywhere
    //    Use ERF's configured surface temperature (erf.most.surf_temp, default 300 K)
    //    This matches one-way coupling behavior
    amrex::Real fallback_sst = 300.0;
    amrex::ParmParse pp_erf("erf");
    pp_erf.query("most.surf_temp", fallback_sst);

    atm_state[coupling::SST]->setVal(fallback_sst);

    // 2. Compute refinement ratio from grid spacings
    const amrex::Real dx_atm = (m_atm_rb.hi(0) - m_atm_rb.lo(0)) / m_atm_n_cell[0];
    const amrex::Real dy_atm = (m_atm_rb.hi(1) - m_atm_rb.lo(1)) / m_atm_n_cell[1];
    const amrex::Real dx_ocn = (m_ocn_rb.hi(0) - m_ocn_rb.lo(0)) / m_ocn_n_cell[0];
    const amrex::Real dy_ocn = (m_ocn_rb.hi(1) - m_ocn_rb.lo(1)) / m_ocn_n_cell[1];

    const amrex::IntVect ratio(
        static_cast<int>(std::round(dx_ocn / dx_atm)),
        static_cast<int>(std::round(dy_ocn / dy_atm)),
        1
    );

    // 3. Create temporary MultiFab at atmosphere resolution on ocean extent
    amrex::BoxArray ba_atm_temp = m_ba_ocn;
    ba_atm_temp.refine(ratio);

    amrex::DistributionMapping dm_atm_temp(ba_atm_temp);
    amrex::MultiFab atm_temp(ba_atm_temp, dm_atm_temp, 1, 0);

    // 4. Interpolate ocean SST to atmosphere resolution using piecewise constant
    //    ParallelCopy will handle MPI communication and simple interpolation
    //    For more sophisticated interpolation, would need InterpFromCoarseLevel with proper BC setup
    atm_temp.ParallelCopy(*ocean_state[coupling::SST], 0, 0, 1, 0, 0);

    // 5. ParallelCopy interpolated SST from ocean extent to full atmosphere extent
    //    ParallelCopy only overwrites overlapping region (ocean extent)
    //    Non-overlapping regions keep fallback SST
    atm_state[coupling::SST]->ParallelCopy(atm_temp, 0, 0, 1, 0, 0);

    // 6. Fill ghost cells for periodicity
    atm_state[coupling::SST]->FillBoundary(m_geom.periodicity());
}

void
ERFRemoraMultiBlockContainer::RemapAtmosToOceanStates(
    const amrex::Vector<amrex::MultiFab*>& atm_states,
    amrex::Vector<amrex::MultiFab*>& ocn_states)
{
    BL_PROFILE("Driver::RemapAtmosToOcean");

    // 1. Compute refinement ratio from grid spacings
    const amrex::Real dx_atm = (m_atm_rb.hi(0) - m_atm_rb.lo(0)) / m_atm_n_cell[0];
    const amrex::Real dy_atm = (m_atm_rb.hi(1) - m_atm_rb.lo(1)) / m_atm_n_cell[1];
    const amrex::Real dx_ocn = (m_ocn_rb.hi(0) - m_ocn_rb.lo(0)) / m_ocn_n_cell[0];
    const amrex::Real dy_ocn = (m_ocn_rb.hi(1) - m_ocn_rb.lo(1)) / m_ocn_n_cell[1];

    const amrex::IntVect ratio(
        static_cast<int>(std::round(dx_ocn / dx_atm)),
        static_cast<int>(std::round(dy_ocn / dy_atm)),
        1
    );

    // 2. Create atmosphere-resolution BoxArray on ocean extent for extraction
    amrex::BoxArray ba_atm_on_ocn = m_ba_ocn;
    ba_atm_on_ocn.refine(ratio);
    amrex::DistributionMapping dm_atm_on_ocn(ba_atm_on_ocn);

    // 3. For each field, extract and remap
    for (int idx = 0; idx < static_cast<int>(ocn_states.size()); ++idx) {
        if (atm_states[idx] == nullptr || ocn_states[idx] == nullptr) {
            continue;
        }

        // Determine field centering (cell vs face)
        const amrex::IntVect ixtype = ocn_states[idx]->boxArray().ixType().toIntVect();
        const bool is_x_face = (ixtype == amrex::IntVect(1, 0, 0));
        const bool is_y_face = (ixtype == amrex::IntVect(0, 1, 0));

        // Create temp BoxArray with appropriate centering
        amrex::BoxArray ba_temp = ba_atm_on_ocn;
        if (is_x_face) {
            ba_temp.convert(amrex::IntVect(1, 0, 0));
        } else if (is_y_face) {
            ba_temp.convert(amrex::IntVect(0, 1, 0));
        }

        // Stack-allocated temp MultiFab - automatically cleaned up at end of iteration
        amrex::MultiFab atm_temp(ba_temp, dm_atm_on_ocn, 1, 0);

        // Extract atmosphere data over ocean extent via ParallelCopy
        atm_temp.ParallelCopy(*atm_states[idx], 0, 0, 1, 0, 0);

        // Conservative averaging down to ocean resolution
        if (is_x_face || is_y_face) {
            // Face-centered field: use average_down_faces
            amrex::average_down_faces(atm_temp, *ocn_states[idx], ratio, 0);
        } else {
            // Cell-centered field: use average_down
            amrex::average_down(atm_temp, *ocn_states[idx], 0, 1, ratio);
        }

        // Fill ghost cells if needed (only flux fields need ghosts)
        if (ocn_states[idx]->nGrow() > 0) {
            ocn_states[idx]->FillBoundary(m_geom_ocn.periodicity());
        }
    }
}

void
ERFRemoraMultiBlockContainer::PrintDiagnostics (const char* label,
                             const amrex::Vector<amrex::MultiFab*>& fields,
                             const amrex::Vector<std::string>& names) const
{
    amrex::Print() << label << "\n";
    for (int idx = 0; idx < static_cast<int>(fields.size()); ++idx) {
        const auto min_val = fields[idx]->min(0);
        const auto max_val = fields[idx]->max(0);
        amrex::Print() << "  " << names[idx] << ": min=" << min_val << " max=" << max_val << "\n";
        if (names[idx] == "rain" && (min_val > 100.0 || max_val > 100.0)) {
            amrex::Print() << "  *** WARNING: implausibly large rain values detected! ***\n";
        }
    }
}

void
ERFRemoraMultiBlockContainer::ValidateFields (const char* label,
                                              const amrex::Vector<amrex::MultiFab*>& fields,
                                              const amrex::Vector<std::string>& names,
                                              amrex::Real time,
                                              long long ocean_step) const
{
    if (fields.size() != names.size()) {
        amrex::Abort("ValidateFields field/name count mismatch.");
    }

    for (int idx = 0; idx < static_cast<int>(fields.size()); ++idx) {
        const auto* mf = fields[idx];
        if (mf == nullptr) {
            amrex::Abort(std::string(label) + ": null field at index " + std::to_string(idx));
        }

        const amrex::Real min_val = mf->min(0);
        const amrex::Real max_val = mf->max(0);

        if (amrex::ParallelDescriptor::IOProcessor()) {
            amrex::Print() << label
                           << " ocean_step=" << ocean_step
                           << " time=" << time
                           << " field=" << names[idx]
                           << " nboxes=" << mf->boxArray().size()
                           << " min=" << min_val
                           << " max=" << max_val
                           << "\n";
        }

        if (!std::isfinite(min_val) || !std::isfinite(max_val)) {
            amrex::Abort(std::string(label) + ": non-finite extrema in field '" +
                         names[idx] + "' at ocean_step=" + std::to_string(ocean_step));
        }
    }
}

void
ERFRemoraMultiBlockContainer::WriteExchangePlotfiles (amrex::Real time, long long ocean_step) const
{
    if (m_use_two_way_coupling) {
        ValidateFields("Exchange validation ocean->atmos",
                       m_ocean_to_atmos, coupling::OceanToAtmosNames(), time, ocean_step);
    }
    ValidateFields("Exchange validation atmos->ocean",
                   m_atmos_to_ocean_active, m_atmos_to_ocean_active_names, time, ocean_step);

    if (m_use_two_way_coupling) {
        WriteExchangePlotfile(m_plot_file_exchange_o2a,
                              m_ocean_to_atmos,
                              coupling::OceanToAtmosNames(),
                              m_geom,
                              time,
                              ocean_step);
    }
    WriteExchangePlotfile(m_plot_file_exchange_a2o,
                          m_atmos_to_ocean_active,
                          m_atmos_to_ocean_active_names,
                          m_geom_ocn,
                          time,
                          ocean_step);
}

void
ERFRemoraMultiBlockContainer::WriteExchangePlotfile (
    const std::string& prefix,
    const amrex::Vector<amrex::MultiFab*>& fields,
    const amrex::Vector<std::string>& names,
    const amrex::Geometry& geom,
    amrex::Real time,
    long long ocean_step) const
{
    if (fields.empty()) {
        return;
    }
    if (fields.size() != names.size()) {
        amrex::Abort("Exchange plotfile field/name count mismatch.");
    }

    const int ncomp = static_cast<int>(fields.size());
    amrex::BoxArray plot_ba = amrex::convert(fields[0]->boxArray(), amrex::IntVect::TheCellVector());
    amrex::MultiFab plot_mf(plot_ba, fields[0]->DistributionMap(), ncomp, 0);
    for (int comp = 0; comp < ncomp; ++comp) {
        if (fields[comp] == nullptr) {
            amrex::Abort("Cannot write exchange plotfile with a null field.");
        }
        if (amrex::convert(fields[comp]->boxArray(), amrex::IntVect::TheCellVector()) != plot_ba) {
            amrex::Abort("Exchange plotfile writer requires a common enclosed-cell slab layout.");
        }
        CopyFieldToCellCenteredPlot(plot_mf, comp, *fields[comp]);
    }

    const std::string plotfilename = amrex::Concatenate(prefix, static_cast<int>(ocean_step), 5);
    amrex::WriteSingleLevelPlotfile(plotfilename, plot_mf, names, geom, time,
                                    static_cast<int>(ocean_step));
}

void
ERFRemoraMultiBlockContainer::AdvanceBlocks (int ocean_steps)
{
    // V1 Hard Invariants: Verify required solver instances are active
    if (!m_erf || !m_remora) {
        amrex::Abort("Coupling loop started without initialized ERF or REMORA solvers.");
    }

    // Verify subcycling ratio is strictly positive
    if (m_atm_substeps_per_ocean_step <= 0) {
        amrex::Abort("Invalid atmosphere subcycling ratio: must be > 0.");
    }

    AdvanceBlocksWithCallbacks(
        [this](amrex::Vector<amrex::MultiFab*>& state, amrex::Real time) {
            BL_PROFILE("Driver::Pack_Surface");
            REMORAPackSurfaceState(m_remora, state, time);

            // Finite check on packed surface state
            for (auto* mf : state) {
                if (mf) {
                    // We use a simple check: if any value is non-finite, abort.
                    // In a production environment, we'd use a reduction.
                    const amrex::Real mn = mf->min(0);
                    const amrex::Real mx = mf->max(0);
                    if (mn != mn || mx != mx) {
                        amrex::Abort("Non-finite value detected in PackSurfaceState output.");
                    }
                }
            }
        },
        [this](const amrex::Vector<amrex::MultiFab*>& state, amrex::Real time) {
            BL_PROFILE("Driver::Apply_Surface");

            if (m_domains_match) {
                // Fast path: domains match, apply directly from ocean-resolution cache
                ERFApplyOceanSurfaceState(m_erf, state, time);
            } else {
                // Slow path: domains don't match, interpolate with fallback for uncovered regions
                BL_PROFILE("Driver::Interp_SST");

                // Create temporary atmosphere-resolution MultiFab
                amrex::Vector<std::unique_ptr<amrex::MultiFab>> atm_temp_storage(1);
                atm_temp_storage[0] = std::make_unique<amrex::MultiFab>(m_ba, m_dm, 1, 0);
                amrex::Vector<amrex::MultiFab*> atm_temp = amrex::GetVecOfPtrs(atm_temp_storage);

                // Interpolate ocean SST to atmosphere resolution with fallback
                InterpOceanToAtmos(state, atm_temp);

                // Apply to ERF
                ERFApplyOceanSurfaceState(m_erf, atm_temp, time);
            }
        },
        [this](amrex::Vector<amrex::MultiFab*>& states, amrex::Real time) {
            BL_PROFILE("Driver::Pack_Atmos");

            if (m_domains_match) {
                // Fast path: domains match, pack directly into ocean-resolution cache
                ERFPackAtmosphericStates(m_erf, states, time);
            } else {
                // Slow path: domains don't match, need spatial remapping
                BL_PROFILE("Driver::Pack_Atmos_Remap");

                // Step 1: Create temporary atmosphere-resolution MultiFabs
                amrex::Vector<std::unique_ptr<amrex::MultiFab>> atm_temp_storage(states.size());
                for (int idx = 0; idx < static_cast<int>(states.size()); ++idx) {
                    if (states[idx] != nullptr) {
                        const amrex::BoxArray ba_temp = AtmosToOceanBoxArray(m_ba, idx);
                        const amrex::IntVect ng = AtmosToOceanNGrow(idx);
                        atm_temp_storage[idx] = std::make_unique<amrex::MultiFab>(ba_temp, m_dm, 1, ng);
                    }
                }
                amrex::Vector<amrex::MultiFab*> atm_temp = amrex::GetVecOfPtrs(atm_temp_storage);

                // Step 2: ERF packs into temp (full atmosphere extent)
                ERFPackAtmosphericStates(m_erf, atm_temp, time);

                // Step 3: Remap from atmosphere to ocean resolution/extent
                RemapAtmosToOceanStates(atm_temp, states);
            }
            static int cloud_warn_count = 0;
            constexpr int cloud_warn_limit = 5;
            constexpr amrex::Real cloud_tol = 1.0e-6;

            // Finite check on packed atmospheric states
            for (int idx = 0; idx < static_cast<int>(states.size()); ++idx) {
                auto* mf = states[idx];
                if (mf) {
                    if (mf->min(0) != mf->min(0) || mf->max(0) != mf->max(0)) {
                        amrex::Abort("Non-finite value detected in PackAtmosphericStates output.");
                    }
                    const bool is_cloud_lane = (idx < static_cast<int>(m_atmos_to_ocean_active_indices.size()) &&
                                                m_atmos_to_ocean_active_indices[idx] == coupling::Cloud);
                    if (is_cloud_lane && cloud_warn_count < cloud_warn_limit) {
                        const amrex::Real cmin = mf->min(0);
                        const amrex::Real cmax = mf->max(0);
                        if (amrex::ParallelDescriptor::IOProcessor()) {
                            if (cmin < -cloud_tol || cmax > amrex::Real(1.0) + cloud_tol) {
                                ++cloud_warn_count;
                                amrex::Print() << "WARNING: packed cloud lane outside tolerance [-1e-6,1+1e-6]: min="
                                                << cmin << " max=" << cmax << "\n";
                            }
                        }
                    }
                }
            }
            const long long ocean_step =
                m_atm_tick / static_cast<long long>(m_atm_substeps_per_ocean_step);
            ValidateFields("Post-pack atmos->ocean validation",
                           states, m_atmos_to_ocean_active_names, time, ocean_step);
        },
        [this](const amrex::Vector<amrex::MultiFab*>& states, amrex::Real time) {
            BL_PROFILE("Driver::Apply_Atmos");
            const long long ocean_step =
                m_atm_tick / static_cast<long long>(m_atm_substeps_per_ocean_step);
            ValidateFields("Pre-apply atmos->ocean validation",
                           states, m_atmos_to_ocean_active_names, time, ocean_step);
            if (ActiveAtmosToOceanContract(ocean_step) == coupling::AtmosToOceanContract::FluxPassing) {
                REMORAApplyAtmosphericFluxes(m_remora, states, time);
            } else {
                REMORAApplyAtmosphericStates(m_remora, states, time);
            }
            ValidateFields("Post-apply atmos->ocean validation",
                           states, m_atmos_to_ocean_active_names, time, ocean_step);
        },
        ocean_steps);
}

void
ERFRemoraMultiBlockContainer::AdvanceBlocksSynthetic (int ocean_steps)
{
    AdvanceBlocksWithCallbacks(
        [] (amrex::Vector<amrex::MultiFab*>& state, amrex::Real time) {
            if (!state.empty() && state[0] != nullptr) {
                state[0]->setVal(amrex::Real(290.0) + amrex::Real(0.01) * time);
            }
        },
        [] (const amrex::Vector<amrex::MultiFab*>&, amrex::Real) {},
        [] (amrex::Vector<amrex::MultiFab*>& states, amrex::Real time) {
            for (int idx = 0; idx < static_cast<int>(states.size()); ++idx) {
                if (states[idx] != nullptr) {
                    states[idx]->setVal(amrex::Real(10.0) + static_cast<amrex::Real>(idx + 1) + time);
                }
            }
        },
        [] (const amrex::Vector<amrex::MultiFab*>& states, amrex::Real time) {
            if (!states.empty() && states[0] != nullptr) {
                const amrex::Real sample = states[0]->min(0);
                if (amrex::ParallelDescriptor::IOProcessor()) {
                    amrex::Print() << "Synthetic ApplyAtmosphericStates time=" << time
                                   << " sample_atmos_field0=" << sample << "\n";
                }
            }
        },
        ocean_steps);
}

amrex::Vector<amrex::MultiFab*>&
ERFRemoraMultiBlockContainer::AtmosToOceanActiveCache ()
{
    return m_atmos_to_ocean_active;
}

const amrex::Vector<std::string>&
ERFRemoraMultiBlockContainer::AtmosToOceanActiveNames () const
{
    return m_atmos_to_ocean_active_names;
}

coupling::AtmosToOceanContract
ERFRemoraMultiBlockContainer::ActiveAtmosToOceanContract (long long ocean_step) const
{
    if (m_atmos_to_ocean_contract == coupling::AtmosToOceanContract::FluxPassing &&
        ocean_step == 0) {
        return coupling::AtmosToOceanContract::StatePassing;
    }
    return m_atmos_to_ocean_contract;
}

void
ERFRemoraMultiBlockContainer::SelectAtmosToOceanActiveView (long long ocean_step)
{
    const auto active_contract = ActiveAtmosToOceanContract(ocean_step);
    if (active_contract == coupling::AtmosToOceanContract::FluxPassing) {
        m_atmos_to_ocean_active = m_atmos_to_ocean_flux_active;
        m_atmos_to_ocean_active_indices = m_atmos_to_ocean_flux_active_indices;
        m_atmos_to_ocean_active_names = m_atmos_to_ocean_flux_active_names;
    } else {
        m_atmos_to_ocean_active = m_atmos_to_ocean_state_active;
        m_atmos_to_ocean_active_indices = m_atmos_to_ocean_state_active_indices;
        m_atmos_to_ocean_active_names = m_atmos_to_ocean_state_active_names;
    }
}

void
ERFRemoraMultiBlockContainer::AdvanceBlocksWithCallbacks (const MutableStateFn& pack_surface_state,
                                                           const ConstStateFn& apply_ocean_surface_state,
                                                           const MutableStateFn& pack_atmos_states,
                                                           const ConstStateFn& apply_atmos_states,
                                                           int ocean_steps)
{
    BL_PROFILE("Driver::Build_Comm_Metadata"); // V1 Baseline marker
    if (m_use_two_way_coupling) {
        PrintDiagnostics("Initial ocean->atmos cache", m_ocean_to_atmos, coupling::OceanToAtmosNames());
    }
    PrintDiagnostics("Initial atmos->ocean cache", m_atmos_to_ocean_active, m_atmos_to_ocean_active_names);


    // Multi-rate clock rationale:
    // - ERF advances a compressible atmosphere where acoustic stability often drives
    //   dt to O(1-10 s).
    // - REMORA advances Boussinesq/split-explicit ocean dynamics with a much larger
    //   macro-step O(10^2-10^3 s).
    // We therefore hold one ocean macro-state while subcycling atmosphere updates.
    const amrex::Real dt_ocn_request = m_erf_dt * static_cast<amrex::Real>(m_atm_substeps_per_ocean_step);
    const long long starting_atm_tick = m_atm_tick;
    const long long expected_atm_ticks = starting_atm_tick +
                                         static_cast<long long>(ocean_steps) *
                                         static_cast<long long>(m_atm_substeps_per_ocean_step);
    const amrex::Real abs_tol = 1.0e-12;
    const amrex::Real rel_tol = 1.0e-8;

    amrex::Print() << "Coupled cadence plan:"
                       << " ocean_steps=" << ocean_steps
                       << " starting_atm_tick=" << starting_atm_tick
                       << " atm_substeps_per_ocean_step=" << m_atm_substeps_per_ocean_step
                       << " target_atm_tick=" << expected_atm_ticks
                       << " expected_run_time=" << (dt_ocn_request * static_cast<amrex::Real>(ocean_steps))
                       << std::endl;

    auto enforce_dt = [&] (const char* solver_name, amrex::Real dt_request, amrex::Real dt_used,
                           int ocean_step, int substep) {
        const amrex::Real err = amrex::Math::abs(dt_used - dt_request);
        const amrex::Real tol = amrex::max(abs_tol, rel_tol * amrex::Math::abs(dt_request));
        if (err > tol) {
            amrex::Print() << "Coupled cadence mismatch before abort:"
                           << " solver=" << solver_name
                           << " ocean_step=" << ocean_step
                           << " substep=" << substep
                           << " dt_request=" << dt_request
                           << " dt_used=" << dt_used
                           << " abs_err=" << err
                           << " tol=" << tol
                           << " m_erf_dt=" << m_erf_dt
                           << " atm_substeps_per_ocean_step=" << m_atm_substeps_per_ocean_step
                           << " expected_ocean_macro_dt=" << dt_ocn_request
                           << std::endl;
            amrex::Abort(std::string("Coupled cadence drift in ") + solver_name +
                         " at ocean_step=" + std::to_string(ocean_step) +
                         " substep=" + std::to_string(substep) +
                         " dt_request=" + std::to_string(dt_request) +
                         " dt_used=" + std::to_string(dt_used) +
                         " abs_err=" + std::to_string(err) +
                         " tol=" + std::to_string(tol));
        }
    };

    for (int ocean_step = 0; ocean_step < ocean_steps; ++ocean_step) {
        bool stop_coupling = false;
        if (m_atm_tick >= expected_atm_ticks) {
            amrex::Print() << "Reached expected coupled end time before ocean_step "
                                << ocean_step << "; stopping cleanly at atmosphere tick "
                                << m_atm_tick << " of " << expected_atm_ticks << std::endl;
            break;
        }
        const amrex::Real ocean_time = m_atm_tick * m_erf_dt;
        const long long global_ocean_step =
            m_atm_tick / static_cast<long long>(m_atm_substeps_per_ocean_step);
        SelectAtmosToOceanActiveView(global_ocean_step);
        const bool use_state_contract =
            (ActiveAtmosToOceanContract(global_ocean_step) == coupling::AtmosToOceanContract::StatePassing);
        if (global_ocean_step == 0 && amrex::ParallelDescriptor::IOProcessor()) {
            amrex::Print() << "Bootstrap atm2ocn contract at ocean_step=0: "
                           << (use_state_contract ? "state" : "flux")
                           << " (configured mode="
                           << coupling::AtmosToOceanContractName(m_atmos_to_ocean_contract)
                           << ")\n";
        }
        DriverPrintLine(
            m_driver_verbose,
            std::string("ocean_step=") + std::to_string(global_ocean_step) +
                " set ERF atm2ocn contract to " + (use_state_contract ? "state" : "flux"));
        SetERFDriverAtmosToOceanStateContract(m_erf, use_state_contract);
        DriverPrintLine(
            m_driver_verbose,
            std::string("ocean_step=") + std::to_string(global_ocean_step) +
                " set REMORA atm2ocn contract to " + (use_state_contract ? "state" : "flux"));
        SetREMORADriverAtmosToOceanStateContract(m_remora, use_state_contract);
        if (m_use_two_way_coupling) {
            DriverPrintLine(
                m_driver_verbose,
                std::string("ocean_step=") + std::to_string(global_ocean_step) +
                    " pack ocean->atmos surface state from REMORA at time=" +
                    std::to_string(ocean_time));
            pack_surface_state(m_ocean_to_atmos, ocean_time);
            FillOceanToAtmosGhosts();
        }
        DriverPrintLine(
            m_driver_verbose,
            std::string("ocean_step=") + std::to_string(global_ocean_step) +
                " pack atmos->ocean " + (use_state_contract ? "state" : "flux") +
                " from ERF at time=" + std::to_string(ocean_time));
        pack_atmos_states(m_atmos_to_ocean_active, ocean_time);
        if (m_plot_int_exchange > 0 &&
            (global_ocean_step % static_cast<long long>(m_plot_int_exchange)) == 0) {
            WriteExchangePlotfiles(ocean_time, global_ocean_step);
        }
        FillAtmosToOceanGhosts();
        DriverPrintLine(
            m_driver_verbose,
            std::string("ocean_step=") + std::to_string(global_ocean_step) +
                " apply atmos->ocean " + (use_state_contract ? "state" : "flux") +
                " to REMORA at time=" + std::to_string(ocean_time));
        apply_atmos_states(m_atmos_to_ocean_active, ocean_time);


        // Initial architecture note: this passes instantaneous slabs at the macro-step boundary.
        // Future update: use old/new cached states and time interpolation during this
        // loop (same temporal pattern as erf-amrwind-driver boundary ingestion in
        // Source/wind_energy/ABLReadERF.cpp) to avoid stepwise forcing artifacts.
        for (int substep = 0; substep < m_atm_substeps_per_ocean_step; ++substep) {
            if (m_atm_tick >= expected_atm_ticks) {
                amrex::Print() << "Reached expected coupled end time before ERF advance at ocean_step="
                                    << ocean_step << " substep=" << substep
                                    << "; skipping remaining substeps and exiting cleanly."
                                    << std::endl;
                stop_coupling = true;
                break;
            }
            const amrex::Real atm_time = m_atm_tick * m_erf_dt;

            BL_PROFILE("Driver::Time_Interpolate"); // V1 Stub region
            if (m_use_two_way_coupling) {
                DriverPrintLine(
                    m_driver_verbose,
                    std::string("atm_tick=") + std::to_string(m_atm_tick) +
                        " apply ocean->atmos surface state to ERF at time=" +
                        std::to_string(atm_time));
                apply_ocean_surface_state(m_ocean_to_atmos, atm_time);
            }
            const amrex::Real dt_used = ERFEvolveOneStep(m_erf, atm_time, m_erf_dt);
            if (dt_used == amrex::Real(0.0)) {
                amrex::Print() << "ERF returned dt=0 before the driver completed the planned"
                                    << " atmosphere subcycling schedule."
                                    << " ocean_step=" << ocean_step
                                    << " substep=" << substep
                                    << " atm_tick=" << m_atm_tick
                                    << " expected_atm_ticks=" << expected_atm_ticks
                                    << std::endl;
                amrex::Abort("ERF terminated before the driver reached a valid coupled macro-step boundary.");
            }
            enforce_dt("ERF", m_erf_dt, dt_used, ocean_step, substep);

            ++m_atm_tick;
        }

        if (stop_coupling) {
            break;
        }

        const amrex::Real ocn_dt_used = REMORAEvolveOneStep(m_remora, ocean_time, dt_ocn_request);
        if (ocn_dt_used == amrex::Real(0.0)) {
            amrex::Print() << "REMORA returned dt=0 before the driver completed the planned"
                                << " ocean macro-step schedule."
                                << " ocean_step=" << ocean_step
                                << " atm_tick=" << m_atm_tick
                                << " expected_atm_ticks=" << expected_atm_ticks
                                << std::endl;
            amrex::Abort("REMORA terminated before the driver reached a valid coupled macro-step boundary.");
        }
        enforce_dt("REMORA", dt_ocn_request, ocn_dt_used, ocean_step, -1);

        amrex::Print() << "Completed ocean step " << ocean_step
                           << " at atmosphere tick " << m_atm_tick << "\n";

        if (m_use_two_way_coupling) {
            PrintDiagnostics("Ocean->atmos cache", m_ocean_to_atmos, coupling::OceanToAtmosNames());
        }
        PrintDiagnostics("Atmos->ocean cache", m_atmos_to_ocean_active, m_atmos_to_ocean_active_names);
    }
    REMORAWriteAtFinalTime(m_remora);
    ERFWriteAtFinalTime(m_erf);

}
} // namespace erf_remora
