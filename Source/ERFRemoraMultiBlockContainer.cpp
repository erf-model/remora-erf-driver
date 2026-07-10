#include "ERFRemoraMultiBlockContainer.H"
#include "StandaloneSolverRuns.H"

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Math.H>
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
                        const char* label)
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
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        fine_len[0] == coarse_len[0] * ratio[0] &&
        fine_len[1] == coarse_len[1] * ratio[1],
        "ValidateCouplingLayout: fine/coarse slab lengths do not match the requested ratio.");

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

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        coarsened_fine == coarse_ba,
        (std::string(label) +
         ": coarsen(fine_ba, ratio) does not match coarse_ba. "
         "Direct average_down is unsafe and helper-required layouts are disallowed.").c_str());

    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "[coupling-layout] " << label
                       << " ratio=" << ratio
                       << " direct_average_down_safe=yes\n";
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
    coupling::AtmosToOceanContract atm2ocn_contract)
    : m_atm_rb(atm_rb), m_atm_n_cell(n_cell_atm),
      m_ocn_rb(ocn_rb), m_ocn_n_cell(n_cell_ocn),
      m_max_level(max_level), m_coord(coord),
      m_ref_ratio(ref_ratio), m_is_periodic(is_periodic),
      m_erf_dt(erf_dt),
      m_atm_substeps_per_ocean_step(atm_substeps_per_ocean_step),
      m_atmos_to_ocean_contract(atm2ocn_contract)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        n_cell_atm[0] % n_cell_ocn[0] == 0 &&
        n_cell_atm[1] % n_cell_ocn[1] == 0,
        "Driver constructor requires integer x/y refinement ratio between atmosphere and ocean slabs.");
    const amrex::IntVect horiz_ratio(n_cell_atm[0] / n_cell_ocn[0],
                                     n_cell_atm[1] / n_cell_ocn[1], 1);

    // Interface slab: horizontal extent of the atmosphere at z = 0.
    amrex::Box domain(amrex::IntVect(0, 0, 0),
                      amrex::IntVect(n_cell_atm[0]-1, n_cell_atm[1]-1, 0));
    m_ba.define(domain);
    m_ba.maxSize(max_grid_size);
    m_dm.define(m_ba);
    PrintBoxArrayLayout("Driver constructor ocean->atmos slab", m_ba, m_dm);

    // Enforce structural compatibility by deriving the ocean slab directly
    // from the atmosphere slab decomposition.
    m_ba_ocn = m_ba;
    m_ba_ocn.coarsen(horiz_ratio);
    m_dm_ocn.define(m_ba_ocn);
    PrintBoxArrayLayout("Driver constructor atmos->ocean slab", m_ba_ocn, m_dm_ocn);
    ValidateCouplingLayout(m_ba, m_dm, m_ba_ocn, m_dm_ocn,
                           horiz_ratio,
                           "ERF->REMORA horizontal coupling");

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
    pp_driver.query("plot_int_exchange", m_plot_int_exchange);
    pp_driver.query("plot_file_exchange_o2a", m_plot_file_exchange_o2a);
    pp_driver.query("plot_file_exchange_a2o", m_plot_file_exchange_a2o);
}

void
ERFRemoraMultiBlockContainer::InitializeSolvers ()
{
    if (m_erf != nullptr) {
        DestroyERFSolver(m_erf);
        m_erf = nullptr;
    }
    if (m_remora != nullptr) {
        DestroyREMORASolver(m_remora);
        m_remora = nullptr;
    }

    m_erf = CreateERFSolver(*this);
    InitERFSolver(m_erf);

    m_remora = CreateREMORASolver(*this);
    InitREMORASolver(m_remora);
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
    m_atm_tick = 0;
    AllocateCaches();
    SelectAtmosToOceanActiveView(/*ocean_step=*/0);
    FillSyntheticData();

    if (m_atmos_to_ocean_contract == coupling::AtmosToOceanContract::FluxPassing) {
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
    // Allocate a superset layer vector once (legacy state + optional flux layers).
    // Runtime coupling modes can then swap active layer views by index without
    // changing memory ownership or reallocating cache slabs.
    m_atmos_to_ocean_storage.resize(coupling::NumAtmosToOceanLayers);
    for (int idx = 0; idx < static_cast<int>(m_atmos_to_ocean_storage.size()); ++idx) {
        auto& mf_ptr = m_atmos_to_ocean_storage[idx];
        mf_ptr = std::make_unique<amrex::MultiFab>(m_ba_ocn, m_dm_ocn, 1, 0);
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

    m_ocean_to_atmos_storage.resize(coupling::NumOceanToAtmosStates);
    const amrex::IntVect ocean_to_atmos_ng(1, 1, 0);
    for (int idx = 0; idx < static_cast<int>(m_ocean_to_atmos_storage.size()); ++idx) {
        auto& mf_ptr = m_ocean_to_atmos_storage[idx];
        mf_ptr = std::make_unique<amrex::MultiFab>(m_ba, m_dm, 1, ocean_to_atmos_ng);
        const std::string label = "Driver AllocateCaches ocean->atmos[" + std::to_string(idx) + "]";
        PrintMultiFabLayout(label.c_str(), *mf_ptr);
    }
    m_ocean_to_atmos = amrex::GetVecOfPtrs(m_ocean_to_atmos_storage);
}

void
ERFRemoraMultiBlockContainer::FillSyntheticData ()
{
    for (int idx = 0; idx < static_cast<int>(m_atmos_to_ocean.size()); ++idx) {
        m_atmos_to_ocean[idx]->setVal(10.0 + static_cast<amrex::Real>(idx));
    }

    m_ocean_to_atmos[coupling::SST]->setVal(300.0);
    FillOceanToAtmosGhosts();
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
    ValidateFields("Exchange validation ocean->atmos",
                   m_ocean_to_atmos, coupling::OceanToAtmosNames(), time, ocean_step);
    ValidateFields("Exchange validation atmos->ocean",
                   m_atmos_to_ocean_active, m_atmos_to_ocean_active_names, time, ocean_step);

    WriteExchangePlotfile(m_plot_file_exchange_o2a,
                          m_ocean_to_atmos,
                          coupling::OceanToAtmosNames(),
                          m_geom,
                          time,
                          ocean_step);
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
    amrex::MultiFab plot_mf(fields[0]->boxArray(), fields[0]->DistributionMap(), ncomp, 0);
    for (int comp = 0; comp < ncomp; ++comp) {
        if (fields[comp] == nullptr) {
            amrex::Abort("Cannot write exchange plotfile with a null field.");
        }
        amrex::MultiFab::Copy(plot_mf, *fields[comp], 0, comp, 1, 0);
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
            ERFApplyOceanSurfaceState(m_erf, state, time);
        },
        [this](amrex::Vector<amrex::MultiFab*>& states, amrex::Real time) {
            BL_PROFILE("Driver::Pack_Atmos");
            ERFPackAtmosphericStates(m_erf, states, time);
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
    PrintDiagnostics("Initial ocean->atmos cache", m_ocean_to_atmos, coupling::OceanToAtmosNames());
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
        pack_surface_state(m_ocean_to_atmos, ocean_time);
        FillOceanToAtmosGhosts();
        pack_atmos_states(m_atmos_to_ocean_active, ocean_time);
        if (m_plot_int_exchange > 0 &&
            (global_ocean_step % static_cast<long long>(m_plot_int_exchange)) == 0) {
            WriteExchangePlotfiles(ocean_time, global_ocean_step);
        }
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
            apply_ocean_surface_state(m_ocean_to_atmos, atm_time);
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

        PrintDiagnostics("Ocean->atmos cache", m_ocean_to_atmos, coupling::OceanToAtmosNames());
        PrintDiagnostics("Atmos->ocean cache", m_atmos_to_ocean_active, m_atmos_to_ocean_active_names);
    }
    REMORAWriteAtFinalTime(m_remora);
    ERFWriteAtFinalTime(m_erf);

}
} // namespace erf_remora
