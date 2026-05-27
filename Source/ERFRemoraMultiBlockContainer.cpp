#include "ERFRemoraMultiBlockContainer.H"
#include "StandaloneSolverRuns.H"

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Math.H>
#include <AMReX_ParmParse.H>

#include <string>

namespace erf_remora
{
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
    int max_grid_size,
    amrex::Real dt_atm,
    int atm_substeps_per_ocean_step,
    coupling::AtmosToOceanContract atm2ocn_contract)
    : m_atm_rb(atm_rb), m_atm_n_cell(n_cell_atm),
      m_ocn_rb(ocn_rb), m_ocn_n_cell(n_cell_ocn),
      m_max_level(max_level), m_coord(coord),
      m_ref_ratio(ref_ratio), m_is_periodic(is_periodic),
      m_dt_atm(dt_atm),
      m_atm_substeps_per_ocean_step(atm_substeps_per_ocean_step),
      m_atmos_to_ocean_contract(atm2ocn_contract)
{
    // Interface slab: horizontal extent of the atmosphere at z = 0.
    amrex::Box domain(amrex::IntVect(0, 0, 0),
                      amrex::IntVect(n_cell_atm[0]-1, n_cell_atm[1]-1, 0));
    m_ba.define(domain);
    m_ba.maxSize(max_grid_size);
    m_dm.define(m_ba);

    // Ocean-resolution slab for atm->ocn fields.
    amrex::Box ocn_domain(amrex::IntVect(0, 0, 0),
                          amrex::IntVect(n_cell_ocn[0]-1, n_cell_ocn[1]-1, 0));
    m_ba_ocn.define(ocn_domain);
    m_ba_ocn.maxSize(max_grid_size);
    m_dm_ocn.define(m_ba_ocn);

    amrex::Real dz_atm = (atm_rb.hi(2) - atm_rb.lo(2)) / n_cell_atm[2];
    amrex::RealBox slab_rb({atm_rb.lo(0), atm_rb.lo(1), 0.0},
                            {atm_rb.hi(0), atm_rb.hi(1), dz_atm});
    m_geom.define(domain, &slab_rb, amrex::CoordSys::cartesian, is_periodic.data());

    m_atmos_to_ocean_active_indices = coupling::AtmosToOceanActiveIndices(m_atmos_to_ocean_contract);
    m_atmos_to_ocean_active_names = coupling::AtmosToOceanNames(m_atmos_to_ocean_contract);
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
    // NOTE: remora.air_humidity is relative humidity [0-1], but the contract slot
    // Qair carries specific humidity [kg/kg] (vec_qair path). A direct copy of 0.776
    // as specific humidity gives ~50× saturation and blows up the bulk flux.
    // Use 0 here; when ERF has prognostic moisture it will overwrite with the
    // correct specific humidity from its state.
    amrex::ParmParse pp_rem("remora");
    amrex::Real cloud_val = 0.0, rain_val = 0.0;
    pp_rem.query("cloud", cloud_val);
    pp_rem.query("rain",  rain_val);

    for (int local = 0; local < static_cast<int>(m_atmos_to_ocean_active_indices.size()); ++local) {
        auto* mf = m_atmos_to_ocean_active[local];
        if (mf == nullptr) {
            continue;
        }

        switch (m_atmos_to_ocean_active_indices[local]) {
        case coupling::Qair:
            mf->setVal(0.0);       // qair [kg/kg]: dry-atmosphere fallback
            break;
        case coupling::Cloud:
            mf->setVal(cloud_val); // cloud fraction [0-1]
            break;
        case coupling::Rain:
            mf->setVal(rain_val);  // rain rate [kg m-2 s-1]
            break;
        case coupling::SWrad:
            mf->setVal(0.0);       // shortwave [W m-2]: no radiation fallback
            break;
        case coupling::LWrad:
            mf->setVal(0.0);       // longwave  [W m-2]: no radiation fallback
            break;
        default:
            break;
        }
    }
}

void
ERFRemoraMultiBlockContainer::InitializeBlocks ()
{
    m_atm_tick = 0;
    AllocateCaches();
    FillSyntheticData();

    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "Atmos->ocean contract mode: "
                       << coupling::AtmosToOceanContractName(m_atmos_to_ocean_contract)
                       << " (" << static_cast<int>(m_atmos_to_ocean_active_indices.size())
                       << " active layers)\n";
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
    for (auto& mf_ptr : m_atmos_to_ocean_storage) {
        mf_ptr = std::make_unique<amrex::MultiFab>(m_ba_ocn, m_dm_ocn, 1, 0);
    }
    m_atmos_to_ocean = amrex::GetVecOfPtrs(m_atmos_to_ocean_storage);
    m_atmos_to_ocean_active.resize(m_atmos_to_ocean_active_indices.size());
    for (int local_idx = 0; local_idx < static_cast<int>(m_atmos_to_ocean_active_indices.size()); ++local_idx) {
        const int global_idx = m_atmos_to_ocean_active_indices[local_idx];
        m_atmos_to_ocean_active[local_idx] = m_atmos_to_ocean[global_idx];
    }

    m_ocean_to_atmos_storage.resize(coupling::NumOceanToAtmosStates);
    for (auto& mf_ptr : m_ocean_to_atmos_storage) {
        mf_ptr = std::make_unique<amrex::MultiFab>(m_ba, m_dm, 1, 0);
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
}

void
ERFRemoraMultiBlockContainer::PrintDiagnostics (const char* label,
                             const amrex::Vector<amrex::MultiFab*>& fields,
                             const amrex::Vector<std::string>& names) const
{
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

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
                    if (mf->min(0) != mf->min(0) || mf->max(0) != mf->max(0)) {
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
                    if (is_cloud_lane && amrex::ParallelDescriptor::IOProcessor() &&
                        cloud_warn_count < cloud_warn_limit) {
                        const amrex::Real cmin = mf->min(0);
                        const amrex::Real cmax = mf->max(0);
                        if (cmin < -cloud_tol || cmax > amrex::Real(1.0) + cloud_tol) {
                            ++cloud_warn_count;
                            amrex::Print() << "WARNING: packed cloud lane outside tolerance [-1e-6,1+1e-6]: min="
                                           << cmin << " max=" << cmax << "\n";
                        }
                    }
                }
            }
        },
        [this](const amrex::Vector<amrex::MultiFab*>& states, amrex::Real time) {
            BL_PROFILE("Driver::Apply_Atmos");
            REMORAApplyAtmosphericStates(m_remora, states, time);
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
            if (amrex::ParallelDescriptor::IOProcessor() && !states.empty() && states[0] != nullptr) {
                const amrex::Real sample = states[0]->min(0);
                amrex::Print() << "Synthetic ApplyAtmosphericStates time=" << time
                               << " sample_atmos_field0=" << sample << "\n";
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
    const amrex::Real dt_ocn_request = m_dt_atm * static_cast<amrex::Real>(m_atm_substeps_per_ocean_step);
    const amrex::Real abs_tol = 1.0e-12;
    const amrex::Real rel_tol = 1.0e-8;

    auto enforce_dt = [&] (const char* solver_name, amrex::Real dt_request, amrex::Real dt_used,
                           int ocean_step, int substep) {
        const amrex::Real err = amrex::Math::abs(dt_used - dt_request);
        const amrex::Real tol = amrex::max(abs_tol, rel_tol * amrex::Math::abs(dt_request));
        if (err > tol) {
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
        const amrex::Real ocean_time = m_atm_tick * m_dt_atm;
        pack_surface_state(m_ocean_to_atmos, ocean_time);

        // Initial architecture note: this passes instantaneous slabs at the macro-step boundary.
        // Future update: use old/new cached states and time interpolation during this
        // loop (same temporal pattern as erf-amrwind-driver boundary ingestion in
        // Source/wind_energy/ABLReadERF.cpp) to avoid stepwise forcing artifacts.
        for (int substep = 0; substep < m_atm_substeps_per_ocean_step; ++substep) {
            const amrex::Real atm_time = m_atm_tick * m_dt_atm;

            BL_PROFILE("Driver::Time_Interpolate"); // V1 Stub region
            apply_ocean_surface_state(m_ocean_to_atmos, atm_time);
            enforce_dt("ERF", m_dt_atm, ERFEvolveOneStep(m_erf, atm_time, m_dt_atm), ocean_step, substep);
            
            BL_PROFILE("Driver::AMReX_ParallelCopy"); // Region around state transfer
            pack_atmos_states(m_atmos_to_ocean_active, atm_time);
            apply_atmos_states(m_atmos_to_ocean_active, atm_time);

            ++m_atm_tick;
        }

        enforce_dt("REMORA", dt_ocn_request, REMORAEvolveOneStep(m_remora, ocean_time, dt_ocn_request),
                   ocean_step, -1);

        if (amrex::ParallelDescriptor::IOProcessor()) {
            amrex::Print() << "Completed ocean step " << ocean_step
                           << " at atmosphere tick " << m_atm_tick << "\n";
        }

        PrintDiagnostics("Ocean->atmos cache", m_ocean_to_atmos, coupling::OceanToAtmosNames());
        PrintDiagnostics("Atmos->ocean cache", m_atmos_to_ocean_active, m_atmos_to_ocean_active_names);
    }
}
} // namespace erf_remora
