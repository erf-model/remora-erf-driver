#include "ERFRemoraMultiBlockContainer.H"

#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

namespace erf_remora
{
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

    amrex::Real dz_atm = (atm_rb.hi(2) - atm_rb.lo(2)) / n_cell_atm[2];
    amrex::RealBox slab_rb({atm_rb.lo(0), atm_rb.lo(1), 0.0},
                            {atm_rb.hi(0), atm_rb.hi(1), dz_atm});
    m_geom.define(domain, &slab_rb, amrex::CoordSys::cartesian, is_periodic.data());

    m_atmos_to_ocean_active_indices = coupling::AtmosToOceanActiveIndices(m_atmos_to_ocean_contract);
    m_atmos_to_ocean_active_names = coupling::AtmosToOceanNames(m_atmos_to_ocean_contract);
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
    // RH carries specific humidity [kg/kg] (vec_qair path). A direct copy of 0.776
    // as specific humidity gives ~50× saturation and blows up the bulk flux.
    // Use 0 here; when ERF has prognostic moisture it will overwrite with the
    // correct specific humidity from its state.
    amrex::ParmParse pp_rem("remora");
    amrex::Real cloud_val = 0.0, rain_val = 0.0;
    pp_rem.query("cloud", cloud_val);
    pp_rem.query("rain",  rain_val);

    auto set_slot = [&] (int global_idx, amrex::Real val) {
        const auto& active = m_atmos_to_ocean_active_indices;
        for (int local = 0; local < static_cast<int>(active.size()); ++local) {
            if (active[local] == global_idx && m_atmos_to_ocean_active[local] != nullptr) {
                m_atmos_to_ocean_active[local]->setVal(val);
                break;
            }
        }
    };

    set_slot(coupling::RH,    0.0);       // qair [kg/kg]: dry-atmosphere fallback
    set_slot(coupling::Cloud, cloud_val); // cloud fraction [0-1]
    set_slot(coupling::Rain,  rain_val);  // rain rate [kg m-2 s-1]
    set_slot(coupling::SWrad, 0.0);       // shortwave [W m-2]: no radiation fallback
    set_slot(coupling::LWrad, 0.0);       // longwave  [W m-2]: no radiation fallback
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
        mf_ptr = std::make_unique<amrex::MultiFab>(m_ba, m_dm, 1, 0);
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
    }
}

void
ERFRemoraMultiBlockContainer::AdvanceBlocks (ERF& erf, REMORA& remora, int ocean_steps)
{
    amrex::ignore_unused(erf, remora, ocean_steps);
    amrex::Abort("Solver-backed AdvanceBlocks() is disabled in initial architecture infrastructure mode. Set driver.synthetic_only = 1.");
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
    PrintDiagnostics("Initial ocean->atmos cache", m_ocean_to_atmos, coupling::OceanToAtmosNames());
    PrintDiagnostics("Initial atmos->ocean cache", m_atmos_to_ocean_active, m_atmos_to_ocean_active_names);

    // Multi-rate clock rationale:
    // - ERF advances a compressible atmosphere where acoustic stability often drives
    //   dt to O(1-10 s).
    // - REMORA advances Boussinesq/split-explicit ocean dynamics with a much larger
    //   macro-step O(10^2-10^3 s).
    // We therefore hold one ocean macro-state while subcycling atmosphere updates.
    for (int ocean_step = 0; ocean_step < ocean_steps; ++ocean_step) {
        const amrex::Real ocean_time = m_atm_tick * m_dt_atm;
        pack_surface_state(m_ocean_to_atmos, ocean_time);

        // Initial architecture note: this passes instantaneous slabs at the macro-step boundary.
        // Future update: use old/new cached states and time interpolation during this
        // loop (same temporal pattern as erf-amrwind-driver boundary ingestion in
        // Source/wind_energy/ABLReadERF.cpp) to avoid stepwise forcing artifacts.
        for (int substep = 0; substep < m_atm_substeps_per_ocean_step; ++substep) {
            const amrex::Real atm_time = m_atm_tick * m_dt_atm;

            apply_ocean_surface_state(m_ocean_to_atmos, atm_time);
            pack_atmos_states(m_atmos_to_ocean_active, atm_time);
            apply_atmos_states(m_atmos_to_ocean_active, atm_time);

            ++m_atm_tick;
        }

        if (amrex::ParallelDescriptor::IOProcessor()) {
            amrex::Print() << "Completed ocean step " << ocean_step
                           << " at atmosphere tick " << m_atm_tick << "\n";
        }

        PrintDiagnostics("Ocean->atmos cache", m_ocean_to_atmos, coupling::OceanToAtmosNames());
        PrintDiagnostics("Atmos->ocean cache", m_atmos_to_ocean_active, m_atmos_to_ocean_active_names);
    }
}
} // namespace erf_remora
