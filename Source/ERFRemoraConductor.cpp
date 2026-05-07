#include "ERFRemoraConductor.H"

#include <AMReX.H>
#include <AMReX_Print.H>

namespace
{
void
fill_constant(amrex::MultiFab& mf, amrex::Real value)
{
    mf.setVal(value);
}
}

namespace erf_remora
{
MultiBlockContainer::MultiBlockContainer (int n_cell, int max_grid_size, amrex::Real dt_atm, int atm_substeps_per_ocean_step)
    : m_dt_atm(dt_atm),
      m_atm_substeps_per_ocean_step(atm_substeps_per_ocean_step)
{
    amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(n_cell - 1, n_cell - 1, 0));
    m_ba.define(domain);
    m_ba.maxSize(max_grid_size);
    m_dm.define(m_ba);

    amrex::RealBox real_box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    amrex::Vector<int> is_periodic(AMREX_SPACEDIM, 0);
    m_geom.define(domain, &real_box, amrex::CoordSys::cartesian, is_periodic.data());
}

void
MultiBlockContainer::InitializeBlocks ()
{
    m_atm_tick = 0;
    AllocateCaches();
    FillSyntheticData();
}

void
MultiBlockContainer::AllocateCaches ()
{
    m_atmos_to_ocean_storage.resize(coupling::NumAtmosToOceanStates);
    m_atmos_to_ocean.resize(coupling::NumAtmosToOceanStates);
    for (int idx = 0; idx < coupling::NumAtmosToOceanStates; ++idx) {
        m_atmos_to_ocean_storage[idx] = std::make_unique<amrex::MultiFab>(m_ba, m_dm, 1, 0);
        m_atmos_to_ocean[idx] = m_atmos_to_ocean_storage[idx].get();
    }

    m_ocean_to_atmos_storage.resize(coupling::NumOceanToAtmosStates);
    m_ocean_to_atmos.resize(coupling::NumOceanToAtmosStates);
    for (int idx = 0; idx < coupling::NumOceanToAtmosStates; ++idx) {
        m_ocean_to_atmos_storage[idx] = std::make_unique<amrex::MultiFab>(m_ba, m_dm, 1, 0);
        m_ocean_to_atmos[idx] = m_ocean_to_atmos_storage[idx].get();
    }
}

void
MultiBlockContainer::FillSyntheticData ()
{
    for (int idx = 0; idx < coupling::NumAtmosToOceanStates; ++idx) {
        fill_constant(*m_atmos_to_ocean[idx], 10.0 + static_cast<amrex::Real>(idx));
    }

    fill_constant(*m_ocean_to_atmos[coupling::SST], 300.0);
}

void
MultiBlockContainer::PrintDiagnostics (const char* label,
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
MultiBlockContainer::AdvanceBlocks (ERF& erf, REMORA& remora, int ocean_steps)
{
    amrex::ignore_unused(erf, remora, ocean_steps);
    amrex::Abort("Solver-backed AdvanceBlocks() is disabled in stage-0 synthetic mode. Set driver.synthetic_only = 1.");
}

void
MultiBlockContainer::AdvanceBlocksSynthetic (int ocean_steps)
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
                               << " sample_Uwind=" << sample << "\n";
            }
        },
        ocean_steps);
}

void
MultiBlockContainer::AdvanceBlocksWithCallbacks (const MutableStateFn& pack_surface_state,
                                                 const ConstStateFn& apply_ocean_surface_state,
                                                 const MutableStateFn& pack_atmos_states,
                                                 const ConstStateFn& apply_atmos_states,
                                                 int ocean_steps)
{
    PrintDiagnostics("Initial ocean->atmos cache", m_ocean_to_atmos, coupling::OceanToAtmosNames());
    PrintDiagnostics("Initial atmos->ocean cache", m_atmos_to_ocean, coupling::AtmosToOceanNames());

    for (int ocean_step = 0; ocean_step < ocean_steps; ++ocean_step) {
        const amrex::Real ocean_time = m_atm_tick * m_dt_atm;
        pack_surface_state(m_ocean_to_atmos, ocean_time);

        for (int substep = 0; substep < m_atm_substeps_per_ocean_step; ++substep) {
            const amrex::Real atm_time = m_atm_tick * m_dt_atm;

            apply_ocean_surface_state(m_ocean_to_atmos, atm_time);
            pack_atmos_states(m_atmos_to_ocean, atm_time);
            apply_atmos_states(m_atmos_to_ocean, atm_time);

            ++m_atm_tick;
        }

        if (amrex::ParallelDescriptor::IOProcessor()) {
            amrex::Print() << "Completed ocean step " << ocean_step
                           << " at atmosphere tick " << m_atm_tick << "\n";
        }

        PrintDiagnostics("Ocean->atmos cache", m_ocean_to_atmos, coupling::OceanToAtmosNames());
        PrintDiagnostics("Atmos->ocean cache", m_atmos_to_ocean, coupling::AtmosToOceanNames());
    }
}
} // namespace erf_remora
