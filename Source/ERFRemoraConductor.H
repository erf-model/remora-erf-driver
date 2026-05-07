#ifndef ERF_REMORA_CONDUCTOR_H_
#define ERF_REMORA_CONDUCTOR_H_

#include "ERFRemoraCouplingContract.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <functional>
#include <memory>

class ERF;
class REMORA;

namespace erf_remora
{
class MultiBlockContainer
{
public:
    MultiBlockContainer(int n_cell, int max_grid_size, amrex::Real dt_atm, int atm_substeps_per_ocean_step);
    ~MultiBlockContainer() = default;

    void InitializeBlocks();
    void AdvanceBlocks(ERF& erf, REMORA& remora, int ocean_steps);
    void AdvanceBlocksSynthetic(int ocean_steps);

private:
    using MutableStateFn = std::function<void(amrex::Vector<amrex::MultiFab*>&, amrex::Real)>;
    using ConstStateFn = std::function<void(const amrex::Vector<amrex::MultiFab*>&, amrex::Real)>;

    void AdvanceBlocksWithCallbacks(const MutableStateFn& pack_surface_state,
                                    const ConstStateFn& apply_ocean_surface_state,
                                    const MutableStateFn& pack_atmos_states,
                                    const ConstStateFn& apply_atmos_states,
                                    int ocean_steps);

    void AllocateCaches();
    void FillSyntheticData();
    void PrintDiagnostics(const char* label,
                          const amrex::Vector<amrex::MultiFab*>& fields,
                          const amrex::Vector<std::string>& names) const;

    amrex::BoxArray m_ba;
    amrex::DistributionMapping m_dm;
    amrex::Geometry m_geom;

    amrex::Real m_dt_atm;
    int m_atm_substeps_per_ocean_step;
    long long m_atm_tick = 0;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_atmos_to_ocean_storage;
    amrex::Vector<amrex::MultiFab*> m_atmos_to_ocean;

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_ocean_to_atmos_storage;
    amrex::Vector<amrex::MultiFab*> m_ocean_to_atmos;
};
} // namespace erf_remora

#endif
