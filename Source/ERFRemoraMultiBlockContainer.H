#ifndef ERF_REMORA_DRIVER_MULTIBLOCK_CONTAINER_H_
#define ERF_REMORA_DRIVER_MULTIBLOCK_CONTAINER_H_

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
/**
 * ERFRemoraMultiBlockContainer is the initial-architecture "conductor" object.
 *
 * Design constraints for the initial architecture:
 * 1. The driver (not ERF/REMORA) owns the interface memory and coupling clock.
 * 2. Interface slabs are native AMReX containers (`Vector<MultiFab*>`), matching
 *    the `erf-amrwind-driver` coupling style instead of introducing wrapper types.
 * 3. Geometry is intentionally restricted to level-0, matched horizontal grids so
 *    we can debug timing and data ownership before adding conservative remapping.
 *    Generalized mismatched-grid remapping is intentionally deferred to a future
 *    spatial-remapping milestone.
 *    We deliberately do not claim or target a specific external interfacing or
 *    remapping scheme at this stage, keeping the architecture flexible.
 *
 * Parallel exchange roadmap:
 * - Initial architecture uses direct local access to shared `MultiFab` objects.
 * - Future distributed-grid updates will move data with AMReX-native `ParallelCopy` plus
 *   `NonLocalBC::MultiBlockCommMetaData` once block/rank ownership diverges.
 */
class ERFRemoraMultiBlockContainer
{
public:
    ERFRemoraMultiBlockContainer(int n_cell,
                                 int max_grid_size,
                                 amrex::Real dt_atm,
                                 int atm_substeps_per_ocean_step,
                                 coupling::AtmosToOceanContract atm2ocn_contract);
    ~ERFRemoraMultiBlockContainer() = default;

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
    coupling::AtmosToOceanContract m_atmos_to_ocean_contract;

    // Atmosphere -> ocean cache slab vector.
    // Initial use is legacy COAWST Block-B state-passing; extra layers are
    // preallocated for optional direct flux-passing so contracts can swap by
    // active index sets without changing ownership/allocation strategy.
    // TODO(initial-temporal-interpolation): add old/new temporal copies for interpolation
    // during fast ERF subcycles
    // (same pattern used by erf-amrwind-driver in boundary-state read paths).
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_atmos_to_ocean_storage;
    amrex::Vector<amrex::MultiFab*> m_atmos_to_ocean;
    amrex::Vector<amrex::MultiFab*> m_atmos_to_ocean_active;
    amrex::Vector<int> m_atmos_to_ocean_active_indices;
    amrex::Vector<std::string> m_atmos_to_ocean_active_names;

    // Ocean -> atmosphere state cache (currently SST only).
    // TODO(initial-temporal-interpolation): same old/new cache split as m_atmos_to_ocean.
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_ocean_to_atmos_storage;
    amrex::Vector<amrex::MultiFab*> m_ocean_to_atmos;
};
} // namespace erf_remora

#endif
