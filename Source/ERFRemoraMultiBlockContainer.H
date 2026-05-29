#ifndef ERF_REMORA_DRIVER_MULTIBLOCK_CONTAINER_H_
#define ERF_REMORA_DRIVER_MULTIBLOCK_CONTAINER_H_

#include "ERFRemoraCouplingContract.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_RealBox.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <array>
#include <functional>
#include <memory>

namespace erf_remora
{
/**
 * ERFRemoraMultiBlockContainer is the driver-owned "conductor" object.
 *
 * Stage-2 design:
 * 1. The driver (not ERF/REMORA) owns the full topology: atmosphere RealBox + n_cell and
 *    ocean RealBox + n_cell are stored here and exposed via getters.
 * 2. The interface slab (cache MultiFabs) is derived from the atmosphere horizontal extent
 *    at z = 0.
 * 3. ERF and REMORA are constructed by caller translation units using GetAtm()/GetOcn() getters
 *    to pass explicit topology into their explicit-construction constructors.
 * 4. Interface memory uses native AMReX containers (Vector<MultiFab*>).
 * 5. Geometry is level-0 only; coupling slabs may remap across integer-ratio
 *    horizontal resolution mismatches where supported by the exchange path.
 *
 * Parallel exchange roadmap:
 * - Current: direct local access to shared MultiFab objects.
 * - Future: AMReX-native ParallelCopy + NonLocalBC::MultiBlockCommMetaData once
 *   block/rank ownership diverges.
 */
class ERFRemoraMultiBlockContainer
{
public:
    ERFRemoraMultiBlockContainer(const amrex::RealBox& atm_rb,
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
                                  coupling::AtmosToOceanContract atm2ocn_contract);
    ~ERFRemoraMultiBlockContainer();

    void InitializeBlocks();
    void InitializeSolvers();
    void PrefillAtmosToOceanFallbacks();
    void AdvanceBlocks(int ocean_steps);
    void AdvanceBlocksSynthetic(int ocean_steps);
    amrex::Vector<amrex::MultiFab*>& AtmosToOceanActiveCache();
    const amrex::Vector<std::string>& AtmosToOceanActiveNames() const;

    // Topology getters — used by RunERFOneStepAndPack / RunREMORAApplyAndOneStep
    // to construct solvers with driver-owned explicit topology.
    const amrex::RealBox& GetAtmRealBox() const { return m_atm_rb; }
    const amrex::Vector<int>& GetAtmNCell() const { return m_atm_n_cell; }
    const amrex::RealBox& GetOcnRealBox() const { return m_ocn_rb; }
    const amrex::Vector<int>& GetOcnNCell() const { return m_ocn_n_cell; }
    int GetMaxLevel() const { return m_max_level; }
    int GetCoord() const { return m_coord; }
    const amrex::Vector<amrex::IntVect>& GetRefRatio() const { return m_ref_ratio; }
    const amrex::Array<int,AMREX_SPACEDIM>& GetIsPeriodic() const { return m_is_periodic; }

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
    void FillOceanToAtmosGhosts();
    void PrintDiagnostics(const char* label,
                          const amrex::Vector<amrex::MultiFab*>& fields,
                          const amrex::Vector<std::string>& names) const;

    // Driver-owned topology
    amrex::RealBox m_atm_rb;
    amrex::Vector<int> m_atm_n_cell;
    amrex::RealBox m_ocn_rb;
    amrex::Vector<int> m_ocn_n_cell;
    int m_max_level;
    int m_coord;
    amrex::Vector<amrex::IntVect> m_ref_ratio;
    amrex::Array<int,AMREX_SPACEDIM> m_is_periodic;

    // Interface slab (z=0 plane, derived from atm horizontal extent)
    amrex::BoxArray m_ba;
    amrex::DistributionMapping m_dm;
    amrex::Geometry m_geom;

    // Ocean-resolution interface slab; may differ from m_ba when grids are non-matching.
    amrex::BoxArray m_ba_ocn;
    amrex::DistributionMapping m_dm_ocn;

    amrex::Real m_dt_atm;
    int m_atm_substeps_per_ocean_step;
    long long m_atm_tick = 0;
    coupling::AtmosToOceanContract m_atmos_to_ocean_contract;

    // Atmosphere -> ocean cache slab vector.
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

    void* m_erf = nullptr;
    void* m_remora = nullptr;
};
} // namespace erf_remora

#endif
