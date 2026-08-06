#ifndef ERF_REMORA_DRIVER_MULTIBLOCK_CONTAINER_H_
#define ERF_REMORA_DRIVER_MULTIBLOCK_CONTAINER_H_

#include "ERFRemoraCouplingContract.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_RealBox.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <array>
#include <functional>
#include <memory>

namespace erf_remora
{
// COAWST-style weight family tag for the A2O remap stencils. The storage key
// is one slot finer than the conceptual three-way family split (TauX/TauY
// need distinct destination BoxArrays even though both are
// STRESS_STAG_MATCHED), so this enum doubles as both the classification tag
// and the m_weight_mf_by_family/m_index_mf_by_family array index.
enum class WeightFamily {
    SCALAR_CENTERED,      // rho-point destination: IntVect(0,0,0)
    TAUX_FACE,             // STRESS_STAG_MATCHED / TauX: u-face destination: IntVect(1,0,0)
    TAUY_FACE,             // STRESS_STAG_MATCHED / TauY: v-face destination: IntVect(0,1,0)
    WIND_FACE_TO_CENTER,  // Uwind/Vwind: face-staggered source -> rho-point destination
    NumWeightFamilies
};

WeightFamily ClassifyAtmosToOceanFieldFamily(int field_idx);

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
                                  const amrex::IntVect& max_grid_size,
                                  amrex::Real erf_dt,
                                  int atm_substeps_per_ocean_step,
                                  coupling::AtmosToOceanContract atm2ocn_contract,
                                  bool use_two_way_coupling);
    ~ERFRemoraMultiBlockContainer();

    void InitializeBlocks();
    void InitializeSolvers();
    void PrefillAtmosToOceanFallbacks();
    void SetAtmosphereTick(long long atm_tick);
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
    void AllocateCachesFromSolverLayouts();
    void AllocateSyntheticCaches();
    void FillSyntheticData();
    void FillOceanToAtmosGhosts();
    void FillAtmosToOceanGhosts();
    void PrintDiagnostics(const char* label,
                          const amrex::Vector<amrex::MultiFab*>& fields,
                          const amrex::Vector<std::string>& names) const;
    void ValidateFields(const char* label,
                        const amrex::Vector<amrex::MultiFab*>& fields,
                        const amrex::Vector<std::string>& names,
                        amrex::Real time,
                        long long ocean_step) const;
    void WriteExchangePlotfiles(amrex::Real time, long long ocean_step) const;
    void WriteExchangePlotfile(const std::string& prefix,
                               const amrex::Vector<amrex::MultiFab*>& fields,
                               const amrex::Vector<std::string>& names,
                               const amrex::Geometry& geom,
                               amrex::Real time,
                               long long ocean_step) const;
    void ValidateCouplingLayoutsFromSolvers() const;
    coupling::AtmosToOceanContract ActiveAtmosToOceanContract(long long ocean_step) const;
    void SelectAtmosToOceanActiveView(long long ocean_step);
    void BuildRemapWeights();

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
    amrex::Geometry m_geom_ocn;

    amrex::Real m_erf_dt;
    int m_atm_substeps_per_ocean_step;
    long long m_atm_tick = 0;
    int m_plot_int_exchange = -1;
    std::string m_plot_file_exchange_o2a = "exchange_o2a";
    std::string m_plot_file_exchange_a2o = "exchange_a2o";
    coupling::AtmosToOceanContract m_atmos_to_ocean_contract;
    bool m_use_two_way_coupling = false;
    int m_driver_verbose = 0;

    // Atmosphere -> ocean cache slab vector.
    // TODO(initial-temporal-interpolation): add old/new temporal copies for interpolation
    // during fast ERF subcycles
    // (same pattern used by erf-amrwind-driver in boundary-state read paths).
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_atmos_to_ocean_storage;
    amrex::Vector<amrex::MultiFab*> m_atmos_to_ocean;
    amrex::Vector<amrex::MultiFab*> m_atmos_to_ocean_state_active;
    amrex::Vector<int> m_atmos_to_ocean_state_active_indices;
    amrex::Vector<std::string> m_atmos_to_ocean_state_active_names;
    amrex::Vector<amrex::MultiFab*> m_atmos_to_ocean_flux_active;
    amrex::Vector<int> m_atmos_to_ocean_flux_active_indices;
    amrex::Vector<std::string> m_atmos_to_ocean_flux_active_names;
    amrex::Vector<amrex::MultiFab*> m_atmos_to_ocean_active;
    amrex::Vector<int> m_atmos_to_ocean_active_indices;
    amrex::Vector<std::string> m_atmos_to_ocean_active_names;

    // Ocean -> atmosphere state cache (currently SST only).
    // TODO(initial-temporal-interpolation): same old/new cache split as m_atmos_to_ocean.
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_ocean_to_atmos_storage;
    amrex::Vector<amrex::MultiFab*> m_ocean_to_atmos;

    static constexpr int MAX_STENCIL_SIZE = 9;
    static constexpr int kNumWeightFamilies = static_cast<int>(WeightFamily::NumWeightFamilies);

    // Atmosphere -> Ocean (A2O) stencils, one per WeightFamily slot (lives on
    // REMORA layout, on that family's destination stagger).
    std::array<std::unique_ptr<amrex::MultiFab>, kNumWeightFamilies> m_weight_mf_by_family;
    std::array<std::unique_ptr<amrex::iMultiFab>, kNumWeightFamilies> m_index_mf_by_family;

    // Ocean -> Atmosphere (O2A) stencils (lives on ERF layout). O2A is
    // SST-only/scalar today; kept as a single pair (SCALAR_CENTERED-equivalent)
    // rather than a full per-family array until O2A gains more lanes.
    std::unique_ptr<amrex::MultiFab> m_weight_o2a_mf;
    std::unique_ptr<amrex::iMultiFab> m_index_o2a_mf;

    void* m_erf = nullptr;
    void* m_remora = nullptr;
};
} // namespace erf_remora

#endif