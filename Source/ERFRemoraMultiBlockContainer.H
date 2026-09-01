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
// Point structure for 2D curvilinear and Cartesian cell geometry
struct Point2D {
    amrex::Real x{0.0};
    amrex::Real y{0.0};
};

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

// What happens to a wet destination cell that no source cell overlaps.
//
// Withhold: leave the cell without a stencil. Its index keeps the -1 sentinel,
//   its weights stay zero, it is reported uncovered, and ApplyConservativeRemap
//   assigns that lane's fallback value. This is what production couplers do:
//   ESMF's unmappedaction defaults to ERROR, CMEPS pre-fills and then maps with
//   zeroregion=REGION_SELECT, YAC does not interpolate unmapped points, and CDO
//   pre-fills with a missing value. It is the default here for that reason.
// Coawst: graft a full-weight stencil pointing at the nearest non-land source
//   cell, reproducing COAWST's check_weights (Lib/SCRIP_COAWST/scrip.f:343).
//   SCRIP and OASIS are the only comparable codes that fill at all, and both
//   fill with zero -- OASIS documents that as a deficiency. ESMF refuses to
//   combine extrapolation with conservative regridding outright, because
//   "doing so would result in non-conservative weights", which is exactly the
//   objection to this mode: a full-weight nearest-neighbour row spliced into a
//   fracarea-normalized conservative matrix is not conservative. Opt-in only.
//
// Both modes are inert in the cartesian frame, where the donor search never runs.
enum class DonorFillMode {
    Withhold,
    Coawst
};

// Cached plan for moving ERF-layout (source) surface geometry onto the REMORA
// (destination) layout. Defined in the implementation file so this header stays
// free of AMReX_NonLocalBC.H.
struct RemapStagingPlan;
struct REMORAStagingPlan;

/**
 * ERFRemoraMultiBlockContainer is the driver-owned "conductor" object.
 *
 * Driver Responsibilities:
 * 1. Owns topology for both atmosphere (ERF) and ocean (REMORA).
 * 2. Allocates and manages driver-level exchange cache slabs.
 * 3. Computes curvilinear polygon-polygon intersection remap weights (A2O and O2A).
 * 4. Owns the source-to-destination staging plan (GetRemapStagingPlan()) that the
 *    solver-side remap apply uses to route source data onto the destination layout.
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
    // Ocean-frame (absolute) stop time for REMORA, or a negative value to leave
    // REMORA's stop_time alone. Must be set before InitializeSolvers().
    void SetOceanFrameStopTime(amrex::Real ocean_frame_stop_time);
    void PrefillAtmosToOceanFallbacks();
    void SetAtmosphereTick(long long atm_tick);
    void AdvanceBlocks(int ocean_steps);
    amrex::Vector<amrex::MultiFab*>& AtmosToOceanActiveCache();
    const amrex::Vector<std::string>& AtmosToOceanActiveNames() const;

    // True when remap weights are built in lon/lat degrees rather than the
    // shared projected frame. Resolved during BuildRemapWeights().
    bool RemapIsSpherical() const { return m_remap_spherical; }

    // Topology getters
    const amrex::RealBox& GetAtmRealBox() const { return m_atm_rb; }
    const amrex::Vector<int>& GetAtmNCell() const { return m_atm_n_cell; }
    const amrex::RealBox& GetOcnRealBox() const { return m_ocn_rb; }
    const amrex::Vector<int>& GetOcnNCell() const { return m_ocn_n_cell; }
    int GetMaxLevel() const { return m_max_level; }
    int GetCoord() const { return m_coord; }
    const amrex::Vector<amrex::IntVect>& GetRefRatio() const { return m_ref_ratio; }
    const amrex::Array<int,AMREX_SPACEDIM>& GetIsPeriodic() const { return m_is_periodic; }

    // Accessors for remap weight and index MultiFabs by family slot
    const amrex::MultiFab* GetA2OWeightMultiFab(WeightFamily family) const {
        return m_weight_mf_by_family[static_cast<int>(family)].get();
    }
    const amrex::iMultiFab* GetA2OIndexMultiFab(WeightFamily family) const {
        return m_index_mf_by_family[static_cast<int>(family)].get();
    }
    const amrex::MultiFab* GetO2AWeightMultiFab() const { return m_weight_o2a_mf.get(); }
    const amrex::iMultiFab* GetO2AIndexMultiFab() const { return m_index_o2a_mf.get(); }

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

    // Weight calculation and stencil safety helpers
    void BuildRemapWeights();
    // wet_mask is REMORA's land/sea mask on this family's stagger (1 = wet), used
    // to report coverage over the cells that actually carry ocean. Pass nullptr
    // when unavailable; coverage is then reported over all cells only.
    void NormalizeWeightsAndCheckConservation(amrex::MultiFab& weight_mf,
                                              const amrex::iMultiFab& index_mf,
                                              const amrex::MultiFab* wet_mask,
                                              const char* label);

    // Decide (and report) which coordinate frame the weights are built in.
    // Aborts rather than silently degrading when exactly one side reports
    // lon/lat, because clipping degrees against meters produces zero overlap
    // for every cell and therefore an all-zero, silently unforced coupling.
    void ResolveRemapFrame();

    // Decide (and report) what an uncovered wet destination cell receives.
    // Called after ResolveRemapFrame() because the mode is inert in the
    // cartesian frame and the report says so. Aborts on an unrecognized
    // request, the same way ResolveRemapFrame() does.
    void ResolveDonorFillMode();

    // Spherical-frame ERF geometry, built on the ERF layout: a nodal lon/lat
    // corner mesh that tiles without gaps, plus the cell-centered land mask.
    // Split from the REMORA-side loop so no MFIter touches two
    // DistributionMappings.
    void BuildERFCornerMesh();

    // Geographic overlap -> per-REMORA-box ERF index regions -> cached
    // NonLocalBC comm plan -> staged copies on the REMORA layout.
    void BuildERFStagingPlan();
    // Geographic overlap -> per-ERF-box REMORA index regions -> cached
    // NonLocalBC comm plan -> staged copies on the ERF layout.
    void BuildREMORAStagingPlan();

    void BuildAtmosToOceanWeights();
    void BuildOceanToAtmosWeights();

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
    amrex::Real m_ocean_time_offset = 0.0;
    // Negative means "leave REMORA's stop_time alone"; see SetOceanFrameStopTime.
    amrex::Real m_ocean_frame_stop_time = -1.0;
    int m_atm_substeps_per_ocean_step;
    long long m_atm_tick = 0;
    int m_plot_int_exchange = -1;
    std::string m_plot_file_exchange_o2a = "exchange_o2a";
    std::string m_plot_file_exchange_a2o = "exchange_a2o";
    coupling::AtmosToOceanContract m_atmos_to_ocean_contract;
    bool m_use_two_way_coupling = false;
    int m_driver_verbose = 0;

    // Remap configuration and resolved state.
    //
    // The frame choice mirrors SCRIP's `spherical` switch (COAWST
    // Lib/SCRIP_COAWST/read_roms.f): either both sides are clipped in lon/lat
    // degrees using their true corner arrays, or both are clipped in the shared
    // projected frame implied by the driver RealBoxes. There is no meaningful
    // mixed mode, so the resolved state is one bool. "auto" | "spherical" |
    // "cartesian" is only ever a request, interpreted once in ResolveRemapFrame().
    int m_remap_verbose = 0;
    std::string m_remap_frame_request = "auto";
    bool m_remap_spherical = false;
    // Extend the ERF corner mesh by one ring (SCRIP create_extra_rho_grid).
    // Off while the REMORA interior is the whole exchange region; needed once a
    // smaller ocean grid nested inside the atmosphere exchanges on its outer ring.
    bool m_remap_extend_corner_ring = false;
    // Target ERF cells per search bin (SCRIP restrict_type=latlon analogue).
    int m_remap_bin_cells = 4;
    // Least ocean area fraction that counts as covering an ERF cell on the
    // ocean->atmos direction. Below this the cell is reported uncovered and ERF
    // keeps its own lower-boundary SST, rather than accepting a value built from
    // a sliver of one ocean cell. Matches CESM's 0.001 destination-fraction cut
    // (OASIS/CERFACS use 1/1000).
    amrex::Real m_remap_min_coverage_frac = amrex::Real(0.001);
    // What a wet destination cell with no source overlap receives. "withhold" |
    // "coawst" is only ever a request, interpreted once in
    // ResolveDonorFillMode(); see DonorFillMode for why withhold is the default.
    std::string m_remap_donor_fill_request = "withhold";
    DonorFillMode m_remap_donor_fill_mode = DonorFillMode::Withhold;

    // Spherical-frame ERF geometry on the ERF layout: nodal corner mesh
    // (comp 0 = lon, comp 1 = lat) and the cell-centered land mask.
    std::unique_ptr<amrex::MultiFab> m_erf_node_coords;
    std::unique_ptr<amrex::MultiFab> m_erf_land_mask;
    // Same data staged onto the REMORA layout, with the comm plan that put it
    // there. Reused every remap apply, not just at weight-build time.
    std::unique_ptr<RemapStagingPlan> m_staging_plan;
    // REMORA psi-corner coordinates and rho wet mask staged onto the ERF layout
    // for spherical ocean->atmosphere SST remapping.
    std::unique_ptr<REMORAStagingPlan> m_remora_staging_plan;

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
    // Same lanes as m_atmos_to_ocean_active, but with any lane ERF cannot fill set
    // to nullptr so REMORA falls back to its own forcing for it. Only this view is
    // handed to REMORA; diagnostics keep using the unmasked one.
    amrex::Vector<amrex::MultiFab*> m_atmos_to_ocean_apply_view;
    bool m_warned_no_erf_radiation = false;
    bool m_warned_no_erf_cloud_water = false;

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
    // Per-cell ocean->atmos coverage: 1 where wet REMORA rho-cell source area
    // exceeds the configured coverage threshold and m_index_o2a_mf holds at
    // least one ocean donor for that ERF cell, 0 otherwise. Handed to ERF
    // alongside the SST slab so ERF can tell "the ocean says 0 K here" from
    // "the wet ocean does not reach here" -- the remapped slab is zero-filled
    // in both cases.
    std::unique_ptr<amrex::iMultiFab> m_valid_o2a_mf;

    void* m_erf = nullptr;
    void* m_remora = nullptr;
};
} // namespace erf_remora

#endif
