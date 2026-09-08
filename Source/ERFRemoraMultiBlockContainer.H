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
    using MutableStateFn = std::function<void(amrex::Vector<amrex::MultiFab*>&,
                                              amrex::Real, long long)>;
    using ConstStateFn = std::function<void(const amrex::Vector<amrex::MultiFab*>&,
                                            amrex::Real, long long)>;

    void AdvanceBlocksWithCallbacks(const MutableStateFn& pack_surface_state,
                                    const ConstStateFn& apply_ocean_surface_state,
                                    const MutableStateFn& pack_atmos_states,
                                    const ConstStateFn& apply_atmos_states,
                                    int ocean_steps);

    void AllocateCaches();
    void AllocateCachesFromSolverLayouts();
    bool AtmosToOceanAveragingActive() const;
    bool TimeAverageAtmosToOceanLane(int global_idx) const;
    void ResetAtmosToOceanAccumulators();
    void AccumulateAtmosToOcean(amrex::Real weight);
    void DrainAtmosToOceanAccumulators();
    void FillOceanToAtmosGhosts();
    void FillAtmosToOceanGhosts();
    void PrintDiagnostics(const char* label,
                          const amrex::Vector<amrex::MultiFab*>& fields,
                          const amrex::Vector<std::string>& names,
                          bool ocean_to_atmos) const;
    void ValidateFields(const char* label,
                        const amrex::Vector<amrex::MultiFab*>& fields,
                        const amrex::Vector<std::string>& names,
                        amrex::Real time,
                        long long ocean_step,
                        bool ocean_to_atmos) const;
    // The two exchange directions are sampled at different instants; each is
    // stamped with the time it actually represents.
    void WriteExchangePlotfiles(amrex::Real o2a_time,
                                amrex::Real a2o_time,
                                amrex::Real dt_ocn,
                                long long ocean_step) const;
    void WriteExchangePlotfile(const std::string& prefix,
                               const amrex::Vector<amrex::MultiFab*>& fields,
                               const amrex::Vector<std::string>& names,
                               const amrex::Geometry& geom,
                               amrex::Real time,
                               amrex::Real ocean_time,
                               amrex::Real dt_ocn,
                               long long ocean_step,
                               bool ocean_to_atmos) const;
    // A plotfile carries one scalar time, but under partial averaging its lanes do
    // not all represent the same instant. Rather than force one number to mean two
    // things, write the per-lane provenance beside the data, the way job_info carries
    // run provenance. Debug metadata only; nothing downstream parses it.
    void WriteExchangeInfo(const std::string& plotfilename,
                           const amrex::Vector<amrex::MultiFab*>& fields,
                           const amrex::Vector<std::string>& names,
                           amrex::Real time,
                           amrex::Real ocean_time,
                           amrex::Real dt_ocn,
                           long long ocean_step,
                           bool ocean_to_atmos) const;
    void ValidateCouplingLayoutsFromSolvers() const;
    coupling::AtmosToOceanContract ActiveAtmosToOceanContract() const;
    void SelectAtmosToOceanActiveView();

    // Decide once who owns each atmos->ocean lane (see LaneSource in
    // ERFRemoraCouplingContract.H) and reconfigure REMORA for the lanes ERF
    // takes over. Must run between ERF InitData, which allocates the radiation
    // fields the capability query reads, and REMORA InitData, which opens
    // NetCDF handles off the selectors this may change.
    void ResolveAtmosToOceanLaneSources();
    void PrintAtmosToOceanLaneProvenance() const;

    // Weight calculation and stencil safety helpers
    void BuildRemapWeights();
    // wet_mask is REMORA's land/sea mask on this family's stagger (1 = wet), used
    // to report coverage over the cells that actually carry ocean. Pass nullptr
    // when unavailable; coverage is then reported over all cells only.
    //
    // donor_fill is that family's per-cell donor flag (1 = this cell's stencil was
    // grafted by the coawst donor search rather than clipped from a real overlap).
    // A donor-filled cell has sum_w > 0 and would otherwise be indistinguishable
    // from a genuinely covered one, so it is reported as its own category. Pass
    // nullptr for a direction that has no donor search (O2A).
    void NormalizeWeightsAndCheckConservation(amrex::MultiFab& weight_mf,
                                              const amrex::iMultiFab& index_mf,
                                              const amrex::MultiFab* wet_mask,
                                              const amrex::iMultiFab* donor_fill,
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
    bool m_time_average_atm2ocn_state = false;
    bool m_time_average_atm2ocn_flux = false;

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
    // Farthest a coawst-mode donor may sit from the destination cell it fills,
    // as a multiple of that donor's own control-volume diagonal, in the same
    // cos(latitude)-corrected metric the search uses. Dimensionless so it
    // survives a change of grid resolution; the search metric is corrected
    // degrees, which no user can pick a threshold in.
    //
    // This bound has NO precedent in any production coupler and is a deliberate
    // departure, not borrowed practice: ESMF's NEAREST_STOD is not
    // distance-bounded at all (when no local candidate exists it queries its
    // process directory with an infinite box, so it always finds the globally
    // nearest source however far away), and creep fill bounds by iteration
    // count rather than by distance. 3.0 is a judgement call carried over from
    // issue #6 and has not been calibrated against a production grid pair.
    // Inert unless the mode is coawst and the frame is spherical.
    amrex::Real m_remap_donor_max_diagonals = amrex::Real(3.0);
    bool m_remap_donor_max_diagonals_set = false;

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
    // Valid-region-only trapezoidal sums. They are reset, filled, and drained
    // within one ocean macro-step and are never checkpoint state.
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_atmos_to_ocean_accumulator_storage;
    amrex::Real m_a2o_accum_weight = 0.0;
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
    // m_atmos_to_ocean_active with every lane the driver does not own nulled, so
    // REMORA keeps its own forcing there. Only this view reaches REMORA;
    // diagnostics use the unmasked one.
    amrex::Vector<amrex::MultiFab*> m_atmos_to_ocean_apply_view;

    // Per-lane source and its reported reason, indexed by AtmosToOceanIndex.
    // Written once by ResolveAtmosToOceanLaneSources(); the per-step view
    // selection only reads them.
    amrex::Vector<coupling::LaneSource> m_atmos_to_ocean_lane_source;
    amrex::Vector<std::string> m_atmos_to_ocean_lane_reason;
    bool m_atmos_to_ocean_lane_sources_resolved = false;

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
    // Per-cell donor flag for each A2O family: 1 where the coawst donor search
    // grafted a full-weight stencil pointing at the nearest non-land source cell,
    // 0 everywhere else. Always 0 under driver.remap_donor_fill=withhold and in
    // the cartesian frame.
    //
    // This is a flag array rather than a counter accumulated in the weight-build
    // loop because faces on a box seam belong to two boxes, so a build-loop
    // counter would double-count seam faces for TAUX_FACE/TAUY_FACE -- and would
    // double-count them a different number of times under a different
    // decomposition, which is a rank-dependent diagnostic inside a change whose
    // whole subject is rank dependence. Carrying the flag on the destination
    // layout lets NormalizeWeightsAndCheckConservation reuse the seam-ownership
    // rule it already applies to every other coverage count.
    //
    // It cannot be folded into m_index_mf_by_family: every component there is a
    // stencil coordinate that ApplyConservativeRemap dereferences.
    //
    // This is also the A2O counterpart of m_valid_o2a_mf, which is the array the
    // O2A direction hands ERF so ERF can tell an uncovered zero from a real ocean
    // value. Nothing yet hands this one to REMORA, so REMORA still cannot tell a
    // coupled 101325 Pa from a placeholder 101325 Pa; closing that asymmetry is a
    // separate change.
    std::array<std::unique_ptr<amrex::iMultiFab>, kNumWeightFamilies> m_donor_fill_mf_by_family;

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
