#include "ERFRemoraMultiBlockContainer.H"
#include "StandaloneSolverRuns.H"

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Math.H>
#include <AMReX_NonLocalBC.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_String.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace erf_remora
{
namespace {

// Compute 2D polygon area using Shoelace formula
amrex::Real PolygonArea (const amrex::Vector<Point2D>& poly) {
    if (poly.size() < 3) return 0.0;
    amrex::Real area = 0.0;
    int n = static_cast<int>(poly.size());
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        area += poly[i].x * poly[j].y - poly[j].x * poly[i].y;
    }
    return 0.5 * std::abs(area);
}

// Ensure polygon vertices are ordered counter-clockwise (CCW)
void EnsureCounterClockwise(amrex::Vector<Point2D>& poly) {
    if (poly.size() < 3) return;
    amrex::Real signed_area = 0.0;
    int n = static_cast<int>(poly.size());
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        signed_area += poly[i].x * poly[j].y - poly[j].x * poly[i].y;
    }
    if (signed_area < 0.0) {
        std::reverse(poly.begin(), poly.end());
    }
}

// Check if point p is on the inside (left) of directed edge A -> B
bool IsInsideEdge(const Point2D& p, const Point2D& A, const Point2D& B) {
    constexpr amrex::Real eps = 1.0e-12;
    return ((B.x - A.x) * (p.y - A.y) - (B.y - A.y) * (p.x - A.x)) >= -eps;
}

// Compute line segment intersection between p1-p2 and edge A-B
Point2D EdgeIntersection(const Point2D& p1, const Point2D& p2, const Point2D& A, const Point2D& B) {
    amrex::Real A1 = p2.y - p1.y;
    amrex::Real B1 = p1.x - p2.x;
    amrex::Real C1 = A1 * p1.x + B1 * p1.y;

    amrex::Real A2 = B.y - A.y;
    amrex::Real B2 = A.x - B.x;
    amrex::Real C2 = A2 * A.x + B2 * A.y;

    amrex::Real det = A1 * B2 - A2 * B1;
    if (std::abs(det) < 1.0e-14) {
        return p2;
    }
    return Point2D{ (B2 * C1 - B1 * C2) / det, (A1 * C2 - A2 * C1) / det };
}

// Sutherland-Hodgman convex polygon clipping algorithm
amrex::Vector<Point2D> ClipPolygon(amrex::Vector<Point2D> subjectPoly, amrex::Vector<Point2D> clipPoly) {
    if (subjectPoly.empty() || clipPoly.empty()) return {};

    EnsureCounterClockwise(subjectPoly);
    EnsureCounterClockwise(clipPoly);

    amrex::Vector<Point2D> outputPoly = subjectPoly;
    int numClipEdges = static_cast<int>(clipPoly.size());

    for (int i = 0; i < numClipEdges; ++i) {
        Point2D A = clipPoly[i];
        Point2D B = clipPoly[(i + 1) % numClipEdges];

        amrex::Vector<Point2D> inputPoly = outputPoly;
        outputPoly.clear();
        if (inputPoly.empty()) break;

        Point2D S = inputPoly.back();
        for (const auto& E : inputPoly) {
            if (IsInsideEdge(E, A, B)) {
                if (IsInsideEdge(S, A, B)) {
                    outputPoly.push_back(E);
                } else {
                    outputPoly.push_back(EdgeIntersection(S, E, A, B));
                    outputPoly.push_back(E);
                }
            } else if (IsInsideEdge(S, A, B)) {
                outputPoly.push_back(EdgeIntersection(S, E, A, B));
            }
            S = E;
        }
    }
    return outputPoly;
}

// Extract 4-corner physical quad for cell (i, j) from 2D coordinate MultiFabs
amrex::Vector<Point2D> ExtractCellQuadFromCoords(int i, int j, const amrex::Array4<const amrex::Real>& x_arr,
                                                 const amrex::Array4<const amrex::Real>& y_arr) {
    amrex::Vector<Point2D> quad = {
        { x_arr(i,   j,   0), y_arr(i,   j,   0) },
        { x_arr(i+1, j,   0), y_arr(i+1, j,   0) },
        { x_arr(i+1, j+1, 0), y_arr(i+1, j+1, 0) },
        { x_arr(i,   j+1, 0), y_arr(i,   j+1, 0) }
    };
    EnsureCounterClockwise(quad);
    return quad;
}

void
DriverPrintHeader (int driver_verbose,
                   const char* title,
                   const char* detail = nullptr)
{
    if (driver_verbose < 2 || !amrex::ParallelDescriptor::IOProcessor()) { return; }
    amrex::Print() << "========================================================\n"
                   << "[driver] " << title << "\n";
    if (detail != nullptr && detail[0] != '\0') {
        amrex::Print() << detail << "\n";
    }
    amrex::Print() << "========================================================\n";
}

void
DriverPrintLine (int driver_verbose, const std::string& line)
{
    if (driver_verbose < 2 || !amrex::ParallelDescriptor::IOProcessor()) { return; }
    amrex::Print() << "[driver] " << line << "\n";
}

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
        const amrex::IntVect ixtype = active_view[idx]->boxArray().ixType().toIntVect();
        if (idx == 0 && ixtype != amrex::IntVect(1,0,0)) {
            amrex::Abort("Flux contract tau_x lane must be x-face centered.");
        }
        if (idx == 1 && ixtype != amrex::IntVect(0,1,0)) {
            amrex::Abort("Flux contract tau_y lane must be y-face centered.");
        }
        if (idx >= 2 && ixtype != amrex::IntVect::TheCellVector()) {
            amrex::Abort("Flux contract scalar flux lanes must remain cell-centered.");
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

amrex::BoxArray
AtmosToOceanBoxArray (const amrex::BoxArray& ba_ocn, int idx)
{
    if (idx == coupling::TauX) {
        return amrex::convert(ba_ocn, amrex::IntVect(1,0,0));
    }
    if (idx == coupling::TauY) {
        return amrex::convert(ba_ocn, amrex::IntVect(0,1,0));
    }
    return ba_ocn;
}

// Field -> expected location, taken verbatim from AtmosToOceanBoxArray() above
// (the trusted development-baseline stagger contract):
//
//   field                        expected location (development)
//   ---------------------------  --------------------------------
//   TauX                         x-face staggered: IntVect(1,0,0)
//   TauY                         y-face staggered: IntVect(0,1,0)
//   Uwind, Vwind, Patm, Qair/RH, cell-centered: IntVect(0,0,0)
//   Tair, Cloud, Rain, SWrad,
//   LWrad, SHflux, LHflux, Evap
//   SST (ocean->atmos)           cell-centered: IntVect(0,0,0)
amrex::IntVect
ExpectedAtmosToOceanFieldStagger (int idx)
{
    if (idx == coupling::TauX) { return amrex::IntVect(1,0,0); }
    if (idx == coupling::TauY) { return amrex::IntVect(0,1,0); }
    return amrex::IntVect(0,0,0);
}

amrex::IntVect
AtmosToOceanNGrow (int idx)
{
    // Flux ingestion in REMORA reads grown tiles around the valid region, so the
    // direct-flux lanes need ghost cells that match REMORA's current 2D stencil width.
    static constexpr amrex::IntVect flux_ngrow(2, 2, 0);
    if (idx >= coupling::TauX) {
        return flux_ngrow;
    }
    return amrex::IntVect(0, 0, 0);
}

void
CopyFieldToCellCenteredPlot (amrex::MultiFab& plot_mf, int comp, const amrex::MultiFab& field)
{
    using namespace amrex;

    const IntVect ixtype = field.boxArray().ixType().toIntVect();
    if (ixtype == IntVect::TheCellVector()) {
        MultiFab::Copy(plot_mf, field, 0, comp, 1, 0);
        return;
    }

    for (MFIter mfi(plot_mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box bx = mfi.validbox();
        bx.makeSlab(2,0);
        auto const& src = field.const_array(mfi);
        auto dst = plot_mf.array(mfi, comp);

        if (ixtype == IntVect(1,0,0)) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                dst(i,j,k) = amrex::Real(0.5) * (src(i,j,k) + src(i+1,j,k));
            });
        } else if (ixtype == IntVect(0,1,0)) {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                dst(i,j,k) = amrex::Real(0.5) * (src(i,j,k) + src(i,j+1,k));
            });
        } else {
            amrex::Abort("Unsupported non-cell-centered field in exchange plotfile writer.");
        }
    }
}
}

// Staging plan for ERF-to-REMORA data routing. Each REMORA box gets the ERF
// index region its geographic footprint needs, and one cached amrex::NonLocalBC
// comm plan moves the data there. Caching matters because the per-step remap
// apply restages every ocean step; MultiBlockCommMetaData construction is
// O(n_remora_boxes * n_erf_boxes) and AMReX's own guidance is to cache it. ERF's
// MultiBlockContainer (ERF<->AMR-Wind) uses the same machinery.
//
// Defined here, ahead of the container destructor, so destroying the owning
// unique_ptr sees a complete type.
struct RemapStagingPlan
{
    // ERF regions per REMORA box, in ERF index space.
    amrex::BoxArray node_ba;   // nodal: corner coordinates
    amrex::BoxArray cell_ba;   // cell-centered: land mask and remapped fields

    // Per REMORA box d, in the active remap frame:
    //   box_bbox_lo/hi[2*d], [2*d+1]  the grown lon/lat bounding box used to
    //                                 select staged ERF boxes.
    //   box_has_erf_overlap[d]        0 when cell_ba[d] is the one-cell corner
    //                                 placeholder rather than a real region.
    amrex::Vector<amrex::Real> box_bbox_lo;
    amrex::Vector<amrex::Real> box_bbox_hi;
    amrex::Vector<char>        box_has_erf_overlap;

    // Staged copies, living on the REMORA DistributionMapping.
    std::unique_ptr<amrex::MultiFab> node_mf;
    std::unique_ptr<amrex::MultiFab> land_mf;

    std::unique_ptr<amrex::NonLocalBC::MultiBlockCommMetaData> node_cmd;
    std::unique_ptr<amrex::NonLocalBC::MultiBlockCommMetaData> cell_cmd;

    amrex::Box erf_node_domain;
    amrex::Box erf_cell_domain;

    // Diagnostics
    long long staged_cells = 0;
    long long erf_domain_cells = 0;
    int max_box_cells = 0;
    int empty_boxes = 0;
};

// Staging plan for REMORA-to-ERF SST data routing. Each ERF box gets the
// REMORA rho index region its geographic footprint needs, and the REMORA
// psi-point lon/lat corners plus rho wet mask are cached on the ERF layout.
//
// This is deliberately separate from RemapStagingPlan so the existing A2O
// staging contract keeps its ERF-source assumptions unchanged.
struct REMORAStagingPlan
{
    // REMORA regions per ERF box, in REMORA index space.
    amrex::BoxArray node_ba;   // nodal: psi-point corner coordinates
    amrex::BoxArray cell_ba;   // cell-centered: rho wet mask and SST source

    // Staged copies, living on the ERF DistributionMapping.
    std::unique_ptr<amrex::MultiFab> node_mf;
    std::unique_ptr<amrex::MultiFab> wet_mf;

    std::unique_ptr<amrex::NonLocalBC::MultiBlockCommMetaData> node_cmd;
    std::unique_ptr<amrex::NonLocalBC::MultiBlockCommMetaData> cell_cmd;

    amrex::Box remora_node_domain;
    amrex::Box remora_cell_domain;

    // Diagnostics
    long long staged_cells = 0;
    long long remora_domain_cells = 0;
    int max_box_cells = 0;
    int empty_boxes = 0;
};

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
    coupling::AtmosToOceanContract atm2ocn_contract,
    bool use_two_way_coupling)
    : m_atm_rb(atm_rb), m_atm_n_cell(n_cell_atm),
      m_ocn_rb(ocn_rb), m_ocn_n_cell(n_cell_ocn),
      m_max_level(max_level), m_coord(coord),
      m_ref_ratio(ref_ratio), m_is_periodic(is_periodic),
      m_erf_dt(erf_dt),
      m_atm_substeps_per_ocean_step(atm_substeps_per_ocean_step),
      m_atmos_to_ocean_contract(atm2ocn_contract),
      m_use_two_way_coupling(use_two_way_coupling)
{
    // Interface slab: horizontal extent of the atmosphere at z = 0.
    amrex::Box domain(amrex::IntVect(0, 0, 0),
                      amrex::IntVect(n_cell_atm[0]-1, n_cell_atm[1]-1, 0));
    m_ba.define(domain);
    m_ba.maxSize(max_grid_size);
    m_dm.define(m_ba);
    PrintBoxArrayLayout("Driver constructor ocean->atmos slab", m_ba, m_dm);

    // Default initialization for m_ba_ocn and m_dm_ocn.
    m_ba_ocn = m_ba;
    m_dm_ocn.define(m_ba_ocn);
    PrintBoxArrayLayout("Driver constructor atmos->ocean slab", m_ba_ocn, m_dm_ocn);

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
    pp_driver.query("v", m_driver_verbose);
    pp_driver.query("plot_int_exchange", m_plot_int_exchange);
    pp_driver.query("ocean_time_offset", m_ocean_time_offset);
    pp_driver.query("plot_file_exchange_o2a", m_plot_file_exchange_o2a);
    pp_driver.query("plot_file_exchange_a2o", m_plot_file_exchange_a2o);
    pp_driver.query("remap_v", m_remap_verbose);
    pp_driver.query("remap_frame", m_remap_frame_request);
    pp_driver.query("remap_extend_corner_ring", m_remap_extend_corner_ring);
    pp_driver.query("remap_bin_cells", m_remap_bin_cells);
    pp_driver.query("remap_min_coverage_frac", m_remap_min_coverage_frac);
    pp_driver.query("remap_donor_fill", m_remap_donor_fill_request);
    pp_driver.query("remap_donor_max_diagonals", m_remap_donor_max_diagonals);
    // contains() rather than the query() return, so "set explicitly to the
    // default" is distinguishable from "not set" when reporting it as inert.
    m_remap_donor_max_diagonals_set = pp_driver.contains("remap_donor_max_diagonals");

    DriverPrintHeader(
        m_driver_verbose,
        "Construct ERFRemoraMultiBlockContainer",
        "Driver topology parsed; base slab geometry and contract metadata are being recorded.");
    DriverPrintLine(
        m_driver_verbose,
        std::string("configured_atm2ocn_mode=") +
            coupling::AtmosToOceanContractName(m_atmos_to_ocean_contract) +
            " two_way=" + (m_use_two_way_coupling ? "1" : "0"));
    DriverPrintLine(
        m_driver_verbose,
        std::string("ocean_time_offset=") + std::to_string(m_ocean_time_offset));
}

void
ERFRemoraMultiBlockContainer::InitializeSolvers ()
{
    DriverPrintHeader(
        m_driver_verbose,
        "Initialize Solvers",
        "Phase 1: construct ERF and REMORA objects. Phase 2: inject driver atm2ocn metadata before InitData. Phase 3: call InitData on both solvers.");

    if (m_erf != nullptr) {
        DriverPrintLine(m_driver_verbose, "destroying previous ERF solver handle");
        DestroyERFSolver(m_erf);
        m_erf = nullptr;
    }
    if (m_remora != nullptr) {
        DriverPrintLine(m_driver_verbose, "destroying previous REMORA solver handle");
        DestroyREMORASolver(m_remora);
        m_remora = nullptr;
    }

    DriverPrintLine(m_driver_verbose, "constructing ERF solver from driver-owned atmosphere topology");
    m_erf = CreateERFSolver(*this);
    DriverPrintLine(m_driver_verbose, "ERF constructor complete");

    // Re-key stop_time for REMORA, after ERF has read it and before REMORA is built.
    //
    // ERF and REMORA both read the unprefixed stop_time, and the deck expresses it as
    // elapsed-from-zero. REMORA, however, runs on the absolute ocean clock once
    // driver.ocean_time_offset is applied, so an elapsed value makes
    // REMORA::EvolveOneStep return dt=0 on the first macro step.
    //
    // The placement matters. ERF captured stop_time into its own state during
    // CreateERFSolver above, so removing the key now cannot un-read it and ERF keeps
    // the dt clamp that makes its last step land on stop_time. Removing it earlier --
    // in main.cpp, before either solver exists -- silently changed ERF's dt and its
    // checkpoint metadata under offset decks.
    //
    // REMORA aborts if both stop_time and remora.stop_time are set
    // (REMORA::ReadParameters), so the unprefixed key must go, not merely be shadowed.
    if (m_ocean_frame_stop_time >= amrex::Real(0.0)) {
        // No warning about a user-supplied remora.stop_time here: main.cpp scrubs that
        // key long before this point, so querying it now would always miss. The warning
        // lives at the scrub site, which is the only place it can still be observed.
        amrex::ParmParse pp_top;
        pp_top.remove("stop_time");
        pp_top.remove("remora.stop_time");
        pp_top.add("remora.stop_time", m_ocean_frame_stop_time);
        DriverPrintLine(m_driver_verbose,
                        std::string("re-keyed REMORA stop_time to ocean frame: ") +
                            std::to_string(m_ocean_frame_stop_time));
    }

    DriverPrintLine(m_driver_verbose, "constructing REMORA solver from driver-owned ocean topology");
    m_remora = CreateREMORASolver(*this);
    DriverPrintLine(m_driver_verbose, "REMORA constructor complete");

    const bool bootstrap_uses_state_contract =
        (ActiveAtmosToOceanContract(/*ocean_step=*/0) == coupling::AtmosToOceanContract::StatePassing);
    DriverPrintLine(
        m_driver_verbose,
        std::string("injecting pre-InitData atm2ocn metadata: bootstrap_contract=") +
            (bootstrap_uses_state_contract ? "state" : "flux") +
            " configured_mode=" + coupling::AtmosToOceanContractName(m_atmos_to_ocean_contract) +
            " two_way=" + (m_use_two_way_coupling ? "1" : "0"));
    ConfigureERFDriverAtmosToOceanCoupling(
        m_erf, true, m_use_two_way_coupling, bootstrap_uses_state_contract);
    ConfigureREMORADriverAtmosToOceanCoupling(
        m_remora, true, m_use_two_way_coupling, bootstrap_uses_state_contract);
    SetREMORAStartTime(m_remora, m_ocean_time_offset);

    DriverPrintLine(m_driver_verbose, "calling ERF InitData (enters ERF InitData_pre -> ERF InitData_post)");
    InitERFSolver(m_erf);
    DriverPrintLine(m_driver_verbose, "ERF InitData complete");
    DriverPrintLine(m_driver_verbose, "calling REMORA InitData");
    InitREMORASolver(m_remora);
    DriverPrintLine(m_driver_verbose, "REMORA InitData complete");
}

void
ERFRemoraMultiBlockContainer::PrefillAtmosToOceanFallbacks ()
{
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
ERFRemoraMultiBlockContainer::SetOceanFrameStopTime (amrex::Real ocean_frame_stop_time)
{
    m_ocean_frame_stop_time = ocean_frame_stop_time;
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
    DriverPrintHeader(
        m_driver_verbose,
        "Initialize Coupling Caches",
        "Phase 4: query finalized post-InitData solver layouts. Phase 5: allocate driver-owned caches from solver layouts. Phase 6: seed bootstrap state and build remap weights and select active atm2ocn view.");
    // Enforce the lifecycle invariant rather than assuming it. Every path below
    // dereferences both solver handles, and the synthetic fallback that used to
    // cover the null case has been removed as unreachable. Assert so a caller who
    // skips InitializeSolvers() gets a clear abort instead of silent bad state.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_erf != nullptr && m_remora != nullptr,
        "InitializeBlocks requires InitializeSolvers to have constructed both solvers");
    m_atm_tick = 0;
    DriverPrintLine(m_driver_verbose, "resetting driver atmosphere tick to 0");
    AllocateCaches();
    BuildRemapWeights();
    DriverPrintLine(m_driver_verbose, "selecting bootstrap active atm2ocn cache view for ocean_step=0");
    SelectAtmosToOceanActiveView(/*ocean_step=*/0);
    if (m_use_two_way_coupling && !m_ocean_to_atmos.empty()) {
        DriverPrintLine(m_driver_verbose, "packing REMORA surface state into ocean->atmos bootstrap cache");
        const amrex::iMultiFab* erf_is_land_bootstrap = nullptr;
        ERFGetLandMask(m_erf, &erf_is_land_bootstrap);
        REMORAPackSurfaceState(m_remora, m_ocean_to_atmos, amrex::Real(0.0),
                               m_weight_o2a_mf.get(), m_index_o2a_mf.get(), MAX_STENCIL_SIZE,
                               erf_is_land_bootstrap);
        FillOceanToAtmosGhosts();
    }

    if (m_atmos_to_ocean_contract == coupling::AtmosToOceanContract::FluxPassing) {
        DriverPrintLine(m_driver_verbose, "validating steady-state flux active view shape and centering");
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
    // Only reachable from InitializeBlocks(), which asserts both solvers exist.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_erf != nullptr && m_remora != nullptr,
        "AllocateCaches requires InitializeSolvers to have constructed both solvers");
    DriverPrintLine(m_driver_verbose, "allocating caches from solver-reported post-InitData layouts");
    AllocateCachesFromSolverLayouts();
}

void
ERFRemoraMultiBlockContainer::AllocateCachesFromSolverLayouts ()
{
    DriverPrintLine(m_driver_verbose, "querying REMORA rho/u/v atm2ocn interface layouts");
    amrex::BoxArray remora_rho_ba;
    amrex::DistributionMapping remora_rho_dm;
    amrex::BoxArray remora_u_ba;
    amrex::DistributionMapping remora_u_dm;
    amrex::BoxArray remora_v_ba;
    amrex::DistributionMapping remora_v_dm;

    REMORAGetAtmosToOceanRhoLayout(m_remora, remora_rho_ba, remora_rho_dm);
    REMORAGetAtmosToOceanUFaceLayout(m_remora, remora_u_ba, remora_u_dm);
    REMORAGetAtmosToOceanVFaceLayout(m_remora, remora_v_ba, remora_v_dm);

    if (m_use_two_way_coupling) {
        DriverPrintLine(m_driver_verbose, "querying ERF ocean-surface layout for ocean->atmos cache allocation");
        ERFGetOceanToAtmosSurfaceLayout(m_erf, m_ba, m_dm);
    }

    m_ba_ocn = remora_rho_ba;
    m_dm_ocn = remora_rho_dm;

    const amrex::Real dz_atm = (m_atm_rb.hi(2) - m_atm_rb.lo(2)) / m_atm_n_cell[2];
    const amrex::RealBox atm_slab_rb({m_atm_rb.lo(0), m_atm_rb.lo(1), 0.0},
                                     {m_atm_rb.hi(0), m_atm_rb.hi(1), dz_atm});
    m_geom.define(m_ba.minimalBox(), &atm_slab_rb, amrex::CoordSys::cartesian, m_is_periodic.data());

    const amrex::Real dz_ocn = (m_ocn_rb.hi(2) - m_ocn_rb.lo(2)) / m_ocn_n_cell[2];
    const amrex::RealBox ocn_slab_rb({m_ocn_rb.lo(0), m_ocn_rb.lo(1), 0.0},
                                     {m_ocn_rb.hi(0), m_ocn_rb.hi(1), dz_ocn});
    m_geom_ocn.define(m_ba_ocn.minimalBox(), &ocn_slab_rb, amrex::CoordSys::cartesian, m_is_periodic.data());

    DriverPrintLine(m_driver_verbose, "validating solver-reported coupling layouts against driver assumptions");
    ValidateCouplingLayoutsFromSolvers();

    m_atmos_to_ocean_storage.resize(coupling::NumAtmosToOceanLayers);
    for (int idx = 0; idx < static_cast<int>(m_atmos_to_ocean_storage.size()); ++idx) {
        auto& mf_ptr = m_atmos_to_ocean_storage[idx];
        const amrex::BoxArray* lane_ba = &remora_rho_ba;
        const amrex::DistributionMapping* lane_dm = &remora_rho_dm;
        if (idx == coupling::TauX) {
            lane_ba = &remora_u_ba;
            lane_dm = &remora_u_dm;
        } else if (idx == coupling::TauY) {
            lane_ba = &remora_v_ba;
            lane_dm = &remora_v_dm;
        }
        mf_ptr = std::make_unique<amrex::MultiFab>(*lane_ba, *lane_dm, 1, AtmosToOceanNGrow(idx));
        mf_ptr->setVal(0.0);
        const std::string label = "Driver AllocateCachesFromSolverLayouts atmos->ocean[" + std::to_string(idx) + "]";
        PrintMultiFabLayout(label.c_str(), *mf_ptr);
    }
    m_atmos_to_ocean = amrex::GetVecOfPtrs(m_atmos_to_ocean_storage);

    // Atmosphere -> Ocean (A2O) stencils, one per WeightFamily slot, each built
    // against its own destination stagger (rho for scalars, u-face for TauX,
    // v-face for TauY, rho for WIND_FACE_TO_CENTER since the face->center
    // averaging happens upstream in the source pack step).
    {
        const std::array<const amrex::BoxArray*, kNumWeightFamilies> family_ba = {
            &remora_rho_ba, &remora_u_ba, &remora_v_ba, &remora_rho_ba
        };
        const std::array<const amrex::DistributionMapping*, kNumWeightFamilies> family_dm = {
            &remora_rho_dm, &remora_u_dm, &remora_v_dm, &remora_rho_dm
        };
        for (int f = 0; f < kNumWeightFamilies; ++f) {
            m_weight_mf_by_family[f] = std::make_unique<amrex::MultiFab>(
                *family_ba[f], *family_dm[f], MAX_STENCIL_SIZE, 0);
            m_index_mf_by_family[f] = std::make_unique<amrex::iMultiFab>(
                *family_ba[f], *family_dm[f], 3 * MAX_STENCIL_SIZE, 0);
        }
    }

    // Ocean -> Atmosphere (O2A) stencils (lives on ERF layout)
    if (m_use_two_way_coupling) {
        m_weight_o2a_mf = std::make_unique<amrex::MultiFab>(m_ba, m_dm, MAX_STENCIL_SIZE, 0);
        m_index_o2a_mf = std::make_unique<amrex::iMultiFab>(m_ba, m_dm, 3 * MAX_STENCIL_SIZE, 0);
        m_valid_o2a_mf = std::make_unique<amrex::iMultiFab>(m_ba, m_dm, 1, 0);
        m_valid_o2a_mf->setVal(0);
    }

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

    if (m_use_two_way_coupling) {
        m_ocean_to_atmos_storage.resize(coupling::NumOceanToAtmosStates);
        const amrex::IntVect ocean_to_atmos_ng(1, 1, 0);
        for (int idx = 0; idx < static_cast<int>(m_ocean_to_atmos_storage.size()); ++idx) {
            auto& mf_ptr = m_ocean_to_atmos_storage[idx];
            mf_ptr = std::make_unique<amrex::MultiFab>(m_ba, m_dm, 1, ocean_to_atmos_ng);
            mf_ptr->setVal(0.0);
            const std::string label = "Driver AllocateCachesFromSolverLayouts ocean->atmos[" + std::to_string(idx) + "]";
            PrintMultiFabLayout(label.c_str(), *mf_ptr);
        }
    } else {
        m_ocean_to_atmos_storage.clear();
    }
    m_ocean_to_atmos = amrex::GetVecOfPtrs(m_ocean_to_atmos_storage);
}

void
ERFRemoraMultiBlockContainer::ValidateCouplingLayoutsFromSolvers () const
{
    // All integer-ratio assertions and layout validations have been removed
    // to support general non-conformal curvilinear grid coupling.
}

// Weight-family classification per field, per the corrected authoritative
// contract for this branch:
//   TAUX_FACE/TAUY_FACE (TauX/TauY): ERF tau_x/tau_y are transferred directly
//     to REMORA vec_sustr/vec_svstr, which already live on the matching
//     x-face/y-face stagger (see AtmosToOceanBoxArray() above and
//     Submodules/REMORA/Source/REMORA.H:423-426). This is a stagger-matched
//     nodal transfer, not a center->face averaging problem.
//   WIND_FACE_TO_CENTER (Uwind/Vwind): ERF U/V live on x-face/y-face (native
//     ERF velocity stagger) but REMORA consumes wind at rho-points
//     (Submodules/REMORA/Source/REMORA.H:428, "Wind ... defined at rho-points").
//     Compatible only when a face->center averaging transform hook exists;
//     none is wired yet, so this family always falls back today.
//   SCALAR_CENTERED (Patm, Qair, Tair, Cloud, Rain, SWrad, LWrad, SHflux,
//     LHflux, Evap, SST): both sides are cell-centered/rho already.
WeightFamily
ClassifyAtmosToOceanFieldFamily (int field_idx) {
    if (field_idx == coupling::TauX) { return WeightFamily::TAUX_FACE; }
    if (field_idx == coupling::TauY) { return WeightFamily::TAUY_FACE; }
    if (field_idx == coupling::Uwind || field_idx == coupling::Vwind) { return WeightFamily::WIND_FACE_TO_CENTER; }
    return WeightFamily::SCALAR_CENTERED;
}

// Normalize each destination cell's stencil by its own covered fraction. Before
// normalization the weight sum is the fraction of the destination cell covered by
// source cells, which is SCRIP's grid2_frac; dividing by it is SCRIP's
// norm_opt_frcarea, the option COAWST uses (COAWST
// Lib/SCRIP_COAWST/remap_conserv.f: "norm_factor = one/grid2_frac ... we use
// this one"). Uncovered cells keep zero weights, so the coverage statistics below
// are the check that the frame and corner meshes actually overlap: an all-zero
// weight set normalizes silently and would otherwise look like a clean run with
// no forcing.
void
ERFRemoraMultiBlockContainer::NormalizeWeightsAndCheckConservation(amrex::MultiFab& weight_mf,
                                                                   const amrex::iMultiFab& index_mf,
                                                                   const amrex::MultiFab* wet_mask,
                                                                   const char* label)
{
    using namespace amrex;
    amrex::Real max_cover_err = 0.0;
    amrex::Real frac_min = std::numeric_limits<amrex::Real>::max();
    amrex::Real frac_max = 0.0;
    amrex::Real frac_sum = 0.0;
    amrex::Long n_cells = 0;
    amrex::Long n_uncovered = 0;
    amrex::Long n_no_stencil = 0;
    // Same counts restricted to REMORA's wet cells. Coverage over dry cells is
    // uninteresting - a land cell has nothing to force - so this is the count that
    // says whether the remap is actually usable.
    amrex::Long n_wet = 0;
    amrex::Long n_wet_uncovered = 0;
    amrex::Real wet_frac_min = std::numeric_limits<amrex::Real>::max();

    // Faces on a box seam belong to both adjacent boxes, so summing validbox over
    // boxes would count them twice. Give each face to the lower box and let the
    // last box in each direction keep its domain-boundary face, which tiles the
    // staggered index space exactly once. Weights are still computed on every
    // face including shared ones; this only makes the counts truthful.
    const Box stagger_domain = weight_mf.boxArray().minimalBox();
    const IntVect stagger = weight_mf.boxArray().ixType().toIntVect();

    for (MFIter mfi(weight_mf); mfi.isValid(); ++mfi) {
        Box bx = mfi.validbox();
        auto w_arr = weight_mf.array(mfi);
        auto const& idx_arr = index_mf.const_array(mfi);
        const bool has_mask = (wet_mask != nullptr);
        auto const& wet_arr = has_mask ? wet_mask->const_array(mfi) : Array4<const Real>{};

        Box count_bx = bx;
        for (int dir = 0; dir < 2; ++dir) {
            if (stagger[dir] == 1 && bx.bigEnd(dir) < stagger_domain.bigEnd(dir)) {
                count_bx.growHi(dir, -1);
            }
        }

        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                amrex::Real sum_w = 0.0;
                for (int m = 0; m < MAX_STENCIL_SIZE; ++m) {
                    sum_w += w_arr(i, j, 0, m);
                }

                const bool counted = count_bx.contains(IntVect(i,j,0));
                const bool is_wet = has_mask ? (wet_arr(i,j,0) > Real(0.5)) : false;
                if (counted) {
                    ++n_cells;
                    if (idx_arr(i, j, 0, 0) < 0) { ++n_no_stencil; }
                    if (has_mask && is_wet) { ++n_wet; }
                }

                if (sum_w > 0.0) {
                    if (counted) {
                        frac_min = std::min(frac_min, sum_w);
                        frac_max = std::max(frac_max, sum_w);
                        frac_sum += sum_w;
                        max_cover_err = std::max(max_cover_err, std::abs(sum_w - 1.0));
                        if (has_mask && is_wet) { wet_frac_min = std::min(wet_frac_min, sum_w); }
                    }
                    for (int m = 0; m < MAX_STENCIL_SIZE; ++m) {
                        w_arr(i, j, 0, m) /= sum_w;
                    }
                } else if (counted) {
                    ++n_uncovered;
                    if (has_mask && is_wet) { ++n_wet_uncovered; }
                }
            }
        }
    }

    ParallelDescriptor::ReduceRealMax(max_cover_err);
    ParallelDescriptor::ReduceRealMin(frac_min);
    ParallelDescriptor::ReduceRealMax(frac_max);
    ParallelDescriptor::ReduceRealSum(frac_sum);
    ParallelDescriptor::ReduceRealMin(wet_frac_min);
    ParallelDescriptor::ReduceLongSum(n_cells);
    ParallelDescriptor::ReduceLongSum(n_uncovered);
    ParallelDescriptor::ReduceLongSum(n_no_stencil);
    ParallelDescriptor::ReduceLongSum(n_wet);
    ParallelDescriptor::ReduceLongSum(n_wet_uncovered);

    // The fallback in ERF_to_REMORA.cpp:174 branches on the index sentinel;
    // n_wet_uncovered above is a weight test. They agree today only because the
    // weight builder never stores an index with zero weight. Assert it, so a
    // future change that breaks the equivalence surfaces here instead of as
    // 0 Pa / 0 K arriving in REMORA's bulk-flux chain.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        n_uncovered == n_no_stencil,
        "NormalizeWeightsAndCheckConservation: a destination cell has a stored "
        "source index but zero total weight. The ERF-side fallback keys on the "
        "index sentinel, so that cell will receive 0 rather than its lane fallback.");

    const amrex::Long n_covered = n_cells - n_uncovered;
    if (n_covered == 0) {
        amrex::Print() << "[driver] WARNING remap (" << label << "): no destination cell overlaps any "
                       << "source cell. Every weight is zero, so this coupling direction carries no "
                       << "forcing. Check driver.remap_frame and the two coordinate ranges reported "
                       << "above.\n";
        return;
    }

    if (ParallelDescriptor::IOProcessor() &&
        (m_remap_verbose >= 1 || m_driver_verbose >= 1 || n_uncovered > 0)) {
        // max|frac-1| is the coverage error, not a conservation error: before
        // normalization the weight sum *is* the covered fraction. A frac_max
        // slightly above 1 is polygon-clipping round-off and is the useful signal
        // hiding in it.
        amrex::Print() << "[driver] remap coverage (" << label << "): cells=" << n_cells
                       << " uncovered=" << n_uncovered
                       << " no_stencil=" << n_no_stencil
                       << " frac min/mean/max=" << frac_min << "/"
                       << (frac_sum / static_cast<amrex::Real>(n_covered)) << "/" << frac_max
                       << " max|frac-1|=" << max_cover_err << "\n";
    }

    // The verdict: an uncovered wet cell is a real hole in the forcing, whereas an
    // uncovered dry cell is expected wherever the atmosphere grid is land.
    if (wet_mask != nullptr) {
        if (n_wet_uncovered > 0) {
            amrex::Print() << "[driver] WARNING remap (" << label << "): " << n_wet_uncovered
                           << " of " << n_wet << " REMORA wet cells have no ERF source overlap. "
                           << "ApplyConservativeRemap in Submodules/ERF/Source/Coupling/"
                           << "ERF_to_REMORA.cpp branches on the -1 index sentinel and assigns a "
                           << "per-lane fallback: Patm=101325 Pa and Tair=288.15 K are the "
                           << "standard-atmosphere values; Uwind/Vwind=0 m/s; Cloud, SWrad, LWrad, "
                           << "and all eight flux lanes=0.0. Qair and Rain are not remapped and "
                           << "retain the PrefillAtmosToOceanFallbacks values from "
                           << "remora.air_humidity and remora.rain. Thus uncovered Patm/Tair cells "
                           << "are forced but not coupled; for wind and flux lanes they are neither "
                           << "forced nor coupled. Expect these along coastlines where the two land "
                           << "masks disagree, or where the ocean grid reaches past the atmosphere.\n";
        } else if (ParallelDescriptor::IOProcessor() && m_remap_verbose >= 1) {
            amrex::Print() << "[driver] remap coverage (" << label << "): all " << n_wet
                           << " wet cells covered, min wet frac=" << wet_frac_min << "\n";
        }
    }
}

namespace
{
// AMReX's Geometry::CellBounds() is not available on this AMReX version;
// ProbLo()/CellSize() are the public accessors for uniform Cartesian geometry.
amrex::RealBox
CellBoundsFromGeom (const amrex::Geometry& g, int i, int j)
{
    const amrex::Real* plo = g.ProbLo();
    const amrex::Real* dx  = g.CellSize();

    const amrex::Real xlo = plo[0] + static_cast<amrex::Real>(i    ) * dx[0];
    const amrex::Real xhi = plo[0] + static_cast<amrex::Real>(i + 1) * dx[0];
    const amrex::Real ylo = plo[1] + static_cast<amrex::Real>(j    ) * dx[1];
    const amrex::Real yhi = plo[1] + static_cast<amrex::Real>(j + 1) * dx[1];

    return amrex::RealBox({xlo, ylo, 0.0}, {xhi, yhi, 0.0});
}

amrex::Long
FlatIndexInBox (const amrex::Box& domain, int i, int j, int k)
{
    const amrex::Long nx = domain.length(0);
    const amrex::Long ny = domain.length(1);
    return ((static_cast<amrex::Long>(k - domain.smallEnd(2)) * ny
             + static_cast<amrex::Long>(j - domain.smallEnd(1))) * nx
            + static_cast<amrex::Long>(i - domain.smallEnd(0)));
}

amrex::Vector<amrex::Real>
BuildReplicatedMaskLookup (const amrex::MultiFab& mask_mf, const amrex::Box& domain)
{
    using namespace amrex;

    Vector<Real> mask_values(domain.numPts(), Real(0.0));
    for (MFIter mfi(mask_mf); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox() & domain;
        if (!bx.ok()) { continue; }
        auto const& mask_arr = mask_mf.const_array(mfi);

        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    mask_values[FlatIndexInBox(domain, i, j, k)] = mask_arr(i, j, k);
                }
            }
        }
    }
    ParallelDescriptor::ReduceRealSum(mask_values.dataPtr(), mask_values.size());
    return mask_values;
}

// Control-volume bounds for a face-staggered destination point on a uniform grid.
//
// The u-cell at face index i is centered on the face, spanning half a cell either
// side; the v-cell likewise in y. Using CellBoundsFromGeom at a face index instead
// shifts the volume half a cell and, worse, sends the outermost face column
// (i = n_cell, which is a valid face index) entirely outside the domain, where it
// finds no source overlap and ends up with zero weights - so tau_x/tau_y carry no
// forcing along that edge. Clamping to the domain gives an edge face the half-width
// control volume it actually has, which is also what the spherical path produces
// from clamped node midpoints.
amrex::RealBox
FaceCellBoundsFromGeom (const amrex::Geometry& g, const amrex::Box& cell_domain,
                        int i, int j, int dir)
{
    const amrex::Real* plo = g.ProbLo();
    const amrex::Real* dx  = g.CellSize();

    const amrex::Real dom_lo_x = plo[0] + static_cast<amrex::Real>(cell_domain.smallEnd(0)    ) * dx[0];
    const amrex::Real dom_hi_x = plo[0] + static_cast<amrex::Real>(cell_domain.bigEnd(0)   + 1) * dx[0];
    const amrex::Real dom_lo_y = plo[1] + static_cast<amrex::Real>(cell_domain.smallEnd(1)    ) * dx[1];
    const amrex::Real dom_hi_y = plo[1] + static_cast<amrex::Real>(cell_domain.bigEnd(1)   + 1) * dx[1];

    amrex::Real xlo, xhi, ylo, yhi;
    if (dir == 0) {
        xlo = plo[0] + (static_cast<amrex::Real>(i) - amrex::Real(0.5)) * dx[0];
        xhi = plo[0] + (static_cast<amrex::Real>(i) + amrex::Real(0.5)) * dx[0];
        ylo = plo[1] + static_cast<amrex::Real>(j    ) * dx[1];
        yhi = plo[1] + static_cast<amrex::Real>(j + 1) * dx[1];
        xlo = std::max(xlo, dom_lo_x); xhi = std::min(xhi, dom_hi_x);
    } else {
        xlo = plo[0] + static_cast<amrex::Real>(i    ) * dx[0];
        xhi = plo[0] + static_cast<amrex::Real>(i + 1) * dx[0];
        ylo = plo[1] + (static_cast<amrex::Real>(j) - amrex::Real(0.5)) * dx[1];
        yhi = plo[1] + (static_cast<amrex::Real>(j) + amrex::Real(0.5)) * dx[1];
        ylo = std::max(ylo, dom_lo_y); yhi = std::min(yhi, dom_hi_y);
    }

    return amrex::RealBox({xlo, ylo, 0.0}, {xhi, yhi, 0.0});
}

// Axis-aligned bounds of a polygon, used both to reject candidate pairs cheaply
// and to size the search structure.
struct PlanarBBox
{
    amrex::Real xlo{ std::numeric_limits<amrex::Real>::max()};
    amrex::Real xhi{-std::numeric_limits<amrex::Real>::max()};
    amrex::Real ylo{ std::numeric_limits<amrex::Real>::max()};
    amrex::Real yhi{-std::numeric_limits<amrex::Real>::max()};

    void include (amrex::Real x, amrex::Real y) {
        xlo = std::min(xlo, x); xhi = std::max(xhi, x);
        ylo = std::min(ylo, y); yhi = std::max(yhi, y);
    }
    void include (const amrex::Vector<Point2D>& poly) {
        for (const auto& p : poly) { include(p.x, p.y); }
    }
    [[nodiscard]] bool ok () const { return xhi >= xlo && yhi >= ylo; }
    [[nodiscard]] bool intersects (const PlanarBBox& o) const {
        return !(o.xlo > xhi || o.xhi < xlo || o.ylo > yhi || o.yhi < ylo);
    }
};

PlanarBBox
BBoxOfPolygon (const amrex::Vector<Point2D>& poly)
{
    PlanarBBox bb;
    bb.include(poly);
    return bb;
}

// Uniform bin index over the source cells of one staged region, used to restrict
// candidate source cells per destination cell. This is SCRIP's search-bin idea
// (COAWST Lib/SCRIP_COAWST/grids.f builds bin_addr from per-cell bounding
// boxes); amrex::DenseBins is the natural utility for binning but assigns each
// item to exactly one bin, which cannot express a source cell whose bounding box
// straddles bins, so this is a small two-pass CSR with multi-bin membership.
class SourceSearchGrid
{
public:
    void define (const amrex::Box& src_box,
                 const amrex::Array4<const amrex::Real>& xn,
                 const amrex::Array4<const amrex::Real>& yn,
                 int cells_per_bin)
    {
        using namespace amrex;

        m_src_box = src_box;
        const int nx_cells = std::max(1, src_box.length(0));
        const int ny_cells = std::max(1, src_box.length(1));
        m_bbox.resize(static_cast<Long>(nx_cells) * static_cast<Long>(ny_cells));

        // Pass 0: per-cell bounds (SCRIP grid1_bound_box) and the global extent.
        int item = 0;
        for (int j = src_box.smallEnd(1); j <= src_box.bigEnd(1); ++j) {
            for (int i = src_box.smallEnd(0); i <= src_box.bigEnd(0); ++i, ++item) {
                PlanarBBox bb;
                bb.include(xn(i  , j  , 0), yn(i  , j  , 0));
                bb.include(xn(i+1, j  , 0), yn(i+1, j  , 0));
                bb.include(xn(i+1, j+1, 0), yn(i+1, j+1, 0));
                bb.include(xn(i  , j+1, 0), yn(i  , j+1, 0));
                m_bbox[item] = bb;
                m_extent.include(bb.xlo, bb.ylo);
                m_extent.include(bb.xhi, bb.yhi);
            }
        }

        const int per_bin = std::max(1, cells_per_bin);
        m_nx = std::max(1, nx_cells / per_bin);
        m_ny = std::max(1, ny_cells / per_bin);

        // Guard against a degenerate extent (single cell, or all-equal coords).
        const Real wx = std::max(m_extent.xhi - m_extent.xlo, Real(1.0e-30));
        const Real wy = std::max(m_extent.yhi - m_extent.ylo, Real(1.0e-30));
        m_inv_dx = static_cast<Real>(m_nx) / wx;
        m_inv_dy = static_cast<Real>(m_ny) / wy;

        // Pass 1: count per bin. Pass 2: prefix sum and scatter.
        const int nbins = m_nx * m_ny;
        m_offsets.assign(nbins + 1, 0);
        for (int n = 0; n < static_cast<int>(m_bbox.size()); ++n) {
            int ilo, ihi, jlo, jhi;
            binRange(m_bbox[n], ilo, ihi, jlo, jhi);
            for (int bj = jlo; bj <= jhi; ++bj) {
                for (int bi = ilo; bi <= ihi; ++bi) { ++m_offsets[bj * m_nx + bi + 1]; }
            }
        }
        for (int b = 0; b < nbins; ++b) { m_offsets[b+1] += m_offsets[b]; }
        m_items.resize(m_offsets[nbins]);
        Vector<int> cursor(m_offsets.begin(), m_offsets.end() - 1);
        for (int n = 0; n < static_cast<int>(m_bbox.size()); ++n) {
            int ilo, ihi, jlo, jhi;
            binRange(m_bbox[n], ilo, ihi, jlo, jhi);
            for (int bj = jlo; bj <= jhi; ++bj) {
                for (int bi = ilo; bi <= ihi; ++bi) { m_items[cursor[bj * m_nx + bi]++] = n; }
            }
        }

        m_stamp.assign(m_bbox.size(), 0);
        m_query = 0;
    }

    // Collect the source cells whose bounds overlap `target`. Multi-bin
    // membership means an item can be reached from several bins, so a query
    // stamp keeps each candidate exactly once: a duplicate would double-count
    // its overlap area into the destination weight sum.
    void query (const PlanarBBox& target, amrex::Vector<int>& out)
    {
        out.clear();
        if (m_bbox.empty() || !m_extent.ok()) { return; }
        ++m_query;
        int ilo, ihi, jlo, jhi;
        binRange(target, ilo, ihi, jlo, jhi);
        for (int bj = jlo; bj <= jhi; ++bj) {
            for (int bi = ilo; bi <= ihi; ++bi) {
                const int b = bj * m_nx + bi;
                for (int k = m_offsets[b]; k < m_offsets[b+1]; ++k) {
                    const int n = m_items[k];
                    if (m_stamp[n] == m_query) { continue; }
                    m_stamp[n] = m_query;
                    if (m_bbox[n].intersects(target)) { out.push_back(n); }
                }
            }
        }
    }

    [[nodiscard]] const PlanarBBox& itemBBox (int item) const { return m_bbox[item]; }
    [[nodiscard]] amrex::IntVect itemIndex (int item) const {
        const int nx_cells = std::max(1, m_src_box.length(0));
        return amrex::IntVect(m_src_box.smallEnd(0) + item % nx_cells,
                              m_src_box.smallEnd(1) + item / nx_cells,
                              0);
    }
    [[nodiscard]] int numBins () const { return m_nx * m_ny; }

private:
    void binRange (const PlanarBBox& bb, int& ilo, int& ihi, int& jlo, int& jhi) const
    {
        ilo = std::clamp(static_cast<int>((bb.xlo - m_extent.xlo) * m_inv_dx), 0, m_nx - 1);
        ihi = std::clamp(static_cast<int>((bb.xhi - m_extent.xlo) * m_inv_dx), 0, m_nx - 1);
        jlo = std::clamp(static_cast<int>((bb.ylo - m_extent.ylo) * m_inv_dy), 0, m_ny - 1);
        jhi = std::clamp(static_cast<int>((bb.yhi - m_extent.ylo) * m_inv_dy), 0, m_ny - 1);
    }

    amrex::Box m_src_box;
    int m_nx = 0, m_ny = 0;
    amrex::Real m_inv_dx = 0.0, m_inv_dy = 0.0;
    PlanarBBox m_extent;
    amrex::Vector<PlanarBBox> m_bbox;
    amrex::Vector<int> m_offsets;
    amrex::Vector<int> m_items;
    amrex::Vector<int> m_stamp;
    int m_query = 0;
};

// Control-volume quad for a face-staggered destination point, built from the
// midpoints of adjacent nodal edges. The u-cell at face index (i,j) spans from
// the center of cell (i-1,j) to the center of cell (i,j), so its corners are the
// midpoints of the node pairs bracketing that span. Without this, a face family
// would reuse the rho-cell polygon at a face index, i.e. a half-cell offset.
//
// `node_domain` clamps the reads: at a domain-edge face one of the bracketing
// nodes lies outside the grid, where FillBoundary has left nothing valid. Clamping
// collapses that side onto the boundary node, which is the correct half-width
// control volume for an edge face rather than a quad built from zeros.
amrex::Vector<Point2D>
FaceCellQuadFromNodes (int i, int j, int dir,
                       const amrex::Array4<const amrex::Real>& xn,
                       const amrex::Array4<const amrex::Real>& yn,
                       const amrex::Box& node_domain)
{
    const int ilo = node_domain.smallEnd(0), ihi = node_domain.bigEnd(0);
    const int jlo = node_domain.smallEnd(1), jhi = node_domain.bigEnd(1);

    auto mid = [&] (int ia, int ja, int ib, int jb) {
        ia = std::clamp(ia, ilo, ihi); ib = std::clamp(ib, ilo, ihi);
        ja = std::clamp(ja, jlo, jhi); jb = std::clamp(jb, jlo, jhi);
        return Point2D{ amrex::Real(0.5) * (xn(ia,ja,0) + xn(ib,jb,0)),
                        amrex::Real(0.5) * (yn(ia,ja,0) + yn(ib,jb,0)) };
    };

    amrex::Vector<Point2D> quad;
    if (dir == 0) {
        quad = { mid(i-1, j  , i, j  ), mid(i, j  , i+1, j  ),
                 mid(i, j+1, i+1, j+1), mid(i-1, j+1, i, j+1) };
    } else {
        quad = { mid(i  , j-1, i  , j), mid(i+1, j-1, i+1, j),
                 mid(i+1, j  , i+1, j+1), mid(i  , j  , i  , j+1) };
    }
    EnsureCounterClockwise(quad);
    return quad;
}

// Extend a nodal corner mesh outward by one ring, so cells on the outer ring of
// the owning grid still have corners. Direct port of SCRIP's
// create_extra_rho_grid (COAWST Lib/SCRIP_COAWST/create_fullgrid.f): linear
// extrapolation along the four edges, then the four corners from the
// already-extended edges. Currently only needed when the exchange region
// includes a grid's outer ring (e.g. an ocean nested inside the atmosphere).
void
ExtendCornerMeshOneRing (amrex::MultiFab& node_mf, const amrex::Box& node_domain)
{
    using namespace amrex;

    for (MFIter mfi(node_mf); mfi.isValid(); ++mfi) {
        const Box gbx = mfi.growntilebox();
        auto arr = node_mf.array(mfi);

        const int ilo = node_domain.smallEnd(0), ihi = node_domain.bigEnd(0);
        const int jlo = node_domain.smallEnd(1), jhi = node_domain.bigEnd(1);

        for (int c = 0; c < node_mf.nComp(); ++c) {
            // Edges
            for (int j = std::max(gbx.smallEnd(1), jlo); j <= std::min(gbx.bigEnd(1), jhi); ++j) {
                if (gbx.contains(IntVect(ilo-1, j, 0))) {
                    arr(ilo-1, j, 0, c) = arr(ilo, j, 0, c) - (arr(ilo+1, j, 0, c) - arr(ilo, j, 0, c));
                }
                if (gbx.contains(IntVect(ihi+1, j, 0))) {
                    arr(ihi+1, j, 0, c) = arr(ihi, j, 0, c) + (arr(ihi, j, 0, c) - arr(ihi-1, j, 0, c));
                }
            }
            for (int i = std::max(gbx.smallEnd(0), ilo); i <= std::min(gbx.bigEnd(0), ihi); ++i) {
                if (gbx.contains(IntVect(i, jlo-1, 0))) {
                    arr(i, jlo-1, 0, c) = arr(i, jlo, 0, c) - (arr(i, jlo+1, 0, c) - arr(i, jlo, 0, c));
                }
                if (gbx.contains(IntVect(i, jhi+1, 0))) {
                    arr(i, jhi+1, 0, c) = arr(i, jhi, 0, c) + (arr(i, jhi, 0, c) - arr(i, jhi-1, 0, c));
                }
            }
            // Corners, from the edges filled just above
            auto fill_corner = [&] (int ic, int jc, int jn, int jn2) {
                if (gbx.contains(IntVect(ic, jc, 0)) && gbx.contains(IntVect(ic, jn2, 0))) {
                    arr(ic, jc, 0, c) = arr(ic, jn, 0, c) + (arr(ic, jn, 0, c) - arr(ic, jn2, 0, c));
                }
            };
            fill_corner(ilo-1, jlo-1, jlo,   jlo+1);
            fill_corner(ilo-1, jhi+1, jhi,   jhi-1);
            fill_corner(ihi+1, jlo-1, jlo,   jlo+1);
            fill_corner(ihi+1, jhi+1, jhi,   jhi-1);
        }
    }
}
}

void
ERFRemoraMultiBlockContainer::BuildRemapWeights()
{
    ResolveRemapFrame();
    ResolveDonorFillMode();

    if (m_remap_spherical) {
        BuildERFCornerMesh();
        BuildERFStagingPlan();
        if (m_use_two_way_coupling) {
            BuildREMORAStagingPlan();
        }
    }

    BuildAtmosToOceanWeights();
    BuildOceanToAtmosWeights();
}

// Decide the frame, and refuse to guess. ERF reports lon/lat degrees
// (XLONG_U/XLAT_V) while REMORA's x_psi/y_psi are meters in a projection the
// grid file does not record, so a half-available pair cannot be reconciled:
// clipping degrees against meters yields zero overlap for every destination
// cell, which normalizes to all-zero weights and a silently unforced ocean.
// SCRIP makes the same either/or choice from the grid file's `spherical` switch
// (COAWST Lib/SCRIP_COAWST/read_roms.f).
void
ERFRemoraMultiBlockContainer::ResolveRemapFrame()
{
    using namespace amrex;

    const MultiFab* remora_lon = nullptr;
    const MultiFab* remora_lat = nullptr;
    const MultiFab* erf_lon = nullptr;
    const MultiFab* erf_lat = nullptr;
    if (m_remora != nullptr) {
        REMORAGetAtmosToOceanPsiLonLat(m_remora, &remora_lon, &remora_lat);
    }
    if (m_erf != nullptr) {
        ERFGetOceanToAtmosCornerCoordinates(m_erf, &erf_lon, &erf_lat);
    }

    const bool remora_has_lonlat = (remora_lon != nullptr && remora_lat != nullptr);
    const bool erf_has_lonlat    = (erf_lon != nullptr && erf_lat != nullptr);

    auto describe = [] (const char* side, const MultiFab* x, const MultiFab* y) {
        std::string s(side);
        if (x == nullptr || y == nullptr) { return s + "=<none>"; }
        return s + "=[" + std::to_string(x->min(0)) + "," + std::to_string(x->max(0)) + "]x[" +
               std::to_string(y->min(0)) + "," + std::to_string(y->max(0)) + "]";
    };

    if (m_remap_frame_request == "cartesian") {
        m_remap_spherical = false;
    } else if (m_remap_frame_request == "spherical") {
        if (!remora_has_lonlat || !erf_has_lonlat) {
            amrex::Abort("driver.remap_frame=spherical requires lon/lat corner coordinates from both "
                         "solvers, but " + describe("remora", remora_lon, remora_lat) + " " +
                         describe("erf", erf_lon, erf_lat) +
                         ". Provide lon_psi/lat_psi in the REMORA grid file and a lat/lon-aware ERF "
                         "initialization, or set driver.remap_frame=cartesian.");
        }
        m_remap_spherical = true;
    } else if (m_remap_frame_request == "auto") {
        if (remora_has_lonlat && erf_has_lonlat) {
            m_remap_spherical = true;
        } else if (!remora_has_lonlat && !erf_has_lonlat) {
            m_remap_spherical = false;
        } else {
            amrex::Abort("driver.remap_frame=auto found lon/lat on only one side of the coupling: " +
                         describe("remora", remora_lon, remora_lat) + " " +
                         describe("erf", erf_lon, erf_lat) +
                         ". These frames cannot be mixed (degrees vs projected meters would produce "
                         "zero overlap everywhere). Set driver.remap_frame=cartesian to remap in the "
                         "shared driver RealBox frame, or supply lon/lat on both sides.");
        }
    } else {
        amrex::Abort("driver.remap_frame must be auto, spherical or cartesian; got '" +
                     m_remap_frame_request + "'");
    }

    amrex::Print() << "[driver] BuildRemapWeights: remap frame = "
                   << (m_remap_spherical ? "spherical" : "cartesian")
                   << " (requested " << m_remap_frame_request << ")\n";
    if (m_remap_verbose >= 1) {
        if (m_remap_spherical) {
            amrex::Print() << "[driver] remap corner ranges: "
                           << describe("remora(lon,lat)", remora_lon, remora_lat) << " "
                           << describe("erf(lon,lat)", erf_lon, erf_lat) << "\n";
        } else {
            amrex::Print() << "[driver] remap Cartesian extents: remora=["
                           << m_ocn_rb.lo(0) << "," << m_ocn_rb.hi(0) << "]x["
                           << m_ocn_rb.lo(1) << "," << m_ocn_rb.hi(1) << "] erf=["
                           << m_atm_rb.lo(0) << "," << m_atm_rb.hi(0) << "]x["
                           << m_atm_rb.lo(1) << "," << m_atm_rb.hi(1) << "]\n";
        }
    }
}

// Decide what an uncovered wet destination cell receives, and refuse to guess:
// like driver.remap_frame, an unrecognized request aborts with the valid values
// rather than silently selecting one of them.
//
// The default is withhold, not COAWST fidelity. A source-level survey of
// production couplers found none that donor-fills an unmapped destination cell:
// ESMF's unmappedaction defaults to ERROR, CMEPS pre-fills the destination and
// then maps with zeroregion=REGION_SELECT so unmapped cells keep the pre-fill,
// YAC does not interpolate unmapped points, and CDO writes a missing value.
// SCRIP and OASIS are the two that fill, and both fill with zero -- OASIS
// documents that as a deficiency, not a feature. ESMF goes further and refuses
// to combine extrapolation with conservative regridding at all, because "doing
// so would result in non-conservative weights". That objection applies exactly
// here: the donor fill grafts a full-weight nearest-neighbour row into a
// fracarea-normalized conservative matrix. So COAWST's check_weights behaviour
// is opt-in and says so loudly, and the shipped default is the one the
// literature supports.
void
ERFRemoraMultiBlockContainer::ResolveDonorFillMode()
{
    if (m_remap_donor_fill_request == "withhold") {
        m_remap_donor_fill_mode = DonorFillMode::Withhold;
    } else if (m_remap_donor_fill_request == "coawst") {
        m_remap_donor_fill_mode = DonorFillMode::Coawst;
    } else {
        amrex::Abort("driver.remap_donor_fill must be withhold or coawst; got '" +
                     m_remap_donor_fill_request + "'");
    }

    if (m_remap_donor_max_diagonals < amrex::Real(0.0)) {
        amrex::Abort("driver.remap_donor_max_diagonals must be >= 0; got " +
                     std::to_string(m_remap_donor_max_diagonals));
    }

    if (m_remap_verbose >= 1) {
        amrex::Print() << "[driver] BuildRemapWeights: donor fill = "
                       << m_remap_donor_fill_request
                       << (m_remap_spherical ? "" : " (inert in the cartesian frame)") << "\n";
    }

    if (m_remap_donor_fill_mode != DonorFillMode::Coawst) {
        if (m_remap_donor_max_diagonals_set) {
            amrex::Print() << "[driver] remap: driver.remap_donor_max_diagonals is inert under "
                           << "driver.remap_donor_fill=withhold; it bounds the coawst donor "
                           << "search, and no donor search runs in withhold mode.\n";
        }
        return;
    }

    // Not verbosity-gated. Opting out of the literature-supported default must
    // never be silent, and the operator needs the word "coawst" in the log to
    // connect a non-conservative A2O operator back to the flag that caused it.
    amrex::Print()
        << "[driver] WARNING remap: driver.remap_donor_fill=coawst grafts a full-weight "
        << "nearest-neighbour stencil into the fracarea-normalized conservative matrix for "
        << "every wet destination cell with no source overlap. This reproduces COAWST's "
        << "check_weights (Lib/SCRIP_COAWST/scrip.f:343), but no production coupler does it: "
        << "ESMF defaults unmappedaction=ERROR and refuses to combine extrapolation with "
        << "conservative regridding because \"doing so would result in non-conservative "
        << "weights\"; CMEPS pre-fills and maps with zeroregion=REGION_SELECT; YAC does not "
        << "interpolate unmapped points; CDO writes a missing value. The resulting "
        << "atmosphere->ocean operator is NOT conservative. The default, "
        << "driver.remap_donor_fill=withhold, instead leaves these cells to the per-lane "
        << "fallbacks in ApplyConservativeRemap and reports them uncovered.\n";

    // Say plainly what the bounds do and do not buy, so nobody reads the clip as
    // a rank-invariance guarantee it does not deliver.
    amrex::Print()
        << "[driver] remap: coawst donors are restricted to source cells whose centroid lies "
        << "inside the destination box's own grown lon/lat bounding box, and to within "
        << "driver.remap_donor_max_diagonals=" << m_remap_donor_max_diagonals
        << " donor diagonals of the destination centroid; destination boxes with no "
        << "geographic source overlap are not donor-filled at all. The geographic clip "
        << "removes the ERF decomposition from the donor choice, because every source cell "
        << "inside that bounding box is staged under every ERF decomposition. It does NOT "
        << "remove the REMORA decomposition: the bounding box is itself per-REMORA-box, so a "
        << "destination cell near its edge can have the same donor accepted under one REMORA "
        << "box layout and rejected under another. The acceptance radius has no precedent in "
        << "any production coupler -- ESMF's nearest-source search is explicitly unbounded -- "
        << "and is a deliberate departure.\n";

    if (!m_remap_spherical) {
        amrex::Print() << "[driver] remap: driver.remap_donor_fill=coawst has no effect in the "
                       << "cartesian frame; the donor search only runs when the frame is "
                       << "spherical.\n";
    }
}

// ERF-side geometry, built entirely on the ERF layout.
//
// ERF stores lon_m from XLONG_U (longitude on the cell's west edge, at the
// cell's mid-latitude) and lat_m from XLAT_V (latitude on the south edge, at the
// mid-longitude). Using those directly as corners is not watertight: cell
// (i,j)'s NE corner would be (lon_m(i+1,j), lat_m(i,j+1)) while cell (i+1,j+1)'s
// SW corner is (lon_m(i+1,j+1), lat_m(i+1,j+1)), and under Lambert rotation
// those differ, leaving gaps and overlaps in a mesh that a conservative remap
// assumes tiles exactly. One staggered average per direction moves each field to
// the node, so every node is computed once and shared by all four adjacent
// cells. This is the same centers-to-corners step as SCRIP's create_psimesh
// (COAWST Lib/SCRIP_COAWST/read_wrf.f), adapted to ERF's edge-staggered storage.
void
ERFRemoraMultiBlockContainer::BuildERFCornerMesh()
{
    using namespace amrex;

    const MultiFab* erf_lon = nullptr;
    const MultiFab* erf_lat = nullptr;
    ERFGetOceanToAtmosCornerCoordinates(m_erf, &erf_lon, &erf_lat);
    AMREX_ALWAYS_ASSERT(erf_lon != nullptr && erf_lat != nullptr);

    const iMultiFab* erf_is_land = nullptr;
    ERFGetLandMask(m_erf, &erf_is_land);

    // Reconcile longitude conventions. ERF carries WRF's XLONG_U, which is
    // signed (-180..180); a ROMS grid file may store lon_psi either signed or as
    // 0..360, and the two never overlap when they disagree - every weight would
    // come out zero. Pick the single 360-degree shift for the ERF longitudes that
    // maximizes overlap with REMORA's range. Choosing from a fixed set of three
    // and scoring by overlap means this is a no-op when the conventions already
    // agree, and it cannot scramble a grid that legitimately straddles the
    // dateline: it either improves the overlap or does nothing.
    amrex::Real lon_shift = 0.0;
    {
        const MultiFab* remora_lon = nullptr;
        const MultiFab* remora_lat = nullptr;
        REMORAGetAtmosToOceanPsiLonLat(m_remora, &remora_lon, &remora_lat);
        AMREX_ALWAYS_ASSERT(remora_lon != nullptr);

        const Real erf_lo = erf_lon->min(0), erf_hi = erf_lon->max(0);
        const Real rem_lo = remora_lon->min(0), rem_hi = remora_lon->max(0);

        Real best_overlap = -1.0;
        for (const Real trial : {Real(-360.0), Real(0.0), Real(360.0)}) {
            const Real overlap = std::min(erf_hi + trial, rem_hi) - std::max(erf_lo + trial, rem_lo);
            if (overlap > best_overlap) { best_overlap = overlap; lon_shift = trial; }
        }

        if (lon_shift != 0.0) {
            amrex::Print() << "[driver] remap: shifting ERF longitudes by " << lon_shift
                           << " deg to match the REMORA convention (erf=[" << erf_lo << "," << erf_hi
                           << "] remora=[" << rem_lo << "," << rem_hi << "])\n";
        } else if (best_overlap <= 0.0) {
            amrex::Print() << "[driver] WARNING remap: ERF longitudes [" << erf_lo << "," << erf_hi
                           << "] and REMORA longitudes [" << rem_lo << "," << rem_hi
                           << "] do not overlap under any 360-degree shift. The two grids do not "
                           << "cover a common region; expect zero coverage below.\n";
        }
    }

    const BoxArray& erf_cell_ba = erf_lon->boxArray();
    const DistributionMapping& erf_dm = erf_lon->DistributionMap();
    const BoxArray erf_node_ba = convert(erf_cell_ba, IntVect(1,1,0));

    // One ghost node so the ring extension (when enabled) has somewhere to write.
    m_erf_node_coords = std::make_unique<MultiFab>(erf_node_ba, erf_dm, 2, IntVect(1,1,0));
    m_erf_node_coords->setVal(0.0);

    const Box node_domain = erf_node_ba.minimalBox();

    for (MFIter mfi(*m_erf_node_coords); mfi.isValid(); ++mfi) {
        const Box nbx = mfi.validbox();
        auto node = m_erf_node_coords->array(mfi);
        auto const& lon = erf_lon->const_array(mfi);
        auto const& lat = erf_lat->const_array(mfi);

        for (int j = nbx.smallEnd(1); j <= nbx.bigEnd(1); ++j) {
            for (int i = nbx.smallEnd(0); i <= nbx.bigEnd(0); ++i) {
                // lon_m is already nodal in i and centered in j; lat_m the reverse.
                node(i,j,0,0) = Real(0.5) * (lon(i,j-1,0) + lon(i,j,0)) + lon_shift;
                node(i,j,0,1) = Real(0.5) * (lat(i-1,j,0) + lat(i,j,0));
            }
        }
    }
    // Interior ghosts only: every valid node above is already correct, and the
    // ring extension needs somewhere outside the domain to write. Deliberately
    // not periodicity-aware, so this cannot depend on the driver slab geometry
    // agreeing with ERF's 2D BoxArray.
    m_erf_node_coords->FillBoundary();

    if (m_remap_extend_corner_ring) {
        ExtendCornerMeshOneRing(*m_erf_node_coords, node_domain);
    }

    // Land mask as Real so it rides the same staging machinery as the coords.
    m_erf_land_mask = std::make_unique<MultiFab>(erf_cell_ba, erf_dm, 1, 0);
    if (erf_is_land != nullptr) {
        for (MFIter mfi(*m_erf_land_mask); mfi.isValid(); ++mfi) {
            const Box bx = mfi.validbox();
            auto mask = m_erf_land_mask->array(mfi);
            auto const& is_land = erf_is_land->const_array(mfi);
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    mask(i,j,0) = static_cast<Real>(is_land(i,j,0));
                }
            }
        }
    } else {
        m_erf_land_mask->setVal(0.0);   // no mask reported: treat everything as water
    }

    if (m_remap_verbose >= 1) {
        // The averaged mesh is watertight by construction (one value per shared
        // node), so measure instead how far it moved the corner relative to
        // using the raw edge-staggered arrays as corners. That displacement is
        // exactly the gap the naive construction would have left between
        // neighboring source cells, and it is the number to watch if the corner
        // source ever changes.
        Real max_shift = 0.0;
        for (MFIter mfi(*m_erf_node_coords); mfi.isValid(); ++mfi) {
            const Box nbx = mfi.validbox() & node_domain;
            auto const& node = m_erf_node_coords->const_array(mfi);
            auto const& lon  = erf_lon->const_array(mfi);
            auto const& lat  = erf_lat->const_array(mfi);
            for (int j = nbx.smallEnd(1); j <= nbx.bigEnd(1); ++j) {
                for (int i = nbx.smallEnd(0); i <= nbx.bigEnd(0); ++i) {
                    // Undo lon_shift so this measures the corner displacement, not
                    // the convention change.
                    const Real dx = (node(i,j,0,0) - lon_shift) - lon(i,j,0);
                    const Real dy = node(i,j,0,1) - lat(i,j,0);
                    max_shift = std::max(max_shift, std::sqrt(dx*dx + dy*dy));
                }
            }
        }
        ParallelDescriptor::ReduceRealMax(max_shift);
        amrex::Print() << "[driver] remap ERF corner mesh: nodes=" << node_domain
                       << " max corner shift vs raw edge arrays=" << max_shift << " deg"
                       << " ring_extension=" << (m_remap_extend_corner_ring ? "on" : "off") << "\n";
    }
}

// Route the ERF corner mesh and land mask onto the REMORA layout.
//
// A shared MFIter is impossible here: mfi.index() is a global box index, so
// reading an ERF-layout MultiFab under a REMORA-layout MFIter requires identical
// BoxArray length and DistributionMapping, which two independently decomposed
// grids never have. So each REMORA box gets told which ERF index region its
// geographic footprint needs, and amrex::NonLocalBC moves exactly that.
//
// The per-box ERF region is derived from bounding boxes only, which is why this
// needs no coordinate data outside its own layout: each rank computes its own
// boxes' geographic bounds, those bounds are reduced to every rank (a few reals
// per box), and then every rank derives the same staging BoxArray. A BoxArray
// must be globally consistent, hence the reduction rather than a local
// computation.
void
ERFRemoraMultiBlockContainer::BuildERFStagingPlan()
{
    using namespace amrex;

    m_staging_plan = std::make_unique<RemapStagingPlan>();
    RemapStagingPlan& plan = *m_staging_plan;

    const MultiFab& erf_node = *m_erf_node_coords;
    const BoxArray& erf_cell_ba = m_erf_land_mask->boxArray();
    plan.erf_cell_domain = erf_cell_ba.minimalBox();
    plan.erf_node_domain = erf_node.boxArray().minimalBox();
    plan.erf_domain_cells = static_cast<long long>(plan.erf_cell_domain.numPts());

    // ERF box bounds in the shared frame, reduced so every rank sees all boxes.
    const int n_erf_boxes = static_cast<int>(erf_cell_ba.size());
    Vector<Real> erf_lo(2 * n_erf_boxes,  std::numeric_limits<Real>::max());
    Vector<Real> erf_hi(2 * n_erf_boxes, -std::numeric_limits<Real>::max());
    for (MFIter mfi(erf_node); mfi.isValid(); ++mfi) {
        const int b = mfi.index();
        const Box nbx = mfi.validbox();
        auto const& node = erf_node.const_array(mfi);
        PlanarBBox bb;
        for (int j = nbx.smallEnd(1); j <= nbx.bigEnd(1); ++j) {
            for (int i = nbx.smallEnd(0); i <= nbx.bigEnd(0); ++i) {
                bb.include(node(i,j,0,0), node(i,j,0,1));
            }
        }
        erf_lo[2*b  ] = bb.xlo; erf_hi[2*b  ] = bb.xhi;
        erf_lo[2*b+1] = bb.ylo; erf_hi[2*b+1] = bb.yhi;
    }
    ParallelDescriptor::ReduceRealMin(erf_lo.dataPtr(), erf_lo.size());
    ParallelDescriptor::ReduceRealMax(erf_hi.dataPtr(), erf_hi.size());

    // REMORA box bounds, same treatment. The rho layout drives this: face
    // families live on staggered versions of the same boxes, so the rho
    // footprint grown by one cell covers them.
    const MultiFab* remora_lon = nullptr;
    const MultiFab* remora_lat = nullptr;
    REMORAGetAtmosToOceanPsiLonLat(m_remora, &remora_lon, &remora_lat);
    AMREX_ALWAYS_ASSERT(remora_lon != nullptr && remora_lat != nullptr);

    const int n_remora_boxes = static_cast<int>(m_ba_ocn.size());
    Vector<Real> remora_lo(2 * n_remora_boxes,  std::numeric_limits<Real>::max());
    Vector<Real> remora_hi(2 * n_remora_boxes, -std::numeric_limits<Real>::max());
    for (MFIter mfi(*remora_lon); mfi.isValid(); ++mfi) {
        const int b = mfi.index();
        // One extra node ring so face-staggered REMORA points and their
        // edge-midpoint control volumes stay inside the requested ERF region.
        const Box nbx = grow(mfi.validbox(), IntVect(1,1,0)) & remora_lon->boxArray().minimalBox();
        auto const& lon = remora_lon->const_array(mfi);
        auto const& lat = remora_lat->const_array(mfi);
        PlanarBBox bb;
        for (int j = nbx.smallEnd(1); j <= nbx.bigEnd(1); ++j) {
            for (int i = nbx.smallEnd(0); i <= nbx.bigEnd(0); ++i) {
                bb.include(lon(i,j,0), lat(i,j,0));
            }
        }
        remora_lo[2*b  ] = bb.xlo; remora_hi[2*b  ] = bb.xhi;
        remora_lo[2*b+1] = bb.ylo; remora_hi[2*b+1] = bb.yhi;
    }
    ParallelDescriptor::ReduceRealMin(remora_lo.dataPtr(), remora_lo.size());
    ParallelDescriptor::ReduceRealMax(remora_hi.dataPtr(), remora_hi.size());

    // Every rank now derives the same staging BoxArray.
    BoxList cell_bl, node_bl;
    plan.box_bbox_lo.resize(2 * n_remora_boxes);
    plan.box_bbox_hi.resize(2 * n_remora_boxes);
    plan.box_has_erf_overlap.resize(n_remora_boxes);
    for (int d = 0; d < n_remora_boxes; ++d) {
        PlanarBBox remora_bb;
        remora_bb.xlo = remora_lo[2*d  ]; remora_bb.xhi = remora_hi[2*d  ];
        remora_bb.ylo = remora_lo[2*d+1]; remora_bb.yhi = remora_hi[2*d+1];
        plan.box_bbox_lo[2*d  ] = remora_bb.xlo;
        plan.box_bbox_hi[2*d  ] = remora_bb.xhi;
        plan.box_bbox_lo[2*d+1] = remora_bb.ylo;
        plan.box_bbox_hi[2*d+1] = remora_bb.yhi;

        Box need;
        for (int s = 0; s < n_erf_boxes; ++s) {
            PlanarBBox erf_bb;
            erf_bb.xlo = erf_lo[2*s  ]; erf_bb.xhi = erf_hi[2*s  ];
            erf_bb.ylo = erf_lo[2*s+1]; erf_bb.yhi = erf_hi[2*s+1];
            if (!remora_bb.ok() || !erf_bb.ok() || !remora_bb.intersects(erf_bb)) { continue; }
            need = need.ok() ? minBox(need, erf_cell_ba[s]) : erf_cell_ba[s];
        }

        plan.box_has_erf_overlap[d] = static_cast<char>(need.ok());
        if (!need.ok()) {
            // No geographic overlap for this REMORA box. A BoxArray cannot hold
            // an empty box, so park a 1-cell placeholder inside the ERF domain;
            // its weights stay zero and the coverage diagnostic counts it.
            ++plan.empty_boxes;
            need = Box(plan.erf_cell_domain.smallEnd(), plan.erf_cell_domain.smallEnd());
        } else {
            need &= plan.erf_cell_domain;
        }

        plan.staged_cells += static_cast<long long>(need.numPts());
        plan.max_box_cells = std::max(plan.max_box_cells, static_cast<int>(need.numPts()));
        cell_bl.push_back(need);
        node_bl.push_back(convert(need, IntVect(1,1,0)));
    }

    plan.cell_ba = BoxArray(std::move(cell_bl));
    plan.node_ba = BoxArray(std::move(node_bl));

    // Staged copies live on the REMORA DistributionMapping, so box index d is
    // owned by the same rank that owns REMORA box d. The staging boxes overlap
    // where neighboring REMORA boxes need overlapping ERF regions; that is fine
    // for a copy, since MultiBlockCommMetaData::define builds tags per
    // destination box independently.
    plan.node_mf = std::make_unique<MultiFab>(plan.node_ba, m_dm_ocn, 2, 0);
    plan.land_mf = std::make_unique<MultiFab>(plan.cell_ba, m_dm_ocn, 1, 0);
    plan.node_mf->setVal(0.0);
    plan.land_mf->setVal(0.0);

    // Identity index mapping: the staging BoxArrays are already expressed in ERF
    // index space, so no shift or permutation is needed.
    NonLocalBC::MultiBlockIndexMapping dtos;
    dtos.permutation = IntVect(AMREX_D_DECL(0,1,2));
    dtos.sign        = IntVect(AMREX_D_DECL(1,1,1));
    dtos.offset      = IntVect(AMREX_D_DECL(0,0,0));

    plan.node_cmd = std::make_unique<NonLocalBC::MultiBlockCommMetaData>(
        plan.node_ba, m_dm_ocn, plan.erf_node_domain,
        erf_node.boxArray(), erf_node.DistributionMap(), IntVect(0), dtos);
    plan.cell_cmd = std::make_unique<NonLocalBC::MultiBlockCommMetaData>(
        plan.cell_ba, m_dm_ocn, plan.erf_cell_domain,
        erf_cell_ba, m_erf_land_mask->DistributionMap(), IntVect(0), dtos);

    NonLocalBC::ParallelCopy(*plan.node_mf, erf_node, *plan.node_cmd, 0, 0, 2, dtos);
    NonLocalBC::ParallelCopy(*plan.land_mf, *m_erf_land_mask, *plan.cell_cmd, 0, 0, 1, dtos);

    if (m_remap_verbose >= 1) {
        amrex::Print() << "[driver] remap staging plan: remora_boxes=" << n_remora_boxes
                       << " erf_boxes=" << n_erf_boxes
                       << " staged_cells=" << plan.staged_cells
                       << " (erf domain " << plan.erf_domain_cells << ")"
                       << " max_box_cells=" << plan.max_box_cells
                       << " remora_boxes_without_overlap=" << plan.empty_boxes << "\n";
    }
    if (m_remap_verbose >= 2) {
        for (int d = 0; d < n_remora_boxes; ++d) {
            amrex::Print() << "[driver] remap staging box " << d
                           << " bbox=[(" << remora_lo[2*d] << ", "
                           << remora_lo[2*d+1] << "), ("
                           << remora_hi[2*d] << ", "
                           << remora_hi[2*d+1] << ")]\n";
        }
    }
}

// Route the REMORA psi corner mesh and rho wet mask onto the ERF layout.
//
// This mirrors BuildERFStagingPlan() with the source and destination reversed:
// each ERF box gets told which REMORA rho-cell index region its geographic
// footprint needs, and amrex::NonLocalBC moves that region onto the ERF
// DistributionMapping. O2A is SST-only today, so there is no weight-family loop
// and only the rho mask is staged.
void
ERFRemoraMultiBlockContainer::BuildREMORAStagingPlan()
{
    using namespace amrex;

    m_remora_staging_plan = std::make_unique<REMORAStagingPlan>();
    REMORAStagingPlan& plan = *m_remora_staging_plan;

    const MultiFab* remora_lon = nullptr;
    const MultiFab* remora_lat = nullptr;
    REMORAGetAtmosToOceanPsiLonLat(m_remora, &remora_lon, &remora_lat);
    AMREX_ALWAYS_ASSERT(remora_lon != nullptr && remora_lat != nullptr);

    const MultiFab* remora_mskr = nullptr;
    REMORAGetMasks(m_remora, &remora_mskr, nullptr, nullptr);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        remora_mskr != nullptr,
        "BuildREMORAStagingPlan requires REMORA's rho wet mask for spherical O2A SST remap.");

    const BoxArray& remora_cell_ba = remora_mskr->boxArray();
    const BoxArray remora_node_ba = convert(remora_cell_ba, IntVect(1,1,0));
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        remora_lon->boxArray() == remora_node_ba,
        "BuildREMORAStagingPlan requires REMORA psi coordinates to be nodal rho boxes.");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        remora_lat->boxArray() == remora_lon->boxArray(),
        "BuildREMORAStagingPlan requires REMORA lon_psi and lat_psi to share box indexing.");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        remora_lat->DistributionMap().ProcessorMap() == remora_lon->DistributionMap().ProcessorMap(),
        "BuildREMORAStagingPlan requires REMORA lon_psi and lat_psi to share a DistributionMapping.");

    const MultiFab& erf_node = *m_erf_node_coords;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_weight_o2a_mf != nullptr,
        "BuildREMORAStagingPlan requires the ERF-layout O2A SST weights to be allocated.");
    const BoxArray& erf_cell_ba = m_weight_o2a_mf->boxArray();
    const DistributionMapping& erf_dm = m_weight_o2a_mf->DistributionMap();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        erf_node.boxArray() == convert(erf_cell_ba, IntVect(1,1,0)),
        "BuildREMORAStagingPlan requires ERF corner coordinates to be nodal O2A weight boxes.");

    plan.remora_cell_domain = remora_cell_ba.minimalBox();
    plan.remora_node_domain = remora_lon->boxArray().minimalBox();
    plan.remora_domain_cells = static_cast<long long>(plan.remora_cell_domain.numPts());

    // REMORA box bounds in the shared frame, reduced so every rank sees all
    // source boxes. REMORA lon_psi/lat_psi are already psi-point cell corners,
    // so unlike the ERF mesh construction there is no averaging step here.
    const int n_remora_boxes = static_cast<int>(remora_cell_ba.size());
    Vector<Real> remora_lo(2 * n_remora_boxes,  std::numeric_limits<Real>::max());
    Vector<Real> remora_hi(2 * n_remora_boxes, -std::numeric_limits<Real>::max());
    for (MFIter mfi(*remora_lon); mfi.isValid(); ++mfi) {
        const int b = mfi.index();
        const Box nbx = mfi.validbox();
        auto const& lon = remora_lon->const_array(mfi);
        auto const& lat = remora_lat->const_array(mfi);
        PlanarBBox bb;
        for (int j = nbx.smallEnd(1); j <= nbx.bigEnd(1); ++j) {
            for (int i = nbx.smallEnd(0); i <= nbx.bigEnd(0); ++i) {
                bb.include(lon(i,j,0), lat(i,j,0));
            }
        }
        remora_lo[2*b  ] = bb.xlo; remora_hi[2*b  ] = bb.xhi;
        remora_lo[2*b+1] = bb.ylo; remora_hi[2*b+1] = bb.yhi;
    }
    ParallelDescriptor::ReduceRealMin(remora_lo.dataPtr(), remora_lo.size());
    ParallelDescriptor::ReduceRealMax(remora_hi.dataPtr(), remora_hi.size());

    // ERF box bounds, same treatment. The destination is SST on ERF cells, so
    // the cell footprint from the nodal corner mesh is sufficient.
    const int n_erf_boxes = static_cast<int>(erf_cell_ba.size());
    Vector<Real> erf_lo(2 * n_erf_boxes,  std::numeric_limits<Real>::max());
    Vector<Real> erf_hi(2 * n_erf_boxes, -std::numeric_limits<Real>::max());
    for (MFIter mfi(erf_node); mfi.isValid(); ++mfi) {
        const int b = mfi.index();
        const Box nbx = mfi.validbox();
        auto const& node = erf_node.const_array(mfi);
        PlanarBBox bb;
        for (int j = nbx.smallEnd(1); j <= nbx.bigEnd(1); ++j) {
            for (int i = nbx.smallEnd(0); i <= nbx.bigEnd(0); ++i) {
                bb.include(node(i,j,0,0), node(i,j,0,1));
            }
        }
        erf_lo[2*b  ] = bb.xlo; erf_hi[2*b  ] = bb.xhi;
        erf_lo[2*b+1] = bb.ylo; erf_hi[2*b+1] = bb.yhi;
    }
    ParallelDescriptor::ReduceRealMin(erf_lo.dataPtr(), erf_lo.size());
    ParallelDescriptor::ReduceRealMax(erf_hi.dataPtr(), erf_hi.size());

    // Every rank now derives the same staging BoxArray.
    BoxList cell_bl, node_bl;
    for (int d = 0; d < n_erf_boxes; ++d) {
        PlanarBBox erf_bb;
        erf_bb.xlo = erf_lo[2*d  ]; erf_bb.xhi = erf_hi[2*d  ];
        erf_bb.ylo = erf_lo[2*d+1]; erf_bb.yhi = erf_hi[2*d+1];

        Box need;
        for (int s = 0; s < n_remora_boxes; ++s) {
            PlanarBBox remora_bb;
            remora_bb.xlo = remora_lo[2*s  ]; remora_bb.xhi = remora_hi[2*s  ];
            remora_bb.ylo = remora_lo[2*s+1]; remora_bb.yhi = remora_hi[2*s+1];
            if (!erf_bb.ok() || !remora_bb.ok() || !erf_bb.intersects(remora_bb)) { continue; }
            need = need.ok() ? minBox(need, remora_cell_ba[s]) : remora_cell_ba[s];
        }

        if (!need.ok()) {
            // No geographic overlap for this ERF box. A BoxArray cannot hold an
            // empty box, so park a 1-cell placeholder inside the REMORA domain;
            // the later O2A weight loop will leave coverage empty for it.
            ++plan.empty_boxes;
            need = Box(plan.remora_cell_domain.smallEnd(), plan.remora_cell_domain.smallEnd());
        } else {
            need &= plan.remora_cell_domain;
        }

        plan.staged_cells += static_cast<long long>(need.numPts());
        plan.max_box_cells = std::max(plan.max_box_cells, static_cast<int>(need.numPts()));
        cell_bl.push_back(need);
        node_bl.push_back(convert(need, IntVect(1,1,0)));
    }

    plan.cell_ba = BoxArray(std::move(cell_bl));
    plan.node_ba = BoxArray(std::move(node_bl));

    // Staged copies live on the ERF DistributionMapping, so box index d is
    // owned by the same rank that owns ERF/O2A box d. The staging boxes overlap
    // where neighboring ERF boxes need overlapping REMORA regions; that is fine
    // for a copy, since MultiBlockCommMetaData::define builds tags per
    // destination box independently.
    plan.node_mf = std::make_unique<MultiFab>(plan.node_ba, erf_dm, 2, 0);
    plan.wet_mf = std::make_unique<MultiFab>(plan.cell_ba, erf_dm, 1, 0);
    plan.node_mf->setVal(0.0);
    plan.wet_mf->setVal(0.0);

    // Identity index mapping: the staging BoxArrays are already expressed in
    // REMORA index space, so no shift or permutation is needed.
    NonLocalBC::MultiBlockIndexMapping dtos;
    dtos.permutation = IntVect(AMREX_D_DECL(0,1,2));
    dtos.sign        = IntVect(AMREX_D_DECL(1,1,1));
    dtos.offset      = IntVect(AMREX_D_DECL(0,0,0));

    plan.node_cmd = std::make_unique<NonLocalBC::MultiBlockCommMetaData>(
        plan.node_ba, erf_dm, plan.remora_node_domain,
        remora_lon->boxArray(), remora_lon->DistributionMap(), IntVect(0), dtos);
    plan.cell_cmd = std::make_unique<NonLocalBC::MultiBlockCommMetaData>(
        plan.cell_ba, erf_dm, plan.remora_cell_domain,
        remora_cell_ba, remora_mskr->DistributionMap(), IntVect(0), dtos);

    NonLocalBC::ParallelCopy(*plan.node_mf, *remora_lon, *plan.node_cmd, 0, 0, 1, dtos);
    NonLocalBC::ParallelCopy(*plan.node_mf, *remora_lat, *plan.node_cmd, 0, 1, 1, dtos);
    NonLocalBC::ParallelCopy(*plan.wet_mf, *remora_mskr, *plan.cell_cmd, 0, 0, 1, dtos);

    if (m_remap_verbose >= 1) {
        amrex::Print() << "[driver] O2A remap staging plan: erf_boxes=" << n_erf_boxes
                       << " remora_boxes=" << n_remora_boxes
                       << " staged_cells=" << plan.staged_cells
                       << " (remora domain " << plan.remora_domain_cells << ")"
                       << " max_box_cells=" << plan.max_box_cells
                       << " erf_boxes_without_overlap=" << plan.empty_boxes << "\n";
    }
}

void
ERFRemoraMultiBlockContainer::BuildAtmosToOceanWeights()
{
    using namespace amrex;

    const MultiFab* remora_lon = nullptr;
    const MultiFab* remora_lat = nullptr;
    if (m_remap_spherical) {
        REMORAGetAtmosToOceanPsiLonLat(m_remora, &remora_lon, &remora_lat);
    }

    // REMORA wet/dry masks, so coverage can be reported over the cells that carry
    // ocean rather than over land the atmosphere grid was never going to force.
    const MultiFab* remora_mskr = nullptr;
    const MultiFab* remora_msku = nullptr;
    const MultiFab* remora_mskv = nullptr;
    if (m_remora != nullptr) {
        REMORAGetMasks(m_remora, &remora_mskr, &remora_msku, &remora_mskv);
    }

    // Reset A2O weight and index caches, one family slot at a time. Each slot
    // was allocated (AllocateCachesFromSolverLayouts) against its own
    // destination BoxArray/stagger (rho for SCALAR_CENTERED/WIND_FACE_TO_CENTER,
    // u-face for TAUX_FACE, v-face for TAUY_FACE), so building weights
    // per-slot means TauX/TauY finally get stencils computed against their own
    // stagger instead of borrowing the rho-point stencil.
    static const std::array<amrex::IntVect, kNumWeightFamilies> kExpectedFamilyStagger = {
        amrex::IntVect(0,0,0), amrex::IntVect(1,0,0), amrex::IntVect(0,1,0), amrex::IntVect(0,0,0)
    };
    // Stencil-capacity accounting: a destination cell that overlaps more than
    // MAX_STENCIL_SIZE source cells keeps the largest-area overlaps and reports
    // the shortfall. Aborting here would kill a run mid-init over a geometry
    // detail; SCRIP itself has no fixed per-cell cap.
    int max_stencil_used = 0;
    amrex::Long dropped_overlaps = 0;
    amrex::Long donor_filled = 0;
    Vector<int> candidates;
    Vector<std::pair<Real,IntVect>> overlaps;

    const Box remora_cell_domain = m_geom_ocn.Domain();

    for (int f = 0; f < kNumWeightFamilies; ++f) {
        MultiFab* weight_mf = m_weight_mf_by_family[f].get();
        iMultiFab* index_mf = m_index_mf_by_family[f].get();
        if (!weight_mf || !index_mf) { continue; }

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            weight_mf->boxArray().ixType().toIntVect() == kExpectedFamilyStagger[f],
            "BuildRemapWeights: A2O weight MultiFab stagger does not match its WeightFamily slot.");

        weight_mf->setVal(0.0);
        index_mf->setVal(-1);

        const IntVect stagger = kExpectedFamilyStagger[f];
        // TauX/TauY carry a face-staggered ERF source as well (tau13/tau23 on
        // x-/y-faces), so their source control volumes and stored indices must
        // be face-staggered too.
        const int src_face_dir =
            (stagger == IntVect(1,0,0)) ? 0 : ((stagger == IntVect(0,1,0)) ? 1 : -1);

        // Index space this family's source cells live in, on the ERF level-0
        // domain rather than on any one rank's staged region. Used only to key
        // the donor-search tie-break below. Empty in the cartesian frame, where
        // there is no staging plan and no donor search.
        const Box donor_src_domain =
            (!m_remap_spherical)
                ? Box()
                : ((src_face_dir < 0)
                       ? m_staging_plan->erf_cell_domain
                       : amrex::convert(m_staging_plan->erf_cell_domain,
                                        amrex::IntVect::TheDimensionVector(src_face_dir)));
        // Global source ID: row-major linear index over donor_src_domain. This
        // is deliberately a *domain* quantity. It contains no rank, no box
        // index and no reference to erf_staged_box, so two ranks that both hold
        // source cell (si,sj) agree on its ID under any decomposition. That is
        // the whole point: an MPI_MINLOC-style tie-break is only reproducible
        // when the key it minimizes over is a global cell ID -- feed it a rank
        // or a local index and the answer becomes rank-count dependent.
        auto donor_global_id = [&] (int si, int sj) -> long long {
            return static_cast<long long>(sj - donor_src_domain.smallEnd(1))
                     * static_cast<long long>(donor_src_domain.length(0))
                 + static_cast<long long>(si - donor_src_domain.smallEnd(0));
        };

        // Each family's mask must sit on that family's stagger, or const_array()
        // would be reading a differently shaped array.
        const MultiFab* family_mask = remora_mskr;
        if (f == static_cast<int>(WeightFamily::TAUX_FACE))      { family_mask = remora_msku; }
        else if (f == static_cast<int>(WeightFamily::TAUY_FACE)) { family_mask = remora_mskv; }
        if (family_mask != nullptr &&
            family_mask->boxArray().ixType() != weight_mf->boxArray().ixType()) {
            family_mask = nullptr;   // stagger mismatch: report unmasked rather than lie
        }

        for (MFIter mfi(*weight_mf); mfi.isValid(); ++mfi) {
            Box bx = mfi.validbox();
            auto weight_arr = weight_mf->array(mfi);
            auto index_arr = index_mf->array(mfi);

            // Spherical frame: REMORA corners from its psi lon/lat, and ERF
            // corners from the staged mesh. All of these arrays are on the REMORA
            // DistributionMapping, so one MFIter index is valid for all of them.
            SourceSearchGrid search;
            Array4<const Real> erf_node{};
            Array4<const Real> erf_land{};
            Box erf_staged_box;
            Box erf_staged_node_box;
            // Geographic bound on this destination box's donor search. G(B) is
            // the grown lon/lat bounding box BuildERFStagingPlan used to pick
            // the ERF boxes to stage, and donor_box_has_erf_overlap is 0 when it
            // intersected none of them.
            bool donor_box_has_erf_overlap = false;
            PlanarBBox donor_box_bbox;
            if (m_remap_spherical) {
                erf_node = m_staging_plan->node_mf->const_array(mfi);
                erf_land = m_staging_plan->land_mf->const_array(mfi);
                erf_staged_box = m_staging_plan->land_mf->boxArray()[mfi.index()];
                erf_staged_node_box = m_staging_plan->node_mf->boxArray()[mfi.index()];
                search.define(erf_staged_box, Array4<const Real>(erf_node, 0),
                              Array4<const Real>(erf_node, 1), m_remap_bin_cells);
                // Same index space the erf_staged_box lookup just above uses:
                // the staging plan's per-box vectors are indexed by REMORA box,
                // and this family's BoxArray is m_ba_ocn on that family's
                // stagger, so mfi.index() is that box index.
                const int d = mfi.index();
                AMREX_ASSERT(d < static_cast<int>(m_staging_plan->box_has_erf_overlap.size()));
                donor_box_has_erf_overlap =
                    (m_staging_plan->box_has_erf_overlap[d] != 0);
                donor_box_bbox.xlo = m_staging_plan->box_bbox_lo[2*d  ];
                donor_box_bbox.xhi = m_staging_plan->box_bbox_hi[2*d  ];
                donor_box_bbox.ylo = m_staging_plan->box_bbox_lo[2*d+1];
                donor_box_bbox.yhi = m_staging_plan->box_bbox_hi[2*d+1];
            }
            // Spherical ERF source geometry on the active source stagger. A
            // face is treated as land when either bracketing ERF cell is land,
            // matching REMORA's own msku/mskv convention.
            auto erf_src_poly = [&] (int si, int sj) {
                if (src_face_dir >= 0) {
                    return FaceCellQuadFromNodes(si, sj, src_face_dir,
                                                 Array4<const Real>(erf_node, 0),
                                                 Array4<const Real>(erf_node, 1),
                                                 erf_staged_node_box);
                }
                return ExtractCellQuadFromCoords(si, sj,
                                                 Array4<const Real>(erf_node, 0),
                                                 Array4<const Real>(erf_node, 1));
            };
            auto erf_src_is_land = [&] (int si, int sj) -> bool {
                if (src_face_dir == 0) {
                    const int im = std::max(si - 1, erf_staged_box.smallEnd(0));
                    const int ip = std::min(si,     erf_staged_box.bigEnd(0));
                    return erf_land(im, sj, 0) >= Real(0.5) || erf_land(ip, sj, 0) >= Real(0.5);
                }
                if (src_face_dir == 1) {
                    const int jm = std::max(sj - 1, erf_staged_box.smallEnd(1));
                    const int jp = std::min(sj,     erf_staged_box.bigEnd(1));
                    return erf_land(si, jm, 0) >= Real(0.5) || erf_land(si, jp, 0) >= Real(0.5);
                }
                return erf_land(si, sj, 0) >= Real(0.5);
            };
            // REMORA wet/dry mask on this family's stagger, used to decide which
            // empty destination cells deserve a donor search.
            const Array4<const Real> wet_arr =
                (family_mask != nullptr) ? family_mask->const_array(mfi) : Array4<const Real>{};
            const Array4<const Real> remora_lon_arr =
                m_remap_spherical ? remora_lon->const_array(mfi) : Array4<const Real>{};
            const Array4<const Real> remora_lat_arr =
                m_remap_spherical ? remora_lat->const_array(mfi) : Array4<const Real>{};
            const Box remora_node_domain =
                m_remap_spherical ? remora_lon->boxArray().minimalBox() : Box();

            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    Vector<Point2D> remora_poly;

                    if (m_remap_spherical) {
                        if (stagger == IntVect(1,0,0)) {
                            remora_poly = FaceCellQuadFromNodes(i, j, 0, remora_lon_arr, remora_lat_arr,
                                                                remora_node_domain);
                        } else if (stagger == IntVect(0,1,0)) {
                            remora_poly = FaceCellQuadFromNodes(i, j, 1, remora_lon_arr, remora_lat_arr,
                                                                remora_node_domain);
                        } else {
                            remora_poly = ExtractCellQuadFromCoords(i, j, remora_lon_arr, remora_lat_arr);
                        }
                    } else {
                        // Face families get their own control volume here too, not
                        // the rho-cell bounds evaluated at a face index.
                        RealBox bounds =
                            (stagger == IntVect(1,0,0))
                                ? FaceCellBoundsFromGeom(m_geom_ocn, remora_cell_domain, i, j, 0)
                                : (stagger == IntVect(0,1,0))
                                    ? FaceCellBoundsFromGeom(m_geom_ocn, remora_cell_domain, i, j, 1)
                                    : CellBoundsFromGeom(m_geom_ocn, i, j);
                        remora_poly = {
                            {bounds.lo(0), bounds.lo(1)},
                            {bounds.hi(0), bounds.lo(1)},
                            {bounds.hi(0), bounds.hi(1)},
                            {bounds.lo(0), bounds.hi(1)}
                        };
                    }

                    // Areas are planar in whatever frame is active. In the
                    // spherical frame that means degrees squared, not true
                    // spherical area as SCRIP computes: the cos(latitude) factor
                    // cancels to first order because only the ratio
                    // intersect/destination is used, and the residual bias across
                    // one ~2 km cell is far below the interpolation error. Revisit
                    // if this is ever used at basin-scale cell sizes.
                    Real remora_area = PolygonArea(remora_poly);
                    if (remora_area <= 0.0) continue;

                    const PlanarBBox remora_bbox = BBoxOfPolygon(remora_poly);

                    overlaps.clear();
                    if (m_remap_spherical) {
                        search.query(remora_bbox, candidates);
                        if (src_face_dir >= 0 && !candidates.empty()) {
                            // SourceSearchGrid indexes source cells. For face
                            // families, scan the candidate-cell hull and include
                            // the high face in the staggered direction.
                            int ci_lo = std::numeric_limits<int>::max();
                            int ci_hi = std::numeric_limits<int>::min();
                            int cj_lo = std::numeric_limits<int>::max();
                            int cj_hi = std::numeric_limits<int>::min();
                            for (const int item : candidates) {
                                const IntVect cij = search.itemIndex(item);
                                ci_lo = std::min(ci_lo, cij[0]); ci_hi = std::max(ci_hi, cij[0]);
                                cj_lo = std::min(cj_lo, cij[1]); cj_hi = std::max(cj_hi, cij[1]);
                            }
                            ci_hi += (src_face_dir == 0) ? 1 : 0;
                            cj_hi += (src_face_dir == 1) ? 1 : 0;
                            for (int erf_j = cj_lo; erf_j <= cj_hi; ++erf_j) {
                                for (int erf_i = ci_lo; erf_i <= ci_hi; ++erf_i) {
                                    if (erf_src_is_land(erf_i, erf_j)) { continue; }
                                    const Real area = PolygonArea(
                                        ClipPolygon(remora_poly, erf_src_poly(erf_i, erf_j)));
                                    if (area > 0.0) {
                                        overlaps.emplace_back(area, IntVect(erf_i, erf_j, 0));
                                    }
                                }
                            }
                        } else {
                            for (const int item : candidates) {
                                const IntVect erf_ij = search.itemIndex(item);
                                const int erf_i = erf_ij[0], erf_j = erf_ij[1];
                                // Wet/dry masking only (no vector rotation): skip
                                // land ERF cells so atmospheric land fluxes cannot
                                // bleed into coastal ocean cells via the remap
                                // stencil.
                                if (erf_land(erf_i, erf_j, 0) >= Real(0.5)) { continue; }

                                const Vector<Point2D> erf_poly =
                                    ExtractCellQuadFromCoords(erf_i, erf_j,
                                                              Array4<const Real>(erf_node, 0),
                                                              Array4<const Real>(erf_node, 1));
                                const Real area = PolygonArea(ClipPolygon(remora_poly, erf_poly));
                                if (area > 0.0) { overlaps.emplace_back(area, erf_ij); }
                            }
                        }
                    } else {
                        // Cartesian frame: both sides are uniform in the shared
                        // driver RealBox, so the candidate window comes straight
                        // from an affine index estimate plus a distortion pad.
                        Box erf_candidate_box = m_geom.Domain();
                        int erf_i_lo = std::floor((remora_bbox.xlo - m_atm_rb.lo(0)) / (m_atm_rb.hi(0) - m_atm_rb.lo(0)) * m_atm_n_cell[0]);
                        int erf_i_hi = std::ceil((remora_bbox.xhi - m_atm_rb.lo(0)) / (m_atm_rb.hi(0) - m_atm_rb.lo(0)) * m_atm_n_cell[0]);
                        int erf_j_lo = std::floor((remora_bbox.ylo - m_atm_rb.lo(1)) / (m_atm_rb.hi(1) - m_atm_rb.lo(1)) * m_atm_n_cell[1]);
                        int erf_j_hi = std::ceil((remora_bbox.yhi - m_atm_rb.lo(1)) / (m_atm_rb.hi(1) - m_atm_rb.lo(1)) * m_atm_n_cell[1]);

                        erf_i_lo = std::max(erf_candidate_box.smallEnd(0), erf_i_lo - 2);
                        erf_i_hi = std::min(erf_candidate_box.bigEnd(0) + ((src_face_dir == 0) ? 1 : 0),
                                            erf_i_hi + 2);
                        erf_j_lo = std::max(erf_candidate_box.smallEnd(1), erf_j_lo - 2);
                        erf_j_hi = std::min(erf_candidate_box.bigEnd(1) + ((src_face_dir == 1) ? 1 : 0),
                                            erf_j_hi + 2);

                        for (int erf_j = erf_j_lo; erf_j <= erf_j_hi; ++erf_j) {
                            for (int erf_i = erf_i_lo; erf_i <= erf_i_hi; ++erf_i) {
                                RealBox erf_bounds =
                                    (src_face_dir >= 0)
                                        ? FaceCellBoundsFromGeom(m_geom, erf_candidate_box,
                                                                 erf_i, erf_j, src_face_dir)
                                        : CellBoundsFromGeom(m_geom, erf_i, erf_j);
                                const Vector<Point2D> erf_poly = {
                                    {erf_bounds.lo(0), erf_bounds.lo(1)},
                                    {erf_bounds.hi(0), erf_bounds.lo(1)},
                                    {erf_bounds.hi(0), erf_bounds.hi(1)},
                                    {erf_bounds.lo(0), erf_bounds.hi(1)}
                                };
                                const Real area = PolygonArea(ClipPolygon(remora_poly, erf_poly));
                                if (area > 0.0) {
                                    overlaps.emplace_back(area, IntVect(erf_i, erf_j, 0));
                                }
                            }
                        }
                    }

                    // Nearest-valid-donor fill, ported from COAWST's check_weights
                    // (Lib/SCRIP_COAWST/scrip.f:343). A wet destination cell whose
                    // overlapping source cells are all land gets no stencil at all;
                    // COAWST detects exactly this case as
                    // (dst_mask_unlim - dst_mask) == 1 and hands the cell the
                    // nearest valid donor rather than leaving it empty. Searching
                    // the staged source cells keeps this rank-local: the staged
                    // region already covers this destination box's footprint plus
                    // margin, so the donor is in hand without a cross-rank search.
                    //
                    // Only wet holes are searched. Dry holes are both uninteresting
                    // and numerous (they are the land-shadowed majority), and
                    // scanning for them would cost orders of magnitude more.
                    //
                    // Off by default (driver.remap_donor_fill=withhold), because no
                    // production coupler donor-fills an unmapped cell and because the
                    // graft below is not conservative -- ESMF refuses to combine
                    // extrapolation with conservative regridding for exactly this
                    // reason. See ResolveDonorFillMode() for the full argument and
                    // for the warning the coawst mode prints. Under withhold this
                    // block is simply not entered: `overlaps` stays empty, the
                    // `continue` below leaves the -1 index sentinel and zero weights
                    // in place, NormalizeWeightsAndCheckConservation counts the cell
                    // uncovered, and ApplyConservativeRemap assigns that lane's
                    // fallback. Nothing else is bypassed.
                    //
                    // donor_box_has_erf_overlap excludes destination boxes whose
                    // grown geographic footprint intersected no ERF box at all.
                    // For those, erf_staged_box is not a real region: it is the
                    // one-cell placeholder BuildERFStagingPlan parks at
                    // erf_cell_domain.smallEnd(), because a BoxArray cannot hold
                    // an empty box. Without this test the search happily selects
                    // that corner and grafts a full-weight stencil, so every wet
                    // cell in such a box takes the ERF domain corner's
                    // atmospheric state -- and is then counted as covered rather
                    // than reported uncovered. A distance threshold cannot
                    // reliably reject it: two nearly-touching grids could place
                    // the corner within any radius. The recorded flag rejects it
                    // by construction instead.
                    if (m_remap_donor_fill_mode == DonorFillMode::Coawst &&
                        overlaps.empty() && m_remap_spherical && wet_arr &&
                        wet_arr(i,j,0) > Real(0.5) && donor_box_has_erf_overlap) {
                        Real dst_cx = 0.0, dst_cy = 0.0;
                        for (const auto& p : remora_poly) { dst_cx += p.x; dst_cy += p.y; }
                        dst_cx /= static_cast<Real>(remora_poly.size());
                        dst_cy /= static_cast<Real>(remora_poly.size());
                        // Longitude degrees are shorter than latitude degrees away
                        // from the equator; without this the nearest donor would be
                        // biased along the meridian.
                        const Real cos_lat = std::max(std::cos(dst_cy * Real(M_PI) / Real(180.0)),
                                                      Real(1.0e-6));

                        // Ordering is lexicographic on (distance^2, global source
                        // ID), not on distance alone. Exact ties are common rather
                        // than exotic on a structured grid -- a destination cell
                        // centered between two source cells produces bitwise-equal
                        // d2 -- and a bare `d2 < best_d2` resolves them by scan
                        // order, which starts at erf_staged_box.smallEnd() and so
                        // moves with the decomposition. ESMF and YAC both break
                        // nearest-neighbour ties by smallest global ID for exactly
                        // this reason; ESMF documents it as "if two points are
                        // equally close, then the point with the smallest sequence
                        // index is arbitrarily used".
                        //
                        // The `d2 == best_d2` test is exact on purpose. The staged
                        // node coordinates are exact NonLocalBC::ParallelCopy copies
                        // of ERF's node values, so d2 for a given (si,sj) is
                        // bitwise identical under any decomposition; ties here are
                        // genuine ties, and a tolerance would introduce a new
                        // ordering ambiguity rather than remove one.
                        Real best_d2 = std::numeric_limits<Real>::max();
                        long long best_gid = std::numeric_limits<long long>::max();
                        IntVect best_ij(-1,-1,0);
                        // The winning donor's own control volume. Kept from the
                        // loop rather than re-derived from the node array
                        // afterwards, because the search runs one index past
                        // erf_staged_box in the staggered direction for the face
                        // families: a raw erf_node(bi+1,bj+1) read at the winning
                        // index would step one node past erf_staged_node_box.
                        Vector<Point2D> best_poly;
                        // The donor index is stored in index_arr and dereferenced by
                        // ApplyConservativeRemap against this family's source layout,
                        // so it has to be selected in the same index space the main
                        // overlap loop uses. For the face families that is face
                        // space: one extra index in the active direction, face
                        // centroids rather than cell centroids, and the bracketing
                        // land rule. A cell index here would be legal but would
                        // reintroduce the half-cell displacement on donor-filled
                        // coastal cells, which is the defect this change removes.
                        const int donor_i_hi =
                            erf_staged_box.bigEnd(0) + ((src_face_dir == 0) ? 1 : 0);
                        const int donor_j_hi =
                            erf_staged_box.bigEnd(1) + ((src_face_dir == 1) ? 1 : 0);
                        for (int sj = erf_staged_box.smallEnd(1); sj <= donor_j_hi; ++sj) {
                            for (int si = erf_staged_box.smallEnd(0); si <= donor_i_hi; ++si) {
                                if (erf_src_is_land(si, sj)) { continue; }
                                const Vector<Point2D> src_poly = erf_src_poly(si, sj);
                                if (src_poly.empty()) { continue; }
                                Real sx = 0.0, sy = 0.0;
                                for (const auto& p : src_poly) { sx += p.x; sy += p.y; }
                                sx /= static_cast<Real>(src_poly.size());
                                sy /= static_cast<Real>(src_poly.size());
                                // Clip the candidate set to G(B). This is what
                                // takes ERF's decomposition out of the answer,
                                // and it is a proof rather than a heuristic: any
                                // ERF cell whose lon/lat lies inside G(B) belongs
                                // to some ERF box whose own geographic bounding
                                // box intersects G(B), and BuildERFStagingPlan
                                // minBox-unions every such ERF box into
                                // erf_staged_box. So the surviving candidates are
                                // exactly the in-G(B) source cells under every ERF
                                // box layout, whatever NProcs() and
                                // amr.max_grid_size* do to ERF. What it does not
                                // remove is dependence on REMORA's decomposition:
                                // G(B) is itself a per-REMORA-box quantity, so a
                                // destination cell near the edge of G(B) can see a
                                // donor accepted under one REMORA box layout and
                                // rejected under another. This mode is bounded,
                                // not proven rank-invariant.
                                if (sx < donor_box_bbox.xlo || sx > donor_box_bbox.xhi ||
                                    sy < donor_box_bbox.ylo || sy > donor_box_bbox.yhi) {
                                    continue;
                                }
                                const Real dx = (sx - dst_cx) * cos_lat;
                                const Real dy = (sy - dst_cy);
                                const Real d2 = dx*dx + dy*dy;
                                const long long gid = donor_global_id(si, sj);
                                if (d2 < best_d2 || (d2 == best_d2 && gid < best_gid)) {
                                    best_d2 = d2;
                                    best_gid = gid;
                                    best_ij = IntVect(si,sj,0);
                                    best_poly = src_poly;
                                }
                            }
                        }

                        if (best_ij[0] >= 0) {
                            // Acceptance radius, measured against the diagonal of
                            // the donor's own control volume in the same
                            // cos(latitude)-corrected metric the search used.
                            // Taken from the polygon the loop already built, not
                            // from a raw erf_node(bi+1,bj+1) read, which would
                            // step one node past erf_staged_node_box for the face
                            // families.
                            //
                            // Unlike the geographic clip above this is a
                            // heuristic, and it is a deliberate departure from
                            // precedent rather than borrowed practice: no
                            // production coupler bounds a nearest-source fill by
                            // distance. ESMF's NEAREST_STOD is explicitly
                            // unbounded -- it queries its process directory with
                            // an infinite box when no local candidate exists, so
                            // it always finds the globally nearest source -- and
                            // creep fill bounds by iteration count instead. The
                            // default of 3 diagonals is uncalibrated.
                            //
                            // Applied uniformly, not conditioned on whether this
                            // destination box was partially or fully uncovered,
                            // because the rank dependence exists in both cases.
                            // (Zero-overlap boxes are already excluded above, so
                            // in practice this governs the partially covered
                            // case.)
                            const PlanarBBox donor_bb = BBoxOfPolygon(best_poly);
                            const Real ddx = (donor_bb.xhi - donor_bb.xlo) * cos_lat;
                            const Real ddy = (donor_bb.yhi - donor_bb.ylo);
                            const Real donor_diag2 = ddx*ddx + ddy*ddy;
                            const Real r = m_remap_donor_max_diagonals;
                            if (best_d2 > r * r * donor_diag2) {
                                // Rejected. `overlaps` stays empty, the continue
                                // below leaves the -1 index sentinel, and the cell
                                // takes the same per-lane ERF-side fallback it
                                // would under withhold. There is no third outcome.
                                best_ij = IntVect(-1,-1,0);
                            }
                        }

                        if (best_ij[0] >= 0) {
                            // Full weight to the donor: after fracarea normalization
                            // this cell takes the donor's value exactly.
                            overlaps.emplace_back(remora_area, best_ij);
                            ++donor_filled;
                        }
                    }

                    if (overlaps.empty()) { continue; }

                    max_stencil_used = std::max(max_stencil_used, static_cast<int>(overlaps.size()));
                    if (static_cast<int>(overlaps.size()) > MAX_STENCIL_SIZE) {
                        dropped_overlaps += static_cast<amrex::Long>(overlaps.size()) - MAX_STENCIL_SIZE;
                        std::partial_sort(overlaps.begin(), overlaps.begin() + MAX_STENCIL_SIZE,
                                          overlaps.end(),
                                          [] (const std::pair<Real,IntVect>& a,
                                              const std::pair<Real,IntVect>& b) {
                                              return a.first > b.first;
                                          });
                        overlaps.resize(MAX_STENCIL_SIZE);
                    }

                    for (int m = 0; m < static_cast<int>(overlaps.size()); ++m) {
                        weight_arr(i, j, 0, m) = overlaps[m].first / remora_area;
                        index_arr(i, j, 0, m * 3    ) = overlaps[m].second[0];
                        index_arr(i, j, 0, m * 3 + 1) = overlaps[m].second[1];
                        index_arr(i, j, 0, m * 3 + 2) = 0;
                    }
                }
            }
        }
        static const std::array<const char*, kNumWeightFamilies> kFamilyLabel = {
            "A2O scalar/rho", "A2O tau_x/u-face", "A2O tau_y/v-face", "A2O wind->rho"
        };
        NormalizeWeightsAndCheckConservation(*weight_mf, *index_mf, family_mask, kFamilyLabel[f]);
    }

    ParallelDescriptor::ReduceIntMax(max_stencil_used);
    ParallelDescriptor::ReduceLongSum(dropped_overlaps);
    ParallelDescriptor::ReduceLongSum(donor_filled);
    if (donor_filled > 0) {
        amrex::Print() << "[driver] remap: " << donor_filled
                       << " wet destination cells had no valid source overlap and were filled from "
                       << "their nearest non-land source cell (COAWST check_weights behaviour). "
                       << "These carry a neighbour's atmospheric state rather than their own.\n";
    }
    if (dropped_overlaps > 0) {
        amrex::Print() << "[driver] WARNING BuildRemapWeights: " << dropped_overlaps
                       << " source overlaps dropped because a destination cell exceeded "
                       << "MAX_STENCIL_SIZE=" << MAX_STENCIL_SIZE
                       << " (max overlaps seen = " << max_stencil_used
                       << "); weights renormalized over the largest-area overlaps kept.\n";
    } else if (m_remap_verbose >= 1) {
        amrex::Print() << "[driver] remap stencils: max overlaps per destination cell = "
                       << max_stencil_used << " of " << MAX_STENCIL_SIZE << "\n";
    }
}

void
ERFRemoraMultiBlockContainer::BuildOceanToAtmosWeights()
{
    using namespace amrex;

    // Reset and compute O2A weight and index caches if active
    if (m_use_two_way_coupling && m_weight_o2a_mf && m_index_o2a_mf) {
        m_weight_o2a_mf->setVal(0.0);
        m_index_o2a_mf->setVal(-1);
        m_valid_o2a_mf->setVal(0);

        int max_stencil_used = 0;
        amrex::Long dropped_overlaps = 0;
        Vector<int> candidates;
        Vector<std::pair<Real,IntVect>> overlaps;

        // Why each uncovered destination cell is uncovered. All three end as
        // valid = 0 and all three correctly leave ERF's own SST in place, but
        // they mean different things and only one of them is ever suspicious:
        //
        //   outside   no ocean cell overlaps at all -- the ocean footprint does
        //             not reach here. Expected around the edge of the ocean grid.
        //   masked    ocean cells do overlap, but every one of them is dry. That
        //             is REMORA calling land what ERF calls water, so it is a
        //             mask disagreement: legitimate along a real coastline,
        //             but also exactly what a grid misalignment or a stale mask
        //             looks like. Collapsing it into the other two is how such a
        //             bug hides, since a silent withhold reads as correct.
        //   thin      wet overlap exists but falls under remap_min_coverage_frac.
        //
        // ESMF draws the same first distinction in ESMF_REGRIDSTATUS, which
        // separates OUTSIDE from SRCMASKED for this reason.
        amrex::Long uncovered_outside = 0;
        amrex::Long uncovered_masked  = 0;
        amrex::Long uncovered_thin    = 0;

        // O2A deliberately does not donor-fill uncovered ERF cells. They keep
        // ERF's own wrflowinp SST, which is better than a distant ocean donor.
        if (m_remap_spherical) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                m_erf_node_coords != nullptr,
                "BuildOceanToAtmosWeights spherical path requires ERF corner coordinates.");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                m_remora_staging_plan != nullptr &&
                m_remora_staging_plan->node_mf != nullptr &&
                m_remora_staging_plan->wet_mf != nullptr,
                "BuildOceanToAtmosWeights spherical path requires the staged REMORA corner mesh and wet mask.");

            for (MFIter mfi(*m_weight_o2a_mf); mfi.isValid(); ++mfi) {
                Box bx = mfi.validbox();
                auto weight_arr = m_weight_o2a_mf->array(mfi);
                auto index_arr = m_index_o2a_mf->array(mfi);
                auto valid_arr = m_valid_o2a_mf->array(mfi);
                auto const& erf_node = m_erf_node_coords->const_array(mfi);
                auto const& remora_node = m_remora_staging_plan->node_mf->const_array(mfi);
                auto const& remora_wet = m_remora_staging_plan->wet_mf->const_array(mfi);
                const Box remora_staged_box = m_remora_staging_plan->wet_mf->boxArray()[mfi.index()];

                SourceSearchGrid search;
                search.define(remora_staged_box,
                              Array4<const Real>(remora_node, 0),
                              Array4<const Real>(remora_node, 1),
                              m_remap_bin_cells);

                for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                    for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                        const Vector<Point2D> erf_poly =
                            ExtractCellQuadFromCoords(i, j,
                                                      Array4<const Real>(erf_node, 0),
                                                      Array4<const Real>(erf_node, 1));
                        const Real dst_area = PolygonArea(erf_poly);
                        if (dst_area <= 0.0) continue;

                        const PlanarBBox erf_bbox = BBoxOfPolygon(erf_poly);
                        search.query(erf_bbox, candidates);

                        overlaps.clear();
                        // Fraction of this ERF cell covered by wet REMORA rho cells.
                        // The raw intersect_area/dst_area weights are the destination
                        // fraction; after normalization their sum is one and cannot
                        // drive the uncovered-cell threshold.
                        Real covered_frac = 0.0;
                        // Dry cells are clipped as well as wet ones, but only to
                        // learn whether the ocean grid reaches this cell at all.
                        // They contribute nothing to coverage or to the stencil.
                        // Skipping them before the clip would be cheaper but would
                        // make "outside the ocean footprint" and "ocean says land"
                        // indistinguishable, which is the pair worth telling apart.
                        bool any_geometric_overlap = false;
                        for (const int item : candidates) {
                            const IntVect remora_ij = search.itemIndex(item);
                            const int src_i = remora_ij[0];
                            const int src_j = remora_ij[1];

                            const Vector<Point2D> remora_poly =
                                ExtractCellQuadFromCoords(src_i, src_j,
                                                          Array4<const Real>(remora_node, 0),
                                                          Array4<const Real>(remora_node, 1));
                            const Vector<Point2D> clipped_poly = ClipPolygon(erf_poly, remora_poly);
                            const Real intersect_area = PolygonArea(clipped_poly);
                            if (intersect_area <= 0.0) { continue; }

                            any_geometric_overlap = true;
                            if (remora_wet(src_i, src_j, 0) <= Real(0.5)) { continue; }

                            covered_frac += intersect_area / dst_area;
                            overlaps.emplace_back(intersect_area, IntVect(src_i, src_j, 0));
                        }

                        if (overlaps.empty()) {
                            if (any_geometric_overlap) { ++uncovered_masked; }
                            else                       { ++uncovered_outside; }
                            continue;
                        }
                        if (covered_frac < m_remap_min_coverage_frac) { ++uncovered_thin; }

                        // Coverage for this destination cell, thresholded on area
                        // rather than on whether any overlap exists at all.
                        //
                        // A sliver overlap is worse than no overlap. The weights are
                        // normalized below before REMORA's raw weighted sum, but
                        // normalization cannot turn 1e-6 destination coverage into
                        // a trustworthy ocean answer. Treating it as uncovered keeps
                        // ERF's own lower-boundary SST.
                        // Thresholds of this kind are standard: CESM drops mapping
                        // below a destination fraction of 0.001, OASIS/CERFACS at
                        // 1/1000.
                        //
                        // Wet REMORA coverage only: dry rho cells do not contribute
                        // to coverage or to the SST stencil. ERF's lmask is still
                        // applied on the destination side after the remap.
                        // The overlap count is tested as well as the fraction: with
                        // remap_min_coverage_frac set to 0, a cell no ocean cell
                        // touches has covered_frac == 0, and 0 >= 0 would otherwise
                        // report every uncovered cell as covered.
                        valid_arr(i, j, 0) =
                            (!overlaps.empty() && covered_frac >= m_remap_min_coverage_frac) ? 1 : 0;

                        max_stencil_used = std::max(max_stencil_used, static_cast<int>(overlaps.size()));
                        if (static_cast<int>(overlaps.size()) > MAX_STENCIL_SIZE) {
                            dropped_overlaps += static_cast<amrex::Long>(overlaps.size()) - MAX_STENCIL_SIZE;
                            std::partial_sort(overlaps.begin(), overlaps.begin() + MAX_STENCIL_SIZE,
                                              overlaps.end(),
                                              [] (const std::pair<Real,IntVect>& a,
                                                  const std::pair<Real,IntVect>& b) {
                                                  return a.first > b.first;
                                              });
                            overlaps.resize(MAX_STENCIL_SIZE);
                        }

                        for (int m = 0; m < static_cast<int>(overlaps.size()); ++m) {
                            weight_arr(i, j, 0, m) = overlaps[m].first / dst_area;
                            index_arr(i, j, 0, m * 3    ) = overlaps[m].second[0];
                            index_arr(i, j, 0, m * 3 + 1) = overlaps[m].second[1];
                            index_arr(i, j, 0, m * 3 + 2) = 0;
                        }
                    }
                }
            }
        } else {
            // Cartesian O2A still uses a replicated rho-mask lookup because it has
            // no staged REMORA source layout. Spherical O2A consumes
            // m_remora_staging_plan->wet_mf instead.
            const MultiFab* remora_mskr = nullptr;
            if (m_remora != nullptr) {
                REMORAGetMasks(m_remora, &remora_mskr, nullptr, nullptr);
            }

            const Box remora_cell_domain = m_geom_ocn.Domain();
            const bool use_remora_mskr =
                remora_mskr != nullptr &&
                remora_mskr->boxArray().ixType() == m_ba_ocn.ixType();
            const Vector<Real> remora_mskr_lookup =
                use_remora_mskr ? BuildReplicatedMaskLookup(*remora_mskr, remora_cell_domain)
                                : Vector<Real>{};
            auto remora_rho_is_wet = [&] (int src_i, int src_j) -> bool {
                if (!use_remora_mskr) { return true; }
                const int src_k = remora_cell_domain.smallEnd(2);
                return remora_mskr_lookup[FlatIndexInBox(remora_cell_domain, src_i, src_j, src_k)]
                       > Real(0.5);
            };

            for (MFIter mfi(*m_weight_o2a_mf); mfi.isValid(); ++mfi) {
                Box bx = mfi.validbox();
                auto weight_arr = m_weight_o2a_mf->array(mfi);
                auto index_arr = m_index_o2a_mf->array(mfi);
                auto valid_arr = m_valid_o2a_mf->array(mfi);

                for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                    for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                        RealBox erf_bounds = CellBoundsFromGeom(m_geom, i, j);
                        Vector<Point2D> erf_poly = {
                            {erf_bounds.lo(0), erf_bounds.lo(1)},
                            {erf_bounds.hi(0), erf_bounds.lo(1)},
                            {erf_bounds.hi(0), erf_bounds.hi(1)},
                            {erf_bounds.lo(0), erf_bounds.hi(1)}
                        };

                        Real dst_area = PolygonArea(erf_poly);
                        if (dst_area <= 0.0) continue;

                        RealBox erf_bbox = erf_bounds;
                        Box remora_candidate_box = m_geom_ocn.Domain();
                        int ocn_i_lo = std::floor((erf_bbox.lo(0) - m_ocn_rb.lo(0)) / (m_ocn_rb.hi(0) - m_ocn_rb.lo(0)) * m_ocn_n_cell[0]);
                        int ocn_i_hi = std::ceil((erf_bbox.hi(0) - m_ocn_rb.lo(0)) / (m_ocn_rb.hi(0) - m_ocn_rb.lo(0)) * m_ocn_n_cell[0]);
                        int ocn_j_lo = std::floor((erf_bbox.lo(1) - m_ocn_rb.lo(1)) / (m_ocn_rb.hi(1) - m_ocn_rb.lo(1)) * m_ocn_n_cell[1]);
                        int ocn_j_hi = std::ceil((erf_bbox.hi(1) - m_ocn_rb.lo(1)) / (m_ocn_rb.hi(1) - m_ocn_rb.lo(1)) * m_ocn_n_cell[1]);

                        ocn_i_lo = std::max(remora_candidate_box.smallEnd(0), ocn_i_lo - 1);
                        ocn_i_hi = std::min(remora_candidate_box.bigEnd(0), ocn_i_hi + 1);
                        ocn_j_lo = std::max(remora_candidate_box.smallEnd(1), ocn_j_lo - 1);
                        ocn_j_hi = std::min(remora_candidate_box.bigEnd(1), ocn_j_hi + 1);

                        overlaps.clear();
                        // Fraction of this ERF cell covered by wet REMORA rho cells.
                        // The raw intersect_area/dst_area weights are the destination
                        // fraction; after normalization their sum is one and cannot
                        // drive the uncovered-cell threshold.
                        Real covered_frac = 0.0;
                        // See the spherical branch: dry cells are clipped too,
                        // solely to distinguish "the ocean does not reach here"
                        // from "the ocean calls this land".
                        bool any_geometric_overlap = false;
                        for (int src_j = ocn_j_lo; src_j <= ocn_j_hi; ++src_j) {
                            for (int src_i = ocn_i_lo; src_i <= ocn_i_hi; ++src_i) {
                                RealBox remora_cell_bounds = CellBoundsFromGeom(m_geom_ocn, src_i, src_j);
                                Vector<Point2D> remora_poly = {
                                    {remora_cell_bounds.lo(0), remora_cell_bounds.lo(1)},
                                    {remora_cell_bounds.hi(0), remora_cell_bounds.lo(1)},
                                    {remora_cell_bounds.hi(0), remora_cell_bounds.hi(1)},
                                    {remora_cell_bounds.lo(0), remora_cell_bounds.hi(1)}
                                };

                                Vector<Point2D> clipped_poly = ClipPolygon(erf_poly, remora_poly);
                                Real intersect_area = PolygonArea(clipped_poly);
                                if (intersect_area <= 0.0) { continue; }

                                any_geometric_overlap = true;
                                if (!remora_rho_is_wet(src_i, src_j)) { continue; }

                                covered_frac += intersect_area / dst_area;
                                overlaps.emplace_back(intersect_area, IntVect(src_i, src_j, 0));
                            }
                        }

                        if (overlaps.empty()) {
                            if (any_geometric_overlap) { ++uncovered_masked; }
                            else                       { ++uncovered_outside; }
                            continue;
                        }
                        if (covered_frac < m_remap_min_coverage_frac) { ++uncovered_thin; }

                        // Coverage for this destination cell, thresholded on area
                        // rather than on whether any overlap exists at all.
                        //
                        // A sliver overlap is worse than no overlap. The weights are
                        // normalized below before REMORA's raw weighted sum, but
                        // normalization cannot turn 1e-6 destination coverage into
                        // a trustworthy ocean answer. Treating it as uncovered keeps
                        // ERF's own lower-boundary SST.
                        // Thresholds of this kind are standard: CESM drops mapping
                        // below a destination fraction of 0.001, OASIS/CERFACS at
                        // 1/1000.
                        //
                        // Wet REMORA coverage only: dry rho cells do not contribute
                        // to coverage or to the SST stencil. ERF's lmask is still
                        // applied on the destination side after the remap.
                        // The overlap count is tested as well as the fraction: with
                        // remap_min_coverage_frac set to 0, a cell no ocean cell
                        // touches has covered_frac == 0, and 0 >= 0 would otherwise
                        // report every uncovered cell as covered.
                        valid_arr(i, j, 0) =
                            (!overlaps.empty() && covered_frac >= m_remap_min_coverage_frac) ? 1 : 0;

                        max_stencil_used = std::max(max_stencil_used, static_cast<int>(overlaps.size()));
                        if (static_cast<int>(overlaps.size()) > MAX_STENCIL_SIZE) {
                            dropped_overlaps += static_cast<amrex::Long>(overlaps.size()) - MAX_STENCIL_SIZE;
                            std::partial_sort(overlaps.begin(), overlaps.begin() + MAX_STENCIL_SIZE,
                                              overlaps.end(),
                                              [] (const std::pair<Real,IntVect>& a,
                                                  const std::pair<Real,IntVect>& b) {
                                                  return a.first > b.first;
                                              });
                            overlaps.resize(MAX_STENCIL_SIZE);
                        }

                        for (int m = 0; m < static_cast<int>(overlaps.size()); ++m) {
                            weight_arr(i, j, 0, m) = overlaps[m].first / dst_area;
                            index_arr(i, j, 0, m * 3    ) = overlaps[m].second[0];
                            index_arr(i, j, 0, m * 3 + 1) = overlaps[m].second[1];
                            index_arr(i, j, 0, m * 3 + 2) = 0;
                        }
                    }
                }
            }
        }

        NormalizeWeightsAndCheckConservation(*m_weight_o2a_mf, *m_index_o2a_mf, nullptr, "O2A SST");

        ParallelDescriptor::ReduceLongSum(uncovered_outside);
        ParallelDescriptor::ReduceLongSum(uncovered_masked);
        ParallelDescriptor::ReduceLongSum(uncovered_thin);

        // Reported unconditionally: an uncovered cell is invisible in the SST
        // field, because ERF simply keeps its own value there, so the only way a
        // mask disagreement surfaces is if the split is written down. A large
        // masked count against a short coastline is the signature of a grid
        // misalignment or a stale ocean mask rather than of correct withholding.
        amrex::Print() << "[driver] O2A coverage: uncovered "
                       << uncovered_outside << " outside ocean footprint, "
                       << uncovered_masked  << " source-masked, "
                       << uncovered_thin    << " below coverage threshold ("
                       << m_remap_min_coverage_frac << ")\n";

        ParallelDescriptor::ReduceIntMax(max_stencil_used);
        ParallelDescriptor::ReduceLongSum(dropped_overlaps);
        if (dropped_overlaps > 0) {
            amrex::Print() << "[driver] WARNING BuildOceanToAtmosWeights: " << dropped_overlaps
                           << " ocean overlaps dropped because a destination cell exceeded "
                           << "MAX_STENCIL_SIZE=" << MAX_STENCIL_SIZE
                           << " (max overlaps seen = " << max_stencil_used
                           << "); weights normalized over the largest-area overlaps kept.\n";
        } else if (m_remap_verbose >= 1) {
            amrex::Print() << "[driver] O2A remap stencils: max overlaps per destination cell = "
                           << max_stencil_used << " of " << MAX_STENCIL_SIZE << "\n";
        }
    }
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
ERFRemoraMultiBlockContainer::FillAtmosToOceanGhosts ()
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_geom_ocn.isPeriodic(0) == static_cast<bool>(m_is_periodic[0]) &&
        m_geom_ocn.isPeriodic(1) == static_cast<bool>(m_is_periodic[1]) &&
        m_geom_ocn.isPeriodic(2) == static_cast<bool>(m_is_periodic[2]),
        "Driver ocean-slab geometry lost the requested periodic flags.");

    for (auto* mf : m_atmos_to_ocean) {
        if (mf != nullptr && mf->nGrowVect().max() > 0) {
            mf->FillBoundary(m_geom_ocn.periodicity());
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
        if (fields[idx] == nullptr) {
            amrex::Print() << "  " << names[idx] << ": withheld (not supplied by ERF)\n";
            continue;
        }
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
            // A null lane is deliberate, not a defect: it marks a field ERF cannot
            // fill, withheld so the receiving solver keeps its own forcing for it
            // (see SelectAtmosToOceanActiveView). Report it and move on.
            if (amrex::ParallelDescriptor::IOProcessor()) {
                amrex::Print() << label
                               << " ocean_step=" << ocean_step
                               << " time=" << time
                               << " field=" << names[idx]
                               << " withheld (not supplied by ERF)\n";
            }
            continue;
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
    if (m_use_two_way_coupling) {
        ValidateFields("Exchange validation ocean->atmos",
                       m_ocean_to_atmos, coupling::OceanToAtmosNames(), time, ocean_step);
    }
    ValidateFields("Exchange validation atmos->ocean",
                   m_atmos_to_ocean_active, m_atmos_to_ocean_active_names, time, ocean_step);

    if (m_use_two_way_coupling) {
        WriteExchangePlotfile(m_plot_file_exchange_o2a,
                              m_ocean_to_atmos,
                              coupling::OceanToAtmosNames(),
                              m_geom,
                              time,
                              ocean_step);
    }
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
    amrex::BoxArray plot_ba = amrex::convert(fields[0]->boxArray(), amrex::IntVect::TheCellVector());
    amrex::MultiFab plot_mf(plot_ba, fields[0]->DistributionMap(), ncomp, 0);
    for (int comp = 0; comp < ncomp; ++comp) {
        if (fields[comp] == nullptr) {
            amrex::Abort("Cannot write exchange plotfile with a null field.");
        }
        if (amrex::convert(fields[comp]->boxArray(), amrex::IntVect::TheCellVector()) != plot_ba) {
            amrex::Abort("Exchange plotfile writer requires a common enclosed-cell slab layout.");
        }
        CopyFieldToCellCenteredPlot(plot_mf, comp, *fields[comp]);
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
            // ERF land mask on the destination (ERF) layout: zeroes out
            // land cells in the remapped ocean->atmos SST. Wet/dry masking
            // only, no vector rotation.
            const amrex::iMultiFab* erf_is_land = nullptr;
            ERFGetLandMask(m_erf, &erf_is_land);
            REMORAPackSurfaceState(m_remora, state, time,
                                  m_weight_o2a_mf.get(), m_index_o2a_mf.get(), MAX_STENCIL_SIZE,
                                  erf_is_land);

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
            // Hand the coverage flag over with the SST. Without it ERF cannot
            // tell an uncovered cell's zero fill from a genuine ocean value,
            // and would overwrite its own wrflowinp SST with 0 K wherever the
            // ocean grid does not reach the atmosphere footprint.
            ERFApplyOceanSurfaceState(m_erf, state, time, m_valid_o2a_mf.get());
        },
        [this](amrex::Vector<amrex::MultiFab*>& states, amrex::Real time) {
            BL_PROFILE("Driver::Pack_Atmos");
            amrex::Vector<const amrex::MultiFab*> weight_mf_by_family(kNumWeightFamilies);
            amrex::Vector<const amrex::iMultiFab*> index_mf_by_family(kNumWeightFamilies);
            for (int f = 0; f < kNumWeightFamilies; ++f) {
                weight_mf_by_family[f] = m_weight_mf_by_family[f].get();
                index_mf_by_family[f] = m_index_mf_by_family[f].get();
            }
            // REMORA wet/dry masks: mskr masks the scalar destination lanes,
            // msku/mskv mask TauX/TauY. Wet/dry masking only, no vector rotation.
            const amrex::MultiFab* mskr = nullptr;
            const amrex::MultiFab* msku = nullptr;
            const amrex::MultiFab* mskv = nullptr;
            REMORAGetMasks(m_remora, &mskr, &msku, &mskv);
            ERFPackAtmosphericStates(m_erf, states, time,
                                    weight_mf_by_family, index_mf_by_family, MAX_STENCIL_SIZE,
                                    mskr, msku, mskv);
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
                    if (is_cloud_lane) {
                        // min/max are collective (MPI_Allreduce), so they are issued
                        // unconditionally here, exactly as the finite check above does.
                        // cloud_warn_count is a process-local static; gating the
                        // reductions on it is what let the IO rank stop entering
                        // collectives that every other rank still executed. The counter
                        // now limits printing only.
                        const amrex::Real cmin = mf->min(0);
                        const amrex::Real cmax = mf->max(0);
                        if ((cmin < -cloud_tol || cmax > amrex::Real(1.0) + cloud_tol) &&
                            cloud_warn_count < cloud_warn_limit) {
                            ++cloud_warn_count;
                            if (amrex::ParallelDescriptor::IOProcessor()) {
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

    // The view actually handed to REMORA. It is the active view with any lane
    // nulled that ERF cannot fill, because REMORA decides per lane whether the
    // driver is authoritative purely from the pointer being non-null
    // (REMORA.cpp:1220 reads NetCDF shortwave only when the driver did not supply
    // it). Handing over a lane of zeros therefore does not mean "no contribution",
    // it means "the atmosphere says zero" - and silently disables REMORA's own
    // NetCDF radiation. Nulling the lane instead lets that fallback stand.
    //
    // State contract only. The flux contract has no such fallback: REMORA's
    // ApplyAtmosphericFluxes requires all eight lanes non-null and skips the whole
    // application if any is missing, so withholding there would disable every flux
    // rather than one. In flux mode a zero radiative flux is at least a truthful
    // statement that the atmosphere computed none.
    m_atmos_to_ocean_apply_view = m_atmos_to_ocean_active;
    if (active_contract == coupling::AtmosToOceanContract::StatePassing &&
        m_erf != nullptr) {
        const bool has_radiation = ERFHasRadiation(m_erf);
        // ERF fills the Cloud lane only when moisture is active and carries cloud
        // liquid or ice; otherwise the lane holds just the driver prefill, and
        // passing it would mark it authoritative in REMORA and suppress REMORA's
        // own cloud forcing (remora.cloud_type).
        const bool has_cloud_water = ERFHasCloudWater(m_erf);
        for (int local_idx = 0;
             local_idx < static_cast<int>(m_atmos_to_ocean_active_indices.size()); ++local_idx) {
            const int global_idx = m_atmos_to_ocean_active_indices[local_idx];
            if (!has_radiation &&
                (global_idx == coupling::SWrad || global_idx == coupling::LWrad)) {
                m_atmos_to_ocean_apply_view[local_idx] = nullptr;
            }
            if (!has_cloud_water && global_idx == coupling::Cloud) {
                m_atmos_to_ocean_apply_view[local_idx] = nullptr;
            }
        }
        if (!has_cloud_water && !m_warned_no_erf_cloud_water) {
            m_warned_no_erf_cloud_water = true;
            amrex::Print() << "[driver] ERF has no cloud condensate (no moisture model, or one "
                              "without cloud liquid/ice), so the Cloud atmos->ocean lane is "
                              "withheld rather than passed as the constant driver prefill. REMORA "
                              "keeps whatever cloud forcing its own configuration provides "
                              "(remora.cloud_type).\n";
        }
        if (!has_radiation && !m_warned_no_erf_radiation) {
            m_warned_no_erf_radiation = true;
            amrex::Print() << "[driver] ERF reports no active radiation model, so the SWrad/LWrad "
                              "atmos->ocean lanes are withheld rather than passed as zeros. REMORA "
                              "keeps whatever radiative forcing its own configuration provides "
                              "(remora.swrad_type / remora.lwrad_type). Enable ERF radiation to "
                              "couple radiation from the atmosphere instead.\n";
        }
    }
}

void
ERFRemoraMultiBlockContainer::AdvanceBlocksWithCallbacks (const MutableStateFn& pack_surface_state,
                                                           const ConstStateFn& apply_ocean_surface_state,
                                                           const MutableStateFn& pack_atmos_states,
                                                           const ConstStateFn& apply_atmos_states,
                                                           int ocean_steps)
{
    BL_PROFILE("Driver::Build_Comm_Metadata");
    if (m_use_two_way_coupling) {
        PrintDiagnostics("Initial ocean->atmos cache", m_ocean_to_atmos, coupling::OceanToAtmosNames());
    }
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
        const amrex::Real ocean_time = m_ocean_time_offset + m_atm_tick * m_erf_dt;
        const long long global_ocean_step =
            m_atm_tick / static_cast<long long>(m_atm_substeps_per_ocean_step);
        SelectAtmosToOceanActiveView(global_ocean_step);
        const bool use_state_contract =
            (ActiveAtmosToOceanContract(global_ocean_step) == coupling::AtmosToOceanContract::StatePassing);
        if (global_ocean_step == 0 && amrex::ParallelDescriptor::IOProcessor()) {
            amrex::Print() << "Bootstrap atm2ocn contract at ocean_step=0: "
                           << (use_state_contract ? "state" : "flux")
                           << " (configured mode="
                           << coupling::AtmosToOceanContractName(m_atmos_to_ocean_contract)
                           << ")\n";
        }
        DriverPrintLine(
            m_driver_verbose,
            std::string("ocean_step=") + std::to_string(global_ocean_step) +
                " set ERF atm2ocn contract to " + (use_state_contract ? "state" : "flux"));
        SetERFDriverAtmosToOceanStateContract(m_erf, use_state_contract);
        DriverPrintLine(
            m_driver_verbose,
            std::string("ocean_step=") + std::to_string(global_ocean_step) +
                " set REMORA atm2ocn contract to " + (use_state_contract ? "state" : "flux"));
        SetREMORADriverAtmosToOceanStateContract(m_remora, use_state_contract);
        if (m_use_two_way_coupling) {
            DriverPrintLine(
                m_driver_verbose,
                std::string("ocean_step=") + std::to_string(global_ocean_step) +
                    " pack ocean->atmos surface state from REMORA at time=" +
                    std::to_string(ocean_time));
            pack_surface_state(m_ocean_to_atmos, ocean_time);
            FillOceanToAtmosGhosts();
        }
        DriverPrintLine(
            m_driver_verbose,
            std::string("ocean_step=") + std::to_string(global_ocean_step) +
                " pack atmos->ocean " + (use_state_contract ? "state" : "flux") +
                " from ERF at time=" + std::to_string(ocean_time));
        pack_atmos_states(m_atmos_to_ocean_active, ocean_time);
        if (m_plot_int_exchange > 0 &&
            (global_ocean_step % static_cast<long long>(m_plot_int_exchange)) == 0) {
            WriteExchangePlotfiles(ocean_time, global_ocean_step);
        }
        FillAtmosToOceanGhosts();
        DriverPrintLine(
            m_driver_verbose,
            std::string("ocean_step=") + std::to_string(global_ocean_step) +
                " apply atmos->ocean " + (use_state_contract ? "state" : "flux") +
                " to REMORA at time=" + std::to_string(ocean_time));
        // Masked view: lanes ERF cannot fill are nulled so REMORA keeps its own
        // forcing for them (see SelectAtmosToOceanActiveView).
        apply_atmos_states(m_atmos_to_ocean_apply_view, ocean_time);


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

            BL_PROFILE("Driver::Time_Interpolate");
            if (m_use_two_way_coupling) {
                DriverPrintLine(
                    m_driver_verbose,
                    std::string("atm_tick=") + std::to_string(m_atm_tick) +
                        " apply ocean->atmos surface state to ERF at time=" +
                        std::to_string(atm_time));
                apply_ocean_surface_state(m_ocean_to_atmos, atm_time);
            }
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

        if (m_use_two_way_coupling) {
            PrintDiagnostics("Ocean->atmos cache", m_ocean_to_atmos, coupling::OceanToAtmosNames());
        }
        PrintDiagnostics("Atmos->ocean cache", m_atmos_to_ocean_active, m_atmos_to_ocean_active_names);
    }
    REMORAWriteAtFinalTime(m_remora);
    ERFWriteAtFinalTime(m_erf);

}
} // namespace erf_remora
