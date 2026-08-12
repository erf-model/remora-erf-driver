#include "ERFRemoraMultiBlockContainer.H"
#include "StandaloneSolverRuns.H"

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Math.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_String.H>

#include <algorithm>
#include <cmath>
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
    m_atm_tick = 0;
    DriverPrintLine(m_driver_verbose, "resetting driver atmosphere tick to 0");
    AllocateCaches();
    BuildRemapWeights();
    DriverPrintLine(m_driver_verbose, "selecting bootstrap active atm2ocn cache view for ocean_step=0");
    SelectAtmosToOceanActiveView(/*ocean_step=*/0);
    if (m_erf != nullptr && m_remora != nullptr) {
        if (m_use_two_way_coupling && !m_ocean_to_atmos.empty()) {
            DriverPrintLine(m_driver_verbose, "packing REMORA surface state into ocean->atmos bootstrap cache");
            const amrex::iMultiFab* erf_is_land_bootstrap = nullptr;
            ERFGetLandMask(m_erf, &erf_is_land_bootstrap);
            REMORAPackSurfaceState(m_remora, m_ocean_to_atmos, amrex::Real(0.0),
                                  m_weight_o2a_mf.get(), m_index_o2a_mf.get(), MAX_STENCIL_SIZE,
                                  erf_is_land_bootstrap);
            FillOceanToAtmosGhosts();
        }
    } else {
        DriverPrintLine(m_driver_verbose, "synthetic-only path: filling synthetic coupling cache values");
        FillSyntheticData();
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
    if (m_erf != nullptr && m_remora != nullptr) {
        DriverPrintLine(m_driver_verbose, "allocating caches from solver-reported post-InitData layouts");
        AllocateCachesFromSolverLayouts();
    } else {
        DriverPrintLine(m_driver_verbose, "allocating synthetic fallback caches from driver-derived layouts");
        AllocateSyntheticCaches();
    }
}

void
ERFRemoraMultiBlockContainer::AllocateSyntheticCaches ()
{
    // Allocate a superset layer vector once (legacy state + optional flux layers).
    // Runtime coupling modes can then swap active layer views by index without
    // changing memory ownership or reallocating cache slabs.
    m_atmos_to_ocean_storage.resize(coupling::NumAtmosToOceanLayers);
    for (int idx = 0; idx < static_cast<int>(m_atmos_to_ocean_storage.size()); ++idx) {
        auto& mf_ptr = m_atmos_to_ocean_storage[idx];
        mf_ptr = std::make_unique<amrex::MultiFab>(
            AtmosToOceanBoxArray(m_ba_ocn, idx), m_dm_ocn, 1, AtmosToOceanNGrow(idx));
        mf_ptr->setVal(0.0);
        const std::string label = "Driver AllocateCaches atmos->ocean[" + std::to_string(idx) + "]";
        PrintMultiFabLayout(label.c_str(), *mf_ptr);
    }
    m_atmos_to_ocean = amrex::GetVecOfPtrs(m_atmos_to_ocean_storage);
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
            const std::string label = "Driver AllocateCaches ocean->atmos[" + std::to_string(idx) + "]";
            PrintMultiFabLayout(label.c_str(), *mf_ptr);
        }
    } else {
        m_ocean_to_atmos_storage.clear();
    }
    m_ocean_to_atmos = amrex::GetVecOfPtrs(m_ocean_to_atmos_storage);
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

void
ERFRemoraMultiBlockContainer::NormalizeWeightsAndCheckConservation(amrex::MultiFab& weight_mf,
                                                                   const amrex::iMultiFab& index_mf,
                                                                   const char* label)
{
    using namespace amrex;
    amrex::Real max_sum_err = 0.0;

    for (MFIter mfi(weight_mf); mfi.isValid(); ++mfi) {
        Box bx = mfi.validbox();
        auto w_arr = weight_mf.array(mfi);

        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                amrex::Real sum_w = 0.0;
                for (int m = 0; m < MAX_STENCIL_SIZE; ++m) {
                    sum_w += w_arr(i, j, 0, m);
                }
                if (sum_w > 0.0) {
                    max_sum_err = std::max(max_sum_err, std::abs(sum_w - 1.0));
                    for (int m = 0; m < MAX_STENCIL_SIZE; ++m) {
                        w_arr(i, j, 0, m) /= sum_w;
                    }
                }
            }
        }
    }

    ParallelDescriptor::ReduceRealMax(max_sum_err);
    if (ParallelDescriptor::IOProcessor() && m_driver_verbose >= 1) {
        amrex::Print() << "[driver] Remap weight normalization (" << label
                       << "): max weight sum drift before normalization = " << max_sum_err << "\n";
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
}

void
ERFRemoraMultiBlockContainer::BuildRemapWeights()
{
    using namespace amrex;

    // Query coordinate-aware sources (psi points from REMORA, corner points
    // from ERF) via the C-bridge. This is plumbing/guard-condition only: the
    // legacy m_geom/m_geom_ocn ProbLo/CellSize-based polygon-clipping loops below
    // are unchanged and remain the only path exercised per-field until a
    // coordinate-aware remap algorithm consumes the per-field compatibility
    // result computed here.
    const amrex::MultiFab* remora_x_psi = nullptr;
    const amrex::MultiFab* remora_y_psi = nullptr;
    const amrex::MultiFab* erf_x_corner = nullptr;
    const amrex::MultiFab* erf_y_corner = nullptr;
    if (m_remora != nullptr) {
        REMORAGetAtmosToOceanPsiCoordinates(m_remora, &remora_x_psi, &remora_y_psi);
    }
    if (m_erf != nullptr) {
        ERFGetOceanToAtmosCornerCoordinates(m_erf, &erf_x_corner, &erf_y_corner);
    }

    // ERF land mask (1 = land, 0 = water). Wet/dry masking only: vector
    // rotation (cos_angler/sin_angler) is intentionally not handled here
    // because REMORA does not currently expose those arrays.
    const amrex::iMultiFab* erf_is_land = nullptr;
    if (m_erf != nullptr) {
        ERFGetLandMask(m_erf, &erf_is_land);
    }
    // Coordinate-aware source is available at all only if both solvers report
    // non-null nodal coordinate MultiFabs (development baseline: coordinate
    // pointers are all-or-nothing per the null-safe getter contract).
    const bool coord_source_available =
        (remora_x_psi != nullptr && remora_y_psi != nullptr &&
         erf_x_corner != nullptr && erf_y_corner != nullptr);

    amrex::Print() << "[driver] BuildRemapWeights: curvilinear coordinate mode = "
                       << (coord_source_available ? "active (true physical corner polygons)" : "fallback (ProbLo/CellSize AABB)")
                       << "\n";

    // Reset A2O weight and index caches, one family slot at a time. Each slot
    // was allocated (AllocateCachesFromSolverLayouts) against its own
    // destination BoxArray/stagger (rho for SCALAR_CENTERED/WIND_FACE_TO_CENTER,
    // u-face for TAUX_FACE, v-face for TAUY_FACE), so building weights
    // per-slot means TauX/TauY finally get stencils computed against their own
    // stagger instead of borrowing the rho-point stencil. The polygon-clipping
    // math itself (ProbLo/CellSize-based AABB overlap) is unchanged; only the
    // destination MultiFab/BoxArray driving the (i,j) loop and dst_area varies.
    static const std::array<amrex::IntVect, kNumWeightFamilies> kExpectedFamilyStagger = {
        amrex::IntVect(0,0,0), amrex::IntVect(1,0,0), amrex::IntVect(0,1,0), amrex::IntVect(0,0,0)
    };
    for (int f = 0; f < kNumWeightFamilies; ++f) {
        MultiFab* weight_mf = m_weight_mf_by_family[f].get();
        iMultiFab* index_mf = m_index_mf_by_family[f].get();
        if (!weight_mf || !index_mf) { continue; }

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            weight_mf->boxArray().ixType().toIntVect() == kExpectedFamilyStagger[f],
            "BuildRemapWeights: A2O weight MultiFab stagger does not match its WeightFamily slot.");

        weight_mf->setVal(0.0);
        index_mf->setVal(-1);

        for (MFIter mfi(*weight_mf); mfi.isValid(); ++mfi) {
            Box bx = mfi.validbox();
            auto weight_arr = weight_mf->array(mfi);
            auto index_arr = index_mf->array(mfi);

            const Array4<const Real> rx_arr = coord_source_available ? remora_x_psi->const_array(mfi) : Array4<const Real>{};
            const Array4<const Real> ry_arr = coord_source_available ? remora_y_psi->const_array(mfi) : Array4<const Real>{};

            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    Vector<Point2D> remora_poly;
                    RealBox remora_bbox;

                    if (coord_source_available) {
                        remora_poly = ExtractCellQuadFromCoords(i, j, rx_arr, ry_arr);
                    } else {
                        RealBox bounds = CellBoundsFromGeom(m_geom_ocn, i, j);
                        remora_poly = {
                            {bounds.lo(0), bounds.lo(1)},
                            {bounds.hi(0), bounds.lo(1)},
                            {bounds.hi(0), bounds.hi(1)},
                            {bounds.lo(0), bounds.hi(1)}
                        };
                    }

                    Real dst_area = PolygonArea(remora_poly);
                    if (dst_area <= 0.0) continue;

                    Real x_min = remora_poly[0].x, x_max = remora_poly[0].x;
                    Real y_min = remora_poly[0].y, y_max = remora_poly[0].y;
                    for (const auto& p : remora_poly) {
                        x_min = std::min(x_min, p.x); x_max = std::max(x_max, p.x);
                        y_min = std::min(y_min, p.y); y_max = std::max(y_max, p.y);
                    }
                    remora_bbox = RealBox({x_min, y_min, 0.0}, {x_max, y_max, 0.0});

                    Box erf_candidate_box = m_geom.Domain();
                    int erf_i_lo = std::floor((remora_bbox.lo(0) - m_atm_rb.lo(0)) / (m_atm_rb.hi(0) - m_atm_rb.lo(0)) * m_atm_n_cell[0]);
                    int erf_i_hi = std::ceil((remora_bbox.hi(0) - m_atm_rb.lo(0)) / (m_atm_rb.hi(0) - m_atm_rb.lo(0)) * m_atm_n_cell[0]);
                    int erf_j_lo = std::floor((remora_bbox.lo(1) - m_atm_rb.lo(1)) / (m_atm_rb.hi(1) - m_atm_rb.lo(1)) * m_atm_n_cell[1]);
                    int erf_j_hi = std::ceil((remora_bbox.hi(1) - m_atm_rb.lo(1)) / (m_atm_rb.hi(1) - m_atm_rb.lo(1)) * m_atm_n_cell[1]);

                    // Add padded safety halo around search bounds for curvilinear distortion
                    erf_i_lo = std::max(erf_candidate_box.smallEnd(0), erf_i_lo - 2);
                    erf_i_hi = std::min(erf_candidate_box.bigEnd(0), erf_i_hi + 2);
                    erf_j_lo = std::max(erf_candidate_box.smallEnd(1), erf_j_lo - 2);
                    erf_j_hi = std::min(erf_candidate_box.bigEnd(1), erf_j_hi + 2);

                    int stencil_idx = 0;
                    for (int src_j = erf_j_lo; src_j <= erf_j_hi; ++src_j) {
                        for (int src_i = erf_i_lo; src_i <= erf_i_hi; ++src_i) {
                            // Wet/dry masking only (no vector rotation): skip land
                            // ERF cells so atmospheric land fluxes cannot bleed
                            // into coastal ocean cells via the remap stencil.
                            if (erf_is_land != nullptr) {
                                auto const& land_arr = erf_is_land->const_array(mfi);
                                if (land_arr(src_i, src_j, 0) == 1) { continue; }
                            }

                            Vector<Point2D> erf_poly;
                            if (coord_source_available && erf_x_corner && erf_y_corner) {
                                auto const& ex_arr = erf_x_corner->const_array(mfi);
                                auto const& ey_arr = erf_y_corner->const_array(mfi);
                                erf_poly = ExtractCellQuadFromCoords(src_i, src_j, ex_arr, ey_arr);
                            } else {
                                RealBox erf_bounds = CellBoundsFromGeom(m_geom, src_i, src_j);
                                erf_poly = {
                                    {erf_bounds.lo(0), erf_bounds.lo(1)},
                                    {erf_bounds.hi(0), erf_bounds.lo(1)},
                                    {erf_bounds.hi(0), erf_bounds.hi(1)},
                                    {erf_bounds.lo(0), erf_bounds.hi(1)}
                                };
                            }

                            Vector<Point2D> clipped_poly = ClipPolygon(remora_poly, erf_poly);
                            Real intersect_area = PolygonArea(clipped_poly);

                            if (intersect_area > 0.0) {
                                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                                    stencil_idx < MAX_STENCIL_SIZE,
                                    "BuildRemapWeights: A2O stencil size exceeded MAX_STENCIL_SIZE=9.");

                                Real weight = intersect_area / dst_area;
                                weight_arr(i, j, 0, stencil_idx) = weight;
                                index_arr(i, j, 0, stencil_idx * 3) = src_i;
                                index_arr(i, j, 0, stencil_idx * 3 + 1) = src_j;
                                index_arr(i, j, 0, stencil_idx * 3 + 2) = 0;
                                stencil_idx++;
                            }
                        }
                    }
                }
            }
        }
        NormalizeWeightsAndCheckConservation(*weight_mf, *index_mf, "A2O Family Slot");
    }

    // Reset and compute O2A weight and index caches if active
    if (m_use_two_way_coupling && m_weight_o2a_mf && m_index_o2a_mf) {
        m_weight_o2a_mf->setVal(0.0);
        m_index_o2a_mf->setVal(-1);

        for (MFIter mfi(*m_weight_o2a_mf); mfi.isValid(); ++mfi) {
            Box bx = mfi.validbox();
            auto weight_arr = m_weight_o2a_mf->array(mfi);
            auto index_arr = m_index_o2a_mf->array(mfi);

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

                    int stencil_idx = 0;
                    for (int src_j = ocn_j_lo; src_j <= ocn_j_hi; ++src_j) {
                        for (int src_i = ocn_i_lo; src_i <= ocn_i_hi; ++src_i) {
                            if (stencil_idx >= MAX_STENCIL_SIZE) break;

                            RealBox remora_cell_bounds = CellBoundsFromGeom(m_geom_ocn, src_i, src_j);
                            Vector<Point2D> remora_poly = {
                                {remora_cell_bounds.lo(0), remora_cell_bounds.lo(1)},
                                {remora_cell_bounds.hi(0), remora_cell_bounds.lo(1)},
                                {remora_cell_bounds.hi(0), remora_cell_bounds.hi(1)},
                                {remora_cell_bounds.lo(0), remora_cell_bounds.hi(1)}
                            };

                            Vector<Point2D> clipped_poly = ClipPolygon(erf_poly, remora_poly);
                            Real intersect_area = PolygonArea(clipped_poly);

                            if (intersect_area > 0.0) {
                                Real weight = intersect_area / dst_area;
                                weight_arr(i, j, 0, stencil_idx) = weight;
                                index_arr(i, j, 0, stencil_idx * 3) = src_i;
                                index_arr(i, j, 0, stencil_idx * 3 + 1) = src_j;
                                index_arr(i, j, 0, stencil_idx * 3 + 2) = 0;
                                stencil_idx++;
                            }
                        }
                    }
                }
            }
        }
    }
}

void
ERFRemoraMultiBlockContainer::FillSyntheticData ()
{
    for (int idx = 0; idx < static_cast<int>(m_atmos_to_ocean.size()); ++idx) {
        m_atmos_to_ocean[idx]->setVal(10.0 + static_cast<amrex::Real>(idx));
    }

    if (m_use_two_way_coupling && !m_ocean_to_atmos.empty()) {
        m_ocean_to_atmos[coupling::SST]->setVal(300.0);
        FillOceanToAtmosGhosts();
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
            amrex::Abort(std::string(label) + ": null field at index " + std::to_string(idx));
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
            ERFApplyOceanSurfaceState(m_erf, state, time);
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
                    if (is_cloud_lane && cloud_warn_count < cloud_warn_limit) {
                        const amrex::Real cmin = mf->min(0);
                        const amrex::Real cmax = mf->max(0);
                        if (amrex::ParallelDescriptor::IOProcessor()) {
                            if (cmin < -cloud_tol || cmax > amrex::Real(1.0) + cloud_tol) {
                                ++cloud_warn_count;
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
            if (!states.empty() && states[0] != nullptr) {
                const amrex::Real sample = states[0]->min(0);
                if (amrex::ParallelDescriptor::IOProcessor()) {
                    amrex::Print() << "Synthetic ApplyAtmosphericStates time=" << time
                                   << " sample_atmos_field0=" << sample << "\n";
                }
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
        apply_atmos_states(m_atmos_to_ocean_active, ocean_time);


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
