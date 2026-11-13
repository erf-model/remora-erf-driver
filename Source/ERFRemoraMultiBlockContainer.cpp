#include "ERFRemoraMultiBlockContainer.H"
#include "StandaloneSolverRuns.H"

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Math.H>
#include <AMReX_NonLocalBC.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Reduce.H>
#include <AMReX_String.H>

#ifdef ERF_USE_NETCDF
#include <netcdf.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>

namespace erf_remora
{
namespace {

struct MaskedMinMaxResult
{
    amrex::Real min_value = 0.0;
    amrex::Real max_value = 0.0;
    amrex::Long included_cells = 0;
    amrex::Long nonfinite_cells = 0;

    // With no finite included cell the extrema are still the reduction
    // sentinels, and must not be reported as values.
    [[nodiscard]] bool has_finite () const { return included_cells > nonfinite_cells; }
};

// Returns extrema over water cells when a mask is supplied, or over every cell
// when it is not. The included-cell count makes an all-masked result explicit,
// so the reduction sentinels are never reported as physical values.
//
// Non-finite values are counted rather than folded into the extrema. Min/max
// cannot be trusted to carry a NaN: amrex::min(d, r) is (r < d) ? r : d, so a
// NaN survives only if it is visited first, which depends on rank count, box
// layout and GPU block layout. Keeping them out of the comparisons also keeps
// amrex.fpe_trap_invalid from firing here instead of the caller's own check.
// The extrema cover the finite included cells only; see has_finite().
// Collective; every caller must invoke this on every rank.
MaskedMinMaxResult
MaskedMinMax (const amrex::MultiFab& mf,
              const amrex::MultiFab* mask,
              bool mask_is_land)
{
    using namespace amrex;

    constexpr Real min_sentinel = std::numeric_limits<Real>::max();
    constexpr Real max_sentinel = -std::numeric_limits<Real>::max();

    ReduceOps<ReduceOpMin, ReduceOpMax, ReduceOpSum, ReduceOpSum> reduce_ops;
    ReduceData<Real, Real, Long, Long> reduce_data(reduce_ops);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    const bool use_mask = (mask != nullptr);
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const auto field_arr = mf.const_array(mfi);
        const auto mask_arr = use_mask ? mask->const_array(mfi) : Array4<const Real>{};

        reduce_ops.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                const bool is_water = !use_mask ||
                    (mask_is_land ? mask_arr(i,j,k) == Real(0.0)
                                  : mask_arr(i,j,k) > Real(0.5));
                const Real value = field_arr(i,j,k);
                const bool is_finite = amrex::isfinite(value);
                const bool use_value = is_water && is_finite;
                return {use_value ? value : min_sentinel,
                        use_value ? value : max_sentinel,
                        is_water ? Long(1) : Long(0),
                        (is_water && !is_finite) ? Long(1) : Long(0)};
            });
    }

    const auto reduced = reduce_data.value(reduce_ops);
    MaskedMinMaxResult result;
    result.min_value = amrex::get<0>(reduced);
    result.max_value = amrex::get<1>(reduced);
    result.included_cells = amrex::get<2>(reduced);
    result.nonfinite_cells = amrex::get<3>(reduced);

    ParallelDescriptor::ReduceRealMin(result.min_value);
    ParallelDescriptor::ReduceRealMax(result.max_value);
    ParallelDescriptor::ReduceLongSum(result.included_cells);
    ParallelDescriptor::ReduceLongSum(result.nonfinite_cells);
    return result;
}

// True if any valid cell of any component holds a NaN or an infinity. Tests
// every value rather than inferring it from min/max, which can drop a NaN (see
// MaskedMinMax). Unmasked, so land cells count too. Collective.
bool
ContainsNonFinite (const amrex::MultiFab& mf)
{
    // One pass and one MPI reduction, where contains_nan || contains_inf
    // would take two of each.
    const auto ma = mf.const_arrays();
    bool r = amrex::ParReduce(amrex::TypeList<amrex::ReduceOpLogicalOr>{},
                              amrex::TypeList<bool>{}, mf, amrex::IntVect(0), mf.nComp(),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
            -> amrex::GpuTuple<bool>
        {
            return !amrex::isfinite(ma[box_no](i,j,k,n));
        });
    amrex::ParallelDescriptor::ReduceBoolOr(r);
    return r;
}

// Overwrite component 0 of one cell with a quiet NaN, for exercising the
// non-finite guards. The cell is the center of the lane's index space at its
// lowest k. A single cell is the case the old min/max guards could miss.
//
// Aborts if the cell lies in no valid box, or, when a wet mask is given (laid
// out like mf; mskr/msku/mskv convention, wet > 0.5), if it is not wet: either
// way the guard under test would have nothing to find, and the run would fail
// or pass for the wrong reason. Collective.
void
PoisonOneCellWithNaN (amrex::MultiFab& mf, const amrex::MultiFab* wet_mask)
{
    const amrex::Box domain = mf.boxArray().minimalBox();
    amrex::IntVect iv = (domain.smallEnd() + domain.bigEnd()) / 2;
    iv[2] = domain.smallEnd(2);
    const amrex::Box target(iv, iv, mf.ixType());
    const amrex::Real nan = std::numeric_limits<amrex::Real>::quiet_NaN();

    amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum> reduce_ops;
    amrex::ReduceData<amrex::Long, amrex::Long> reduce_data(reduce_ops);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    const bool use_mask = (wet_mask != nullptr);
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox() & target;
        if (!bx.ok()) { continue; }
        const auto arr = mf.array(mfi);
        const auto mask_arr = use_mask ? wet_mask->const_array(mfi)
                                       : amrex::Array4<const amrex::Real>{};
        reduce_ops.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                return {amrex::Long(1),
                        (!use_mask || mask_arr(i,j,k) > amrex::Real(0.5)) ? amrex::Long(1)
                                                                          : amrex::Long(0)};
            });
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            arr(i,j,k,0) = nan;
        });
    }
    const auto reduced = reduce_data.value(reduce_ops);
    amrex::Long counts[2] = {amrex::get<0>(reduced), amrex::get<1>(reduced)};
    amrex::ParallelDescriptor::ReduceLongSum(counts, 2);

    std::ostringstream where;
    where << iv;
    if (counts[0] == 0) {
        amrex::Abort("driver.debug_inject_nan_lane: cell " + where.str() +
                     " is in no valid box of the lane; nothing was poisoned.");
    }
    if (counts[1] == 0) {
        amrex::Abort("driver.debug_inject_nan_lane: cell " + where.str() +
                     " is land; the wet-cell validation cannot see it.");
    }
    amrex::Print() << "driver.debug_inject_nan_lane: NaN written at " << where.str() << "\n";
}

// Gather a remap array onto a single box on the IO rank and write its data raw
// (<dir>/<name>.bin) with a text description (<dir>/<name>.txt), so two runs'
// arrays compare exactly with cmp. The bytes do not depend on rank count or
// DistributionMapping. Across different BoxArrays a cell's stencil can list the
// same overlaps in a different order, since each box builds its own. Collective.
template <class FAB>
void
DumpFabArrayRaw (const amrex::FabArray<FAB>& src, const std::string& dir, const std::string& name)
{
    using namespace amrex;
    using value_type = typename FAB::value_type;

    const int io_rank = ParallelDescriptor::IOProcessorNumber();
    const Box domain = src.boxArray().minimalBox();
    const int ncomp = src.nComp();

    // A face-staggered array's boxes overlap on the faces between them, and each
    // box builds its own stencil there. ParallelCopy would keep whichever copy
    // lands last. Give each such face to its lower box (the rule the coverage
    // counts in NormalizeWeightsAndCheckConservation use), zero it in the other,
    // and add: every point then comes from exactly one box, whatever the order.
    FabArray<FAB> owned(src.boxArray(), src.DistributionMap(), ncomp, 0, MFInfo().SetArena(The_Pinned_Arena()));
    amrex::Copy(owned, src, 0, 0, ncomp, 0);
    const IntVect stagger = src.ixType().toIntVect();
    Gpu::streamSynchronize();   // the Copy above may still be running; setVal below is on the host
    for (MFIter mfi(owned); mfi.isValid(); ++mfi) {
        const Box& vbx = mfi.validbox();
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            if (stagger[dir] == 1 && vbx.bigEnd(dir) < domain.bigEnd(dir)) {
                Box hi_face = vbx;
                hi_face.setSmall(dir, vbx.bigEnd(dir));
                owned[mfi].template setVal<RunOn::Host>(value_type(0), hi_face, 0, ncomp);
            }
        }
    }

    // Pinned so the IO rank can read the gathered data from the host on GPU builds.
    FabArray<FAB> gathered(BoxArray(domain), DistributionMapping(Vector<int>{io_rank}),
                           ncomp, 0, MFInfo().SetArena(The_Pinned_Arena()));
    gathered.setVal(value_type(0));
    gathered.ParallelAdd(owned, 0, 0, ncomp);
    Gpu::streamSynchronize();

    for (MFIter mfi(gathered); mfi.isValid(); ++mfi) {
        const FAB& fab = gathered[mfi];
        std::ofstream bin(dir + "/" + name + ".bin", std::ios::binary);
        bin.write(reinterpret_cast<const char*>(fab.dataPtr()),
                  static_cast<std::streamsize>(sizeof(value_type) * fab.box().numPts() * fab.nComp()));
        std::ofstream txt(dir + "/" + name + ".txt");
        txt << "box " << fab.box() << "\n"
            << "ncomp " << fab.nComp() << "\n"
            << "type " << (std::is_integral_v<value_type> ? "int" : "float") << sizeof(value_type) * 8 << "\n"
            << "order i fastest, then j, k, component\n";
        if (!bin || !txt) { amrex::Abort("driver.remap_dump_weights: failed to write " + dir + "/" + name); }
    }
}

// The remap weights are built in host loops over FabArray data (polygon
// clipping with dynamic overlap lists does not map onto ParallelFor), so every
// array those loops touch must be host accessible even when the default arena
// is device memory. Driver-owned setup arrays are allocated with HostArenaInfo();
// solver-owned inputs are read through HostView(); the weight, index and valid
// arrays the exchanges read are built on the pinned arena and moved to the
// default arena by MoveToDefaultArena() once the build is done. Pinned memory is
// also device accessible, so the ParallelCopy, FillBoundary and ParallelFor
// calls between the host loops work on it unchanged.
amrex::MFInfo
HostArenaInfo ()
{
    return amrex::MFInfo().SetArena(amrex::The_Pinned_Arena());
}

// Every host loop over remap arrays starts here. Aborts with a message instead
// of segfaulting when an array is in device memory, and waits for device work
// still writing the arrays: setVal, Copy, ParallelCopy and FillBoundary return
// before their kernels finish on GPU builds.
template <class... FAs>
void
RequireHostAccess (const char* where, const FAs&... fas)
{
    if (!(... && fas.arena()->isHostAccessible())) {
        amrex::Abort(std::string(where) + ": a remap array is in device memory, but this loop reads "
                     "it on the host. Allocate it with HostArenaInfo() or read it through HostView().");
    }
    amrex::Gpu::streamSynchronize();
}

// Same, for arrays that may be absent: null pointers are skipped.
template <class... FAs>
void
RequireHostAccessIfPresent (const char* where, const FAs*... fas)
{
    if (!(... && (fas == nullptr || fas->arena()->isHostAccessible()))) {
        amrex::Abort(std::string(where) + ": a remap array is in device memory, but this loop reads "
                     "it on the host. Allocate it with HostArenaInfo() or read it through HostView().");
    }
    amrex::Gpu::streamSynchronize();
}

// Host-readable view of a solver-owned array: src itself when its arena is host
// accessible (CPU builds, or amrex.the_arena_is_managed=1), otherwise a pinned
// copy in `copy`, ghost cells included, since the staging plans read REMORA's
// coordinates one cell past the valid box. Null in, null out.
template <class FA>
const FA*
HostView (const FA* src, FA& copy)
{
    if (src == nullptr || src->arena()->isHostAccessible()) { return src; }
    copy.define(src->boxArray(), src->DistributionMap(), src->nComp(), src->nGrowVect(), HostArenaInfo());
    amrex::Copy(copy, *src, 0, 0, src->nComp(), src->nGrowVect());
    amrex::Gpu::streamSynchronize();
    return &copy;
}

// Replace a remap array built on the pinned arena with a default-arena copy, so
// the per-exchange remap kernels read device (or managed) memory rather than
// host memory over the bus. On CPU builds both arenas are host memory and this
// does nothing. Null is left alone.
template <class FA>
void
MoveToDefaultArena (std::unique_ptr<FA>& fa)
{
#ifndef AMREX_USE_GPU
    amrex::ignore_unused(fa);
    return;
#endif
    if (!fa) { return; }
    auto moved = std::make_unique<FA>(fa->boxArray(), fa->DistributionMap(), fa->nComp(), fa->nGrowVect());
    amrex::Copy(*moved, *fa, 0, 0, fa->nComp(), fa->nGrowVect());
    amrex::Gpu::streamSynchronize();
    fa = std::move(moved);
}

// Which model owns an O2A destination cell, and whether ERF's land mask and
// REMORA's wet mask agree about it. Every ERF surface cell is exactly one of these.
enum class OwnerMatchStatus : int {
    Coupled = 0,        // ERF water, wet REMORA overlap >= threshold: takes the coupled SST
    WaterOutside,       // ERF water, ocean grid does not reach: keeps ERF's SST
    WaterOceanDry,      // ERF water, every overlapping REMORA cell dry: keeps ERF's SST
    WaterThin,          // ERF water, wet overlap below remap_min_coverage_frac
    LandOceanWet,       // ERF land, wet REMORA overlap >= threshold: SST dropped
    Land,               // ERF land, no or thin wet overlap: both models agree
    Degenerate,         // zero-area ERF cell in the corner mesh: skipped
    NumStatus
};

// Cell count for one OwnerMatchStatus (or one A2O verdict), one example cell and
// the i/j bounding box. The example is the smallest flat (i,j) index over the
// destination domain, so it names the same cell under any decomposition; the box
// tells a coastline sliver from a band of misalignment. reduce() is collective.
struct OwnerMatchTally
{
    amrex::Long count = 0;
    amrex::Long first = std::numeric_limits<amrex::Long>::max();
    int ilo = std::numeric_limits<int>::max();
    int jlo = std::numeric_limits<int>::max();
    int ihi = std::numeric_limits<int>::lowest();
    int jhi = std::numeric_limits<int>::lowest();

    void add (const amrex::Box& domain, int i, int j) {
        ++count;
        const amrex::Long flat =
            static_cast<amrex::Long>(j - domain.smallEnd(1)) * domain.length(0)
            + static_cast<amrex::Long>(i - domain.smallEnd(0));
        first = std::min(first, flat);
        ilo = std::min(ilo, i);  ihi = std::max(ihi, i);
        jlo = std::min(jlo, j);  jhi = std::max(jhi, j);
    }

    void reduce () {
        amrex::ParallelDescriptor::ReduceLongSum(count);
        amrex::ParallelDescriptor::ReduceLongMin(first);
        amrex::ParallelDescriptor::ReduceIntMin(ilo);
        amrex::ParallelDescriptor::ReduceIntMin(jlo);
        amrex::ParallelDescriptor::ReduceIntMax(ihi);
        amrex::ParallelDescriptor::ReduceIntMax(jhi);
    }

    // "i[a..b] j[c..d]"; only meaningful when count > 0.
    std::string extent () const {
        return "i[" + std::to_string(ilo) + ".." + std::to_string(ihi) + "] j["
             + std::to_string(jlo) + ".." + std::to_string(jhi) + "]";
    }

    amrex::IntVect example (const amrex::Box& domain) const {
        return amrex::IntVect(domain.smallEnd(0) + static_cast<int>(first % domain.length(0)),
                              domain.smallEnd(1) + static_cast<int>(first / domain.length(0)),
                              0);
    }
};

// Twice the signed shoelace area, CCW positive. Each term is taken relative to
// the first vertex: absolute lon/lat products are ~1e7 times the area of a
// small cell, so the result would carry only ~9 digits and would depend on
// compiler contraction choices (issue #34). The area is translation-invariant.
amrex::Real SignedArea2 (const amrex::Vector<Point2D>& poly) {
    amrex::Real area = 0.0;
    const int n = static_cast<int>(poly.size());
    const Point2D o = poly[0];
    for (int i = 1; i + 1 < n; ++i) {
        const amrex::Real xi = poly[i].x - o.x, yi = poly[i].y - o.y;
        const amrex::Real xj = poly[i+1].x - o.x, yj = poly[i+1].y - o.y;
        area += xi * yj - xj * yi;
    }
    return area;
}

// Compute 2D polygon area using Shoelace formula
amrex::Real PolygonArea (const amrex::Vector<Point2D>& poly) {
    if (poly.size() < 3) return 0.0;
    return 0.5 * std::abs(SignedArea2(poly));
}

// Ensure polygon vertices are ordered counter-clockwise (CCW)
void EnsureCounterClockwise(amrex::Vector<Point2D>& poly) {
    if (poly.size() < 3) return;
    if (SignedArea2(poly) < 0.0) {
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

// Area of subjectPoly clipped by clipPoly. Both polygons are first shifted so
// subjectPoly's first vertex is the origin; otherwise EdgeIntersection's
// line constants are built from absolute coordinates and the clipped vertices
// lose the same ~7 digits as an absolute shoelace sum (issue #34). Neighbouring
// coordinates are within a factor of two of each other, so the shift itself is
// exact (Sterbenz) away from lon/lat = 0.
amrex::Real ClippedArea (amrex::Vector<Point2D> subjectPoly, amrex::Vector<Point2D> clipPoly) {
    if (subjectPoly.empty() || clipPoly.empty()) return 0.0;
    const Point2D o = subjectPoly[0];
    for (auto& p : subjectPoly) { p.x -= o.x; p.y -= o.y; }
    for (auto& p : clipPoly)    { p.x -= o.x; p.y -= o.y; }
    return PolygonArea(ClipPolygon(std::move(subjectPoly), std::move(clipPoly)));
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
    amrex::Real ocean_dt,
    int atm_steps_per_exchange,
    int ocean_steps_per_exchange,
    coupling::AtmosToOceanContract atm2ocn_contract,
    bool use_two_way_coupling)
    : m_atm_rb(atm_rb), m_atm_n_cell(n_cell_atm),
      m_ocn_rb(ocn_rb), m_ocn_n_cell(n_cell_ocn),
      m_max_level(max_level), m_coord(coord),
      m_ref_ratio(ref_ratio), m_is_periodic(is_periodic),
      m_erf_dt(erf_dt),
      m_ocean_dt(ocean_dt),
      m_atm_steps_per_exchange(atm_steps_per_exchange),
      m_ocean_steps_per_exchange(ocean_steps_per_exchange),
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
    pp_driver.query("scrip_weights_file", m_scrip_weights_file);
    pp_driver.query("scrip_a2o_group", m_scrip_a2o_group);
    pp_driver.query("scrip_o2a_group", m_scrip_o2a_group);
    pp_driver.query("remap_v", m_remap_verbose);
    pp_driver.query("remap_frame", m_remap_frame_request);
    pp_driver.query("remap_extend_corner_ring", m_remap_extend_corner_ring);
    pp_driver.query("remap_bin_cells", m_remap_bin_cells);
    pp_driver.query("remap_min_coverage_frac", m_remap_min_coverage_frac);
    pp_driver.query("remap_donor_fill", m_remap_donor_fill_request);
    pp_driver.query("remap_donor_max_diagonals", m_remap_donor_max_diagonals);
    pp_driver.query("remap_dump_weights", m_remap_dump_weights);
    pp_driver.query("debug_inject_nan_lane", m_debug_inject_nan_lane);
    pp_driver.query("debug_inject_nan_after_guard", m_debug_inject_nan_after_guard);
    // contains() rather than the query() return, so "set explicitly to the
    // default" is distinguishable from "not set" when reporting it as inert.
    m_remap_donor_max_diagonals_set = pp_driver.contains("remap_donor_max_diagonals");
    pp_driver.query("time_average_atm2ocn_state", m_time_average_atm2ocn_state);
    pp_driver.query("time_average_atm2ocn_flux", m_time_average_atm2ocn_flux);

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

    const bool use_state_contract =
        (m_atmos_to_ocean_contract == coupling::AtmosToOceanContract::StatePassing);
    DriverPrintLine(
        m_driver_verbose,
        std::string("injecting pre-InitData atm2ocn metadata: configured_contract=") +
            (use_state_contract ? "state" : "flux") +
            " two_way=" + (m_use_two_way_coupling ? "1" : "0"));
    ConfigureERFDriverAtmosToOceanCoupling(
        m_erf, true, m_use_two_way_coupling, use_state_contract);
    ConfigureREMORADriverAtmosToOceanCoupling(
        m_remora, true, m_use_two_way_coupling, use_state_contract);
    SetREMORAStartTime(m_remora, m_ocean_time_offset);

    DriverPrintLine(m_driver_verbose, "calling ERF InitData (enters ERF InitData_pre -> ERF InitData_post)");
    InitERFSolver(m_erf);
    DriverPrintLine(m_driver_verbose, "ERF InitData complete");

    // Between the two InitData calls on purpose, and the only window that works:
    // the ERF capability query reads fields ERF allocates above, and the REMORA
    // longwave reconfiguration this may perform has to land before REMORA
    // InitData opens NetCDF forcing handles off the same selector.
    DriverPrintLine(m_driver_verbose, "resolving atm2ocn lane ownership from ERF capability and REMORA deck");
    ResolveAtmosToOceanLaneSources();

    DriverPrintLine(m_driver_verbose, "calling REMORA InitData");
    InitREMORASolver(m_remora);
    DriverPrintLine(m_driver_verbose, "REMORA InitData complete");
}

void
ERFRemoraMultiBlockContainer::SetLandMaskFromDomainOverlap(const amrex::iMultiFab* erf_is_land)
{
    if (erf_is_land == nullptr) {
        if (m_driver_verbose > 0) {
            DriverPrintLine(m_driver_verbose,
                "ERF land mask not available; skipping domain overlap adjustment");
        }
        return;
    }

    DriverPrintLine(m_driver_verbose, "adjusting land mask for atmosphere cells outside ocean domain");

    // Cast away const for modification (same pattern as buildings update in ERF_MakeNewArrays.cpp:882)
    auto* lmask = const_cast<amrex::iMultiFab*>(erf_is_land);

    // Get ocean domain bounds
    const amrex::Real ocn_xlo = m_ocn_rb.lo(0);
    const amrex::Real ocn_xhi = m_ocn_rb.hi(0);
    const amrex::Real ocn_ylo = m_ocn_rb.lo(1);
    const amrex::Real ocn_yhi = m_ocn_rb.hi(1);

    // By value: ProbLo()/CellSize() point into the host-side Geometry, which a
    // device lambda cannot dereference.
    const auto plo = m_geom.ProbLoArray();
    const auto dx  = m_geom.CellSizeArray();

    // Iterate over land mask and mark cells outside ocean domain as land
    for (amrex::MFIter mfi(*lmask); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto lmask_arr = lmask->array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // Only modify cells currently marked as water (lmask = 0)
            // Preserve land (lmask = 1) and buildings (lmask = 2)
            if (lmask_arr(i, j, k) != 0) {
                return;
            }

            // Compute cell center physical position
            const amrex::Real x_center = plo[0] + (static_cast<amrex::Real>(i) + 0.5) * dx[0];
            const amrex::Real y_center = plo[1] + (static_cast<amrex::Real>(j) + 0.5) * dx[1];

            // Check if cell center lies outside REMORA domain
            const bool outside_ocean = (x_center < ocn_xlo || x_center > ocn_xhi ||
                                        y_center < ocn_ylo || y_center > ocn_yhi);

            if (outside_ocean) {
                lmask_arr(i, j, k) = 1;  // Mark as land
            }
        });
    }

    DriverPrintLine(m_driver_verbose, "land mask adjustment complete");
}

void
ERFRemoraMultiBlockContainer::PrefillAtmosToOceanFallbacks ()
{
    // Only the flux lanes need seeding, and only with zero.
    //
    // The state lanes used to be prefilled with REMORA's own constants read back
    // out of ParmParse. That was redundant - a lane the driver does not own is
    // withheld and never read, so REMORA keeps its value directly - and fragile,
    // since it depended on REMORA's queryAdd having already written its defaults
    // into the table. Resolving lane ownership up front removes the need for the
    // driver to guess at REMORA's defaults at all.
    for (int global_idx = 0; global_idx < static_cast<int>(m_atmos_to_ocean.size()); ++global_idx) {
        auto* mf = m_atmos_to_ocean[global_idx];
        if (mf == nullptr) {
            continue;
        }

        switch (global_idx) {
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

    if (m_atm_extends_beyond_ocn) {
        const amrex::iMultiFab* erf_is_land = nullptr;
        ERFGetLandMask(m_erf, &erf_is_land);
        SetLandMaskFromDomainOverlap(erf_is_land);
    }

    BuildRemapWeights();
    DriverPrintLine(m_driver_verbose, "selecting configured active atm2ocn cache view");
    SelectAtmosToOceanActiveView();
    if (m_use_two_way_coupling && !m_ocean_to_atmos.empty()) {
        DriverPrintLine(m_driver_verbose, "packing REMORA surface state into ocean->atmos bootstrap cache");
        const amrex::iMultiFab* erf_is_land_bootstrap = nullptr;
        ERFGetLandMask(m_erf, &erf_is_land_bootstrap);
        REMORAPackSurfaceState(m_remora, m_ocean_to_atmos, amrex::Real(0.0),
                               m_weight_o2a_mf.get(), m_index_o2a_mf.get(), MAX_STENCIL_SIZE,
                               erf_is_land_bootstrap);
        AbortOnNonFinite("PackSurfaceState (bootstrap)",
                         m_ocean_to_atmos, coupling::OceanToAtmosNames());
        FillOceanToAtmosGhosts();
    }
    CheckDebugInjectNaNLane();

    if (m_atmos_to_ocean_contract == coupling::AtmosToOceanContract::FluxPassing) {
        DriverPrintLine(m_driver_verbose, "validating flux active view shape and centering");
        ValidateFluxActiveViewOrAbort(m_atmos_to_ocean_flux_active, m_atmos_to_ocean_flux_active_names);
    }

    amrex::Print() << "Atmos->ocean contract mode: "
                       << coupling::AtmosToOceanContractName(m_atmos_to_ocean_contract)
                       << " (" << coupling::NumAtmosToOceanActiveLayers(m_atmos_to_ocean_contract)
                       << " active layers)\n";
    amrex::Print() << "Atmos->ocean active lanes: "
                   << JoinLaneNames(m_atmos_to_ocean_active_names) << "\n";

    // Every active lane in the flux contract is flux-character (tau_x, tau_y, SHflux,
    // LHflux, SWrad, LWrad, rain, evap), so driver.time_average_atm2ocn_state selects
    // no lane at all there. It cannot average anything, yet it is enough to make the
    // report below say ENABLED and to suppress the instantaneous-forcing warning.
    // Reject the combination rather than announce averaging that does not happen.
    if (m_atmos_to_ocean_contract == coupling::AtmosToOceanContract::FluxPassing &&
        m_time_average_atm2ocn_state) {
        amrex::Abort("driver.time_average_atm2ocn_state=1 has no effect when "
                     "driver.atm2ocn_mode=flux: every lane in the flux contract is "
                     "flux-character. Set driver.time_average_atm2ocn_flux=1 to average "
                     "the flux lanes, or drop the state flag.");
    }

    amrex::Print() << "[driver] Atmos->ocean time averaging for state-character lanes is "
                   << (m_time_average_atm2ocn_state ? "ENABLED" : "disabled")
                   << " (driver.time_average_atm2ocn_state).\n"
                   << "[driver] Atmos->ocean time averaging for flux-character lanes is "
                   << (m_time_average_atm2ocn_flux ? "ENABLED" : "disabled")
                   << " (driver.time_average_atm2ocn_flux).\n";
    // Key the warning on whether any lane in the ACTIVE VIEW is actually averaged, not
    // on the flags. A flag set for a lane character absent from the view averages
    // nothing, and inspecting the flags alone let exactly that case suppress this
    // warning while reporting ENABLED.
    if (m_atm_steps_per_exchange > 1 && !AtmosToOceanAveragingActive()) {
        amrex::Print()
            << "[driver] WARNING: atmosphere->ocean forcing uses the instantaneous "
               "interval-end sample. Enable driver.time_average_atm2ocn_state=1 and/or "
               "driver.time_average_atm2ocn_flux=1 to force REMORA with the atmosphere "
               "mean over each ocean interval.\n";
    }
    PrintAtmosToOceanLaneProvenance();

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

    // Compute domain extent comparison once for reuse
    constexpr amrex::Real eps = 1.0e-6;
    m_atm_extends_beyond_ocn = (m_atm_rb.lo(0) < m_ocn_rb.lo(0) - eps) ||
                               (m_atm_rb.hi(0) > m_ocn_rb.hi(0) + eps) ||
                               (m_atm_rb.lo(1) < m_ocn_rb.lo(1) - eps) ||
                               (m_atm_rb.hi(1) > m_ocn_rb.hi(1) + eps);

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

    if (m_time_average_atm2ocn_state || m_time_average_atm2ocn_flux) {
        m_atmos_to_ocean_accumulator_storage.resize(coupling::NumAtmosToOceanLayers);
        for (int idx = 0;
             idx < static_cast<int>(m_atmos_to_ocean_accumulator_storage.size()); ++idx) {
            if (!TimeAverageAtmosToOceanLane(idx)) { continue; }
            const amrex::MultiFab& lane = *m_atmos_to_ocean[idx];
            auto& accumulator = m_atmos_to_ocean_accumulator_storage[idx];
            accumulator = std::make_unique<amrex::MultiFab>(
                lane.boxArray(), lane.DistributionMap(), 1, 0);
            accumulator->setVal(0.0);
        }
    }

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
        // Built on the host; BuildRemapWeights moves them to the default arena.
        for (int f = 0; f < kNumWeightFamilies; ++f) {
            m_weight_mf_by_family[f] = std::make_unique<amrex::MultiFab>(
                *family_ba[f], *family_dm[f], MAX_STENCIL_SIZE, 0, HostArenaInfo());
            m_index_mf_by_family[f] = std::make_unique<amrex::iMultiFab>(
                *family_ba[f], *family_dm[f], 3 * MAX_STENCIL_SIZE, 0, HostArenaInfo());
            // The donor-fill flag array is NOT allocated here. It is only ever
            // written by the coawst donor search, which runs solely in the
            // spherical frame, and neither the mode nor the frame is resolved
            // yet: this function is reached from AllocateCaches() one line before
            // BuildRemapWeights(), which is what calls ResolveRemapFrame() and
            // ResolveDonorFillMode(). Allocating it unconditionally here gave four
            // full-layout iMultiFabs that are structurally all-zero under the
            // shipped default, read back per cell only to add zero. It is
            // allocated lazily in BuildAtmosToOceanWeights instead, where both are
            // known.
        }
    }

    // Ocean -> Atmosphere (O2A) stencils (lives on ERF layout). Built on the host
    // like the A2O ones above.
    if (m_use_two_way_coupling) {
        m_weight_o2a_mf = std::make_unique<amrex::MultiFab>(m_ba, m_dm, MAX_STENCIL_SIZE, 0, HostArenaInfo());
        m_index_o2a_mf = std::make_unique<amrex::iMultiFab>(m_ba, m_dm, 3 * MAX_STENCIL_SIZE, 0, HostArenaInfo());
        m_valid_o2a_mf = std::make_unique<amrex::iMultiFab>(m_ba, m_dm, 1, 0, HostArenaInfo());
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

bool
ERFRemoraMultiBlockContainer::TimeAverageAtmosToOceanLane (int global_idx) const
{
    const bool flux_character =
        global_idx == coupling::Rain || global_idx == coupling::SWrad ||
        global_idx == coupling::LWrad || global_idx == coupling::TauX ||
        global_idx == coupling::TauY || global_idx == coupling::SHflux ||
        global_idx == coupling::LHflux || global_idx == coupling::Evap;
    return flux_character ? m_time_average_atm2ocn_flux : m_time_average_atm2ocn_state;
}

bool
ERFRemoraMultiBlockContainer::AtmosToOceanAveragingActive () const
{
    for (int global_idx : m_atmos_to_ocean_active_indices) {
        if (TimeAverageAtmosToOceanLane(global_idx)) {
            return true;
        }
    }
    return false;
}

void
ERFRemoraMultiBlockContainer::ResetAtmosToOceanAccumulators ()
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_atmos_to_ocean_accumulator_storage.empty(),
        "Atmos->ocean averaging requested without allocated accumulators");
    for (auto& accumulator : m_atmos_to_ocean_accumulator_storage) {
        if (accumulator) { accumulator->setVal(0.0); }
    }
    m_a2o_accum_weight = amrex::Real(0.0);
}

void
ERFRemoraMultiBlockContainer::AccumulateAtmosToOcean (amrex::Real weight)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        weight >= amrex::Real(0.0),
        "Atmos->ocean averaging sample weight must not be negative");
    if (weight == amrex::Real(0.0)) { return; }
    for (int local_idx = 0;
         local_idx < static_cast<int>(m_atmos_to_ocean_active_indices.size()); ++local_idx) {
        const int global_idx = m_atmos_to_ocean_active_indices[local_idx];
        auto* accumulator = m_atmos_to_ocean_accumulator_storage[global_idx].get();
        if (accumulator != nullptr) {
            amrex::MultiFab::Saxpy(*accumulator,
                                   weight, *m_atmos_to_ocean_active[local_idx],
                                   0, 0, 1, 0);
        }
    }
    m_a2o_accum_weight += weight;
}

void
ERFRemoraMultiBlockContainer::DrainAtmosToOceanAccumulators ()
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_a2o_accum_weight > amrex::Real(0.0),
        "Atmos->ocean averaging drained an empty accumulator");
    const amrex::Real inverse_weight = amrex::Real(1.0) / m_a2o_accum_weight;
    for (int local_idx = 0;
         local_idx < static_cast<int>(m_atmos_to_ocean_active_indices.size()); ++local_idx) {
        const int global_idx = m_atmos_to_ocean_active_indices[local_idx];
        const auto* accumulator = m_atmos_to_ocean_accumulator_storage[global_idx].get();
        if (accumulator != nullptr) {
            amrex::MultiFab::Copy(*m_atmos_to_ocean_active[local_idx],
                                  *accumulator,
                                  0, 0, 1, 0);
            m_atmos_to_ocean_active[local_idx]->mult(inverse_weight, 0);
        }
    }
    m_a2o_accum_weight = amrex::Real(0.0);
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
                                                                   const amrex::iMultiFab* donor_fill,
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
    // Cells whose stencil is a donor graft rather than a real overlap. These have
    // sum_w > 0, so every count above already treats them as covered; without a
    // separate tally the verdict below would report a clean sheet for a run in
    // which some cells took a neighbour's atmospheric state instead of their own.
    // The donor search only fires for cells that are wet under this same mask, so
    // this is a wet-cell count whenever it is nonzero.
    amrex::Long n_donor_filled = 0;
    // Example cells for the two wet-cell verdicts, in this family's own index
    // space (face indices for the tau families).
    OwnerMatchTally donor_example;
    OwnerMatchTally wet_uncovered_example;

    // Faces on a box seam belong to both adjacent boxes, so summing validbox over
    // boxes would count them twice. Give each face to the lower box and let the
    // last box in each direction keep its domain-boundary face, which tiles the
    // staggered index space exactly once. Weights are still computed on every
    // face including shared ones; this only makes the counts truthful.
    const Box stagger_domain = weight_mf.boxArray().minimalBox();
    const IntVect stagger = weight_mf.boxArray().ixType().toIntVect();

    RequireHostAccess("NormalizeWeightsAndCheckConservation", weight_mf, index_mf);
    RequireHostAccessIfPresent("NormalizeWeightsAndCheckConservation", wet_mask, donor_fill);

    for (MFIter mfi(weight_mf); mfi.isValid(); ++mfi) {
        Box bx = mfi.validbox();
        auto w_arr = weight_mf.array(mfi);
        auto const& idx_arr = index_mf.const_array(mfi);
        const bool has_mask = (wet_mask != nullptr);
        auto const& wet_arr = has_mask ? wet_mask->const_array(mfi) : Array4<const Real>{};
        const bool has_donor_flag = (donor_fill != nullptr);
        auto const& donor_arr = has_donor_flag ? donor_fill->const_array(mfi)
                                               : Array4<const int>{};

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
                const bool is_donor_filled = has_donor_flag && (donor_arr(i,j,0) != 0);
                if (counted) {
                    ++n_cells;
                    if (idx_arr(i, j, 0, 0) < 0) { ++n_no_stencil; }
                    if (has_mask && is_wet) { ++n_wet; }
                    if (is_donor_filled) {
                        ++n_donor_filled;
                        donor_example.add(stagger_domain, i, j);
                    }
                }

                if (sum_w > 0.0) {
                    if (counted) {
                        frac_min = std::min(frac_min, sum_w);
                        frac_max = std::max(frac_max, sum_w);
                        frac_sum += sum_w;
                        max_cover_err = std::max(max_cover_err, std::abs(sum_w - 1.0));
                        // A donor graft is a single full-area overlap, so its
                        // sum_w is exactly the destination area and normalizes to
                        // a flawless 1.0. Letting that into wet_frac_min would
                        // report the worst real coverage as better than it is, so
                        // the statistic is taken over genuinely covered wet cells.
                        if (has_mask && is_wet && !is_donor_filled) {
                            wet_frac_min = std::min(wet_frac_min, sum_w);
                        }
                    }
                    for (int m = 0; m < MAX_STENCIL_SIZE; ++m) {
                        w_arr(i, j, 0, m) /= sum_w;
                    }
                } else if (counted) {
                    ++n_uncovered;
                    if (has_mask && is_wet) {
                        ++n_wet_uncovered;
                        wet_uncovered_example.add(stagger_domain, i, j);
                    }
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
    ParallelDescriptor::ReduceLongSum(n_donor_filled);
    donor_example.reduce();
    wet_uncovered_example.reduce();

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
                       << " no_stencil=" << n_no_stencil;
        // Only for directions that can donor-fill at all, so the O2A line keeps
        // its existing shape rather than carrying a field that is structurally
        // zero there.
        if (donor_fill != nullptr) {
            if (wet_mask != nullptr) {
                amrex::Print() << " donor_filled=" << n_donor_filled;
            } else {
                amrex::Print() << " donor_filled=n/a (no wet mask)";
            }
        }
        amrex::Print() << " frac min/mean/max=" << frac_min << "/"
                       << (frac_sum / static_cast<amrex::Real>(n_covered)) << "/" << frac_max
                       << " max|frac-1|=" << max_cover_err << "\n";
    }

    // The verdict: an uncovered wet cell is a real hole in the forcing, whereas an
    // uncovered dry cell is expected wherever the atmosphere grid is land.
    if (wet_mask != nullptr) {
        // Reported before, and independently of, the uncovered warning. A
        // donor-filled cell is not uncovered -- it has a stencil and sum_w > 0, so
        // n_uncovered and n_wet_uncovered deliberately exclude it, and those two
        // counters keep exactly the meaning they had before this category existed.
        // But it is not covered either, in the sense the coverage numbers are read
        // for: no source cell of its own overlaps it. Without this line a run in
        // which every uncovered cell was donor-filled prints the clean-sheet
        // verdict below, and a reader cannot tell a genuinely covered cell from
        // one holding a neighbour's value.
        if (n_donor_filled > 0) {
            const IntVect ex = donor_example.example(stagger_domain);
            amrex::Print() << "[driver] WARNING remap (" << label << "): " << n_donor_filled
                           << " of " << n_wet << " REMORA wet cells (e.g. REMORA (" << ex[0]
                           << "," << ex[1] << "); " << donor_example.extent()
                           << ") have no ERF source cell of "
                           << "their own and were donor-filled from the nearest non-land ERF "
                           << "cell (driver.remap_donor_fill=coawst). Each carries a "
                           << "neighbour's atmospheric state, not its own, and the full-weight "
                           << "graft makes this atmosphere->ocean operator non-conservative. "
                           << "They are excluded from the min wet frac below and are NOT part "
                           << "of the uncovered counts. Under the default "
                           << "driver.remap_donor_fill=withhold these cells would instead be "
                           << "reported uncovered and take the per-lane ApplyConservativeRemap "
                           << "fallbacks.\n";
        }
        if (n_wet_uncovered > 0) {
            const IntVect ex = wet_uncovered_example.example(stagger_domain);
            amrex::Print() << "[driver] WARNING remap (" << label << "): " << n_wet_uncovered
                           << " of " << n_wet << " REMORA wet cells (e.g. REMORA (" << ex[0]
                           << "," << ex[1] << "); " << wet_uncovered_example.extent()
                           << ") have no ERF water source overlap. "
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
        } else if (n_donor_filled == 0 &&
                   ParallelDescriptor::IOProcessor() && m_remap_verbose >= 1) {
            // The clean sheet, and it says donor_filled=0 out loud so the reader
            // does not have to infer it from the line's absence. When any cell was
            // donor-filled this message is suppressed and the warning above stands
            // as the verdict instead.
            amrex::Print() << "[driver] remap coverage (" << label << "): all " << n_wet
                           << " wet cells covered by real source overlap (donor_filled=0)"
                           << ", min wet frac=" << wet_frac_min << "\n";
        }
    } else if (donor_fill != nullptr && ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "[driver] WARNING remap (" << label
                       << "): wet-cell verdict skipped because no wet mask was available; "
                       << "uncovered= counts all cells rather than wet cells.\n";
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

    RequireHostAccess("BuildReplicatedMaskLookup", mask_mf);
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

    RequireHostAccess("ExtendCornerMeshOneRing", node_mf);
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
    CompareWeightsWithSCRIP();

    // Built on the pinned arena for the host loops above; the per-exchange remap
    // kernels read these, so move them to device memory. Must finish before
    // InitializeBlocks packs the bootstrap ocean->atmos state with the O2A set.
    for (int f = 0; f < kNumWeightFamilies; ++f) {
        MoveToDefaultArena(m_weight_mf_by_family[f]);
        MoveToDefaultArena(m_index_mf_by_family[f]);
    }
    MoveToDefaultArena(m_weight_o2a_mf);
    MoveToDefaultArena(m_index_o2a_mf);
    MoveToDefaultArena(m_valid_o2a_mf);

    if (m_remap_dump_weights) {
        DumpRemapWeights();
    }

    // The corner mesh, the staged coordinates and masks, and the donor flags are
    // read only by the weight build above. They live on the pinned arena, so free
    // them rather than hold page-locked host memory for the rest of the run. The
    // staging plans keep their layouts and comm metadata.
    m_erf_node_coords.reset();
    if (m_staging_plan) {
        m_staging_plan->node_mf.reset();
        m_staging_plan->land_mf.reset();
    }
    if (m_remora_staging_plan) {
        m_remora_staging_plan->node_mf.reset();
        m_remora_staging_plan->wet_mf.reset();
    }
    for (auto& donor_fill : m_donor_fill_mf_by_family) { donor_fill.reset(); }
}

// Write every remap array the exchanges read, for exact comparison of two runs
// (driver.remap_dump_weights; see DumpFabArrayRaw). Collective.
void
ERFRemoraMultiBlockContainer::DumpRemapWeights () const
{
    using namespace amrex;

    const std::string dir = "remap_weights";
    UtilCreateCleanDirectory(dir, true);

    static constexpr std::array<const char*, kNumWeightFamilies> family_names = {
        "scalar", "taux", "tauy", "wind"
    };
    for (int f = 0; f < kNumWeightFamilies; ++f) {
        const std::string family = family_names[f];
        DumpFabArrayRaw(*m_weight_mf_by_family[f], dir, "a2o_weight_" + family);
        DumpFabArrayRaw(*m_index_mf_by_family[f], dir, "a2o_index_" + family);
        if (m_donor_fill_mf_by_family[f]) {
            DumpFabArrayRaw(*m_donor_fill_mf_by_family[f], dir, "a2o_donor_fill_" + family);
        }
    }
    if (m_weight_o2a_mf) {
        DumpFabArrayRaw(*m_weight_o2a_mf, dir, "o2a_weight");
        DumpFabArrayRaw(*m_index_o2a_mf, dir, "o2a_index");
        DumpFabArrayRaw(*m_valid_o2a_mf, dir, "o2a_valid");
    }
    if (m_o2a_uncoupled_mask) {
        DumpFabArrayRaw(*m_o2a_uncoupled_mask, dir, "o2a_uncoupled_mask");
    }
    amrex::Print() << "[driver] remap weights written to " << dir << "/\n";
}

// Compare our weights with an external SCRIP file (driver.scrip_weights_file,
// e.g. the scrip_coawst.nc a COAWST run used), once, after they are built. Per
// direction: which destination cells each side covers, and where both do, the L1
// difference of the first-order weights. Both are normalized per destination cell
// (SCRIP fracarea), so L1 runs 0..2; it is only meaningful when the SCRIP map is
// also conservative, and map_method is printed so the reader can tell. SCRIP
// addresses are 1-based, i fastest; the grids may carry a boundary ring (ROMS rho
// points), so the index offset is half the dims difference. Collective.
void
ERFRemoraMultiBlockContainer::CompareWeightsWithSCRIP ()
{
    using namespace amrex;
    if (m_scrip_weights_file.empty()) { return; }
#ifndef ERF_USE_NETCDF
    amrex::Print() << "[driver] driver.scrip_weights_file needs a NetCDF build; SCRIP check skipped.\n";
#else

    auto compare = [&] (const std::string& group, const MultiFab& weight_mf,
                        const iMultiFab& index_mf, const Box& src_domain,
                        const char* label, const char* dst_name)
    {
        // Read on the IO rank and broadcast: dims, addresses, first weight.
        int dims[4] = {0, 0, 0, 0};   // src nx, ny, dst nx, ny
        Long nlinks = 0;
        Vector<int> src_addr, dst_addr;
        Vector<Real> w;
        std::string method;
        const int io = ParallelDescriptor::IOProcessorNumber();
        if (ParallelDescriptor::IOProcessor()) {
            auto ok = [&] (int status, const char* what) {
                if (status != NC_NOERR) {
                    amrex::Abort("driver.scrip_weights_file " + m_scrip_weights_file + ": " +
                                 what + ": " + nc_strerror(status));
                }
            };
            int file = 0, grp = 0, id = 0, nw_dim = 0;
            size_t len = 0, nw = 1;
            ok(nc_open(m_scrip_weights_file.c_str(), NC_NOWRITE, &file), "open");
            // COAWST files hold one group per direction; a plain SCRIP file has none.
            if (nc_inq_grp_ncid(file, group.c_str(), &grp) != NC_NOERR) { grp = file; }
            // Only a 2-D logically rectangular grid has (nx, ny); any other rank
            // leaves zeros, which the size check below reports as skipped.
            auto get_dims = [&] (const char* name, int* out) {
                int dimid = 0;
                size_t rank = 0;
                ok(nc_inq_varid(grp, name, &id), name);
                ok(nc_inq_vardimid(grp, id, &dimid), name);
                ok(nc_inq_dimlen(grp, dimid, &rank), name);
                if (rank == 2) { ok(nc_get_var_int(grp, id, out), name); }
            };
            get_dims("src_grid_dims", dims);
            get_dims("dst_grid_dims", dims + 2);
            ok(nc_inq_dimid(grp, "num_links", &id), "num_links");
            ok(nc_inq_dimlen(grp, id, &len), "num_links");
            if (nc_inq_dimid(grp, "num_wgts", &nw_dim) == NC_NOERR) { nc_inq_dimlen(grp, nw_dim, &nw); }
            nlinks = static_cast<Long>(len);
            src_addr.resize(len);  dst_addr.resize(len);  w.resize(len);
            Vector<double> all(len * nw);
            ok(nc_inq_varid(grp, "src_address", &id), "src_address");
            ok(nc_get_var_int(grp, id, src_addr.data()), "src_address");
            ok(nc_inq_varid(grp, "dst_address", &id), "dst_address");
            ok(nc_get_var_int(grp, id, dst_addr.data()), "dst_address");
            ok(nc_inq_varid(grp, "remap_matrix", &id), "remap_matrix");
            ok(nc_get_var_double(grp, id, all.data()), "remap_matrix");
            for (size_t l = 0; l < len; ++l) { w[l] = all[l * nw]; }   // first-order weight
            size_t alen = 0;
            if (nc_inq_attlen(grp, NC_GLOBAL, "map_method", &alen) == NC_NOERR) {
                method.resize(alen);
                nc_get_att_text(grp, NC_GLOBAL, "map_method", method.data());
            }
            nc_close(file);
        }
        ParallelDescriptor::Bcast(dims, 4, io);
        ParallelDescriptor::Bcast(&nlinks, 1, io);
        src_addr.resize(nlinks);  dst_addr.resize(nlinks);  w.resize(nlinks);
        ParallelDescriptor::Bcast(src_addr.data(), nlinks, io);
        ParallelDescriptor::Bcast(dst_addr.data(), nlinks, io);
        ParallelDescriptor::Bcast(w.data(), nlinks, io);

        const Box dst_domain = weight_mf.boxArray().minimalBox();
        const int off[4] = {(dims[0] - src_domain.length(0)) / 2, (dims[1] - src_domain.length(1)) / 2,
                            (dims[2] - dst_domain.length(0)) / 2, (dims[3] - dst_domain.length(1)) / 2};
        if (*std::min_element(off, off + 4) < 0) {
            amrex::Print() << "[driver] SCRIP check " << label << ": skipped, " << group
                           << " grids (" << dims[0] << "x" << dims[1] << " -> " << dims[2] << "x"
                           << dims[3] << ") are smaller than ours.\n";
            return;
        }

        // Links grouped by destination: links of 0-based dst address a are
        // [start[a], start[a+1]). Weights > 0 only.
        const Long ndst = Long(dims[2]) * dims[3];
        const Long nsrc = Long(dims[0]) * dims[1];
        for (Long l = 0; l < nlinks; ++l) {
            if (src_addr[l] < 1 || src_addr[l] > nsrc || dst_addr[l] < 1 || dst_addr[l] > ndst) {
                amrex::Abort("driver.scrip_weights_file " + m_scrip_weights_file + ": " + group +
                             " link " + std::to_string(l) + " addresses a cell outside its grid_dims.");
            }
        }
        Vector<Long> start(ndst + 1, 0);
        for (Long l = 0; l < nlinks; ++l) { if (w[l] > 0) { ++start[dst_addr[l]]; } }   // 1-based
        for (Long a = 0; a < ndst; ++a) { start[a + 1] += start[a]; }
        Vector<Long> cursor(start.begin(), start.end() - 1);
        Vector<int> link_src(start[ndst]);
        Vector<Real> link_w(start[ndst]);
        for (Long l = 0; l < nlinks; ++l) {
            if (w[l] > 0) {
                const Long k = cursor[dst_addr[l] - 1]++;
                link_src[k] = src_addr[l] - 1;
                link_w[k] = w[l];
            }
        }

        OwnerMatchTally both, ours_only, scrip_only;
        Real l1_sum = 0, l1_max = -1;
        Long l1_first = std::numeric_limits<Long>::max();
        RequireHostAccess("CompareWeightsWithSCRIP", weight_mf, index_mf);
        for (MFIter mfi(weight_mf); mfi.isValid(); ++mfi) {
            auto const& wa = weight_mf.const_array(mfi);
            auto const& ia = index_mf.const_array(mfi);
            const Box& bx = mfi.validbox();
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    const Long a = Long(j - dst_domain.smallEnd(1) + off[3]) * dims[2]
                                 + (i - dst_domain.smallEnd(0) + off[2]);
                    const bool scrip_has = start[a + 1] > start[a];
                    const bool ours_has = ia(i,j,0,0) >= 0;
                    if (!ours_has && !scrip_has) { continue; }
                    if (!scrip_has) { ours_only.add(dst_domain, i, j);  continue; }
                    if (!ours_has)  { scrip_only.add(dst_domain, i, j); continue; }
                    both.add(dst_domain, i, j);
                    // L1 over the union of the two stencils.
                    Real l1 = 0, matched = 0;
                    for (int m = 0; m < MAX_STENCIL_SIZE && ia(i,j,0,3*m) >= 0; ++m) {
                        const int s = (ia(i,j,0,3*m+1) - src_domain.smallEnd(1) + off[1]) * dims[0]
                                    + (ia(i,j,0,3*m) - src_domain.smallEnd(0) + off[0]);
                        Real ws = 0;
                        for (Long k = start[a]; k < start[a + 1]; ++k) { if (link_src[k] == s) { ws += link_w[k]; } }
                        l1 += std::abs(wa(i,j,0,m) - ws);
                        matched += ws;
                    }
                    for (Long k = start[a]; k < start[a + 1]; ++k) { l1 += link_w[k]; }
                    l1 -= matched;   // SCRIP weights on sources we do not use
                    l1_sum += l1;
                    const Long flat = Long(j - dst_domain.smallEnd(1)) * dst_domain.length(0)
                                    + (i - dst_domain.smallEnd(0));
                    if (l1 > l1_max || (l1 == l1_max && flat < l1_first)) { l1_max = l1; l1_first = flat; }
                }
            }
        }
        both.reduce();  ours_only.reduce();  scrip_only.reduce();
        ParallelDescriptor::ReduceRealSum(l1_sum);
        Real global_max = l1_max;
        ParallelDescriptor::ReduceRealMax(global_max);
        if (l1_max != global_max) { l1_first = std::numeric_limits<Long>::max(); }
        ParallelDescriptor::ReduceLongMin(l1_first);

        auto where = [&] (const OwnerMatchTally& t) {
            if (t.count == 0) { return std::string{}; }
            const IntVect e = t.example(dst_domain);
            return " (e.g. " + std::string(dst_name) + " (" + std::to_string(e[0]) + "," +
                   std::to_string(e[1]) + "); " + t.extent() + ")";
        };
        amrex::Print() << "[driver] SCRIP check " << label << " vs " << group << " ("
                       << (method.empty() ? "map_method unknown" : method) << "): "
                       << both.count << " cells covered by both, " << ours_only.count
                       << " only by ours" << where(ours_only) << ", " << scrip_only.count
                       << " only by SCRIP" << where(scrip_only) << ".\n";
        if (both.count > 0) {
            const Long nx = dst_domain.length(0);
            amrex::Print() << "[driver] SCRIP check " << label << ": weight L1 over cells covered by both"
                           << " mean " << l1_sum / Real(both.count) << " max " << global_max
                           << " at " << dst_name << " (" << dst_domain.smallEnd(0) + l1_first % nx
                           << "," << dst_domain.smallEnd(1) + l1_first / nx << ")"
                           << " (0 = identical, 2 = disjoint stencils).\n";
        }
    };

    const int f = static_cast<int>(WeightFamily::SCALAR_CENTERED);
    compare(m_scrip_a2o_group, *m_weight_mf_by_family[f], *m_index_mf_by_family[f],
            m_geom.Domain(), "A2O scalar/rho", "REMORA");
    if (m_use_two_way_coupling && m_weight_o2a_mf) {
        compare(m_scrip_o2a_group, *m_weight_o2a_mf, *m_index_o2a_mf,
                m_geom_ocn.Domain(), "O2A SST", "ERF");
    }
#endif
}

// Decide the frame, and refuse to guess. ERF reports lon/lat degrees
// (XLONG/XLAT) while REMORA's x_psi/y_psi are meters in a projection the
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
        ERFGetOceanToAtmosCellCenterLonLat(m_erf, &erf_lon, &erf_lat);
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
            // ERF reports cell centres; its corner mesh, built from them, reaches
            // half a cell further on each side.
            amrex::Print() << "[driver] remap lon/lat ranges: "
                           << describe("remora corners(lon,lat)", remora_lon, remora_lat) << " "
                           << describe("erf centers(lon,lat)", erf_lon, erf_lat) << "\n";
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
// ERF reports lon/lat at cell centers (WRF's XLONG/XLAT or metgrid's
// XLONG_M/XLAT_M), so every corner has to be derived. Each node is the average
// of the four cells around it, computed once and shared by all four, so the mesh
// tiles exactly as a conservative remap assumes. This is SCRIP's create_psimesh
// (COAWST Lib/SCRIP_COAWST/read_wrf.f). ERF's ghost cells at the domain boundary
// are clamped copies of the edge cell, and averaging with those would put the
// outer ring of nodes on the edge-cell centers, half a cell inward. So the
// centers are first extended one ring by linear extrapolation, which is SCRIP's
// create_extra_rho_grid, before averaging.
//
// Not handled, as before: a grid whose cells straddle a longitude wrap. The
// average would be wrong there, and so would the clipping downstream.
void
ERFRemoraMultiBlockContainer::BuildERFCornerMesh()
{
    using namespace amrex;

    const MultiFab* erf_lon = nullptr;
    const MultiFab* erf_lat = nullptr;
    ERFGetOceanToAtmosCellCenterLonLat(m_erf, &erf_lon, &erf_lat);
    AMREX_ALWAYS_ASSERT(erf_lon != nullptr && erf_lat != nullptr);

    const iMultiFab* erf_is_land = nullptr;
    ERFGetLandMask(m_erf, &erf_is_land);

    // Reconcile longitude conventions. ERF carries WRF's XLONG, which is
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
    // Setup only: the staging plans and the O2A weight build read it on the host.
    m_erf_node_coords = std::make_unique<MultiFab>(erf_node_ba, erf_dm, 2, IntVect(1,1,0), HostArenaInfo());
    m_erf_node_coords->setVal(0.0);

    const Box node_domain = erf_node_ba.minimalBox();
    const Box cell_domain = erf_cell_ba.minimalBox();

    // Cell centers with one ghost ring: interior ghosts from neighbors, domain
    // ghosts extrapolated (ExtendCornerMeshOneRing is index-generic).
    MultiFab centers(erf_cell_ba, erf_dm, 2, IntVect(1,1,0), HostArenaInfo());
    MultiFab::Copy(centers, *erf_lon, 0, 0, 1, 0);
    MultiFab::Copy(centers, *erf_lat, 0, 1, 1, 0);
    centers.FillBoundary();
    ExtendCornerMeshOneRing(centers, cell_domain);

    RequireHostAccess("BuildERFCornerMesh", *m_erf_node_coords, centers);
    for (MFIter mfi(*m_erf_node_coords); mfi.isValid(); ++mfi) {
        const Box nbx = mfi.validbox();
        auto node = m_erf_node_coords->array(mfi);
        auto const& cc = centers.const_array(mfi);

        for (int j = nbx.smallEnd(1); j <= nbx.bigEnd(1); ++j) {
            for (int i = nbx.smallEnd(0); i <= nbx.bigEnd(0); ++i) {
                for (int c = 0; c < 2; ++c) {
                    node(i,j,0,c) = Real(0.25) * (cc(i-1,j-1,0,c) + cc(i,j-1,0,c) +
                                                  cc(i-1,j  ,0,c) + cc(i,j  ,0,c));
                }
                node(i,j,0,0) += lon_shift;
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
    // Default arena: the per-exchange diagnostics read it on the device, and
    // nothing reads it on the host (the staging plan copies it).
    m_erf_land_mask = std::make_unique<MultiFab>(erf_cell_ba, erf_dm, 1, 0);
    if (erf_is_land != nullptr) {
        amrex::Copy(*m_erf_land_mask, *erf_is_land, 0, 0, 1, 0);
    } else {
        m_erf_land_mask->setVal(0.0);   // no mask reported: treat everything as water
    }

    if (m_remap_verbose >= 1) {
        // The mesh is watertight by construction, so measure whether it sits
        // where the cells are: how far each cell's corner mean lies from the
        // cell center ERF reported. On a smooth grid that is a small curvature
        // term; a mesh misplaced by a staggering mistake shows up here as a
        // sizable fraction of a cell.
        Real max_offset = 0.0;
        MultiFab lon_copy, lat_copy;
        const MultiFab* erf_lon_host = HostView(erf_lon, lon_copy);
        const MultiFab* erf_lat_host = HostView(erf_lat, lat_copy);
        RequireHostAccess("BuildERFCornerMesh", *m_erf_node_coords, *erf_lon_host, *erf_lat_host);
        for (MFIter mfi(*m_erf_node_coords); mfi.isValid(); ++mfi) {
            // The node MultiFab's validbox is nodal; iterate its cells.
            const Box bx = amrex::enclosedCells(mfi.validbox()) & cell_domain;
            auto const& node = m_erf_node_coords->const_array(mfi);
            auto const& lon  = erf_lon_host->const_array(mfi);
            auto const& lat  = erf_lat_host->const_array(mfi);
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    Real mean[2];
                    for (int c = 0; c < 2; ++c) {
                        mean[c] = Real(0.25) * (node(i,j,0,c) + node(i+1,j,0,c) +
                                                node(i,j+1,0,c) + node(i+1,j+1,0,c));
                    }
                    // Undo lon_shift so this measures placement, not the
                    // convention change.
                    const Real dx = (mean[0] - lon_shift) - lon(i,j,0);
                    const Real dy = mean[1] - lat(i,j,0);
                    max_offset = std::max(max_offset, std::sqrt(dx*dx + dy*dy));
                }
            }
        }
        ParallelDescriptor::ReduceRealMax(max_offset);
        amrex::Print() << "[driver] remap ERF corner mesh: nodes=" << node_domain
                       << " max corner-mean offset from cell centers=" << max_offset << " deg"
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
    RequireHostAccess("BuildERFStagingPlan", erf_node);
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
    const MultiFab* remora_lon_dev = nullptr;
    const MultiFab* remora_lat_dev = nullptr;
    REMORAGetAtmosToOceanPsiLonLat(m_remora, &remora_lon_dev, &remora_lat_dev);
    AMREX_ALWAYS_ASSERT(remora_lon_dev != nullptr && remora_lat_dev != nullptr);
    MultiFab lon_copy, lat_copy;
    const MultiFab* remora_lon = HostView(remora_lon_dev, lon_copy);
    const MultiFab* remora_lat = HostView(remora_lat_dev, lat_copy);

    const int n_remora_boxes = static_cast<int>(m_ba_ocn.size());
    Vector<Real> remora_lo(2 * n_remora_boxes,  std::numeric_limits<Real>::max());
    Vector<Real> remora_hi(2 * n_remora_boxes, -std::numeric_limits<Real>::max());
    RequireHostAccess("BuildERFStagingPlan", *remora_lon, *remora_lat);
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
    // Read by the A2O weight build's host loops.
    plan.node_mf = std::make_unique<MultiFab>(plan.node_ba, m_dm_ocn, 2, 0, HostArenaInfo());
    plan.land_mf = std::make_unique<MultiFab>(plan.cell_ba, m_dm_ocn, 1, 0, HostArenaInfo());
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
    MultiFab lon_copy, lat_copy;
    const MultiFab* remora_lon_host = HostView(remora_lon, lon_copy);
    const MultiFab* remora_lat_host = HostView(remora_lat, lat_copy);
    RequireHostAccess("BuildREMORAStagingPlan", *remora_lon_host, *remora_lat_host);
    for (MFIter mfi(*remora_lon_host); mfi.isValid(); ++mfi) {
        const int b = mfi.index();
        const Box nbx = mfi.validbox();
        auto const& lon = remora_lon_host->const_array(mfi);
        auto const& lat = remora_lat_host->const_array(mfi);
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
    RequireHostAccess("BuildREMORAStagingPlan", erf_node);
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
    // Read by the O2A weight build's host loops.
    plan.node_mf = std::make_unique<MultiFab>(plan.node_ba, erf_dm, 2, 0, HostArenaInfo());
    plan.wet_mf = std::make_unique<MultiFab>(plan.cell_ba, erf_dm, 1, 0, HostArenaInfo());
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

    // Everything read from REMORA below goes through a host view: the loops that
    // read it run on the host.
    const MultiFab* remora_lon_dev = nullptr;
    const MultiFab* remora_lat_dev = nullptr;
    if (m_remap_spherical) {
        REMORAGetAtmosToOceanPsiLonLat(m_remora, &remora_lon_dev, &remora_lat_dev);
    }
    MultiFab lon_copy, lat_copy;
    const MultiFab* remora_lon = HostView(remora_lon_dev, lon_copy);
    const MultiFab* remora_lat = HostView(remora_lat_dev, lat_copy);

    // REMORA wet/dry masks, so coverage can be reported over the cells that carry
    // ocean rather than over land the atmosphere grid was never going to force.
    const MultiFab* remora_mskr_dev = nullptr;
    const MultiFab* remora_msku_dev = nullptr;
    const MultiFab* remora_mskv_dev = nullptr;
    if (m_remora != nullptr) {
        REMORAGetMasks(m_remora, &remora_mskr_dev, &remora_msku_dev, &remora_mskv_dev);
    }
    MultiFab mskr_copy, msku_copy, mskv_copy;
    const MultiFab* remora_mskr = HostView(remora_mskr_dev, mskr_copy);
    const MultiFab* remora_msku = HostView(remora_msku_dev, msku_copy);
    const MultiFab* remora_mskv = HostView(remora_mskv_dev, mskv_copy);

    // Reset A2O weight and index caches, one family slot at a time. Each slot
    // was allocated (AllocateCachesFromSolverLayouts) against its own
    // destination BoxArray/stagger (rho for SCALAR_CENTERED/WIND_FACE_TO_CENTER,
    // u-face for TAUX_FACE, v-face for TAUY_FACE), so building weights
    // per-slot means TauX/TauY finally get stencils computed against their own
    // stagger instead of borrowing the rho-point stencil.
    static const std::array<amrex::IntVect, kNumWeightFamilies> kExpectedFamilyStagger = {
        amrex::IntVect(0,0,0), amrex::IntVect(1,0,0), amrex::IntVect(0,1,0), amrex::IntVect(0,0,0)
    };
    static const std::array<const char*, kNumWeightFamilies> kFamilyLabel = {
        "A2O scalar/rho", "A2O tau_x/u-face", "A2O tau_y/v-face", "A2O wind->rho"
    };
    static const std::array<const char*, kNumWeightFamilies> kExpectedRemoraMask = {
        "mskr", "msku", "mskv", "mskr"
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

        // Allocate the donor-fill flag array only when the donor search can
        // actually run. Under the default driver.remap_donor_fill=withhold, or in
        // the cartesian frame, nothing ever writes it. A null array here is
        // therefore normal, not a failure, and must NOT skip the family the way
        // a missing weight or index array would -- the remap itself is unaffected.
        // NormalizeWeightsAndCheckConservation already treats a null donor_fill as
        // "this direction cannot donor-fill" and omits the donor_filled= field,
        // which is exactly right here; ResolveDonorFillMode has already reported
        // the mode, so the reader is not left guessing why it is absent.
        const bool donor_fill_can_run =
            (m_remap_donor_fill_mode == DonorFillMode::Coawst) && m_remap_spherical;
        if (donor_fill_can_run && !m_donor_fill_mf_by_family[f]) {
            m_donor_fill_mf_by_family[f] = std::make_unique<iMultiFab>(
                weight_mf->boxArray(), weight_mf->DistributionMap(), 1, 0, HostArenaInfo());
        }
        iMultiFab* donor_fill_mf = m_donor_fill_mf_by_family[f].get();

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            weight_mf->boxArray().ixType().toIntVect() == kExpectedFamilyStagger[f],
            "BuildRemapWeights: A2O weight MultiFab stagger does not match its WeightFamily slot.");

        weight_mf->setVal(0.0);
        index_mf->setVal(-1);
        if (donor_fill_mf) { donor_fill_mf->setVal(0); }

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
        // Each family's mask must sit on that family's stagger, or const_array()
        // would be reading a differently shaped array.
        const MultiFab* family_mask = remora_mskr;
        if (f == static_cast<int>(WeightFamily::TAUX_FACE))      { family_mask = remora_msku; }
        else if (f == static_cast<int>(WeightFamily::TAUY_FACE)) { family_mask = remora_mskv; }
        const bool donor_fill_needs_mask =
            m_remap_donor_fill_mode == DonorFillMode::Coawst && m_remap_spherical;
        if (family_mask == nullptr) {
            if (donor_fill_needs_mask && ParallelDescriptor::IOProcessor()) {
                amrex::Print() << "[driver] WARNING remap (" << kFamilyLabel[f]
                               << "): REMORA supplied no " << kExpectedRemoraMask[f]
                               << " wet mask. Donor fill is inactive for this family; wet cells "
                               << "without usable source overlap will take the same per-lane "
                               << "fallback they would under driver.remap_donor_fill=withhold.\n";
            }
        } else if (family_mask->boxArray().ixType() != weight_mf->boxArray().ixType()) {
            if (donor_fill_needs_mask && ParallelDescriptor::IOProcessor()) {
                amrex::Print() << "[driver] WARNING remap (" << kFamilyLabel[f] << "): REMORA "
                               << kExpectedRemoraMask[f] << " wet mask ixType="
                               << family_mask->boxArray().ixType().toIntVect()
                               << " does not match expected ixType="
                               << weight_mf->boxArray().ixType().toIntVect()
                               << ". Donor fill is inactive for this family; wet cells without "
                               << "usable source overlap will take the same per-lane fallback "
                               << "they would under driver.remap_donor_fill=withhold.\n";
            }
            family_mask = nullptr;   // stagger mismatch: report unmasked rather than lie
        }

        RequireHostAccess("BuildAtmosToOceanWeights", *weight_mf, *index_mf);
        RequireHostAccessIfPresent("BuildAtmosToOceanWeights", donor_fill_mf, family_mask);
        if (m_remap_spherical) {
            RequireHostAccess("BuildAtmosToOceanWeights", *remora_lon, *remora_lat,
                              *m_staging_plan->node_mf, *m_staging_plan->land_mf);
        }
        for (MFIter mfi(*weight_mf); mfi.isValid(); ++mfi) {
            Box bx = mfi.validbox();
            auto weight_arr = weight_mf->array(mfi);
            auto index_arr = index_mf->array(mfi);
            // Null under withhold / cartesian; the donor search that writes it is
            // gated on exactly the condition that allocated it, so an empty Array4
            // here is never dereferenced.
            auto donor_fill_arr = donor_fill_mf ? donor_fill_mf->array(mfi)
                                                : amrex::Array4<int>{};

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
                                    const Real area =
                                        ClippedArea(remora_poly, erf_src_poly(erf_i, erf_j));
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
                                const Real area = ClippedArea(remora_poly, erf_poly);
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
                                const Real area = ClippedArea(remora_poly, erf_poly);
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
                            // Mark the cell so the coverage report can separate it
                            // from a genuinely covered one. Both this flag and
                            // donor_filled count the same event; they differ in
                            // that donor_filled is a raw per-box tally used only
                            // for the one-line summary, whereas this flag is
                            // de-duplicated across box seams downstream.
                            donor_fill_arr(i, j, 0) = 1;
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
        NormalizeWeightsAndCheckConservation(*weight_mf, *index_mf, family_mask,
                                             donor_fill_mf, kFamilyLabel[f]);
    }

    ParallelDescriptor::ReduceIntMax(max_stencil_used);
    ParallelDescriptor::ReduceLongSum(dropped_overlaps);
    ParallelDescriptor::ReduceLongSum(donor_filled);
    if (donor_filled > 0) {
        // Summed over all four families and incremented once per box visit, so for
        // TAUX_FACE/TAUY_FACE a donor-filled face on a box seam is counted once
        // per adjacent box that fills it and this total moves with the
        // decomposition. The per-family donor_filled= field on the coverage lines
        // above applies the same seam-ownership rule as every other coverage
        // count, so that is the number to compare between runs.
        amrex::Print() << "[driver] remap: " << donor_filled
                       << " wet destination cells had no valid source overlap and were filled from "
                       << "their nearest non-land source cell (COAWST check_weights behaviour). "
                       << "These carry a neighbour's atmospheric state rather than their own. "
                       << "Counted per box visit, so a donor-filled face on a box seam can be "
                       << "counted once per adjacent box; the per-family donor_filled= counts "
                       << "above are seam-de-duplicated.\n";
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

        // Every destination cell gets one OwnerMatchStatus, and only Coupled cells
        // are valid. Water cells that are not Coupled keep ERF's own SST, but the
        // reasons differ and only some are suspicious: WaterOutside is expected
        // around the edge of the ocean grid; WaterOceanDry and LandOceanWet are the
        // two directions of mask disagreement, legitimate along a real coastline
        // but also exactly what a grid misalignment or a stale mask looks like.
        // ESMF draws the same first distinction in ESMF_REGRIDSTATUS, which
        // separates OUTSIDE from SRCMASKED for this reason.
        //
        // LandOceanWet cells keep their weights, since REMORA's pack already zeroes
        // ERF land, but are not valid: the coverage flag then means "ERF water with
        // a wet ocean donor" for every consumer, including ERF's range check.
        const Box dst_domain = m_valid_o2a_mf->boxArray().minimalBox();
        std::array<OwnerMatchTally, static_cast<int>(OwnerMatchStatus::NumStatus)> owner_match;

        // ERF's land mask, read from ERF directly: m_erf_land_mask exists only in
        // the spherical frame. Copied onto the O2A layout and checked, since a copy
        // between non-intersecting layouts would silently report every cell as water.
        const iMultiFab* erf_is_land = nullptr;
        if (m_erf != nullptr) { ERFGetLandMask(m_erf, &erf_is_land); }
        iMultiFab land_dst(m_ba, m_dm, 1, 0, HostArenaInfo());
        land_dst.setVal(0);
        if (erf_is_land != nullptr) {
            land_dst.ParallelCopy(*erf_is_land, 0, 0, 1);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                land_dst.sum(0) == erf_is_land->sum(0),
                "BuildOceanToAtmosWeights: ERF land mask did not copy onto the O2A layout.");
        }

        auto classify = [&] (bool is_land, bool any_wet, bool any_geometric, Real covered_frac)
        {
            const bool enough = any_wet && covered_frac >= m_remap_min_coverage_frac;
            if (is_land)       { return enough ? OwnerMatchStatus::LandOceanWet : OwnerMatchStatus::Land; }
            if (!any_wet)      { return any_geometric ? OwnerMatchStatus::WaterOceanDry
                                                      : OwnerMatchStatus::WaterOutside; }
            return enough ? OwnerMatchStatus::Coupled : OwnerMatchStatus::WaterThin;
        };

        // Shared by the spherical and Cartesian loops once `overlaps` holds this
        // destination cell's wet REMORA overlaps: classify, flag, store the stencil.
        //
        // Coverage is thresholded on area rather than on whether any overlap exists
        // at all. A sliver overlap is worse than no overlap: normalization cannot
        // turn 1e-6 destination coverage into a trustworthy ocean answer, so it
        // keeps ERF's own lower-boundary SST. Thresholds of this kind are standard:
        // CESM drops mapping below a destination fraction of 0.001, OASIS/CERFACS
        // at 1/1000. classify() tests the overlap count as well as the fraction:
        // with remap_min_coverage_frac set to 0, a cell no ocean cell touches has
        // covered_frac == 0, and 0 >= 0 would otherwise report it as covered.
        // ERF land is never valid, even under a wet ocean cell.
        auto finish_cell = [&] (int i, int j, Real dst_area, Real covered_frac,
                                bool any_geometric_overlap, bool is_land,
                                Array4<Real> const& weight_arr, Array4<int> const& index_arr,
                                Array4<int> const& valid_arr)
        {
            const OwnerMatchStatus status =
                classify(is_land, !overlaps.empty(), any_geometric_overlap, covered_frac);
            owner_match[static_cast<int>(status)].add(dst_domain, i, j);
            if (overlaps.empty()) { return; }

            valid_arr(i, j, 0) = (status == OwnerMatchStatus::Coupled) ? 1 : 0;

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
        };

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

            RequireHostAccess("BuildOceanToAtmosWeights", *m_weight_o2a_mf, *m_index_o2a_mf,
                              *m_valid_o2a_mf, land_dst, *m_erf_node_coords,
                              *m_remora_staging_plan->node_mf, *m_remora_staging_plan->wet_mf);
            for (MFIter mfi(*m_weight_o2a_mf); mfi.isValid(); ++mfi) {
                Box bx = mfi.validbox();
                auto weight_arr = m_weight_o2a_mf->array(mfi);
                auto index_arr = m_index_o2a_mf->array(mfi);
                auto valid_arr = m_valid_o2a_mf->array(mfi);
                auto const& land_arr = land_dst.const_array(mfi);
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
                        if (dst_area <= 0.0) {
                            owner_match[static_cast<int>(OwnerMatchStatus::Degenerate)]
                                .add(dst_domain, i, j);
                            continue;
                        }

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
                            const Real intersect_area = ClippedArea(erf_poly, remora_poly);
                            if (intersect_area <= 0.0) { continue; }

                            any_geometric_overlap = true;
                            if (remora_wet(src_i, src_j, 0) <= Real(0.5)) { continue; }

                            covered_frac += intersect_area / dst_area;
                            overlaps.emplace_back(intersect_area, IntVect(src_i, src_j, 0));
                        }

                        finish_cell(i, j, dst_area, covered_frac, any_geometric_overlap,
                                    land_arr(i,j,0) != 0, weight_arr, index_arr, valid_arr);
                    }
                }
            }
        } else {
            // Cartesian O2A still uses a replicated rho-mask lookup because it has
            // no staged REMORA source layout. Spherical O2A consumes
            // m_remora_staging_plan->wet_mf instead.
            const MultiFab* remora_mskr_dev = nullptr;
            if (m_remora != nullptr) {
                REMORAGetMasks(m_remora, &remora_mskr_dev, nullptr, nullptr);
            }
            MultiFab mskr_copy;
            const MultiFab* remora_mskr = HostView(remora_mskr_dev, mskr_copy);

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

            RequireHostAccess("BuildOceanToAtmosWeights", *m_weight_o2a_mf, *m_index_o2a_mf,
                              *m_valid_o2a_mf, land_dst);
            for (MFIter mfi(*m_weight_o2a_mf); mfi.isValid(); ++mfi) {
                Box bx = mfi.validbox();
                auto weight_arr = m_weight_o2a_mf->array(mfi);
                auto index_arr = m_index_o2a_mf->array(mfi);
                auto valid_arr = m_valid_o2a_mf->array(mfi);
                auto const& land_arr = land_dst.const_array(mfi);

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
                        if (dst_area <= 0.0) {
                            owner_match[static_cast<int>(OwnerMatchStatus::Degenerate)]
                                .add(dst_domain, i, j);
                            continue;
                        }

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

                                Real intersect_area = ClippedArea(erf_poly, remora_poly);
                                if (intersect_area <= 0.0) { continue; }

                                any_geometric_overlap = true;
                                if (!remora_rho_is_wet(src_i, src_j)) { continue; }

                                covered_frac += intersect_area / dst_area;
                                overlaps.emplace_back(intersect_area, IntVect(src_i, src_j, 0));
                            }
                        }

                        finish_cell(i, j, dst_area, covered_frac, any_geometric_overlap,
                                    land_arr(i,j,0) != 0, weight_arr, index_arr, valid_arr);
                    }
                }
            }
        }

        // No wet mask and no donor flag: O2A reports coverage over all ERF cells,
        // and the donor search is an atmosphere->ocean-only mechanism, so there is
        // no donor category to separate out here.
        // Its counts are stencils, which ERF-land cells under wet ocean keep; the
        // owner-match table below is the count of cells that take the coupled SST.
        NormalizeWeightsAndCheckConservation(*m_weight_o2a_mf, *m_index_o2a_mf,
                                             nullptr, nullptr,
                                             "O2A SST stencils; coupled cells in owner match below");

        for (auto& t : owner_match) { t.reduce(); }

        // Where an example cell is, for the log. Collective: the owning rank
        // supplies the value and the others a sentinel.
        auto locate = [&] (const IntVect& ij) -> std::string {
            Real a = std::numeric_limits<Real>::max();
            Real b = std::numeric_limits<Real>::max();
            for (MFIter mfi(*m_weight_o2a_mf); mfi.isValid(); ++mfi) {
                if (!mfi.validbox().contains(ij)) { continue; }
                if (m_remap_spherical) {
                    auto const& node = m_erf_node_coords->const_array(mfi);
                    const int i = ij[0], j = ij[1];
                    a = Real(0.25) * (node(i,j,0,0) + node(i+1,j,0,0) + node(i,j+1,0,0) + node(i+1,j+1,0,0));
                    b = Real(0.25) * (node(i,j,0,1) + node(i+1,j,0,1) + node(i,j+1,0,1) + node(i+1,j+1,0,1));
                } else {
                    a = m_geom.ProbLo(0) + (ij[0] + Real(0.5)) * m_geom.CellSize(0);
                    b = m_geom.ProbLo(1) + (ij[1] + Real(0.5)) * m_geom.CellSize(1);
                }
            }
            ParallelDescriptor::ReduceRealMin(a);
            ParallelDescriptor::ReduceRealMin(b);
            std::ostringstream os;
            os << "ERF (" << ij[0] << "," << ij[1] << ") "
               << (m_remap_spherical ? "lon " : "x ") << a
               << (m_remap_spherical ? " lat " : " y ") << b;
            return os.str();
        };

        // Reported unconditionally: a cell that is not Coupled is invisible in the
        // SST field, because ERF simply keeps its own value there, so the only way a
        // mask disagreement surfaces is if the split is written down. A large
        // disagreement count against a short coastline is the signature of a grid
        // misalignment or a stale ocean mask rather than of correct withholding.
        const std::string thin_effect = "keeps ERF SST; wet fraction < remap_min_coverage_frac="
                                      + std::to_string(m_remap_min_coverage_frac);
        struct Line { OwnerMatchStatus s; const char* what; std::string effect; };
        const Line lines[] = {
            {OwnerMatchStatus::Coupled,       "ERF water, wet ocean donor",    "coupled SST"},
            {OwnerMatchStatus::WaterOutside,  "ERF water, outside ocean grid", "keeps ERF SST"},
            {OwnerMatchStatus::WaterOceanDry, "ERF water, REMORA dry",         "keeps ERF SST; mask disagreement"},
            {OwnerMatchStatus::WaterThin,     "ERF water, thin wet overlap",   thin_effect},
            {OwnerMatchStatus::LandOceanWet,  "ERF land, REMORA wet",          "SST dropped; mask disagreement"},
            {OwnerMatchStatus::Land,          "ERF land, no/thin wet ocean",   "masks agree"},
            {OwnerMatchStatus::Degenerate,    "zero-area ERF cell",            "skipped; degenerate corner mesh"},
        };
        const amrex::Long n_cells = m_valid_o2a_mf->boxArray().numPts();
        amrex::Long n_classified = 0;
        for (const auto& t : owner_match) { n_classified += t.count; }
        amrex::Print() << "[driver] O2A SST owner match over " << n_cells << " ERF surface cells"
                       << (m_remap_spherical ? " (lon/lat in REMORA's convention)" : " (x/y in m)")
                       << ":\n";
        for (const auto& l : lines) {
            const auto& t = owner_match[static_cast<int>(l.s)];
            // Rank-uniform: counts are reduced, so every rank takes the same branch.
            const bool show_example = t.count > 0 && l.s != OwnerMatchStatus::Coupled
                                                  && l.s != OwnerMatchStatus::Land;
            const std::string where = show_example ? locate(t.example(dst_domain)) : std::string{};
            amrex::Print() << "[driver]   " << std::left << std::setw(30) << l.what << std::right
                           << std::setw(9) << t.count << "  " << l.effect;
            if (show_example) { amrex::Print() << " (e.g. " << where << "; " << t.extent() << ")"; }
            amrex::Print() << "\n";
        }
        if (n_classified != n_cells) {
            amrex::Print() << "[driver] WARNING O2A owner match: " << (n_cells - n_classified)
                           << " ERF surface cells were not classified.\n";
        }

        // 1 where the cell is not Coupled, for the O2A min/max diagnostics, which
        // otherwise count the zero-filled cells as ocean.
        m_o2a_uncoupled_mask = std::make_unique<MultiFab>(m_ba, m_dm, 1, 0);
        for (MFIter mfi(*m_o2a_uncoupled_mask); mfi.isValid(); ++mfi) {
            auto const& valid = m_valid_o2a_mf->const_array(mfi);
            auto const& mask = m_o2a_uncoupled_mask->array(mfi);
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                mask(i,j,k) = (valid(i,j,k) != 0) ? Real(0.0) : Real(1.0);
            });
        }

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
                             const amrex::Vector<std::string>& names,
                             bool ocean_to_atmos) const
{
    // m_driver_verbose comes from ParmParse and is identical on every rank, so
    // this is a rank-uniform gate around the collective reductions below.
    if (m_driver_verbose < 1) { return; }

    const amrex::MultiFab* remora_mskr = nullptr;
    const amrex::MultiFab* remora_msku = nullptr;
    const amrex::MultiFab* remora_mskv = nullptr;
    if (!ocean_to_atmos && m_remora != nullptr) {
        REMORAGetMasks(m_remora, &remora_mskr, &remora_msku, &remora_mskv);
    }

    auto mask_for = [&] (const amrex::MultiFab& field) -> const amrex::MultiFab* {
        const amrex::MultiFab* candidate = nullptr;
        if (ocean_to_atmos) {
            // Not ERF land alone: water cells with no ocean donor are zero-filled too.
            candidate = m_o2a_uncoupled_mask ? m_o2a_uncoupled_mask.get()
                                             : m_erf_land_mask.get();
        } else {
            const amrex::IntVect stagger = field.boxArray().ixType().toIntVect();
            if (stagger == amrex::IntVect(1,0,0)) {
                candidate = remora_msku;
            } else if (stagger == amrex::IntVect(0,1,0)) {
                candidate = remora_mskv;
            } else if (stagger == amrex::IntVect::TheCellVector()) {
                candidate = remora_mskr;
            }
        }
        // MaskedMinMax indexes the mask with an MFIter built on the field, which is
        // only defined when the two share a BoxArray *and* a DistributionMapping --
        // matching ixType alone is not enough. The A2O masks come off the same
        // REMORA layout as their fields, but m_erf_land_mask is built on the ERF
        // solver's cell layout while the O2A cache lives on the driver's, so the
        // O2A pairing is the one that can differ. Fall back to an unmasked
        // reduction (reported as such) rather than reading across layouts.
        if (candidate != nullptr &&
            (candidate->boxArray() != field.boxArray() ||
             candidate->DistributionMap() != field.DistributionMap())) {
            candidate = nullptr;
        }
        return candidate;
    };

    amrex::Print() << label << "\n";
    for (int idx = 0; idx < static_cast<int>(fields.size()); ++idx) {
        if (fields[idx] == nullptr) {
            amrex::Print() << "  " << names[idx] << ": withheld (not supplied by ERF)\n";
            continue;
        }
        const amrex::MultiFab* mask = mask_for(*fields[idx]);
        const auto extrema = MaskedMinMax(*fields[idx], mask, ocean_to_atmos);
        const char* marker = (mask != nullptr) ? " (wet only)" : " (unmasked)";
        if (extrema.included_cells == 0) {
            amrex::Print() << "  " << names[idx] << ": no wet cells" << marker << "\n";
            continue;
        }
        if (!extrema.has_finite()) {
            amrex::Print() << "  " << names[idx] << ": all " << extrema.included_cells
                           << " cells non-finite" << marker << "\n";
            continue;
        }
        amrex::Print() << "  " << names[idx] << ": min=" << extrema.min_value
                       << " max=" << extrema.max_value;
        if (extrema.nonfinite_cells != 0) {
            amrex::Print() << " non-finite=" << extrema.nonfinite_cells;
        }
        amrex::Print() << marker << "\n";
        if (names[idx] == "rain" &&
            (extrema.min_value > 100.0 || extrema.max_value > 100.0)) {
            amrex::Print() << "  *** WARNING: implausibly large rain values detected! ***\n";
        }
    }
}

void
ERFRemoraMultiBlockContainer::ValidateFields (const char* label,
                                              const amrex::Vector<amrex::MultiFab*>& fields,
                                              const amrex::Vector<std::string>& names,
                                              amrex::Real time,
                                              long long ocean_step,
                                              bool ocean_to_atmos) const
{
    if (fields.size() != names.size()) {
        amrex::Abort("ValidateFields field/name count mismatch.");
    }
    // Rank-uniform: m_driver_verbose is a ParmParse value shared by all ranks.
    if (m_driver_verbose < 1) { return; }

    const amrex::MultiFab* remora_mskr = nullptr;
    const amrex::MultiFab* remora_msku = nullptr;
    const amrex::MultiFab* remora_mskv = nullptr;
    if (!ocean_to_atmos && m_remora != nullptr) {
        REMORAGetMasks(m_remora, &remora_mskr, &remora_msku, &remora_mskv);
    }

    auto mask_for = [&] (const amrex::MultiFab& field) -> const amrex::MultiFab* {
        const amrex::MultiFab* candidate = nullptr;
        if (ocean_to_atmos) {
            // Not ERF land alone: water cells with no ocean donor are zero-filled too.
            candidate = m_o2a_uncoupled_mask ? m_o2a_uncoupled_mask.get()
                                             : m_erf_land_mask.get();
        } else {
            const amrex::IntVect stagger = field.boxArray().ixType().toIntVect();
            if (stagger == amrex::IntVect(1,0,0)) {
                candidate = remora_msku;
            } else if (stagger == amrex::IntVect(0,1,0)) {
                candidate = remora_mskv;
            } else if (stagger == amrex::IntVect::TheCellVector()) {
                candidate = remora_mskr;
            }
        }
        // MaskedMinMax indexes the mask with an MFIter built on the field, which is
        // only defined when the two share a BoxArray *and* a DistributionMapping --
        // matching ixType alone is not enough. The A2O masks come off the same
        // REMORA layout as their fields, but m_erf_land_mask is built on the ERF
        // solver's cell layout while the O2A cache lives on the driver's, so the
        // O2A pairing is the one that can differ. Fall back to an unmasked
        // reduction (reported as such) rather than reading across layouts.
        if (candidate != nullptr &&
            (candidate->boxArray() != field.boxArray() ||
             candidate->DistributionMap() != field.DistributionMap())) {
            candidate = nullptr;
        }
        return candidate;
    };

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

        const amrex::MultiFab* mask = mask_for(*mf);
        const auto extrema = MaskedMinMax(*mf, mask, ocean_to_atmos);
        const char* marker = (mask != nullptr) ? " (wet only)" : " (unmasked)";

        if (amrex::ParallelDescriptor::IOProcessor()) {
            amrex::Print() << label
                           << " ocean_step=" << ocean_step
                           << " time=" << time
                           << " field=" << names[idx]
                           << " nboxes=" << mf->boxArray().size();
            if (extrema.included_cells == 0) {
                amrex::Print() << " no wet cells" << marker << "\n";
            } else if (!extrema.has_finite()) {
                amrex::Print() << " all " << extrema.included_cells
                               << " cells non-finite" << marker << "\n";
            } else {
                amrex::Print() << " min=" << extrema.min_value
                               << " max=" << extrema.max_value;
                if (extrema.nonfinite_cells != 0) {
                    amrex::Print() << " non-finite=" << extrema.nonfinite_cells;
                }
                amrex::Print() << marker << "\n";
            }
        }

        if (extrema.nonfinite_cells != 0) {
            amrex::Abort(std::string(label) + ": " + std::to_string(extrema.nonfinite_cells) +
                         " non-finite value(s) in field '" + names[idx] + "'" + marker +
                         " at ocean_step=" + std::to_string(ocean_step));
        }
    }
}

void
ERFRemoraMultiBlockContainer::AbortOnNonFinite (const char* what,
                                                const amrex::Vector<amrex::MultiFab*>& lanes,
                                                const amrex::Vector<std::string>& names) const
{
    AMREX_ALWAYS_ASSERT(lanes.size() == names.size());
    // Every rank holds the same lane vector, so every rank enters the same
    // collectives in the same order.
    for (int idx = 0; idx < static_cast<int>(lanes.size()); ++idx) {
        if (lanes[idx] != nullptr && ContainsNonFinite(*lanes[idx])) {
            amrex::Abort(std::string("Non-finite value detected in ") + what +
                         " output, lane '" + names[idx] + "'.");
        }
    }
}

void
ERFRemoraMultiBlockContainer::CheckDebugInjectNaNLane () const
{
    if (m_debug_inject_nan_lane.empty()) { return; }
    const auto& o2a = coupling::OceanToAtmosNames();
    const bool is_o2a = m_use_two_way_coupling &&
        std::find(o2a.begin(), o2a.end(), m_debug_inject_nan_lane) != o2a.end();
    const bool is_a2o =
        std::find(m_atmos_to_ocean_active_names.begin(), m_atmos_to_ocean_active_names.end(),
                  m_debug_inject_nan_lane) != m_atmos_to_ocean_active_names.end();
    if (!is_o2a && !is_a2o) {
        amrex::Abort("driver.debug_inject_nan_lane='" + m_debug_inject_nan_lane +
                     "' is not an exchange lane of this run.");
    }
    if (is_o2a && m_debug_inject_nan_after_guard) {
        amrex::Abort("driver.debug_inject_nan_after_guard applies to atmos->ocean lanes only; "
                     "nothing checks an ocean->atmos lane between its pack and its apply.");
    }
    amrex::Print() << "WARNING: driver.debug_inject_nan_lane=" << m_debug_inject_nan_lane
                   << " will poison every exchange; this run is expected to abort.\n";
}

void
ERFRemoraMultiBlockContainer::MaybeInjectDebugNaN (const amrex::Vector<amrex::MultiFab*>& lanes,
                                                   const amrex::Vector<std::string>& names,
                                                   bool after_guard) const
{
    if (m_debug_inject_nan_lane.empty() || after_guard != m_debug_inject_nan_after_guard) {
        return;
    }
    for (int idx = 0; idx < static_cast<int>(names.size()); ++idx) {
        if (names[idx] != m_debug_inject_nan_lane) { continue; }
        if (lanes[idx] == nullptr) {
            amrex::Abort("driver.debug_inject_nan_lane='" + m_debug_inject_nan_lane +
                         "' names a withheld lane.");
        }
        // Behind the post-pack guard only ValidateFields' wet-cell count can
        // see the NaN, so the cell has to be wet under the mask it will use.
        // Ahead of it the unmasked guard sees any cell.
        const amrex::MultiFab* wet_mask = nullptr;
        if (after_guard && m_remora != nullptr) {
            const amrex::MultiFab* mskr = nullptr;
            const amrex::MultiFab* msku = nullptr;
            const amrex::MultiFab* mskv = nullptr;
            REMORAGetMasks(m_remora, &mskr, &msku, &mskv);
            const amrex::IntVect stagger = lanes[idx]->boxArray().ixType().toIntVect();
            wet_mask = (stagger == amrex::IntVect(1,0,0)) ? msku
                     : (stagger == amrex::IntVect(0,1,0)) ? mskv
                     : (stagger == amrex::IntVect::TheCellVector()) ? mskr : nullptr;
            // Same rule as ValidateFields: a mask on another layout is not used.
            if (wet_mask != nullptr &&
                (wet_mask->boxArray() != lanes[idx]->boxArray() ||
                 wet_mask->DistributionMap() != lanes[idx]->DistributionMap())) {
                wet_mask = nullptr;
            }
        }
        PoisonOneCellWithNaN(*lanes[idx], wet_mask);
        return;
    }
}

void
ERFRemoraMultiBlockContainer::WriteExchangePlotfiles (amrex::Real o2a_time,
                                                      amrex::Real a2o_time,
                                                      amrex::Real dt_ocn,
                                                      long long ocean_step) const
{
    // The two directions are sampled at different instants once the atmosphere
    // substeps run before the A2O pack: ocean->atmos is packed from REMORA at the
    // start of the ocean interval, atmos->ocean after ERF has advanced across it.
    // Stamping both with one time silently mislabels the O2A output by a full ocean
    // step -- invisible to every test here, because each fcompare is a self-comparison
    // in which both sides move together.
    if (m_use_two_way_coupling) {
        ValidateFields("Exchange validation ocean->atmos",
                       m_ocean_to_atmos, coupling::OceanToAtmosNames(), o2a_time, ocean_step,
                       true);
    }
    ValidateFields("Exchange validation atmos->ocean",
                   m_atmos_to_ocean_active, m_atmos_to_ocean_active_names, a2o_time, ocean_step,
                   false);

    if (m_use_two_way_coupling) {
        WriteExchangePlotfile(m_plot_file_exchange_o2a,
                              m_ocean_to_atmos,
                              coupling::OceanToAtmosNames(),
                              m_geom,
                              o2a_time,
                              o2a_time,
                              dt_ocn,
                              ocean_step,
                              true);
    }
    WriteExchangePlotfile(m_plot_file_exchange_a2o,
                          m_atmos_to_ocean_active,
                          m_atmos_to_ocean_active_names,
                          m_geom_ocn,
                          a2o_time,
                          o2a_time,
                          dt_ocn,
                          ocean_step,
                          false);
}

void
ERFRemoraMultiBlockContainer::WriteExchangePlotfile (
    const std::string& prefix,
    const amrex::Vector<amrex::MultiFab*>& fields,
    const amrex::Vector<std::string>& names,
    const amrex::Geometry& geom,
    amrex::Real time,
    amrex::Real ocean_time,
    amrex::Real dt_ocn,
    long long ocean_step,
    bool ocean_to_atmos) const
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
    // Written after the plotfile so the directory exists.
    WriteExchangeInfo(plotfilename, fields, names, time, ocean_time, dt_ocn, ocean_step,
                      ocean_to_atmos);
}

void
ERFRemoraMultiBlockContainer::WriteExchangeInfo (const std::string& plotfilename,
                                                 const amrex::Vector<amrex::MultiFab*>& fields,
                                                 const amrex::Vector<std::string>& names,
                                                 amrex::Real time,
                                                 amrex::Real ocean_time,
                                                 amrex::Real dt_ocn,
                                                 long long ocean_step,
                                                 bool ocean_to_atmos) const
{
    if (!amrex::ParallelDescriptor::IOProcessor()) { return; }

    std::ofstream os(plotfilename + "/exchange_info", std::ios::trunc);
    if (!os.good()) { return; }

    const amrex::Real t_end = ocean_time + dt_ocn;
    os << "=== ERF/REMORA driver exchange ===\n"
       << "direction          : " << (ocean_to_atmos ? "ocean->atmos" : "atmos->ocean") << "\n"
       << "ocean_step         : " << ocean_step << "\n"
       << "plotfile time      : " << time << "\n"
       << "ocean interval     : [" << ocean_time << ", " << t_end << "]\n\n";

    os << "--- coupling ---\n"
       << "atm2ocn contract   : "
       << coupling::AtmosToOceanContractName(m_atmos_to_ocean_contract) << "\n"
       << "two-way coupling   : " << (m_use_two_way_coupling ? "yes" : "no") << "\n"
       << "remap frame        : " << (m_remap_spherical ? "spherical" : "cartesian") << "\n"
       << "donor fill mode    : "
       << (m_remap_donor_fill_mode == DonorFillMode::Coawst ? "coawst" : "withhold") << "\n"
       << "atm steps/exchange : " << m_atm_steps_per_exchange
       << "  (erf dt=" << m_erf_dt << ")\n"
       << "ocn steps/exchange : " << m_ocean_steps_per_exchange
       << "  (ocean dt=" << m_ocean_dt << ")\n\n";

    if (ocean_to_atmos) {
        os << "--- lanes ---\n"
           << "The ocean->atmos direction is packed from REMORA at the start of the ocean\n"
           << "interval and is never time-averaged, so every lane below represents t="
           << ocean_time << ".\n\n";
        for (int i = 0; i < static_cast<int>(names.size()); ++i) {
            const bool present = (i < static_cast<int>(fields.size()) && fields[i] != nullptr);
            os << "  " << names[i] << (present ? "" : "   [withheld]") << "\n";
        }
        return;
    }

    os << "--- time averaging ---\n"
       << "driver.time_average_atm2ocn_state : " << (m_time_average_atm2ocn_state ? 1 : 0) << "\n"
       << "driver.time_average_atm2ocn_flux  : " << (m_time_average_atm2ocn_flux ? 1 : 0) << "\n"
       << "any lane averaged                 : "
       << (AtmosToOceanAveragingActive() ? "yes" : "no") << "\n\n";

    os << "--- lanes ---\n"
       << "A lane is averaged according to its CHARACTER, so a view holding both\n"
       << "characters can carry lanes representing different instants under a single\n"
       << "plotfile time. 'represents' below is authoritative; 'plotfile time' is not.\n\n"
       << "  " << std::left << std::setw(10) << "name"
       << std::setw(12) << "character"
       << std::setw(11) << "averaged" << "represents\n";

    for (int local = 0; local < static_cast<int>(names.size()); ++local) {
        const bool present = (local < static_cast<int>(fields.size()) && fields[local] != nullptr);
        const int global_idx =
            (local < static_cast<int>(m_atmos_to_ocean_active_indices.size()))
                ? m_atmos_to_ocean_active_indices[local] : -1;
        const bool averaged = (global_idx >= 0) && TimeAverageAtmosToOceanLane(global_idx);
        const bool flux_character = (global_idx >= 0) &&
            (global_idx == coupling::Rain   || global_idx == coupling::SWrad  ||
             global_idx == coupling::LWrad  || global_idx == coupling::TauX   ||
             global_idx == coupling::TauY   || global_idx == coupling::SHflux ||
             global_idx == coupling::LHflux || global_idx == coupling::Evap);

        std::string represents;
        if (!present) {
            represents = "withheld -- not supplied by ERF";
        } else if (averaged) {
            represents = "interval mean over [" + std::to_string(ocean_time) + ", " +
                         std::to_string(t_end) + "]";
        } else {
            represents = "instantaneous at t=" + std::to_string(t_end);
        }

        os << "  " << std::left << std::setw(10) << names[local]
           << std::setw(12) << (flux_character ? "flux" : "state")
           << std::setw(11) << (present ? (averaged ? "yes" : "no") : "n/a")
           << represents << "\n";
    }
}

void
ERFRemoraMultiBlockContainer::AdvanceBlocks (int ocean_steps)
{
    // V1 Hard Invariants: Verify required solver instances are active
    if (!m_erf || !m_remora) {
        amrex::Abort("Coupling loop started without initialized ERF or REMORA solvers.");
    }

    // Verify subcycling ratio is strictly positive
    if (m_atm_steps_per_exchange <= 0 || m_ocean_steps_per_exchange <= 0 ||
        ocean_steps < 0 || ocean_steps % m_ocean_steps_per_exchange != 0 ||
        m_atm_tick % m_atm_steps_per_exchange != 0) {
        amrex::Abort("Invalid exchange schedule or restart boundary: ocean_steps=" +
                     std::to_string(ocean_steps) + " atm_tick=" + std::to_string(m_atm_tick) +
                     " atm_steps_per_exchange=" + std::to_string(m_atm_steps_per_exchange) +
                     " ocean_steps_per_exchange=" + std::to_string(m_ocean_steps_per_exchange));
    }

    AdvanceBlocksWithCallbacks(
        [this](amrex::Vector<amrex::MultiFab*>& state, amrex::Real time,
               long long /*ocean_step*/) {
            BL_PROFILE("Driver::Pack_Surface");
            // ERF land mask on the destination (ERF) layout: zeroes out
            // land cells in the remapped ocean->atmos SST. Wet/dry masking
            // only, no vector rotation.
            const amrex::iMultiFab* erf_is_land = nullptr;
            ERFGetLandMask(m_erf, &erf_is_land);
            REMORAPackSurfaceState(m_remora, state, time,
                                  m_weight_o2a_mf.get(), m_index_o2a_mf.get(), MAX_STENCIL_SIZE,
                                  erf_is_land);

            MaybeInjectDebugNaN(state, coupling::OceanToAtmosNames(), false);
            AbortOnNonFinite("PackSurfaceState", state, coupling::OceanToAtmosNames());
        },
        [this](const amrex::Vector<amrex::MultiFab*>& state, amrex::Real time,
               long long /*ocean_step*/) {
            BL_PROFILE("Driver::Apply_Surface");
            // Hand the coverage flag over with the SST. Without it ERF cannot
            // tell an uncovered cell's zero fill from a genuine ocean value,
            // and would overwrite its own wrflowinp SST with 0 K wherever the
            // ocean grid does not reach the atmosphere footprint.
            ERFApplyOceanSurfaceState(m_erf, state, time, m_valid_o2a_mf.get());
        },
        [this](amrex::Vector<amrex::MultiFab*>& states, amrex::Real time,
               long long ocean_step) {
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

            MaybeInjectDebugNaN(states, m_atmos_to_ocean_active_names, false);
            AbortOnNonFinite("PackAtmosphericStates", states, m_atmos_to_ocean_active_names);
            MaybeInjectDebugNaN(states, m_atmos_to_ocean_active_names, true);

            for (int idx = 0; idx < static_cast<int>(states.size()); ++idx) {
                auto* mf = states[idx];
                if (mf) {
                    const bool is_cloud_lane = (idx < static_cast<int>(m_atmos_to_ocean_active_indices.size()) &&
                                                m_atmos_to_ocean_active_indices[idx] == coupling::Cloud);
                    if (is_cloud_lane) {
                        // MaskedMinMax is collective (MPI_Allreduce), so it is issued
                        // unconditionally here, exactly as the finite check above is.
                        // cloud_warn_count is a process-local static; gating the
                        // reductions on it is what let the IO rank stop entering
                        // collectives that every other rank still executed. The counter
                        // now limits printing only.
                        const amrex::MultiFab* cloud_mask = mskr;
                        if (cloud_mask != nullptr &&
                            cloud_mask->boxArray().ixType() != mf->boxArray().ixType()) {
                            cloud_mask = nullptr;
                        }
                        const auto cloud_extrema = MaskedMinMax(*mf, cloud_mask, false);
                        const char* marker =
                            (cloud_mask != nullptr) ? " (wet only)" : " (unmasked)";
                        const bool no_wet_cells = (cloud_extrema.included_cells == 0);
                        const bool outside_tolerance =
                            !no_wet_cells &&
                            (cloud_extrema.min_value < -cloud_tol ||
                             cloud_extrema.max_value > amrex::Real(1.0) + cloud_tol);
                        if ((no_wet_cells || outside_tolerance) &&
                            cloud_warn_count < cloud_warn_limit) {
                            ++cloud_warn_count;
                            if (amrex::ParallelDescriptor::IOProcessor()) {
                                if (no_wet_cells) {
                                    amrex::Print() << "WARNING: packed cloud lane has no wet cells; "
                                                     "range check skipped" << marker << "\n";
                                } else {
                                    amrex::Print() << "WARNING: packed cloud lane outside tolerance "
                                                     "[-1e-6,1+1e-6]: min="
                                                   << cloud_extrema.min_value
                                                   << " max=" << cloud_extrema.max_value
                                                   << marker << "\n";
                                }
                            }
                        }
                    }
                }
            }
            ValidateFields("Post-pack atmos->ocean validation",
                           states, m_atmos_to_ocean_active_names, time, ocean_step, false);
        },
        [this](const amrex::Vector<amrex::MultiFab*>& states, amrex::Real time,
               long long ocean_step) {
            BL_PROFILE("Driver::Apply_Atmos");
            ValidateFields("Pre-apply atmos->ocean validation",
                           states, m_atmos_to_ocean_active_names, time, ocean_step, false);
            if (ActiveAtmosToOceanContract() == coupling::AtmosToOceanContract::FluxPassing) {
                REMORAApplyAtmosphericFluxes(m_remora, states, time);
            } else {
                REMORAApplyAtmosphericStates(m_remora, states, time);
            }
            ValidateFields("Post-apply atmos->ocean validation",
                           states, m_atmos_to_ocean_active_names, time, ocean_step, false);
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
ERFRemoraMultiBlockContainer::ActiveAtmosToOceanContract () const
{
    return m_atmos_to_ocean_contract;
}

namespace {
// Fixed-column padding for the provenance table, without pulling in <iomanip>.
std::string
PadTo (const std::string& text, std::size_t width)
{
    std::string padded = text;
    if (padded.size() < width) { padded.append(width - padded.size(), ' '); }
    return padded;
}
} // namespace

void
ERFRemoraMultiBlockContainer::ResolveAtmosToOceanLaneSources ()
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_erf != nullptr && m_remora != nullptr,
        "ResolveAtmosToOceanLaneSources requires both solvers to exist");

    const int n_lanes = coupling::NumAtmosToOceanLayers;
    m_atmos_to_ocean_lane_source.assign(n_lanes, coupling::LaneSource::RemoraDefault);
    m_atmos_to_ocean_lane_reason.assign(n_lanes, std::string("not in either active view"));

    // Both signals arrive in active-view slot order, the same order the pack and
    // apply vectors use, so map back to lane identity through *_active_indices.
    amrex::Vector<int> erf_can_fill_state;
    amrex::Vector<int> erf_can_fill_flux;
    ERFGetFillableAtmosToOceanLanes(m_erf, /*state_contract=*/true,  erf_can_fill_state);
    ERFGetFillableAtmosToOceanLanes(m_erf, /*state_contract=*/false, erf_can_fill_flux);

    amrex::Vector<int> deck_configured;
    REMORAGetDeckConfiguredAtmosLanes(m_remora, deck_configured);

    // driver.atm2ocn_lane_source.<lane> = erf | remora | auto
    const std::string lane_source_prefix = "driver.atm2ocn_lane_source";
    amrex::ParmParse pp_lane(lane_source_prefix);

    // Reject a misspelled lane rather than ignoring it: these are read with
    // query(), so "...lane_source.Uwnd = erf" would otherwise do nothing and
    // leave the lane wherever the default rule put it.
    for (const auto& entry : amrex::ParmParse::getEntries(lane_source_prefix)) {
        const std::string lane_name = entry.substr(lane_source_prefix.size() + 1);
        bool known = false;
        for (const auto& candidate : coupling::AtmosToOceanLaneSpecs()) {
            if (lane_name == candidate.name) { known = true; break; }
        }
        if (!known) {
            std::string valid;
            for (const auto& candidate : coupling::AtmosToOceanLaneSpecs()) {
                if (candidate.index >= coupling::NumAtmosToOceanStateLayers) { break; }
                if (!valid.empty()) { valid += ", "; }
                valid += candidate.name;
            }
            amrex::Abort(entry + " names no atmos->ocean lane. Valid lanes: " + valid + ".");
        }
    }

    const auto& state_indices = m_atmos_to_ocean_state_active_indices;
    for (int local = 0; local < static_cast<int>(state_indices.size()); ++local) {
        const int global_idx = state_indices[local];
        const auto& spec = coupling::AtmosToOceanLane(global_idx);

        const bool erf_can_fill =
            (local < static_cast<int>(erf_can_fill_state.size())) && erf_can_fill_state[local] != 0;
        const bool deck_owns =
            (local < static_cast<int>(deck_configured.size())) && deck_configured[local] != 0;

        std::string requested = "auto";
        const bool has_override = pp_lane.query(spec.name, requested);
        requested = amrex::toLower(requested);
        if (has_override && requested != "erf" && requested != "remora" && requested != "auto") {
            amrex::Abort("driver.atm2ocn_lane_source." + std::string(spec.name) + " = " + requested +
                         " is not a source. Valid values: erf, remora, auto.");
        }

        coupling::LaneSource source = coupling::LaneSource::RemoraDefault;
        std::string reason;

        if (has_override && requested == "erf") {
            // Unsatisfiable, not a preference: the lane would reach REMORA
            // holding the driver's prefill. Say so rather than quietly obeying.
            if (!erf_can_fill) {
                amrex::Abort("driver.atm2ocn_lane_source." + std::string(spec.name) +
                             " = erf, but ERF's pack does not fill the " + spec.name +
                             " lane in this configuration, so the lane would carry only the "
                             "driver's prefill. Configure ERF to produce it, or drop the override.");
            }
            source = coupling::LaneSource::Erf;
            reason = "forced by driver.atm2ocn_lane_source." + std::string(spec.name) + "=erf";
        } else if (has_override && requested == "remora") {
            source = deck_owns ? coupling::LaneSource::RemoraDeck
                               : coupling::LaneSource::RemoraDefault;
            reason = "withheld by driver.atm2ocn_lane_source." + std::string(spec.name) + "=remora";
        } else if (deck_owns) {
            // The one case where ERF could supply a lane and deliberately does not.
            source = coupling::LaneSource::RemoraDeck;
            reason = spec.remora_type_key != nullptr
                         ? "remora." + std::string(spec.remora_type_key) + " or its value set in deck"
                         : std::string("set in deck");
        } else if (erf_can_fill) {
            source = coupling::LaneSource::Erf;
            reason = "packed by ERF";
        } else {
            source = coupling::LaneSource::RemoraDefault;
            reason = "ERF does not fill it; REMORA default stands";
        }

        m_atmos_to_ocean_lane_source[global_idx] = source;
        m_atmos_to_ocean_lane_reason[global_idx] = reason;
    }

    // Flux lanes are reported but never withheld: ApplyAtmosphericFluxes skips
    // the whole application if any of the eight is null. Recorded so the report
    // can still say which ones ERF computes.
    const auto& flux_indices = m_atmos_to_ocean_flux_active_indices;
    for (int local = 0; local < static_cast<int>(flux_indices.size()); ++local) {
        const int global_idx = flux_indices[local];
        if (global_idx < coupling::TauX) {
            continue; // SWrad/LWrad/Rain sit in both views; the state pass owns them.
        }
        const bool erf_can_fill =
            (local < static_cast<int>(erf_can_fill_flux.size())) && erf_can_fill_flux[local] != 0;
        m_atmos_to_ocean_lane_source[global_idx] =
            erf_can_fill ? coupling::LaneSource::Erf : coupling::LaneSource::RemoraDefault;
        m_atmos_to_ocean_lane_reason[global_idx] =
            erf_can_fill ? "packed by ERF"
                         : "ERF does not compute it; passed as zero (flux mode passes every lane)";
    }

    // The one lane ERF fills that REMORA would otherwise ignore: lwrad_type
    // defaults to "computed", under which bulk_fluxes gets a null longwave
    // pointer and uses Berliand, so copying ERF's longwave in accomplishes
    // nothing. A deck that set lwrad_type resolved to RemoraDeck above and never
    // reaches here, so this cannot override an explicit choice.
    if (m_atmos_to_ocean_lane_source[coupling::LWrad] == coupling::LaneSource::Erf) {
        SetREMORALongwaveFromDriver(m_remora);
        m_atmos_to_ocean_lane_reason[coupling::LWrad] +=
            "; REMORA longwave set to external downward";
    }

    m_atmos_to_ocean_lane_sources_resolved = true;
}

void
ERFRemoraMultiBlockContainer::PrintAtmosToOceanLaneProvenance () const
{
    amrex::Print() << "Atmos->ocean lane provenance (configured mode="
                   << coupling::AtmosToOceanContractName(m_atmos_to_ocean_contract) << ")\n"
                   << "  " << PadTo("lane", 8) << PadTo("source", 8) << "reason\n";

    auto print_lane = [&](int global_idx) {
        const auto& spec = coupling::AtmosToOceanLane(global_idx);
        amrex::Print() << "  " << PadTo(spec.name, 8)
                       << PadTo(coupling::LaneSourceName(m_atmos_to_ocean_lane_source[global_idx]), 8)
                       << m_atmos_to_ocean_lane_reason[global_idx] << "\n";
    };

    for (int global_idx : m_atmos_to_ocean_state_active_indices) {
        print_lane(global_idx);
    }
    if (m_atmos_to_ocean_contract == coupling::AtmosToOceanContract::FluxPassing) {
        amrex::Print() << "  (above is the ocean_step=0 bootstrap; steady-state flux lanes follow)\n";
        for (int global_idx : m_atmos_to_ocean_flux_active_indices) {
            if (global_idx >= coupling::TauX) { print_lane(global_idx); }
        }
    }
    amrex::Print() << "  A lane sourced from REMORA is withheld from the apply view, so REMORA's own\n"
                      "  forcing for it stands. Override per lane with "
                      "driver.atm2ocn_lane_source.<lane>.\n";
}

void
ERFRemoraMultiBlockContainer::SelectAtmosToOceanActiveView ()
{
    const auto active_contract = ActiveAtmosToOceanContract();
    if (active_contract == coupling::AtmosToOceanContract::FluxPassing) {
        m_atmos_to_ocean_active = m_atmos_to_ocean_flux_active;
        m_atmos_to_ocean_active_indices = m_atmos_to_ocean_flux_active_indices;
        m_atmos_to_ocean_active_names = m_atmos_to_ocean_flux_active_names;
    } else {
        m_atmos_to_ocean_active = m_atmos_to_ocean_state_active;
        m_atmos_to_ocean_active_indices = m_atmos_to_ocean_state_active_indices;
        m_atmos_to_ocean_active_names = m_atmos_to_ocean_state_active_names;
    }

    // The view handed to REMORA: the active view with every lane the driver does
    // not own nulled. REMORA reads authority purely from the pointer being
    // non-null, so a lane the atmosphere never wrote would read as "the
    // atmosphere says this" and displace REMORA's own forcing. Withholding is
    // how the driver declines a lane.
    //
    // ResolveAtmosToOceanLaneSources() decided which; this only applies it.
    //
    // State contract only. The flux contract has no per-lane fallback: REMORA's
    // ApplyAtmosphericFluxes requires all eight lanes non-null and skips the whole
    // application if any is missing, so withholding there would disable every flux
    // rather than one. In flux mode a zero radiative flux is at least a truthful
    // statement that the atmosphere computed none.
    m_atmos_to_ocean_apply_view = m_atmos_to_ocean_active;
    if (active_contract == coupling::AtmosToOceanContract::StatePassing &&
        m_atmos_to_ocean_lane_sources_resolved) {
        for (int local_idx = 0;
             local_idx < static_cast<int>(m_atmos_to_ocean_active_indices.size()); ++local_idx) {
            const int global_idx = m_atmos_to_ocean_active_indices[local_idx];
            if (m_atmos_to_ocean_lane_source[global_idx] != coupling::LaneSource::Erf) {
                m_atmos_to_ocean_apply_view[local_idx] = nullptr;
            }
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
        PrintDiagnostics("Initial ocean->atmos cache", m_ocean_to_atmos,
                         coupling::OceanToAtmosNames(), true);
    }
    PrintDiagnostics("Initial atmos->ocean cache", m_atmos_to_ocean_active,
                     m_atmos_to_ocean_active_names, false);


    // Multi-rate clock rationale:
    // - ERF advances a compressible atmosphere where acoustic stability often drives
    //   dt to O(1-10 s).
    // - REMORA advances Boussinesq/split-explicit ocean dynamics with a much larger
    //   macro-step O(10^2-10^3 s).
    // Hold interval-start SST while advancing the atmosphere, then hold its
    // interval-mean flux throughout the ocean substeps.
    const amrex::Real exchange_interval = m_erf_dt * static_cast<amrex::Real>(m_atm_steps_per_exchange);
    const long long starting_atm_tick = m_atm_tick;
    const long long expected_atm_ticks = starting_atm_tick +
                                         static_cast<long long>(ocean_steps / m_ocean_steps_per_exchange) *
                                         static_cast<long long>(m_atm_steps_per_exchange);
    const amrex::Real abs_tol = 1.0e-12;
    const amrex::Real rel_tol = 1.0e-8;

    amrex::Print() << "Coupled cadence plan:"
                       << " ocean_steps=" << ocean_steps
                       << " starting_atm_tick=" << starting_atm_tick
                       << " atm_steps_per_exchange=" << m_atm_steps_per_exchange
                       << " target_atm_tick=" << expected_atm_ticks
                       << " expected_run_time=" << (m_ocean_dt * static_cast<amrex::Real>(ocean_steps))
                       << std::endl;

    auto enforce_dt = [&] (const char* solver_name, amrex::Real dt_request, amrex::Real dt_used,
                           int ocean_step, int substep) {
        const amrex::Real err = amrex::Math::abs(dt_used - dt_request);
        const amrex::Real tol = amrex::max(abs_tol, rel_tol * amrex::Math::abs(dt_request));
        if (!std::isfinite(dt_used) || err > tol) {
            amrex::Print() << "Coupled cadence mismatch before abort:"
                           << " solver=" << solver_name
                           << " ocean_step=" << ocean_step
                           << " substep=" << substep
                           << " dt_request=" << dt_request
                           << " dt_used=" << dt_used
                           << " abs_err=" << err
                           << " tol=" << tol
                           << " m_erf_dt=" << m_erf_dt
                           << " atm_steps_per_exchange=" << m_atm_steps_per_exchange
                           << " expected_ocean_macro_dt=" << exchange_interval
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

    SelectAtmosToOceanActiveView();
    const bool averaging_active = AtmosToOceanAveragingActive();
    for (int ocean_step = 0; ocean_step < ocean_steps; ocean_step += m_ocean_steps_per_exchange) {
        bool stop_coupling = false;
        if (m_atm_tick >= expected_atm_ticks) {
            amrex::Print() << "Reached expected coupled end time before ocean_step "
                                << ocean_step << "; stopping cleanly at atmosphere tick "
                                << m_atm_tick << " of " << expected_atm_ticks << std::endl;
            break;
        }
        const amrex::Real ocean_time = m_ocean_time_offset + m_atm_tick * m_erf_dt;
        const long long global_ocean_step =
            (m_atm_tick / m_atm_steps_per_exchange) * m_ocean_steps_per_exchange;
        const bool use_state_contract =
            (ActiveAtmosToOceanContract() == coupling::AtmosToOceanContract::StatePassing);
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
            pack_surface_state(m_ocean_to_atmos, ocean_time, global_ocean_step);
            FillOceanToAtmosGhosts();
        }

        if (averaging_active) {
            ResetAtmosToOceanAccumulators();
        }

        // At global_ocean_step 0 on a cold start the atmosphere has not advanced
        // yet, so ERF's diagnostic flux arrays (Tau, SFS_*) still hold the
        // setVal(zero) from allocation in ERF_MakeNewArrays.cpp rather than a
        // physical state at t. Sampling them as the trapezoid's left endpoint
        // would fold an exact zero into the interval mean at weight 0.5/n. So
        // skip that sample and reconstruct f(t) from the first two post-advance
        // samples instead -- see the weight correction below. On a restart Tau is
        // checkpointed and m_atm_tick is restored, so global_ocean_step > 0 and
        // the ordinary trapezoid runs untouched.
        const bool cold_start_extrapolate = (global_ocean_step == 0);

        for (int substep = 0; substep < m_atm_steps_per_exchange; ++substep) {
            if (m_atm_tick >= expected_atm_ticks) {
                amrex::Print() << "Reached expected coupled end time before ERF advance at ocean_step="
                                    << ocean_step << " substep=" << substep
                                    << "; skipping remaining substeps and exiting cleanly."
                                    << std::endl;
                stop_coupling = true;
                break;
            }

            if (averaging_active && substep == 0 && !cold_start_extrapolate) {
                pack_atmos_states(m_atmos_to_ocean_active, ocean_time, global_ocean_step);
                AccumulateAtmosToOcean(amrex::Real(0.5));
            }
            const amrex::Real atm_time = m_atm_tick * m_erf_dt;

            BL_PROFILE("Driver::Time_Interpolate");
            if (m_use_two_way_coupling) {
                DriverPrintLine(
                    m_driver_verbose,
                    std::string("atm_tick=") + std::to_string(m_atm_tick) +
                        " apply ocean->atmos surface state to ERF at time=" +
                        std::to_string(atm_time));
                apply_ocean_surface_state(m_ocean_to_atmos, atm_time, global_ocean_step);
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

            if (averaging_active) {
                const amrex::Real substep_end_time =
                    ocean_time + static_cast<amrex::Real>(substep + 1) * m_erf_dt;
                pack_atmos_states(m_atmos_to_ocean_active,
                                  substep_end_time, global_ocean_step);
                amrex::Real sample_weight =
                    (substep == m_atm_steps_per_exchange - 1)
                        ? amrex::Real(0.5) : amrex::Real(1.0);
                // Cold-start reconstruction of the left endpoint skipped above.
                // Linear extrapolation f(t) ~= 2*f1 - f2 is exact for linear
                // forcing, which is the property the trapezoid was chosen for, so
                // folding 0.5*f(t) into the sum adds +1.0 to f1's weight and
                // -0.5 to f2's and leaves the total at n:
                //   n=10 -> 2.0 + 0.5 + 7*1.0 + 0.5 = 10.
                // At n=2 the last-substep 0.5 and this -0.5 cancel, so f2 drops
                // out and the mean collapses to f1 -- the interval midpoint,
                // which is the exact linear answer over [t, t+2*dt]. f2 then
                // carries zero weight and AccumulateAtmosToOcean adds nothing,
                // but it is still packed: non-averaged lanes keep that
                // interval-end sample. At n=1 there
                // is no f2 to extrapolate from, so the correction is skipped
                // entirely; DrainAtmosToOceanAccumulators divides by the total
                // accumulated weight, so that interval degrades cleanly to the
                // single sample f1 rather than being mis-normalized.
                if (cold_start_extrapolate && m_atm_steps_per_exchange >= 2) {
                    if (substep == 0) {
                        sample_weight += amrex::Real(1.0);
                    } else if (substep == 1) {
                        sample_weight -= amrex::Real(0.5);
                    }
                }
                AccumulateAtmosToOcean(sample_weight);
            }
        }

        if (stop_coupling) {
            break;
        }

        // An interval mean is stamped at its midpoint. With averaging disabled,
        // preserve commit 1's interval-end sample and label exactly.
        const amrex::Real a2o_sample_time =
            averaging_active ? ocean_time + amrex::Real(0.5) * exchange_interval
                             : ocean_time + exchange_interval;
        if (averaging_active) {
            DriverPrintLine(
                m_driver_verbose,
                std::string("ocean_step=") + std::to_string(global_ocean_step) +
                    " drain atmos->ocean time average at time=" +
                    std::to_string(a2o_sample_time));
            DrainAtmosToOceanAccumulators();
        } else {
            DriverPrintLine(
                m_driver_verbose,
                std::string("ocean_step=") + std::to_string(global_ocean_step) +
                    " pack atmos->ocean " + (use_state_contract ? "state" : "flux") +
                    " from ERF at time=" + std::to_string(a2o_sample_time));
            pack_atmos_states(m_atmos_to_ocean_active, a2o_sample_time, global_ocean_step);
        }
        if (m_plot_int_exchange > 0 &&
            (global_ocean_step % static_cast<long long>(m_plot_int_exchange)) == 0) {
            // O2A was packed at ocean_time, A2O at a2o_sample_time; stamp each
            // with the instant it actually represents.
            WriteExchangePlotfiles(ocean_time, a2o_sample_time, exchange_interval, global_ocean_step);
        }
        FillAtmosToOceanGhosts();
        DriverPrintLine(
            m_driver_verbose,
            std::string("ocean_step=") + std::to_string(global_ocean_step) +
                " apply atmos->ocean " + (use_state_contract ? "state" : "flux") +
                " to REMORA at time=" + std::to_string(a2o_sample_time));
        // Masked view: lanes ERF cannot fill are nulled so REMORA keeps its own
        // forcing for them (see SelectAtmosToOceanActiveView).
        for (int ocn_substep = 0; ocn_substep < m_ocean_steps_per_exchange; ++ocn_substep) {
            const amrex::Real step_time = ocean_time + ocn_substep * m_ocean_dt;
            // Reapply the same interval-mean cache before every ocean advance;
            // REMORA may refresh forcing arrays inside an individual timestep.
            apply_atmos_states(m_atmos_to_ocean_apply_view, a2o_sample_time,
                               global_ocean_step + ocn_substep);
            const amrex::Real used = REMORAEvolveOneStep(m_remora, step_time, m_ocean_dt);
            if (!std::isfinite(used) || used <= amrex::Real(0.0)) {
                amrex::Abort("REMORA terminated before completing the exchange interval.");
            }
            enforce_dt("REMORA", m_ocean_dt, used, ocean_step + ocn_substep, ocn_substep);
        }
        amrex::Print() << "Completed exchange " << (m_atm_tick / m_atm_steps_per_exchange)
                       << " elapsed_time=" << m_atm_tick * m_erf_dt
                       << " atmosphere_tick=" << m_atm_tick
                       << " completed_ocean_steps=" << global_ocean_step + m_ocean_steps_per_exchange
                       << "\n";

        if (m_use_two_way_coupling) {
            PrintDiagnostics("Ocean->atmos cache", m_ocean_to_atmos,
                             coupling::OceanToAtmosNames(), true);
        }
        PrintDiagnostics("Atmos->ocean cache", m_atmos_to_ocean_active,
                         m_atmos_to_ocean_active_names, false);
    }
    REMORAWriteAtFinalTime(m_remora);
    ERFWriteAtFinalTime(m_erf);

}
} // namespace erf_remora
