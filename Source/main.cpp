#include "ERFRemoraMultiBlockContainer.H"
#include "CoupledRestartUtils.H"
#include "StandaloneSolverRuns.H"

#include <AMReX.H>
#include <AMReX_Math.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_RealBox.H>

#ifdef ERF_USE_KOKKOS
#include <Kokkos_Core.hpp>
#include "ERF_Kokkos.H"
#endif

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <string>

std::string inputs_name;

namespace
{
struct RestartState
{
    bool active = false;
    long long atm_tick = 0;
    int completed_ocean_steps = 0;
    amrex::Real restored_time = amrex::Real(0.0);
};

bool
NearlyEqual (amrex::Real lhs, amrex::Real rhs,
             amrex::Real abs_tol = amrex::Real(1.0e-12),
             amrex::Real rel_tol = amrex::Real(1.0e-8))
{
    const amrex::Real err = amrex::Math::abs(lhs - rhs);
    const amrex::Real scale = amrex::max(amrex::Math::abs(lhs), amrex::Math::abs(rhs));
    return err <= amrex::max(abs_tol, rel_tol * scale);
}

void
WarnIfDtChanged (const char* solver_name,
                 amrex::Real checkpoint_dt,
                 amrex::Real configured_dt)
{
    if (!NearlyEqual(checkpoint_dt, configured_dt)) {
        amrex::Print() << "WARNING: " << solver_name
                        << " checkpoint dt=" << checkpoint_dt
                        << " differs from configured dt=" << configured_dt
                        << "; continuing because checkpoint time aligns with the configured coupled cadence.\n";
    }
}

long long
RequireIntegralTick (const char* label, amrex::Real time, amrex::Real dt)
{
    if (dt <= amrex::Real(0.0)) {
        amrex::Abort(std::string("Non-positive dt while restoring coupled restart state for ") + label);
    }

    const amrex::Real raw_tick = time / dt;
    const long long tick = static_cast<long long>(std::llround(raw_tick));
    if (!NearlyEqual(static_cast<amrex::Real>(tick) * dt, time)) {
        amrex::Abort(std::string("Restored ") + label +
                     " time does not align with the configured coupled cadence:"
                     " time=" + std::to_string(time) +
                     " dt=" + std::to_string(dt));
    }

    return tick;
}

void
RequireCheckpointSelfConsistency (const erf_remora::restart::CheckpointHeaderMetadata& meta,
                                  amrex::Real time_offset = amrex::Real(0.0))
{
    const amrex::Real reconstructed_time =
        time_offset + static_cast<amrex::Real>(meta.step) * meta.dt;
    if (!NearlyEqual(reconstructed_time, meta.time)) {
        amrex::Abort(std::string("Restart header for ") + erf_remora::restart::SolverKindName(meta.solver) +
                     " is internally inconsistent:"
                     " istep[0]=" + std::to_string(meta.step) +
                     " dt[0]=" + std::to_string(meta.dt) +
                     " t_new[0]=" + std::to_string(meta.time));
    }
}

RestartState
BuildRestartState (const erf_remora::restart::CheckpointHeaderMetadata& erf_meta,
                   const erf_remora::restart::CheckpointHeaderMetadata& rem_meta,
                   amrex::Real erf_dt,
                   amrex::Real ocean_macro_dt,
                   amrex::Real ocean_time_offset)
{
    RequireCheckpointSelfConsistency(erf_meta);
    RequireCheckpointSelfConsistency(rem_meta, ocean_time_offset);

    // REMORA checkpoint time includes driver.ocean_time_offset; ERF uses elapsed time.
    const amrex::Real rem_elapsed_time = rem_meta.time - ocean_time_offset;

    if (!NearlyEqual(erf_meta.time, rem_elapsed_time)) {
        amrex::Abort("ERF and REMORA restart checkpoints do not represent the same coupled time.");
    }

    const long long restored_atm_tick =
        RequireIntegralTick("ERF", erf_meta.time, erf_dt);
    const long long restored_ocean_tick =
        RequireIntegralTick("REMORA", rem_elapsed_time, ocean_macro_dt);

    WarnIfDtChanged("ERF", erf_meta.dt, erf_dt);
    WarnIfDtChanged("REMORA", rem_meta.dt, ocean_macro_dt);

    RestartState state;
    state.active = true;
    state.atm_tick = restored_atm_tick;
    state.completed_ocean_steps = static_cast<int>(restored_ocean_tick);
    state.restored_time = erf_meta.time;
    return state;
}
}

int
main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    if (argc < 2) {
        amrex::Abort("No input file provided."
                     " Usage: erf-remora-driver <inputs-file> [<key>=<value> ...]");
    }

    if (!strchr(argv[1], '=')) {
      inputs_name = argv[1];
    }

#ifdef ERF_USE_KOKKOS
    // Initialize kokkos
    if (!Kokkos::is_initialized()) {
        Kokkos::initialize(Kokkos::InitializationSettings()
                           .set_device_id(amrex::Gpu::Device::deviceId()));
    }
#endif

    {
        amrex::ParmParse pp_amr("amr");
        if (!pp_amr.contains("refine_grid_layout_x")) {
            pp_amr.add("refine_grid_layout_x", 1);
        }
        if (!pp_amr.contains("refine_grid_layout_y")) {
            pp_amr.add("refine_grid_layout_y", 1);
        }
        if (!pp_amr.contains("refine_grid_layout_z")) {
            pp_amr.add("refine_grid_layout_z", 0);
        }

        // Match the standalone ERF/REMORA entry points: keep the default AMR
        // max-grid size large so any additional splitting for rank count happens
        // in x/y only, unless the user explicitly overrides AMR grid sizes.
        int default_max_grid_size = 2048;
        pp_amr.queryAdd("max_grid_size", default_max_grid_size);

        int default_blocking_factor = 1;
        pp_amr.queryAdd("blocking_factor", default_blocking_factor);

        int default_n_proper = 2;
        pp_amr.queryAdd("n_proper", default_n_proper);

        // The driver owns the combined ParmParse environment, so drop
        // the prefixed alias; the canonical max_step/stop_time keys remain.
        // Capture remora.stop_time before removing it: under a nonzero
        // ocean_time_offset the driver later re-keys an ocean-frame value, and a
        // user who set this explicitly deserves to be told it was discarded.
        // Querying after the scrub would never see it.
        amrex::Real user_remora_stop_time = 0.0;
        const bool user_set_remora_stop_time =
            amrex::ParmParse().query("remora.stop_time", user_remora_stop_time);
        amrex::ParmParse().remove("remora.max_step");
        amrex::ParmParse().remove("remora.stop_time");

        amrex::ParmParse pp_driver("driver");

        int max_grid_size = default_max_grid_size;
        int ocean_steps = 2;
        std::string atm2ocn_mode = "state";
        std::string coupling_direction;

        // The coupled interface layout is a horizontal z=0 slab, so the
        // divisibility constraint that motivated the AMR fallback is really an
        // x/y concern. We still inherit all three AMR max-grid knobs here so
        // driver.max_grid_size[_x/_y/_z] behaves like a normal override of the
        // canonical AMR defaults.
        pp_amr.query("max_grid_size", max_grid_size);
        int max_grid_size_x = max_grid_size;
        int max_grid_size_y = max_grid_size;
        int max_grid_size_z = max_grid_size;
        pp_amr.query("max_grid_size_x", max_grid_size_x);
        pp_amr.query("max_grid_size_y", max_grid_size_y);
        pp_amr.query("max_grid_size_z", max_grid_size_z);
        const bool has_driver_max_grid_size = pp_driver.query("max_grid_size", max_grid_size);
        if (has_driver_max_grid_size) {
            max_grid_size_x = max_grid_size;
            max_grid_size_y = max_grid_size;
            max_grid_size_z = max_grid_size;
        }
        pp_driver.query("max_grid_size_x", max_grid_size_x);
        pp_driver.query("max_grid_size_y", max_grid_size_y);
        pp_driver.query("max_grid_size_z", max_grid_size_z);
        const amrex::IntVect driver_max_grid_size(max_grid_size_x, max_grid_size_y, max_grid_size_z);

        pp_driver.query("atm2ocn_mode", atm2ocn_mode);
        pp_driver.get("coupling_direction", coupling_direction);

        amrex::ParmParse pp_erf("erf");
        amrex::ParmParse pp_rem("remora");
        int driver_verbose = 0;
        pp_driver.query("v", driver_verbose);

        amrex::Real ocean_time_offset = 0.0;
        pp_driver.query("ocean_time_offset", ocean_time_offset);

        amrex::Real erf_dt;
        amrex::Real ocean_macro_dt;
        pp_erf.get("fixed_dt", erf_dt);
        pp_rem.get("fixed_dt", ocean_macro_dt);

        // Both dts divide or scale the cadence arithmetic below, so reject
        // non-positive values before dividing rather than letting a non-finite
        // quotient reach std::llround.
        if (erf_dt <= amrex::Real(0.0)) {
            amrex::Abort("erf.fixed_dt must be positive; got " + std::to_string(erf_dt));
        }
        if (ocean_macro_dt <= amrex::Real(0.0)) {
            amrex::Abort("remora.fixed_dt must be positive; got " + std::to_string(ocean_macro_dt));
        }

        // Round the ratio and validate within tolerance rather than truncating and
        // testing exact equality: decimal dt pairs with mathematically integral
        // ratios (e.g. 0.1 and 0.3) are not exactly representable in binary, and
        // both the truncation and the multiply-back can land off by 1 ulp.
        const long long substep_ratio = std::llround(ocean_macro_dt / erf_dt);
        int atm_substeps_per_ocean_step = static_cast<int>(substep_ratio);
        if (substep_ratio < 1 ||
            !NearlyEqual(static_cast<amrex::Real>(substep_ratio) * erf_dt, ocean_macro_dt)) {
            amrex::Print() << "ERF    dt " << erf_dt << std::endl;
            amrex::Print() << "REMORA dt " << ocean_macro_dt << std::endl;
            amrex::Abort("Ocean dt is not a multiple of atmosphere dt");
        }

        amrex::Real stop_time = -1.0;
        const bool has_stop_time = amrex::ParmParse().query("stop_time", stop_time);
        // Round when stop_time sits (within tolerance) on a macro-step boundary;
        // otherwise keep truncating to the previous boundary, since stop_time is
        // legitimately allowed to fall between boundaries.
        const amrex::Real stop_time_ratio = stop_time / ocean_macro_dt;
        const long long stop_time_ratio_rounded = std::llround(stop_time_ratio);
        const int stop_time_macro_steps =
            NearlyEqual(static_cast<amrex::Real>(stop_time_ratio_rounded) * ocean_macro_dt, stop_time)
                ? static_cast<int>(stop_time_ratio_rounded)
                : static_cast<int>(stop_time_ratio);
        // REMORA advances on the absolute ocean clock once driver.ocean_time_offset
        // is applied, but it reads the unprefixed stop_time, which the deck and the
        // macro-step plan above express as elapsed-from-zero. Compute the ocean-frame
        // value here, where the macro-step plan lives, and hand it to the container:
        // the re-key must happen between the ERF and REMORA constructions so ERF
        // still reads the elapsed value and keeps its own dt clamp.
        amrex::Real ocean_frame_stop_time = -1.0;
        if (has_stop_time && ocean_time_offset != amrex::Real(0.0)) {
            ocean_frame_stop_time =
                ocean_time_offset + ocean_macro_dt * static_cast<amrex::Real>(stop_time_macro_steps);
            if (user_set_remora_stop_time) {
                amrex::Print() << "WARNING: remora.stop_time=" << user_remora_stop_time
                               << " is discarded; the driver derives REMORA's ocean-frame stop time "
                               << ocean_frame_stop_time
                               << " from stop_time and driver.ocean_time_offset.\n";
            }
        }
        int top_level_max_step = -1;
        const bool has_top_level_max_step = amrex::ParmParse().query("max_step", top_level_max_step);

        // ERF consumes amr.restart inside InitializeSolvers(), before the scrub below.
        std::string amr_restart_path;
        if (pp_amr.query("restart", amr_restart_path)) {
            amrex::Abort("amr.restart is not supported by the coupled driver; "
                         "set erf.restart and remora.restart together instead.");
        }

        std::string erf_restart_path;
        std::string remora_restart_path;
        const bool has_erf_restart = pp_erf.query("restart", erf_restart_path);
        const bool has_remora_restart = pp_rem.query("restart", remora_restart_path);
        if (has_erf_restart != has_remora_restart) {
            amrex::Abort("Coupled restart requires both erf.restart and remora.restart to be set together.");
        }

        RestartState restart_state;
        if (has_erf_restart) {
            const auto erf_meta =
                erf_remora::restart::ReadCheckpointHeaderMetadata(
                    erf_restart_path, erf_remora::restart::SolverKind::ERF);
            const auto rem_meta =
                erf_remora::restart::ReadCheckpointHeaderMetadata(
                    remora_restart_path, erf_remora::restart::SolverKind::REMORA);
            restart_state = BuildRestartState(erf_meta, rem_meta, erf_dt,
                                              ocean_macro_dt, ocean_time_offset);
        }

        if (pp_driver.contains("coupled_one_step")) {
            amrex::Print() << "WARNING: driver.coupled_one_step is deprecated and ignored; "
                           << "use deck naming plus driver.coupling_direction instead.\n";
        }
        if (pp_driver.contains("synthetic_only")) {
            amrex::Print() << "WARNING: driver.synthetic_only is retired and ignored; "
                           << "synthetic-only mode has been removed.\n";
        }

        auto lowercase = [](std::string value) {
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            return value;
        };

        if (has_stop_time) {
            amrex::Print() << "stop_time = " << stop_time << "s, corresponding to " << stop_time_macro_steps << " steps" << std::endl;
        }
        if (has_top_level_max_step) {
            amrex::Print() << "max_step = " << top_level_max_step << std::endl;
        }
        if (has_stop_time && has_top_level_max_step) {
            ocean_steps = std::min(stop_time_macro_steps, top_level_max_step);
            amrex::Print() << "Clipping to minimum of stop_time steps and max_step = " << ocean_steps << std::endl;
        } else if (has_stop_time && !has_top_level_max_step) {
            ocean_steps = stop_time_macro_steps;
        } else if (!has_stop_time && has_top_level_max_step) {
            ocean_steps = top_level_max_step;
        } else {
            amrex::Abort("Either stop_time or max_step must be specified");
        }

        if (has_stop_time) {
            const amrex::Real clipped_stop_time =
                ocean_macro_dt * static_cast<amrex::Real>(stop_time_macro_steps);
            if (stop_time_macro_steps >= 0 &&
                stop_time > clipped_stop_time) {
                amrex::Print() << "WARNING: coupled driver only stops on ocean macro-step boundaries;"
                               << " stop_time=" << stop_time
                               << " is being clipped to " << clipped_stop_time
                               << " for driver-side termination.\n";
            }
        }

        if (restart_state.active) {
            amrex::Print() << "Restored coupled restart state:"
                                << " time=" << restart_state.restored_time
                                << " atm_tick=" << restart_state.atm_tick
                                << " completed_ocean_steps=" << restart_state.completed_ocean_steps
                                << std::endl;
            const int requested_ocean_steps = ocean_steps;
            ocean_steps = requested_ocean_steps - restart_state.completed_ocean_steps;
            if (ocean_steps < 0) {
                amrex::Abort("Restarted coupled run has already advanced beyond the requested max_step/stop_time target.");
            }
        } else {
            // Cold-start path keeps the original max_step/stop_time logic above.
        }

        const amrex::Real expected_run_time =
            ocean_macro_dt * static_cast<amrex::Real>(ocean_steps);
        amrex::Print() << "Coupled run plan:"
                       << " expected_run_time=" << expected_run_time
                       << " erf_dt=" << erf_dt
                       << " atm_substeps_per_ocean_step=" << atm_substeps_per_ocean_step
                       << " ocean_steps=" << ocean_steps
                       << " driver_max_grid_size=" << driver_max_grid_size;
        if (restart_state.active) {
            amrex::Print() << " restored_time=" << restart_state.restored_time
                           << " completed_ocean_steps=" << restart_state.completed_ocean_steps
                           << " final_target_ocean_steps="
                           << (ocean_steps + restart_state.completed_ocean_steps);
        }
        if (has_stop_time) {
            amrex::Print() << " stop_time=" << stop_time
                           << " (either stop_time or max_step may terminate the run first)";
        } else {
            amrex::Print() << " stop_time=<unset>";
        }
        amrex::Print() << "\n";

        erf_remora::coupling::AtmosToOceanContract atm2ocn_contract = erf_remora::coupling::AtmosToOceanContract::StatePassing;
        if (atm2ocn_mode == "state") {
            atm2ocn_contract = erf_remora::coupling::AtmosToOceanContract::StatePassing;
        } else if (atm2ocn_mode == "flux") {
            atm2ocn_contract = erf_remora::coupling::AtmosToOceanContract::FluxPassing;
        } else {
            amrex::Abort("Unknown driver.atm2ocn_mode. Valid values: state, flux.");
        }

        bool use_two_way_coupling = false;
        const std::string coupling_direction_lc = lowercase(coupling_direction);

        if (coupling_direction_lc == "two_way") {
            use_two_way_coupling = true;
            // Tell ERF to expect coupled SST. This used to inject
            // land_surface_model=OceanSurf, which consumed ERF's one land
            // surface model slot and so made a two-way run mutually exclusive
            // with Noah-MP. Coupled SST now arrives through the lower-boundary
            // path instead, leaving erf.land_surface_model free for whatever
            // land surface the deck actually wants.
            pp_erf.add("use_coupled_sst", 1);
        } else if (coupling_direction_lc == "erftoremora") {
            use_two_way_coupling = false;
        } else {
            amrex::Abort("Unknown driver.coupling_direction. Valid values: ERFtoREMORA, two_way.");
        }

        amrex::Vector<amrex::Real> atm_lo;
        amrex::Vector<amrex::Real> atm_hi;
        amrex::Vector<int> n_cell_atm;
        amrex::Array<int,AMREX_SPACEDIM> atm_is_periodic{1, 1, 0};
        amrex::Vector<amrex::Real> ocn_lo;
        amrex::Vector<amrex::Real> ocn_hi;
        amrex::Vector<int> n_cell_ocn;
        amrex::Array<int,AMREX_SPACEDIM> ocn_is_periodic{1, 1, 0};
        if (!pp_erf.queryarr("prob_lo", atm_lo, 0, AMREX_SPACEDIM)) {
            amrex::Abort("Missing required input: erf.prob_lo");
        }
        const bool has_erf_prob_hi = pp_erf.queryarr("prob_hi", atm_hi, 0, AMREX_SPACEDIM);
        amrex::Vector<amrex::Real> atm_extent;
        const bool has_erf_prob_extent = pp_erf.queryarr("prob_extent", atm_extent, 0, AMREX_SPACEDIM);
        if (has_erf_prob_hi == has_erf_prob_extent) {
            amrex::Abort("ERF topology requires exactly one of erf.prob_hi or erf.prob_extent.");
        }
        if (has_erf_prob_extent) {
            atm_hi.resize(AMREX_SPACEDIM);
            for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                atm_hi[i] = atm_lo[i] + atm_extent[i];
            }
        }
        if (!pp_erf.queryarr("n_cell", n_cell_atm, 0, AMREX_SPACEDIM)) {
            amrex::Abort("Missing required input: erf.n_cell");
        }
        amrex::Vector<int> atm_periodic_vec;
        if (!pp_erf.queryarr("is_periodic", atm_periodic_vec, 0, AMREX_SPACEDIM)) {
            amrex::Abort("Missing required input: erf.is_periodic");
        }
        atm_is_periodic = {atm_periodic_vec[0], atm_periodic_vec[1], atm_periodic_vec[2]};

        if (!pp_rem.queryarr("prob_lo", ocn_lo, 0, AMREX_SPACEDIM)) {
            amrex::Abort("Missing required input: remora.prob_lo");
        }
        if (!pp_rem.queryarr("prob_hi", ocn_hi, 0, AMREX_SPACEDIM)) {
            amrex::Abort("Missing required input: remora.prob_hi");
        }
        if (!pp_rem.queryarr("n_cell", n_cell_ocn, 0, AMREX_SPACEDIM)) {
            amrex::Abort("Missing required input: remora.n_cell");
        }
        amrex::Vector<int> ocn_periodic_vec;
        if (!pp_rem.queryarr("is_periodic", ocn_periodic_vec, 0, AMREX_SPACEDIM)) {
            amrex::Abort("Missing required input: remora.is_periodic");
        }
        ocn_is_periodic = {ocn_periodic_vec[0], ocn_periodic_vec[1], ocn_periodic_vec[2]};

        amrex::RealBox atm_rb(atm_lo.data(), atm_hi.data());
        amrex::RealBox ocn_rb(ocn_lo.data(), ocn_hi.data());

        // The mapping branch supports non-conformal atmosphere/ocean layouts.
        // Keep these details visible in the logs, but do not abort when the
        // coupled grids have different horizontal extents or incompatible
        // integer ratios. The remap layer is responsible for handling those
        // cases.
        const bool matching_horizontal_bounds =
            atm_lo[0] == ocn_lo[0] && atm_lo[1] == ocn_lo[1] &&
            atm_hi[0] == ocn_hi[0] && atm_hi[1] == ocn_hi[1];
        const bool integer_ratio_horizontal_counts =
            (n_cell_ocn[0] != 0 && n_cell_ocn[1] != 0) &&
            (n_cell_atm[0] % n_cell_ocn[0] == 0) &&
            (n_cell_atm[1] % n_cell_ocn[1] == 0);
        if (amrex::ParallelDescriptor::IOProcessor() &&
            (!matching_horizontal_bounds || !integer_ratio_horizontal_counts))
        {
            amrex::Print()
                << "WARNING: coupled grids are non-conformal."
                << " atm_lo=(" << atm_lo[0] << "," << atm_lo[1] << ")"
                << " atm_hi=(" << atm_hi[0] << "," << atm_hi[1] << ")"
                << " atm_ncell=(" << n_cell_atm[0] << "," << n_cell_atm[1] << ")"
                << " ocn_lo=(" << ocn_lo[0] << "," << ocn_lo[1] << ")"
                << " ocn_hi=(" << ocn_hi[0] << "," << ocn_hi[1] << ")"
                << " ocn_ncell=(" << n_cell_ocn[0] << "," << n_cell_ocn[1] << ")"
                << ". Proceeding with mapping-branch remap support.\n";
        }

        int max_level = 0;
        int coord = 0;
        amrex::Vector<amrex::IntVect> ref_ratio(1, amrex::IntVect(2));
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            atm_is_periodic[0] == ocn_is_periodic[0] &&
            atm_is_periodic[1] == ocn_is_periodic[1],
            "ERF and REMORA horizontal periodicity must match for coupling.");
        amrex::Array<int,AMREX_SPACEDIM> is_periodic = atm_is_periodic;

        amrex::Print() << "Running coupled ERF <-> REMORA\n";
        if (driver_verbose >= 2 && amrex::ParallelDescriptor::IOProcessor()) {
            amrex::Print()
                << "========================================================\n"
                << "[driver] Enter Coupled ERF <-> REMORA Execution\n"
                << "main.cpp is handing control to the coupled driver container.\n"
                << "========================================================\n";
        }
        erf_remora::ERFRemoraMultiBlockContainer container(
            atm_rb, n_cell_atm, ocn_rb, n_cell_ocn,
            max_level, coord, ref_ratio, is_periodic,
            driver_max_grid_size, erf_dt, atm_substeps_per_ocean_step,
            atm2ocn_contract, use_two_way_coupling);
        container.SetOceanFrameStopTime(ocean_frame_stop_time);
        container.InitializeSolvers();
        container.InitializeBlocks();
        if (restart_state.active) {
            container.SetAtmosphereTick(restart_state.atm_tick);
        }
        container.PrefillAtmosToOceanFallbacks();

        amrex::ParmParse().remove("amr.check_int");
        amrex::ParmParse().remove("amr.check_per");
        amrex::ParmParse().remove("amr.check_int_time");
        amrex::ParmParse().remove("amr.restart");
        amrex::ParmParse().remove("amr.regrid_int");
        amrex::ParmParse().remove("amr.do_substep");
        amrex::ParmParse().remove("amr.iterate_grids");

        container.AdvanceBlocks(ocean_steps);
    }
    amrex::Print() << "erf-remora-driver run completed\n";

#ifdef ERF_USE_KOKKOS
    if (Kokkos::is_initialized()) {
        Kokkos::finalize();
    }
#endif

    amrex::Finalize();

    return 0;
}
