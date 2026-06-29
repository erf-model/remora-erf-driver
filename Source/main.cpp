#include "ERFRemoraMultiBlockContainer.H"
#include "StandaloneSolverRuns.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_RealBox.H>

#include <algorithm>
#include <cctype>
#include <string>

std::string inputs_name = "erf-remora-driver.synthetic";

int
main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        // The driver owns the combined ParmParse environment, so drop
        // the prefixed alias; the canonical max_step/stop_time keys remain.
        amrex::ParmParse().remove("remora.max_step");
        amrex::ParmParse().remove("remora.stop_time");

        amrex::ParmParse pp_driver("driver");

        int max_grid_size = 32;
        amrex::Real erf_dt = 10.0; // Replaced dt_atm with erf_dt
        int atm_substeps_per_ocean_step = 6;
        int ocean_steps = 2;
        int synthetic_only = 1;
        std::string atm2ocn_mode = "state";
        std::string coupling_direction;

        pp_driver.query("max_grid_size", max_grid_size);
        pp_driver.query("erf_dt", erf_dt);
        pp_driver.query("atm_substeps_per_ocean_step", atm_substeps_per_ocean_step);
        const bool has_driver_ocean_steps = pp_driver.query("ocean_steps", ocean_steps);
        pp_driver.query("synthetic_only", synthetic_only);
        pp_driver.query("atm2ocn_mode", atm2ocn_mode);
        const bool has_coupling_direction = pp_driver.query("coupling_direction", coupling_direction);

        amrex::Real stop_time = -1.0;
        const bool has_stop_time = amrex::ParmParse().query("stop_time", stop_time);

        if (pp_driver.contains("coupled_one_step")) {
            amrex::Print() << "WARNING: driver.coupled_one_step is deprecated and ignored; "
                           << "use deck naming plus driver.coupling_direction instead.\n";
        }

        auto lowercase = [](std::string value) {
            std::transform(value.begin(), value.end(), value.begin(),
                           [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            return value;
        };

        int top_level_max_step = ocean_steps;
        const bool has_top_level_max_step = amrex::ParmParse().query("max_step", top_level_max_step);
        if (has_top_level_max_step) {
            if (has_driver_ocean_steps && top_level_max_step != ocean_steps) {
                amrex::Print() << "WARNING: both max_step (" << top_level_max_step
                               << ") and driver.ocean_steps (" << ocean_steps
                               << ") are set; using max_step for ocean macro-step count.\n";
            }
            ocean_steps = top_level_max_step;
        }

        const amrex::Real ocean_macro_dt =
            erf_dt * static_cast<amrex::Real>(atm_substeps_per_ocean_step);
        if (has_stop_time) {
            const int stop_time_macro_steps =
                static_cast<int>(stop_time / ocean_macro_dt);
            if (stop_time_macro_steps < ocean_steps) {
                if (amrex::ParallelDescriptor::IOProcessor()) {
                    amrex::Print() << "Clipping coupled ocean_steps from " << ocean_steps
                                   << " to " << stop_time_macro_steps
                                   << " so the driver stops cleanly at stop_time="
                                   << stop_time << ".\n";
                }
                ocean_steps = stop_time_macro_steps;
            }
            const amrex::Real clipped_stop_time =
                ocean_macro_dt * static_cast<amrex::Real>(stop_time_macro_steps);
            if (amrex::ParallelDescriptor::IOProcessor() &&
                stop_time_macro_steps >= 0 &&
                stop_time > clipped_stop_time) {
                amrex::Print() << "WARNING: coupled driver only stops on ocean macro-step boundaries;"
                               << " stop_time=" << stop_time
                               << " is being clipped to " << clipped_stop_time
                               << " for driver-side termination.\n";
            }
        }

        const amrex::Real expected_run_time = erf_dt * static_cast<amrex::Real>(atm_substeps_per_ocean_step)
                                             * static_cast<amrex::Real>(ocean_steps);
        if (amrex::ParallelDescriptor::IOProcessor()) {
            amrex::Print() << "Coupled run plan:"
                           << " expected_run_time=" << expected_run_time
                           << " erf_dt=" << erf_dt
                           << " atm_substeps_per_ocean_step=" << atm_substeps_per_ocean_step
                           << " ocean_steps=" << ocean_steps;
            if (has_stop_time) {
                amrex::Print() << " stop_time=" << stop_time
                               << " (either stop_time or max_step may terminate the run first)";
            } else {
                amrex::Print() << " stop_time=<unset>";
            }
            amrex::Print() << "\n";
        }

        if (pp_driver.contains("ocean_step") && !pp_driver.contains("ocean_steps")) {
            amrex::Print() << "WARNING: found driver.ocean_step, but driver.ocean_steps is the supported key; ignoring driver.ocean_step.\n";
        }

        erf_remora::coupling::AtmosToOceanContract atm2ocn_contract = erf_remora::coupling::AtmosToOceanContract::StatePassing;
        if (atm2ocn_mode == "state") {
            atm2ocn_contract = erf_remora::coupling::AtmosToOceanContract::StatePassing;
        } else if (atm2ocn_mode == "flux") {
            amrex::Abort("driver.atm2ocn_mode=flux is disabled for V1 infrastructure. Only 'state' mode is supported.");
        } else {
            amrex::Abort("Unknown driver.atm2ocn_mode. Valid values: state, flux.");
        }

        if (has_coupling_direction) {
            const std::string coupling_direction_lc = lowercase(coupling_direction);
            if (!synthetic_only) {
                amrex::ParmParse pp_erf("erf");
                pp_erf.get("fixed_dt", erf_dt);
                std::string land_surface_model;
                const bool has_land_surface_model = pp_erf.query("land_surface_model", land_surface_model);
                const std::string land_surface_model_lc = lowercase(land_surface_model);

                if (coupling_direction_lc == "two_way") {
                    if (!has_land_surface_model || land_surface_model_lc != "oceansurf") {
                        amrex::Abort("driver.coupling_direction=two_way requires "
                                     "erf.land_surface_model=OceanSurf so REMORA SST feeds back into ERF.");
                    }
                } else if (coupling_direction_lc == "erftoremora") {
                    if (has_land_surface_model && land_surface_model_lc == "oceansurf") {
                        amrex::Abort("driver.coupling_direction=ERFtoREMORA is inconsistent with "
                                     "erf.land_surface_model=OceanSurf. Use driver.coupling_direction=two_way instead.");
                    }
                } else {
                    amrex::Abort("Unknown driver.coupling_direction. Valid values: ERFtoREMORA, two_way.");
                }
            } else {
                amrex::Print() << "WARNING: driver.coupling_direction is ignored when "
                               << "driver.synthetic_only=1.\n";
            }
        }

        // Parse driver-owned topology from input file.
        // Atmosphere: z >= 0, Ocean: z <= 0.
        amrex::ParmParse pp_atm("driver.atm");
        amrex::ParmParse pp_ocn("driver.ocn");

        amrex::Vector<amrex::Real> atm_lo = {0.0, 0.0, 0.0};
        amrex::Vector<amrex::Real> atm_hi = {1.0, 1.0, 1.0};
        amrex::Vector<int> n_cell_atm = {32, 32, 32};
        amrex::Vector<amrex::Real> ocn_lo = {0.0, 0.0, -1.0};
        amrex::Vector<amrex::Real> ocn_hi = {1.0, 1.0,  0.0};
        amrex::Vector<int> n_cell_ocn = {32, 32, 10};

        pp_atm.queryarr("prob_lo", atm_lo, 0, AMREX_SPACEDIM);
        pp_atm.queryarr("prob_hi", atm_hi, 0, AMREX_SPACEDIM);
        pp_atm.queryarr("n_cell",  n_cell_atm, 0, AMREX_SPACEDIM);
        pp_ocn.queryarr("prob_lo", ocn_lo, 0, AMREX_SPACEDIM);
        pp_ocn.queryarr("prob_hi", ocn_hi, 0, AMREX_SPACEDIM);
        pp_ocn.queryarr("n_cell",  n_cell_ocn, 0, AMREX_SPACEDIM);

        amrex::RealBox atm_rb(atm_lo.data(), atm_hi.data());
        amrex::RealBox ocn_rb(ocn_lo.data(), ocn_hi.data());

        // Current supported topology is atmosphere same-or-finer than ocean on a
        // shared horizontal domain. ERF -> REMORA uses averaging; REMORA -> ERF
        // handles the reverse SST prolongation internally.
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            atm_lo[0] == ocn_lo[0] && atm_lo[1] == ocn_lo[1] &&
            atm_hi[0] == ocn_hi[0] && atm_hi[1] == ocn_hi[1],
            "Atmosphere and ocean horizontal domains must be identical.");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            n_cell_atm[0] % n_cell_ocn[0] == 0 && n_cell_atm[1] % n_cell_ocn[1] == 0,
            "Current coupling requires atmosphere horizontal cell counts to be integer multiples of ocean horizontal cell counts; ocean-finer-than-atmosphere topologies are unsupported.");

        int max_level = 0;
        int coord = 0;
        amrex::Vector<amrex::IntVect> ref_ratio(1, amrex::IntVect(2));
        amrex::Vector<int> is_periodic_vec = {1, 1, 0};
        const bool have_driver_periodic = pp_atm.queryarr("is_periodic", is_periodic_vec, 0, AMREX_SPACEDIM);
        if (!have_driver_periodic) {
            amrex::ParmParse pp_geom("geometry");
            pp_geom.queryarr("is_periodic", is_periodic_vec, 0, AMREX_SPACEDIM);
        }
        amrex::Array<int,AMREX_SPACEDIM> is_periodic{is_periodic_vec[0], is_periodic_vec[1], is_periodic_vec[2]};

        if (synthetic_only) {
            erf_remora::ERFRemoraMultiBlockContainer container(
                atm_rb, n_cell_atm, ocn_rb, n_cell_ocn,
                max_level, coord, ref_ratio, is_periodic,
                max_grid_size, erf_dt, atm_substeps_per_ocean_step,
                atm2ocn_contract);
            container.InitializeBlocks();
            container.AdvanceBlocksSynthetic(ocean_steps);
        } else {
            amrex::Print() << "Running coupled ERF <-> REMORA path\n";
            erf_remora::ERFRemoraMultiBlockContainer container(
                atm_rb, n_cell_atm, ocn_rb, n_cell_ocn,
                max_level, coord, ref_ratio, is_periodic,
                max_grid_size, erf_dt, atm_substeps_per_ocean_step,
                atm2ocn_contract);
            container.InitializeBlocks();
            container.InitializeSolvers();
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
    }
    amrex::Print() << "erf-remora-driver run completed\n";
    amrex::Finalize();

    return 0;
}
