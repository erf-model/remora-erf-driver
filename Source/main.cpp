#include "ERFRemoraMultiBlockContainer.H"
#include "StandaloneSolverRuns.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_RealBox.H>

#include <string>

std::string inputs_name = "erf-remora-driver.synthetic";

int
main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        // When both FILE=CoupleToERF/inputs (remora.max_step) and the ERF shear
        // deck (max_step) are loaded, REMORA::ReadParameters aborts on the
        // conflict. The driver owns the combined ParmParse environment, so drop
        // the prefixed alias; the canonical max_step/stop_time keys remain.
        amrex::ParmParse().remove("remora.max_step");
        amrex::ParmParse().remove("remora.stop_time");

        amrex::ParmParse pp_driver("driver");

        int max_grid_size = 32;
        amrex::Real dt_atm = 10.0;
        int atm_substeps_per_ocean_step = 6;
        int ocean_steps = 2;
        int synthetic_only = 1;
        int coupled_one_step = 0;
        std::string atm2ocn_mode = "state";

        pp_driver.query("max_grid_size", max_grid_size);
        pp_driver.query("dt_atm", dt_atm);
        pp_driver.query("atm_substeps_per_ocean_step", atm_substeps_per_ocean_step);
        const bool has_driver_ocean_steps = pp_driver.query("ocean_steps", ocean_steps);
        pp_driver.query("synthetic_only", synthetic_only);
        pp_driver.query("coupled_one_step", coupled_one_step);
        pp_driver.query("atm2ocn_mode", atm2ocn_mode);

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

        int max_level = 0;
        int coord = 0;
        amrex::Vector<amrex::IntVect> ref_ratio(1, amrex::IntVect(2));
        amrex::Vector<int> is_periodic_vec = {0, 0, 0};
        pp_atm.queryarr("is_periodic", is_periodic_vec, 0, AMREX_SPACEDIM);
        amrex::Array<int,AMREX_SPACEDIM> is_periodic{is_periodic_vec[0], is_periodic_vec[1], is_periodic_vec[2]};

        if (synthetic_only) {
            erf_remora::ERFRemoraMultiBlockContainer container(
                atm_rb, n_cell_atm, ocn_rb, n_cell_ocn,
                max_level, coord, ref_ratio, is_periodic,
                max_grid_size, dt_atm, atm_substeps_per_ocean_step,
                atm2ocn_contract);
            container.InitializeBlocks();
            container.AdvanceBlocksSynthetic(ocean_steps);
        } else {
            amrex::Print() << "Running coupled ERF <-> REMORA path\n";
            erf_remora::ERFRemoraMultiBlockContainer container(
                atm_rb, n_cell_atm, ocn_rb, n_cell_ocn,
                max_level, coord, ref_ratio, is_periodic,
                max_grid_size, dt_atm, atm_substeps_per_ocean_step,
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
