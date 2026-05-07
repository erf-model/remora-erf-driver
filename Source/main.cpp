#include "ERFRemoraMultiBlockContainer.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <string>

std::string inputs_name = "erf-remora-driver.synthetic";

int
main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        amrex::ParmParse pp_driver("driver");

        int n_cell = 32;
        int max_grid_size = 32;
        amrex::Real dt_atm = 10.0;
        int atm_substeps_per_ocean_step = 6;
        int ocean_steps = 2;
        int synthetic_only = 1;

        pp_driver.query("n_cell", n_cell);
        pp_driver.query("max_grid_size", max_grid_size);
        pp_driver.query("dt_atm", dt_atm);
        pp_driver.query("atm_substeps_per_ocean_step", atm_substeps_per_ocean_step);
        pp_driver.query("ocean_steps", ocean_steps);
        pp_driver.query("synthetic_only", synthetic_only);

        erf_remora::ERFRemoraMultiBlockContainer container(
            n_cell,
            max_grid_size,
            dt_atm,
            atm_substeps_per_ocean_step);
        container.InitializeBlocks();

        if (synthetic_only) {
            container.AdvanceBlocksSynthetic(ocean_steps);
        } else {
            amrex::Abort("driver.synthetic_only = 0 is not supported in stage-0.");
        }

        amrex::Print() << "erf-remora-driver synthetic exchange completed\n";
    }
    amrex::Finalize();

    return 0;
}
