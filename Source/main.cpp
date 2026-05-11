#include "ERFRemoraMultiBlockContainer.H"
#include "StandaloneSolverRuns.H"

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
        int coupled_one_step = 0;
        std::string atm2ocn_mode = "state";

        pp_driver.query("n_cell", n_cell);
        pp_driver.query("max_grid_size", max_grid_size);
        pp_driver.query("dt_atm", dt_atm);
        pp_driver.query("atm_substeps_per_ocean_step", atm_substeps_per_ocean_step);
        pp_driver.query("ocean_steps", ocean_steps);
        pp_driver.query("synthetic_only", synthetic_only);
        pp_driver.query("coupled_one_step", coupled_one_step);
        pp_driver.query("atm2ocn_mode", atm2ocn_mode);

        erf_remora::coupling::AtmosToOceanContract atm2ocn_contract = erf_remora::coupling::AtmosToOceanContract::StatePassing;
        if (atm2ocn_mode == "state") {
            atm2ocn_contract = erf_remora::coupling::AtmosToOceanContract::StatePassing;
        } else if (atm2ocn_mode == "flux") {
            atm2ocn_contract = erf_remora::coupling::AtmosToOceanContract::FluxPassing;
        } else {
            amrex::Abort("Unknown driver.atm2ocn_mode. Valid values: state, flux.");
        }

        if (synthetic_only) {
            erf_remora::ERFRemoraMultiBlockContainer container(
                n_cell,
                max_grid_size,
                dt_atm,
                atm_substeps_per_ocean_step,
                atm2ocn_contract);
            container.InitializeBlocks();
            container.AdvanceBlocksSynthetic(ocean_steps);
        } else {
            if (coupled_one_step) {
                amrex::Print() << "Running one-step ERF->cache->REMORA path\n";
                erf_remora::ERFRemoraMultiBlockContainer container(
                    n_cell,
                    max_grid_size,
                    dt_atm,
                    atm_substeps_per_ocean_step,
                    atm2ocn_contract);
                container.InitializeBlocks();
                const amrex::Real coupling_time = amrex::Real(0.0);
                erf_remora::RunERFOneStepAndPack(container.AtmosToOceanActiveCache(), coupling_time);
                erf_remora::RunREMORAApplyAndOneStep(container.AtmosToOceanActiveCache(), coupling_time);
            } else {
                amrex::Print() << "Running uncoupled ERF then REMORA paths\n";
                erf_remora::RunERFStandalone();
                erf_remora::RunREMORAStandalone();
            }
        }

        amrex::Print() << "erf-remora-driver run completed\n";
    }
    amrex::Finalize();

    return 0;
}
