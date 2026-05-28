#!/bin/bash

# --- Configuration ---
DRIVER_EXE="./build/erf-remora-driver"
PROFILE="${1:-ERFtoREMORA}"
ERF_DEP_DIR="external-src/ERF/.Exec_dev/RemoraCoupling"

if [ "$PROFILE" = "two_way" ] || [ "$PROFILE" = "oceansurf" ]; then
    INPUTS_FILE="Exec/CouplingReproducers/inputs_two_way_coupleToERF"
    PROFILE_ARGS=(
        "max_step=600"
        "remora.v=1"
        "erf.v=1"
        "amr.v=1"
        "amrex.tiny_profiler=1"
    )
elif [ "$PROFILE" = "ERFtoREMORA" ] || [ "$PROFILE" = "erftoremora" ] || [ "$PROFILE" = "plumbing" ]; then
    INPUTS_FILE="Exec/CouplingReproducers/inputs_ERFtoREMORA_coupleToERF"
    PROFILE_ARGS=(
        "max_step=10"
        "driver.dt_atm=1"
        "driver.atm_substeps_per_ocean_step=1"
        "remora.fixed_dt=1"
        "erf.plot_int_1=1"
        "remora.plot_int=1"
        "remora.v=1"
        "erf.v=1"
        "amr.v=1"
        "amrex.tiny_profiler=1"
    )
else
    echo "ERROR: Unknown profile '$PROFILE'. Use 'ERFtoREMORA' or 'two_way'."
    exit 1
fi

# --- Pre-flight Checks ---
echo "[1/3] Checking dependencies..."

# Check for the driver executable
if [ ! -f "$DRIVER_EXE" ]; then
    echo "ERROR: Driver executable not found at $DRIVER_EXE. Please build the project first."
    exit 1
fi

# Check for required ERF shear deck files
if [ "$INPUTS_FILE" = "Exec/CouplingReproducers/inputs_two_way_coupleToERF" ]; then
    REQUIRED_FILES=("$ERF_DEP_DIR/inputs_shear_two_way" "$ERF_DEP_DIR/input_sounding")
else
    REQUIRED_FILES=("$ERF_DEP_DIR/inputs_shear_ERFtoREMORA" "$ERF_DEP_DIR/input_sounding")
fi
for f in "${REQUIRED_FILES[@]}"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Missing required dependency: $f"
        echo "Ensure your ERF checkout is up to date and correctly located in external-src/ERF."
        exit 1
    fi
done
echo "All dependencies verified."

# --- Execution ---
echo "[2/3] Executing profile '$PROFILE' with timing verbosity..."

# Run the simulation.
$DRIVER_EXE "$INPUTS_FILE" "${PROFILE_ARGS[@]}" | tee simulation.log

# --- Performance Summary ---
echo "[3/3] Performance Summary (extracted from simulation.log):"

# Extract the final completion message
grep "erf-remora-driver run completed" simulation.log

# Highlight the per-step timing from the log
echo "--------------------------------------------------"
echo "Per-step completion logs:"
grep "Completed ocean step" simulation.log

if [ "$PROFILE" = "two_way" ] || [ "$PROFILE" = "oceansurf" ]; then
    echo "--------------------------------------------------"
    echo "OceanSurf SST-apply diagnostics:"
    grep "OceanSurf apply" simulation.log
fi

echo "--------------------------------------------------"
echo "Tiny Profiler Summary (look for REMORA::Evolve time):"
sed -n '/TinyProfiler/,/End of TinyProfiler/p' simulation.log
