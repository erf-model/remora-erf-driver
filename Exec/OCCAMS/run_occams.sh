#!/usr/bin/env bash
#
# Run the OCCAMS coupled cases and summarize the driver's coupling diagnostics.
#
# Build and environment setup are deliberately NOT here -- they are machine- and
# user-specific. See the Perlmutter recipe in INSTALL.md. This script assumes an
# already-built driver and an allocation you can srun into:
#
#   salloc -N 1 --constraint cpu -A m4304 -t 01:00:00 --qos interactive
#   ./run_occams.sh
#
# Usage:
#   ./run_occams.sh [case ...]
#
#   With no arguments, runs oneway, state, and flux in that order.
#   Cases are the suffixes of inputs_OCCAMS_*: oneway, state, flux, smoke.
#
# Environment overrides:
#   BUILD_DIR   build directory holding erf-remora-driver (default ../../build_ncutil)
#   NRANKS      MPI ranks (default 16)
#   NNODES      nodes (default 1)
#   LOGDIR      where to write logs (default .)
#   EXTRA_ARGS  additional driver key=value overrides

set -euo pipefail

BUILD_DIR="${BUILD_DIR:-../../build_ncutil}"
NRANKS="${NRANKS:-16}"
NNODES="${NNODES:-1}"
LOGDIR="${LOGDIR:-.}"

# Overrides applied to the production cases only.
#
# FPE trapping is disabled because the coupled OCCAMS case currently raises
# invalid operations during startup. This is a workaround, not a setting to
# keep: if the run is clean, drop it and let the trap catch real problems.
# Note the parameter is amrex.fpe_trap_invalid -- amrex.trap_fpe_invalid is not
# a real AMReX parameter and is silently ignored.
#
# The fixed ocean timestep keeps the production cases comparable.
#
# These are NOT applied to inputs_OCCAMS_smoke. AMReX reads the inputs file
# first and then the command line (ppinit in AMReX_ParmParse.cpp), so command
# line entries win -- and the smoke deck deliberately sets
# amrex.fpe_trap_invalid = 1 and remora.fixed_dt = 20.0 to trap FPEs and to get
# two ERF steps per ocean step. Overriding either one defeats the deck.
CASE_ARGS=(amrex.fpe_trap_invalid=0 remora.fixed_dt=10)

EXTRA_ARGS="${EXTRA_ARGS:-}"

cd "$(dirname "${BASH_SOURCE[0]}")"

DRIVER="${BUILD_DIR}/erf-remora-driver"
if [[ ! -x "${DRIVER}" ]]; then
    echo "error: no driver at ${DRIVER}" >&2
    echo "       build first, or set BUILD_DIR" >&2
    exit 1
fi

# Guard against reading a stale binary's output as if it were a new build. A
# failed build leaves the previous executable in place, so without this an
# aborted rebuild silently reruns old code -- exactly the case where you are
# looking for diagnostics that the old binary does not emit.
echo "driver: ${DRIVER}"
echo "built : $(date -r "${DRIVER}" '+%Y-%m-%d %H:%M:%S')"
echo

if ! command -v srun >/dev/null 2>&1; then
    echo "error: srun not found; are you in an allocation?" >&2
    exit 1
fi

CASES=("$@")
if [[ ${#CASES[@]} -eq 0 ]]; then
    CASES=(oneway state flux)
fi

# Pull out the lines that say whether the coupling did what it should. These are
# documented in Exec/README.md under "Reading driver output".
summarize () {
    local log="$1"
    echo
    echo "--- coupling diagnostics from $(basename "${log}") ---"
    grep -E 'remap frame =|O2A coverage:|remap staging plan|O2A remap staging plan' "${log}" || true
    # Warnings are the silent-degradation cases: donor fill, stencil overflow,
    # wet cells with no atmospheric source, non-overlapping longitudes.
    if grep -qE '^\[driver\] (WARNING|remap: [0-9])' "${log}"; then
        echo "  warnings:"
        grep -E '^\[driver\] (WARNING|remap: [0-9])' "${log}" | sed 's/^/    /'
    else
        echo "  warnings: none"
    fi
}

failed=()
for case_name in "${CASES[@]}"; do
    deck="inputs_OCCAMS_${case_name}"
    log="${LOGDIR}/out_${case_name}.txt"

    if [[ ! -f "${deck}" ]]; then
        echo "error: no deck ${deck}" >&2
        echo "       available: $(ls inputs_OCCAMS_* 2>/dev/null | tr '\n' ' ')" >&2
        exit 1
    fi

    # The smoke deck is self-contained and its settings are the point of it,
    # so it runs with no overrides.
    args=()
    if [[ "${case_name}" != "smoke" ]]; then
        args=("${CASE_ARGS[@]}")
    fi

    echo "========================================================"
    echo "case ${case_name}: ${deck} -> ${log}"
    echo "overrides: ${args[*]:-none} ${EXTRA_ARGS}"
    echo "========================================================"

    # Each case gets its own log. Reusing one filename across cases silently
    # discards the earlier run.
    if srun -N "${NNODES}" -n "${NRANKS}" "${DRIVER}" "${deck}" \
            ${args[@]+"${args[@]}"} ${EXTRA_ARGS} 2>&1 | tee "${log}"; then
        summarize "${log}"
    else
        echo "case ${case_name} FAILED; see ${log}" >&2
        summarize "${log}"
        failed+=("${case_name}")
    fi
    echo
done

if [[ ${#failed[@]} -gt 0 ]]; then
    echo "failed cases: ${failed[*]}" >&2
    exit 1
fi

echo "all cases completed: ${CASES[*]}"
