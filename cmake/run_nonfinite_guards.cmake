# Regression check for the non-finite guards on exchanged fields (driver issue #36).
#
# Each case writes a single NaN into one cell of one packed lane
# (driver.debug_inject_nan_lane) and requires the run to abort with the message of
# the guard that is supposed to catch it. A single cell is the case the old
# min/max-based guards could miss: amrex::min(d, r) is (r < d) ? r : d, so a NaN
# survives the reduction only if it happens to be visited first.
#
# amrex.fpe_trap_invalid is switched off, so it is the driver's own guard that
# stops the run and not a floating-point trap. On GPU nothing traps, and these
# guards are all there is.
#
# Optional MPI: pass MPIEXEC_EXECUTABLE, MPIEXEC_NUMPROC_FLAG and
# ERF_REMORA_MPI_NPROCS to run every case under mpiexec.

if(NOT DEFINED ERF_REMORA_DRIVER_EXE)
  message(FATAL_ERROR "ERF_REMORA_DRIVER_EXE is required")
endif()
if(NOT DEFINED ERF_REMORA_SOURCE_DIR)
  message(FATAL_ERROR "ERF_REMORA_SOURCE_DIR is required")
endif()
if(NOT DEFINED ERF_REMORA_INPUT_FILE)
  message(FATAL_ERROR "ERF_REMORA_INPUT_FILE is required")
endif()
if(NOT EXISTS "${ERF_REMORA_DRIVER_EXE}")
  message(FATAL_ERROR "Missing driver executable: ${ERF_REMORA_DRIVER_EXE}")
endif()

set(_input_file "${ERF_REMORA_SOURCE_DIR}/${ERF_REMORA_INPUT_FILE}")
if(NOT EXISTS "${_input_file}")
  message(FATAL_ERROR "Missing non-finite guard deck: ${_input_file}")
endif()
get_filename_component(_input_name "${ERF_REMORA_INPUT_FILE}" NAME)
get_filename_component(_deck_dir "${_input_file}" DIRECTORY)

# The serial run is tagged "serial", not "np1", so that it never shares a work
# directory with an MPI run at one rank.
set(_launcher "")
set(_rank_tag "serial")
if(DEFINED ERF_REMORA_MPI_NPROCS AND NOT "${ERF_REMORA_MPI_NPROCS}" STREQUAL "")
  if(NOT DEFINED MPIEXEC_EXECUTABLE OR NOT DEFINED MPIEXEC_NUMPROC_FLAG)
    message(FATAL_ERROR "ERF_REMORA_MPI_NPROCS needs MPIEXEC_EXECUTABLE and MPIEXEC_NUMPROC_FLAG")
  endif()
  set(_launcher "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${ERF_REMORA_MPI_NPROCS}")
  if(DEFINED MPIEXEC_PREFLAGS AND NOT "${MPIEXEC_PREFLAGS}" STREQUAL "")
    list(APPEND _launcher ${MPIEXEC_PREFLAGS})
  endif()
  set(_rank_tag "np${ERF_REMORA_MPI_NPROCS}")
endif()

set(_work_dir "${CMAKE_CURRENT_BINARY_DIR}/nonfinite_guards_${_rank_tag}")
file(REMOVE_RECURSE "${_work_dir}")
file(MAKE_DIRECTORY "${_work_dir}")

# The deck's FILE= includes resolve against the working directory, so the
# imported base decks have to travel with it.
foreach(_fname "${_input_name}" "inputs_ERF" "inputs_REMORA" "input_sounding_for_ERF")
  set(_src "${_deck_dir}/${_fname}")
  if(NOT EXISTS "${_src}")
    message(FATAL_ERROR "Missing non-finite guard deck dependency: ${_src}")
  endif()
  file(COPY "${_src}" DESTINATION "${_work_dir}")
endforeach()

set(_common_args
  "max_step=1"
  "amrex.fpe_trap_invalid=0"
  "erf.plot_int_1=-1"
  "remora.plot_int=-1"
  "erf.fixed_dt=1"
  "remora.fixed_dt=1"
  "remora.ndtfast=1"
)
if(DEFINED ERF_REMORA_RRTMGP_DIR AND NOT "${ERF_REMORA_RRTMGP_DIR}" STREQUAL "")
  list(APPEND _common_args "erf.rrtmgp_file_path=${ERF_REMORA_RRTMGP_DIR}")
else()
  list(APPEND _common_args "erf.radiation_model=None")
endif()

# Run one poisoned configuration and require it to abort with _expect.
#
#   _tag         log name for this case
#   _extra_args  runtime arguments selecting the lane, stage and contract
#   _expect      regex the abort message must match
function(_check_case _tag _extra_args _expect)
  set(_log "${_work_dir}/${_tag}.log")
  execute_process(
    COMMAND ${_launcher} "${ERF_REMORA_DRIVER_EXE}" "${_work_dir}/${_input_name}"
            ${_common_args} ${_extra_args}
    WORKING_DIRECTORY "${_work_dir}"
    OUTPUT_FILE "${_log}"
    ERROR_FILE "${_log}"
    RESULT_VARIABLE _result
  )
  file(READ "${_log}" _log_contents)
  if(_result EQUAL 0)
    message(FATAL_ERROR
      "Case '${_tag}' (${_rank_tag}) ran to completion with a NaN in an exchanged "
      "lane. The guard that should have stopped it did not fire.\nLog: ${_log}")
  endif()
  if(NOT _log_contents MATCHES "driver.debug_inject_nan_lane: NaN written at")
    message(FATAL_ERROR
      "Case '${_tag}' (${_rank_tag}) failed before the NaN was injected, so it "
      "says nothing about the guard (exit code ${_result}).\nLog: ${_log}\n${_log_contents}")
  endif()
  if(NOT _log_contents MATCHES "${_expect}")
    message(FATAL_ERROR
      "Case '${_tag}' (${_rank_tag}) aborted (exit code ${_result}), but not with "
      "the expected message\n  ${_expect}\nLog: ${_log}\n${_log_contents}")
  endif()
  message(STATUS "non-finite guard '${_tag}' (${_rank_tag}): caught")
endfunction()

# Ocean->atmos: the post-pack guard on the surface state.
_check_case("o2a_pack_sst"
  "driver.atm2ocn_mode=state;driver.debug_inject_nan_lane=SST"
  "Non-finite value detected in PackSurfaceState output, lane 'SST'")

# Atmos->ocean: the post-pack guard on the cell-centred state lanes.
_check_case("a2o_pack_tair"
  "driver.atm2ocn_mode=state;driver.debug_inject_nan_lane=Tair"
  "Non-finite value detected in PackAtmosphericStates output, lane 'Tair'")

# Atmos->ocean, behind the post-pack guard: only ValidateFields' masked count
# stands between the NaN and REMORA. Once on a state lane and once on a stress
# lane; both are on the REMORA rho layout and masked with mskr.
_check_case("a2o_validate_tair"
  "driver.atm2ocn_mode=state;driver.v=1;driver.debug_inject_nan_lane=Tair;driver.debug_inject_nan_after_guard=1"
  "Post-pack atmos->ocean validation: 1 non-finite value\\(s\\) in field 'Tair' \\(wet only\\)")
_check_case("a2o_validate_taux"
  "driver.atm2ocn_mode=flux;driver.v=1;driver.debug_inject_nan_lane=tau_x;driver.debug_inject_nan_after_guard=1"
  "Post-pack atmos->ocean validation: 1 non-finite value\\(s\\) in field 'tau_x' \\(wet only\\)")
