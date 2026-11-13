# Cold-start atmos->ocean averaging over two atmosphere substeps (issue #18).
#
# On the first exchange interval of a cold start the driver skips the
# trapezoid's left endpoint, because ERF's diagnostic flux arrays are still
# zero, and reconstructs it as f(t) ~= 2*f1 - f2. With two substeps per
# exchange the weights are 2 on f1 and 0 on f2, so the interval mean is f1
# exactly. A lane that is not averaged keeps the interval-end sample, f2.
#
# Both references come from non-averaged runs: one substep per exchange writes
# f1, two substeps per exchange write f2. Coupling is forced one-way, so ERF's
# first two steps do not depend on the ocean step and every run sees the same
# f1 and f2. Each lane of the averaged run's first exchange plotfile must equal
# its reference at zero tolerance.
#
# What the check can distinguish depends on the build. Without RRTMGP,
# radiation is off and rain, SWrad, LWrad, Qair and cloud are constant over the
# first two steps of this deck, so for those lanes f1 = f2 and a wrong weight
# would go unnoticed. With RRTMGP (ERF_REMORA_RRTMGP_DIR set) SWrad and LWrad
# vary as well.

if(NOT DEFINED ERF_REMORA_DRIVER_EXE)
  message(FATAL_ERROR "ERF_REMORA_DRIVER_EXE is required")
endif()
if(NOT DEFINED ERF_REMORA_SOURCE_DIR)
  message(FATAL_ERROR "ERF_REMORA_SOURCE_DIR is required")
endif()
if(NOT DEFINED ERF_REMORA_INPUT_FILE)
  message(FATAL_ERROR "ERF_REMORA_INPUT_FILE is required")
endif()
if(NOT DEFINED FCOMPARE_EXE)
  message(FATAL_ERROR "FCOMPARE_EXE is required")
endif()
if(NOT EXISTS "${ERF_REMORA_DRIVER_EXE}")
  message(FATAL_ERROR "Missing driver executable: ${ERF_REMORA_DRIVER_EXE}")
endif()

set(_input_file "${ERF_REMORA_SOURCE_DIR}/${ERF_REMORA_INPUT_FILE}")
if(NOT EXISTS "${_input_file}")
  message(FATAL_ERROR "Missing cold-start averaging deck: ${_input_file}")
endif()
get_filename_component(_input_name "${ERF_REMORA_INPUT_FILE}" NAME)
get_filename_component(_deck_dir "${_input_file}" DIRECTORY)

set(_work_dir "${CMAKE_CURRENT_BINARY_DIR}/cold_start_average")
file(REMOVE_RECURSE "${_work_dir}")
file(MAKE_DIRECTORY "${_work_dir}")

# The deck's FILE= includes resolve against the working directory.
foreach(_fname "${_input_name}" inputs_ERF inputs_REMORA input_sounding_for_ERF)
  set(_src "${_deck_dir}/${_fname}")
  if(NOT EXISTS "${_src}")
    message(FATAL_ERROR "Missing cold-start averaging deck dependency: ${_src}")
  endif()
  file(COPY "${_src}" DESTINATION "${_work_dir}")
endforeach()

separate_arguments(_fcompare_cmd UNIX_COMMAND "${FCOMPARE_EXE}")
if(NOT _fcompare_cmd)
  message(FATAL_ERROR "FCOMPARE_EXE resolved to an empty command")
endif()

set(_common_args
  "max_step=1"
  "driver.coupling_direction=ERFtoREMORA"
  "erf.plot_int_1=-1"
  "erf.check_int=-1"
  "remora.plot_int=-1"
  "remora.check_int=-1"
  "erf.fixed_dt=1"
  "remora.ndtfast=1"
  "driver.plot_int_exchange=1"
)
if(DEFINED ERF_REMORA_RRTMGP_DIR AND NOT "${ERF_REMORA_RRTMGP_DIR}" STREQUAL "")
  list(APPEND _common_args "erf.rrtmgp_file_path=${ERF_REMORA_RRTMGP_DIR}")
else()
  list(APPEND _common_args "erf.radiation_model=None")
endif()

# Lane vocabularies and characters, from Source/ERFRemoraCouplingContract.H.
set(_flux_lanes_flux_contract tau_x tau_y SHflux LHflux SWrad LWrad rain evap)
set(_state_lanes_state_contract Uwind Vwind Patm Qair Tair cloud)
set(_flux_lanes_state_contract rain SWrad LWrad)

function(_run _tag)
  set(_log "${_work_dir}/${_tag}.log")
  execute_process(
    COMMAND "${ERF_REMORA_DRIVER_EXE}" "${_work_dir}/${_input_name}"
            ${_common_args}
            "driver.plot_file_exchange_a2o=${_tag}_a2o"
            "driver.plot_file_exchange_o2a=${_tag}_o2a"
            ${ARGN}
    WORKING_DIRECTORY "${_work_dir}"
    OUTPUT_FILE "${_log}"
    ERROR_FILE "${_log}"
    RESULT_VARIABLE _result
  )
  if(NOT _result EQUAL 0)
    file(READ "${_log}" _log_contents)
    message(FATAL_ERROR
      "Run '${_tag}' failed with exit code ${_result}\nLog: ${_log}\n${_log_contents}")
  endif()
  if(NOT IS_DIRECTORY "${_work_dir}/${_tag}_a2o00000")
    message(FATAL_ERROR "Run '${_tag}' wrote no ${_tag}_a2o00000 exchange plotfile")
  endif()
endfunction()

# Require each named lane to have zero absolute error between the first
# exchange plotfiles of two runs. fcompare cannot select variables, so its
# per-variable table is parsed; its exit status is ignored because lanes
# outside the list are allowed to differ.
function(_expect_lanes_equal _case _test_tag _ref_tag _what)
  execute_process(
    COMMAND ${_fcompare_cmd} "${_test_tag}_a2o00000" "${_ref_tag}_a2o00000"
    WORKING_DIRECTORY "${_work_dir}"
    OUTPUT_VARIABLE _out
    ERROR_VARIABLE _out
  )
  foreach(_lane IN LISTS ARGN)
    if(NOT _out MATCHES "\n *${_lane} +([^ \n]+) +[^ \n]+\n")
      message(FATAL_ERROR
        "${_case}: fcompare reported no ${_lane} for ${_test_tag} vs ${_ref_tag}\n${_out}")
    endif()
    if(NOT CMAKE_MATCH_1 STREQUAL "0")
      message(FATAL_ERROR
        "${_case}: ${_lane} is not ${_what} (absolute error ${CMAKE_MATCH_1} "
        "against ${_ref_tag})\n${_out}")
    endif()
  endforeach()
  message(STATUS "${_case}: ${ARGN} equal ${_what}")
endfunction()

foreach(_contract flux state)
  _run(${_contract}_instant_n1 "driver.atm2ocn_mode=${_contract}" "remora.fixed_dt=1")
endforeach()
_run(state_instant_n2 "driver.atm2ocn_mode=state" "remora.fixed_dt=2")

# Every lane averaged: all of them must be f1.
_run(flux_averaged_n2
  "driver.atm2ocn_mode=flux" "remora.fixed_dt=2"
  "driver.time_average_atm2ocn_flux=1")
_expect_lanes_equal("flux contract, all lanes averaged"
  flux_averaged_n2 flux_instant_n1 "the mean f1" ${_flux_lanes_flux_contract})

_run(state_averaged_n2
  "driver.atm2ocn_mode=state" "remora.fixed_dt=2"
  "driver.time_average_atm2ocn_state=1" "driver.time_average_atm2ocn_flux=1")
_expect_lanes_equal("state contract, all lanes averaged"
  state_averaged_n2 state_instant_n1 "the mean f1"
  ${_state_lanes_state_contract} ${_flux_lanes_state_contract})

# Only the flux-character lanes averaged: those must be f1, and the state lanes
# must keep the interval-end sample f2, which carries zero averaging weight.
_run(state_flux_averaged_n2
  "driver.atm2ocn_mode=state" "remora.fixed_dt=2"
  "driver.time_average_atm2ocn_flux=1")
_expect_lanes_equal("state contract, flux lanes averaged"
  state_flux_averaged_n2 state_instant_n1 "the mean f1" ${_flux_lanes_state_contract})
_expect_lanes_equal("state contract, flux lanes averaged"
  state_flux_averaged_n2 state_instant_n2 "the interval-end sample f2"
  ${_state_lanes_state_contract})
