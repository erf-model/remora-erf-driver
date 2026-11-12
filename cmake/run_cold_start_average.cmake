# Cold-start atmos->ocean averaging over two atmosphere substeps (issue #18).
#
# On the first exchange interval of a cold start the driver skips the
# trapezoid's left endpoint, because ERF's diagnostic flux arrays are still
# zero, and reconstructs it as f(t) ~= 2*f1 - f2. With two substeps per
# exchange the weights are 2 on f1 and 0 on f2, so the interval mean is f1
# exactly.
#
# The reference for f1 comes from a second run: one atmosphere substep per
# exchange with averaging off, which packs the atmosphere state after its first
# step. The deck is one-way, so ERF's first step does not depend on the ocean
# step and both runs see the same f1. The first exchange plotfiles of the two
# runs must then match to zero tolerance.

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

# Each contract is checked with every lane it carries averaged. The flux
# contract carries only flux-character lanes, and the driver rejects the state
# flag there.
foreach(_contract flux state)
  set(_average_args "driver.time_average_atm2ocn_flux=1")
  if(_contract STREQUAL "state")
    list(APPEND _average_args "driver.time_average_atm2ocn_state=1")
  endif()
  _run(${_contract}_averaged_n2
    "driver.atm2ocn_mode=${_contract}"
    "remora.fixed_dt=2"
    ${_average_args})
  _run(${_contract}_instant_n1
    "driver.atm2ocn_mode=${_contract}"
    "remora.fixed_dt=1")

  execute_process(
    COMMAND ${_fcompare_cmd} -r 0 --abs_tol 0
            "${_contract}_averaged_n2_a2o00000"
            "${_contract}_instant_n1_a2o00000"
    WORKING_DIRECTORY "${_work_dir}"
    OUTPUT_VARIABLE _out
    ERROR_VARIABLE _out
    RESULT_VARIABLE _result
  )
  if(NOT _result EQUAL 0)
    message(FATAL_ERROR
      "${_contract} contract: the cold-start two-substep mean is not the "
      "atmosphere state after its first substep.\n${_out}")
  endif()
  message(STATUS "${_contract} contract: cold-start two-substep mean equals f1")
endforeach()
