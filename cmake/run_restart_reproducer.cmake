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
if(NOT DEFINED ERF_REMORA_TEST_FCOMPARE_RTOL)
  set(ERF_REMORA_TEST_FCOMPARE_RTOL "2.0e-10")
endif()
if(NOT DEFINED ERF_REMORA_TEST_FCOMPARE_ATOL)
  set(ERF_REMORA_TEST_FCOMPARE_ATOL "2.0e-10")
endif()
if(NOT DEFINED ERF_REMORA_ARGS)
  set(ERF_REMORA_ARGS "")
endif()

if(NOT EXISTS "${ERF_REMORA_DRIVER_EXE}")
  message(FATAL_ERROR "Missing driver executable: ${ERF_REMORA_DRIVER_EXE}")
endif()

set(_input_file "${ERF_REMORA_SOURCE_DIR}/${ERF_REMORA_INPUT_FILE}")
if(NOT EXISTS "${_input_file}")
  message(FATAL_ERROR "Missing restart deck: ${_input_file}")
endif()

get_filename_component(_input_name "${ERF_REMORA_INPUT_FILE}" NAME)
set(_work_subdir "restart_reproducer")
if(DEFINED ERF_REMORA_WORK_SUBDIR AND NOT "${ERF_REMORA_WORK_SUBDIR}" STREQUAL "")
  set(_work_subdir "${ERF_REMORA_WORK_SUBDIR}")
endif()
set(_work_dir "${CMAKE_CURRENT_BINARY_DIR}/${_work_subdir}")
file(REMOVE_RECURSE "${_work_dir}")
file(MAKE_DIRECTORY "${_work_dir}")

set(_upwelling_dir "${ERF_REMORA_SOURCE_DIR}/Exec/Upwelling")
set(_required_files
  "${_input_name}"
  "inputs_ERF"
  "inputs_REMORA"
  "input_sounding_for_ERF"
)
foreach(_fname IN LISTS _required_files)
  set(_src "${_upwelling_dir}/${_fname}")
  if(NOT EXISTS "${_src}")
    message(FATAL_ERROR "Missing restart input dependency: ${_src}")
  endif()
  file(COPY "${_src}" DESTINATION "${_work_dir}")
endforeach()

set(_common_args
  "max_step=2"
  "erf.plot_int_1=1"
  "erf.check_int=1"
  "remora.plot_int=1"
  "remora.check_int=1"
  "driver.plot_int_exchange=1"
)
foreach(_arg IN LISTS ERF_REMORA_ARGS)
  if(NOT "${_arg}" STREQUAL "")
    list(APPEND _common_args "${_arg}")
  endif()
endforeach()

if(DEFINED ERF_REMORA_RRTMGP_DIR AND NOT "${ERF_REMORA_RRTMGP_DIR}" STREQUAL "")
  if(NOT IS_DIRECTORY "${ERF_REMORA_RRTMGP_DIR}")
    message(FATAL_ERROR "Missing RRTMGP directory: ${ERF_REMORA_RRTMGP_DIR}")
  endif()
  list(APPEND _common_args "erf.rrtmgp_file_path=${ERF_REMORA_RRTMGP_DIR}")
else()
  list(APPEND _common_args "erf.radiation_model=None")
endif()

separate_arguments(_fcompare_cmd UNIX_COMMAND "${FCOMPARE_EXE}")
if(NOT _fcompare_cmd)
  message(FATAL_ERROR "FCOMPARE_EXE resolved to an empty command")
endif()

set(_cold_args
  "erf.plot_file_1=cold_erfplt"
  "erf.check_file=cold_erfchk"
  "remora.plot_file=cold_remoraplt"
  "remora.check_file=cold_remorachk"
  "driver.plot_file_exchange_a2o=cold_exchange_a2o"
  "driver.plot_file_exchange_o2a=cold_exchange_o2a"
)
list(APPEND _cold_args ${_common_args})

execute_process(
  COMMAND "${ERF_REMORA_DRIVER_EXE}" "${_work_dir}/${_input_name}" ${_cold_args}
  WORKING_DIRECTORY "${_work_dir}"
  RESULT_VARIABLE _cold_result
)
if(NOT _cold_result EQUAL 0)
  message(FATAL_ERROR "Cold-start restart reproducer failed with exit code ${_cold_result}")
endif()

set(_restart_args
  "erf.restart=cold_erfchk00001"
  "remora.restart=cold_remorachk00001"
  "erf.plot_file_1=restart_erfplt"
  "erf.check_file=restart_erfchk"
  "remora.plot_file=restart_remoraplt"
  "remora.check_file=restart_remorachk"
  "driver.plot_file_exchange_a2o=restart_exchange_a2o"
  "driver.plot_file_exchange_o2a=restart_exchange_o2a"
)
list(APPEND _restart_args ${_common_args})

execute_process(
  COMMAND "${ERF_REMORA_DRIVER_EXE}" "${_work_dir}/${_input_name}" ${_restart_args}
  WORKING_DIRECTORY "${_work_dir}"
  RESULT_VARIABLE _restart_result
)
if(NOT _restart_result EQUAL 0)
  message(FATAL_ERROR "Restarted reproducer failed with exit code ${_restart_result}")
endif()

function(_erf_remora_compare_output _lhs _rhs)
  if(EXISTS "${_lhs}" AND EXISTS "${_rhs}")
    execute_process(
      COMMAND ${_fcompare_cmd}
              -a
              -r "${ERF_REMORA_TEST_FCOMPARE_RTOL}"
              --abs_tol "${ERF_REMORA_TEST_FCOMPARE_ATOL}"
              "${_lhs}"
              "${_rhs}"
      WORKING_DIRECTORY "${_work_dir}"
      RESULT_VARIABLE _compare_result
    )
    if(NOT _compare_result EQUAL 0)
      message(FATAL_ERROR
        "Restart compare failed for ${_lhs} vs ${_rhs} with exit code ${_compare_result}")
    endif()
  elseif(EXISTS "${_lhs}" OR EXISTS "${_rhs}")
    message(FATAL_ERROR "Restart compare output mismatch for ${_lhs} vs ${_rhs}")
  endif()
endfunction()

foreach(_pair
    "cold_remoraplt00002|restart_remoraplt00002"
    "cold_erfplt00002|restart_erfplt00002"
    "cold_exchange_a2o00002|restart_exchange_a2o00002"
    "cold_exchange_o2a00002|restart_exchange_o2a00002")
  string(REPLACE "|" ";" _pair_fields "${_pair}")
  list(GET _pair_fields 0 _lhs_rel)
  list(GET _pair_fields 1 _rhs_rel)
  _erf_remora_compare_output("${_work_dir}/${_lhs_rel}" "${_work_dir}/${_rhs_rel}")
endforeach()
