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
if(NOT DEFINED ERF_REMORA_ERF_RESTART_SUFFIX)
  set(ERF_REMORA_ERF_RESTART_SUFFIX "00001")
endif()
if(NOT DEFINED ERF_REMORA_REMORA_RESTART_SUFFIX)
  set(ERF_REMORA_REMORA_RESTART_SUFFIX "00001")
endif()
if(NOT DEFINED ERF_REMORA_ERF_COMPARE_SUFFIX)
  set(ERF_REMORA_ERF_COMPARE_SUFFIX "00002")
endif()
if(NOT DEFINED ERF_REMORA_REMORA_COMPARE_SUFFIX)
  set(ERF_REMORA_REMORA_COMPARE_SUFFIX "00002")
endif()
if(NOT DEFINED ERF_REMORA_EXCHANGE_COMPARE_SUFFIX)
  set(ERF_REMORA_EXCHANGE_COMPARE_SUFFIX "${ERF_REMORA_REMORA_COMPARE_SUFFIX}")
endif()
if(NOT DEFINED ERF_REMORA_ERF_PLOT_INT)
  set(ERF_REMORA_ERF_PLOT_INT "1")
endif()
if(NOT DEFINED ERF_REMORA_ERF_CHECK_INT)
  set(ERF_REMORA_ERF_CHECK_INT "1")
endif()
if(NOT DEFINED ERF_REMORA_REMORA_PLOT_INT)
  set(ERF_REMORA_REMORA_PLOT_INT "1")
endif()
if(NOT DEFINED ERF_REMORA_REMORA_CHECK_INT)
  set(ERF_REMORA_REMORA_CHECK_INT "1")
endif()
if(NOT DEFINED ERF_REMORA_EXCHANGE_PLOT_INT)
  set(ERF_REMORA_EXCHANGE_PLOT_INT "1")
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
  "erf.plot_int_1=${ERF_REMORA_ERF_PLOT_INT}"
  "erf.check_int=${ERF_REMORA_ERF_CHECK_INT}"
  "remora.plot_int=${ERF_REMORA_REMORA_PLOT_INT}"
  "remora.check_int=${ERF_REMORA_REMORA_CHECK_INT}"
  "driver.plot_int_exchange=${ERF_REMORA_EXCHANGE_PLOT_INT}"
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
set(_cold_log "${_work_dir}/cold_run.log")
set(_restart_log "${_work_dir}/restart_run.log")

execute_process(
  COMMAND "${ERF_REMORA_DRIVER_EXE}" "${_work_dir}/${_input_name}" ${_cold_args}
  WORKING_DIRECTORY "${_work_dir}"
  OUTPUT_FILE "${_cold_log}"
  ERROR_FILE "${_cold_log}"
  RESULT_VARIABLE _cold_result
)
if(NOT _cold_result EQUAL 0)
  file(READ "${_cold_log}" _cold_log_contents)
  string(REGEX REPLACE "([^\n]*\n){0,200}$" "\\0" _cold_log_tail "${_cold_log_contents}")
  message(FATAL_ERROR
    "Cold-start restart reproducer failed with exit code ${_cold_result}\n"
    "Log: ${_cold_log}\n"
    "${_cold_log_tail}")
endif()

set(_restart_args
  "erf.restart=cold_erfchk${ERF_REMORA_ERF_RESTART_SUFFIX}"
  "remora.restart=cold_remorachk${ERF_REMORA_REMORA_RESTART_SUFFIX}"
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
  OUTPUT_FILE "${_restart_log}"
  ERROR_FILE "${_restart_log}"
  RESULT_VARIABLE _restart_result
)
if(NOT _restart_result EQUAL 0)
  file(READ "${_restart_log}" _restart_log_contents)
  string(REGEX REPLACE "([^\n]*\n){0,200}$" "\\0" _restart_log_tail "${_restart_log_contents}")
  message(FATAL_ERROR
    "Restarted reproducer failed with exit code ${_restart_result}\n"
    "Log: ${_restart_log}\n"
    "${_restart_log_tail}")
endif()

set(_erf_remora_failed_compares)
set(_erf_remora_compare_status_lines)

function(_erf_remora_compare_output _label _lhs _rhs)
  message(STATUS "Comparing ${_label}: ${_lhs} vs ${_rhs}")
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
      set(_failure
        "Restart compare failed for ${_label}: ${_lhs} vs ${_rhs} with exit code ${_compare_result}")
      set(_erf_remora_failed_compares "${_erf_remora_failed_compares}" "${_failure}" PARENT_SCOPE)
      set(_erf_remora_compare_status_lines
        "${_erf_remora_compare_status_lines}"
        "  - ${_label}: FAILS"
        PARENT_SCOPE)
    else()
      set(_erf_remora_compare_status_lines
        "${_erf_remora_compare_status_lines}"
        "  - ${_label}: AGREES"
        PARENT_SCOPE)
    endif()
  elseif(EXISTS "${_lhs}" OR EXISTS "${_rhs}")
    set(_failure "Restart compare output mismatch for ${_label}: ${_lhs} vs ${_rhs}")
    set(_erf_remora_failed_compares "${_erf_remora_failed_compares}" "${_failure}" PARENT_SCOPE)
    set(_erf_remora_compare_status_lines
      "${_erf_remora_compare_status_lines}"
      "  - ${_label}: FAILS"
      PARENT_SCOPE)
  else()
    set(_erf_remora_compare_status_lines
      "${_erf_remora_compare_status_lines}"
      "  - ${_label}: SKIPPED (BOTH OUTPUTS ABSENT)"
      PARENT_SCOPE)
  endif()
endfunction()

foreach(_pair
    "REMORA compare|cold_remoraplt${ERF_REMORA_REMORA_COMPARE_SUFFIX}|restart_remoraplt${ERF_REMORA_REMORA_COMPARE_SUFFIX}"
    "ERF compare|cold_erfplt${ERF_REMORA_ERF_COMPARE_SUFFIX}|restart_erfplt${ERF_REMORA_ERF_COMPARE_SUFFIX}"
    "REMORA packed fields compare|cold_exchange_a2o${ERF_REMORA_EXCHANGE_COMPARE_SUFFIX}|restart_exchange_a2o${ERF_REMORA_EXCHANGE_COMPARE_SUFFIX}"
    "ERF packed fields compare|cold_exchange_o2a${ERF_REMORA_EXCHANGE_COMPARE_SUFFIX}|restart_exchange_o2a${ERF_REMORA_EXCHANGE_COMPARE_SUFFIX}")
  string(REPLACE "|" ";" _pair_fields "${_pair}")
  list(GET _pair_fields 0 _label)
  list(GET _pair_fields 1 _lhs_rel)
  list(GET _pair_fields 2 _rhs_rel)
  _erf_remora_compare_output("${_label}" "${_work_dir}/${_lhs_rel}" "${_work_dir}/${_rhs_rel}")
endforeach()

if(_erf_remora_compare_status_lines)
  string(JOIN "\n" _status_summary ${_erf_remora_compare_status_lines})
  message(STATUS "===== RESTART COMPARE SUMMARY =====\n${_status_summary}")
endif()

if(_erf_remora_failed_compares)
  string(JOIN "\n" _failure_summary ${_erf_remora_failed_compares})
  message(FATAL_ERROR "${_failure_summary}")
endif()
