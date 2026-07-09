if(NOT DEFINED ERF_REMORA_DRIVER_EXE)
  message(FATAL_ERROR "ERF_REMORA_DRIVER_EXE is required")
endif()
if(NOT DEFINED ERF_REMORA_SOURCE_DIR)
  message(FATAL_ERROR "ERF_REMORA_SOURCE_DIR is required")
endif()
if(NOT DEFINED FCOMPARE_EXE)
  message(FATAL_ERROR "FCOMPARE_EXE is required")
endif()
if(NOT DEFINED ERF_REMORA_RRTMGP_DIR)
  message(FATAL_ERROR "ERF_REMORA_RRTMGP_DIR is required")
endif()
if(NOT DEFINED ERF_REMORA_TEST_FCOMPARE_RTOL)
  set(ERF_REMORA_TEST_FCOMPARE_RTOL "2.0e-10")
endif()
if(NOT DEFINED ERF_REMORA_TEST_FCOMPARE_ATOL)
  set(ERF_REMORA_TEST_FCOMPARE_ATOL "2.0e-10")
endif()

if(NOT EXISTS "${ERF_REMORA_DRIVER_EXE}")
  message(FATAL_ERROR "Missing driver executable: ${ERF_REMORA_DRIVER_EXE}")
endif()
if(NOT IS_DIRECTORY "${ERF_REMORA_RRTMGP_DIR}")
  message(FATAL_ERROR "Missing RRTMGP directory: ${ERF_REMORA_RRTMGP_DIR}")
endif()

set(_work_dir "${CMAKE_CURRENT_BINARY_DIR}/restart_substep_work")
file(REMOVE_RECURSE "${_work_dir}")
file(MAKE_DIRECTORY "${_work_dir}")

set(_upwelling_dir "${ERF_REMORA_SOURCE_DIR}/Exec/Upwelling")
set(_required_files
  "inputs_debug_cold_start_substep"
  "inputs_debug_restart_substep"
  "inputs_perlin_ERF_REMORA"
  "inputs_ERF"
  "inputs_REMORA"
  "input_sounding_for_ERF"
)

foreach(_fname IN LISTS _required_files)
  set(_src "${_upwelling_dir}/${_fname}")
  if(NOT EXISTS "${_src}")
    message(FATAL_ERROR "Missing restart reproducer input dependency: ${_src}")
  endif()
  file(COPY "${_src}" DESTINATION "${_work_dir}")
endforeach()

set(_common_args "erf.rrtmgp_file_path=${ERF_REMORA_RRTMGP_DIR}")
separate_arguments(_fcompare_cmd UNIX_COMMAND "${FCOMPARE_EXE}")
if(NOT _fcompare_cmd)
  message(FATAL_ERROR "FCOMPARE_EXE resolved to an empty command")
endif()

execute_process(
  COMMAND "${ERF_REMORA_DRIVER_EXE}" "${_work_dir}/inputs_debug_cold_start_substep" "${_common_args}"
  WORKING_DIRECTORY "${_work_dir}"
  RESULT_VARIABLE _cold_result
)
if(NOT _cold_result EQUAL 0)
  message(FATAL_ERROR "Cold-start restart reproducer failed with exit code ${_cold_result}")
endif()

execute_process(
  COMMAND "${ERF_REMORA_DRIVER_EXE}" "${_work_dir}/inputs_debug_restart_substep" "${_common_args}"
  WORKING_DIRECTORY "${_work_dir}"
  RESULT_VARIABLE _restart_result
)
if(NOT _restart_result EQUAL 0)
  message(FATAL_ERROR "Restarted reproducer failed with exit code ${_restart_result}")
endif()

execute_process(
  COMMAND ${_fcompare_cmd} -a
          "${_work_dir}/cold_remoraplt00002"
          "${_work_dir}/remoraplt00002"
  WORKING_DIRECTORY "${_work_dir}"
  RESULT_VARIABLE _remora_compare_result
)
if(NOT _remora_compare_result EQUAL 0)
  message(FATAL_ERROR "REMORA restart compare failed with exit code ${_remora_compare_result}")
endif()

execute_process(
  COMMAND ${_fcompare_cmd} -a
          -r "${ERF_REMORA_TEST_FCOMPARE_RTOL}"
          --abs_tol "${ERF_REMORA_TEST_FCOMPARE_ATOL}"
          "${_work_dir}/cold_erfplt00020"
          "${_work_dir}/restart_erfplt00020"
  WORKING_DIRECTORY "${_work_dir}"
  RESULT_VARIABLE _erf_compare_result
)
if(NOT _erf_compare_result EQUAL 0)
  message(FATAL_ERROR "ERF restart compare failed with exit code ${_erf_compare_result}")
endif()
