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
if(NOT DEFINED MPIEXEC_EXECUTABLE)
  message(FATAL_ERROR "MPIEXEC_EXECUTABLE is required")
endif()
if(NOT DEFINED MPIEXEC_NUMPROC_FLAG)
  message(FATAL_ERROR "MPIEXEC_NUMPROC_FLAG is required")
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
if(NOT DEFINED ERF_REMORA_REMORA_PLOT_INT)
  set(ERF_REMORA_REMORA_PLOT_INT "1")
endif()
if(NOT DEFINED ERF_REMORA_EXCHANGE_PLOT_INT)
  set(ERF_REMORA_EXCHANGE_PLOT_INT "1")
endif()

if(NOT EXISTS "${ERF_REMORA_DRIVER_EXE}")
  message(FATAL_ERROR "Missing driver executable: ${ERF_REMORA_DRIVER_EXE}")
endif()

set(_input_file "${ERF_REMORA_SOURCE_DIR}/${ERF_REMORA_INPUT_FILE}")
if(NOT EXISTS "${_input_file}")
  message(FATAL_ERROR "Missing upwelling deck: ${_input_file}")
endif()

get_filename_component(_input_name "${ERF_REMORA_INPUT_FILE}" NAME)
set(_work_subdir "mpi_upwelling_compare")
if(DEFINED ERF_REMORA_WORK_SUBDIR AND NOT "${ERF_REMORA_WORK_SUBDIR}" STREQUAL "")
  set(_work_subdir "${ERF_REMORA_WORK_SUBDIR}")
endif()
set(_work_dir "${CMAKE_CURRENT_BINARY_DIR}/${_work_subdir}")
file(REMOVE_RECURSE "${_work_dir}")
file(MAKE_DIRECTORY "${_work_dir}")
file(MAKE_DIRECTORY "${_work_dir}/np1")
file(MAKE_DIRECTORY "${_work_dir}/np2")

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
    message(FATAL_ERROR "Missing MPI compare input dependency: ${_src}")
  endif()
  file(COPY "${_src}" DESTINATION "${_work_dir}")
endforeach()

set(_common_args
  "erf.plot_int_1=${ERF_REMORA_ERF_PLOT_INT}"
  "remora.plot_int=${ERF_REMORA_REMORA_PLOT_INT}"
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

function(_erf_remora_run_case _nranks _tag)
  set(_cmd "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "${_nranks}")
  if(DEFINED MPIEXEC_PREFLAGS AND NOT "${MPIEXEC_PREFLAGS}" STREQUAL "")
    separate_arguments(_mpiexec_preflags UNIX_COMMAND "${MPIEXEC_PREFLAGS}")
    list(APPEND _cmd ${_mpiexec_preflags})
  endif()
  list(APPEND _cmd "${ERF_REMORA_DRIVER_EXE}" "${_work_dir}/${_input_name}")
  list(APPEND _cmd
    "erf.plot_file_1=${_tag}/erfplt"
    "remora.plot_file=${_tag}/remoraplt"
    "driver.plot_file_exchange_a2o=${_tag}/exchange_a2o"
    "driver.plot_file_exchange_o2a=${_tag}/exchange_o2a"
  )
  list(APPEND _cmd ${_common_args})
  execute_process(
    COMMAND ${_cmd}
    WORKING_DIRECTORY "${_work_dir}"
    RESULT_VARIABLE _run_result
  )
  if(NOT _run_result EQUAL 0)
    message(FATAL_ERROR "MPI compare run (${_tag}) failed with exit code ${_run_result}")
  endif()
endfunction()

function(_erf_remora_compare_output _relative_path)
  set(_lhs "${_work_dir}/np1/${_relative_path}")
  set(_rhs "${_work_dir}/np2/${_relative_path}")
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
        "MPI compare failed for ${_relative_path} with exit code ${_compare_result}")
    endif()
  elseif(EXISTS "${_lhs}" OR EXISTS "${_rhs}")
    message(FATAL_ERROR
      "MPI compare output mismatch for ${_relative_path}: one side is missing")
  endif()
endfunction()

_erf_remora_run_case("1" "np1")
_erf_remora_run_case("2" "np2")

foreach(_relative_path
    "remoraplt${ERF_REMORA_REMORA_COMPARE_SUFFIX}"
    "erfplt${ERF_REMORA_ERF_COMPARE_SUFFIX}"
    "exchange_a2o${ERF_REMORA_EXCHANGE_COMPARE_SUFFIX}"
    "exchange_o2a${ERF_REMORA_EXCHANGE_COMPARE_SUFFIX}")
  _erf_remora_compare_output("${_relative_path}")
endforeach()
