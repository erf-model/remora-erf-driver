if(NOT DEFINED ERF_REMORA_DRIVER_EXE)
  message(FATAL_ERROR "ERF_REMORA_DRIVER_EXE is required")
endif()
if(NOT DEFINED ERF_REMORA_SOURCE_DIR)
  message(FATAL_ERROR "ERF_REMORA_SOURCE_DIR is required")
endif()
if(NOT DEFINED ERF_REMORA_INPUT_FILE)
  message(FATAL_ERROR "ERF_REMORA_INPUT_FILE is required")
endif()
if(NOT DEFINED MPIEXEC_EXECUTABLE)
  message(FATAL_ERROR "MPIEXEC_EXECUTABLE is required")
endif()
if(NOT DEFINED MPIEXEC_NUMPROC_FLAG)
  message(FATAL_ERROR "MPIEXEC_NUMPROC_FLAG is required")
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
set(_work_subdir "mpi_upwelling_smoke")
if(DEFINED ERF_REMORA_WORK_SUBDIR AND NOT "${ERF_REMORA_WORK_SUBDIR}" STREQUAL "")
  set(_work_subdir "${ERF_REMORA_WORK_SUBDIR}")
endif()
set(_work_dir "${CMAKE_CURRENT_BINARY_DIR}/${_work_subdir}")
file(REMOVE_RECURSE "${_work_dir}")
file(MAKE_DIRECTORY "${_work_dir}")
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
    message(FATAL_ERROR "Missing MPI smoke input dependency: ${_src}")
  endif()
  file(COPY "${_src}" DESTINATION "${_work_dir}")
endforeach()

set(_runtime_args
  "erf.plot_file_1=np2/erfplt"
  "erf.plot_int_1=${ERF_REMORA_ERF_PLOT_INT}"
  "remora.plot_file=np2/remoraplt"
  "remora.plot_int=${ERF_REMORA_REMORA_PLOT_INT}"
  "driver.plot_file_exchange_a2o=np2/exchange_a2o"
  "driver.plot_file_exchange_o2a=np2/exchange_o2a"
  "driver.plot_int_exchange=${ERF_REMORA_EXCHANGE_PLOT_INT}"
)
foreach(_arg IN LISTS ERF_REMORA_ARGS)
  if(NOT "${_arg}" STREQUAL "")
    list(APPEND _runtime_args "${_arg}")
  endif()
endforeach()

if(DEFINED ERF_REMORA_RRTMGP_DIR AND NOT "${ERF_REMORA_RRTMGP_DIR}" STREQUAL "")
  if(NOT IS_DIRECTORY "${ERF_REMORA_RRTMGP_DIR}")
    message(FATAL_ERROR "Missing RRTMGP directory: ${ERF_REMORA_RRTMGP_DIR}")
  endif()
  list(APPEND _runtime_args "erf.rrtmgp_file_path=${ERF_REMORA_RRTMGP_DIR}")
else()
  list(APPEND _runtime_args "erf.radiation_model=None")
endif()

set(_cmd "${MPIEXEC_EXECUTABLE}" "${MPIEXEC_NUMPROC_FLAG}" "2")
if(DEFINED MPIEXEC_PREFLAGS AND NOT "${MPIEXEC_PREFLAGS}" STREQUAL "")
  separate_arguments(_mpiexec_preflags UNIX_COMMAND "${MPIEXEC_PREFLAGS}")
  list(APPEND _cmd ${_mpiexec_preflags})
endif()
list(APPEND _cmd "${ERF_REMORA_DRIVER_EXE}" "${_work_dir}/${_input_name}")
list(APPEND _cmd ${_runtime_args})

execute_process(
  COMMAND ${_cmd}
  WORKING_DIRECTORY "${_work_dir}"
  RESULT_VARIABLE _result
)
if(NOT _result EQUAL 0)
  message(FATAL_ERROR "MPI smoke run failed with exit code ${_result}")
endif()

foreach(_expected
    "np2/remoraplt${ERF_REMORA_REMORA_COMPARE_SUFFIX}"
    "np2/exchange_a2o${ERF_REMORA_EXCHANGE_COMPARE_SUFFIX}"
    "np2/erfplt${ERF_REMORA_ERF_COMPARE_SUFFIX}")
  if(NOT EXISTS "${_work_dir}/${_expected}")
    message(FATAL_ERROR "MPI smoke run did not produce ${_work_dir}/${_expected}")
  endif()
endforeach()
