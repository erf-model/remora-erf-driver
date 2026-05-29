if(NOT DEFINED ERF_REMORA_DRIVER_EXE)
  message(FATAL_ERROR "ERF_REMORA_DRIVER_EXE is required")
endif()
if(NOT DEFINED ERF_REMORA_SOURCE_DIR)
  message(FATAL_ERROR "ERF_REMORA_SOURCE_DIR is required")
endif()
if(NOT DEFINED ERF_REMORA_INPUT_FILE)
  message(FATAL_ERROR "ERF_REMORA_INPUT_FILE is required")
endif()
if(NOT DEFINED ERF_REMORA_ARGS)
  set(ERF_REMORA_ARGS "")
endif()

set(_input_file "${ERF_REMORA_SOURCE_DIR}/${ERF_REMORA_INPUT_FILE}")
if(NOT EXISTS "${_input_file}")
  message(FATAL_ERROR "Missing reproducer deck: ${_input_file}")
endif()

set(_required_files "")
if(ERF_REMORA_INPUT_FILE STREQUAL "Exec/CouplingReproducers/inputs_two_way_coupleToERF")
  list(APPEND _required_files
    "${ERF_REMORA_SOURCE_DIR}/Submodules/ERF/.Exec_dev/RemoraCoupling/inputs_shear_two_way"
    "${ERF_REMORA_SOURCE_DIR}/Submodules/ERF/.Exec_dev/RemoraCoupling/input_sounding")
elseif(ERF_REMORA_INPUT_FILE STREQUAL "Exec/CouplingReproducers/inputs_ERFtoREMORA_coupleToERF")
  list(APPEND _required_files
    "${ERF_REMORA_SOURCE_DIR}/Submodules/ERF/.Exec_dev/RemoraCoupling/inputs_shear_ERFtoREMORA"
    "${ERF_REMORA_SOURCE_DIR}/Submodules/ERF/.Exec_dev/RemoraCoupling/input_sounding")
elseif(ERF_REMORA_INPUT_FILE STREQUAL "Exec/CouplingReproducers/inputs_two_way_ABL")
  list(APPEND _required_files
    "${ERF_REMORA_SOURCE_DIR}/Submodules/ERF/.Exec_dev/RemoraCoupling/inputs_abl")
elseif(ERF_REMORA_INPUT_FILE STREQUAL "Exec/CouplingReproducers/inputs_two_way_couette")
  list(APPEND _required_files
    "${ERF_REMORA_SOURCE_DIR}/Submodules/ERF/.Exec_dev/RemoraCoupling/inputs_couette")
endif()

foreach(_required_file IN LISTS _required_files)
  if(NOT EXISTS "${_required_file}")
    message(FATAL_ERROR "Missing reproducer dependency: ${_required_file}")
  endif()
endforeach()

set(_cmd "${ERF_REMORA_DRIVER_EXE}" "${_input_file}")
foreach(_arg IN LISTS ERF_REMORA_ARGS)
  if(NOT "${_arg}" STREQUAL "")
    list(APPEND _cmd "${_arg}")
  endif()
endforeach()

execute_process(
  COMMAND ${_cmd}
  RESULT_VARIABLE _result
)
if(NOT _result EQUAL 0)
  message(FATAL_ERROR "Reproducer failed with exit code ${_result}")
endif()
