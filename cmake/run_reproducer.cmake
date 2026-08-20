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

# Per-deck dependency lists go here, keyed on ERF_REMORA_INPUT_FILE. The decks
# that needed one (the CouplingReproducers and IdealizedFlows coupleToERF/ABL/
# couette cases, which imported base decks from ERF's .Exec_dev) were removed,
# so there is nothing to declare at present.
set(_required_files "")

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
