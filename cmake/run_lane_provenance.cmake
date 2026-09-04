# Regression check for atmos->ocean lane provenance (driver issue #4).
#
# Every lane the driver hands REMORA has to be traceable to exactly one source,
# and the record of that source has to agree with the data actually written.
# Three surfaces describe the active view, and they are produced by different
# code paths:
#
#   1. the "Atmos->ocean active lanes:" line InitializeBlocks prints
#   2. the variable list in the exchange plotfile's Header
#   3. the lane table in the exchange_info sidecar beside that plotfile
#
# If they disagree, one of them is lying about what REMORA received: a lane in
# the Header with no row in the table is forcing applied with no provenance, and
# a row with no Header entry is provenance for a lane ERF never wrote.
#
# The table also names each lane's CHARACTER, and character alone decides whether
# a lane is time-averaged. That matters because a lane whose global index fails to
# resolve is reported as state/not-averaged -- indistinguishable from a genuine
# un-averaged state lane. So the character column is checked against this file's
# own statement of the contract rather than against whatever the driver emitted,
# and the averaged column is checked against what the character plus the flags
# imply. A lane that lost its identity shows up as a character mismatch.
#
# The state contract is the case worth spending a run on: it carries six
# state-character lanes and three flux-character ones, so with exactly one
# averaging flag set the table has to split its own view down the middle.

# Script mode starts with no policies set, and this file uses IN_LIST.
cmake_policy(SET CMP0057 NEW)

if(NOT DEFINED ERF_REMORA_DRIVER_EXE)
  message(FATAL_ERROR "ERF_REMORA_DRIVER_EXE is required")
endif()
if(NOT DEFINED ERF_REMORA_SOURCE_DIR)
  message(FATAL_ERROR "ERF_REMORA_SOURCE_DIR is required")
endif()
if(NOT DEFINED ERF_REMORA_INPUT_FILE)
  message(FATAL_ERROR "ERF_REMORA_INPUT_FILE is required")
endif()
if(NOT DEFINED ERF_REMORA_MAX_STEP)
  set(ERF_REMORA_MAX_STEP "1")
endif()

if(NOT EXISTS "${ERF_REMORA_DRIVER_EXE}")
  message(FATAL_ERROR "Missing driver executable: ${ERF_REMORA_DRIVER_EXE}")
endif()

set(_input_file "${ERF_REMORA_SOURCE_DIR}/${ERF_REMORA_INPUT_FILE}")
if(NOT EXISTS "${_input_file}")
  message(FATAL_ERROR "Missing lane provenance deck: ${_input_file}")
endif()

get_filename_component(_input_name "${ERF_REMORA_INPUT_FILE}" NAME)
get_filename_component(_deck_dir "${_input_file}" DIRECTORY)

set(_work_subdir "lane_provenance")
if(DEFINED ERF_REMORA_WORK_SUBDIR AND NOT "${ERF_REMORA_WORK_SUBDIR}" STREQUAL "")
  set(_work_subdir "${ERF_REMORA_WORK_SUBDIR}")
endif()
set(_work_dir "${CMAKE_CURRENT_BINARY_DIR}/${_work_subdir}")
file(REMOVE_RECURSE "${_work_dir}")
file(MAKE_DIRECTORY "${_work_dir}")

# The deck's FILE= includes resolve against the working directory, so the
# imported base decks have to travel with it.
set(_required_files
  "${_input_name}"
  "inputs_ERF"
  "inputs_REMORA"
  "input_sounding_for_ERF"
)
foreach(_fname IN LISTS _required_files)
  set(_src "${_deck_dir}/${_fname}")
  if(NOT EXISTS "${_src}")
    message(FATAL_ERROR "Missing lane provenance deck dependency: ${_src}")
  endif()
  file(COPY "${_src}" DESTINATION "${_work_dir}")
endforeach()

# One atmosphere substep per ocean step keeps the runs short. Averaging over a
# single substep is still reported as an interval mean, and it is the label this
# test reads, not the value.
set(_common_args
  "max_step=${ERF_REMORA_MAX_STEP}"
  "erf.plot_int_1=-1"
  "remora.plot_int=-1"
  "erf.fixed_dt=1"
  "remora.fixed_dt=1"
  "remora.ndtfast=1"
  "driver.plot_int_exchange=1"
)
if(DEFINED ERF_REMORA_RRTMGP_DIR AND NOT "${ERF_REMORA_RRTMGP_DIR}" STREQUAL "")
  if(NOT IS_DIRECTORY "${ERF_REMORA_RRTMGP_DIR}")
    message(FATAL_ERROR "Missing RRTMGP directory: ${ERF_REMORA_RRTMGP_DIR}")
  endif()
  list(APPEND _common_args "erf.rrtmgp_file_path=${ERF_REMORA_RRTMGP_DIR}")
else()
  list(APPEND _common_args "erf.radiation_model=None")
endif()

# This file's own statement of the contract, deliberately not read from the
# driver. Source of truth is Source/ERFRemoraCouplingContract.H:
# AtmosToOceanNames(contract) for the vocabularies, and the flux-character index
# set in WriteExchangeInfo for the characters.
set(_vocab_flux  tau_x tau_y SHflux LHflux SWrad LWrad rain evap)
set(_vocab_state Uwind Vwind Patm Qair Tair cloud rain SWrad LWrad)
set(_flux_character tau_x tau_y SHflux LHflux SWrad LWrad rain evap)

# Run one configuration and assert its three provenance surfaces agree with each
# other and with the contract.
#
#   _tag            work/plotfile prefix for this case
#   _contract       state|flux, the value driver.atm2ocn_mode is given
#   _extra_args     additional runtime arguments
#   _flux_averaged  yes|no, expected "averaged" for every flux-character lane
#   _state_averaged yes|no, expected "averaged" for every state-character lane
function(_check_case _tag _contract _extra_args _flux_averaged _state_averaged)
  if(_contract STREQUAL "flux")
    set(_vocab ${_vocab_flux})
  else()
    set(_vocab ${_vocab_state})
  endif()

  set(_log "${_work_dir}/${_tag}.log")
  set(_prefix "${_tag}_a2o")
  set(_cmd "${ERF_REMORA_DRIVER_EXE}" "${_work_dir}/${_input_name}")
  list(APPEND _cmd
    ${_common_args}
    "driver.atm2ocn_mode=${_contract}"
    "driver.plot_file_exchange_a2o=${_prefix}"
    "driver.plot_file_exchange_o2a=${_tag}_o2a"
    ${_extra_args})

  execute_process(
    COMMAND ${_cmd}
    WORKING_DIRECTORY "${_work_dir}"
    OUTPUT_FILE "${_log}"
    ERROR_FILE "${_log}"
    RESULT_VARIABLE _result
  )
  if(NOT _result EQUAL 0)
    file(READ "${_log}" _log_contents)
    message(FATAL_ERROR
      "Lane provenance run '${_tag}' failed with exit code ${_result}\n"
      "Log: ${_log}\n${_log_contents}")
  endif()

  # ---- surface 1: the active view the driver announced --------------------
  file(READ "${_log}" _log_contents)
  if(NOT _log_contents MATCHES "Atmos->ocean active lanes: ([^\n]*)")
    message(FATAL_ERROR
      "Run '${_tag}' never announced an atmos->ocean active view. Without it "
      "there is nothing to check the written lanes against.\nLog: ${_log}")
  endif()
  string(STRIP "${CMAKE_MATCH_1}" _announced_raw)
  if(_announced_raw STREQUAL "")
    message(FATAL_ERROR
      "Run '${_tag}' announced an empty atmos->ocean active view, so REMORA "
      "received no forcing at all from ERF.\nLog: ${_log}")
  endif()
  string(REPLACE ", " ";" _announced "${_announced_raw}")

  # ---- surface 2: the variables actually written --------------------------
  file(GLOB _plots "${_work_dir}/${_prefix}?????")
  if(NOT _plots)
    message(FATAL_ERROR
      "Run '${_tag}' wrote no ${_prefix} exchange plotfile, so the announced "
      "view was never checked against anything REMORA could receive.\n"
      "Log: ${_log}")
  endif()
  list(SORT _plots)
  list(GET _plots -1 _plot)

  file(STRINGS "${_plot}/Header" _hdr)
  list(GET _hdr 1 _ncomp)
  set(_written "")
  math(EXPR _last_var "${_ncomp} + 1")
  foreach(_i RANGE 2 ${_last_var})
    list(GET _hdr ${_i} _name)
    list(APPEND _written "${_name}")
  endforeach()

  # ---- surface 3: the provenance record beside the data -------------------
  set(_info "${_plot}/exchange_info")
  if(NOT EXISTS "${_info}")
    message(FATAL_ERROR
      "Run '${_tag}' wrote ${_plot} with no exchange_info beside it. The data "
      "reached REMORA with no record of where each lane came from.")
  endif()
  file(READ "${_info}" _info_contents)

  # The coupling block is what makes the record self-describing; a plotfile whose
  # sidecar omits it cannot be interpreted after the fact.
  foreach(_field "atm2ocn contract" "two-way coupling" "remap frame" "donor fill mode")
    if(NOT _info_contents MATCHES "${_field} +: +[^\n]+")
      message(FATAL_ERROR
        "exchange_info for '${_tag}' records no '${_field}'.\nFile: ${_info}")
    endif()
  endforeach()
  if(NOT _info_contents MATCHES "atm2ocn contract +: +([A-Za-z_]+)")
    message(FATAL_ERROR "Could not read the contract from ${_info}")
  endif()
  if(NOT CMAKE_MATCH_1 STREQUAL _contract)
    message(FATAL_ERROR
      "Run '${_tag}' was given driver.atm2ocn_mode=${_contract} but its "
      "exchange_info records the contract as '${CMAKE_MATCH_1}'. The provenance "
      "record does not describe the run that produced it.\nFile: ${_info}")
  endif()

  string(FIND "${_info_contents}" "--- lanes ---" _lanes_at)
  if(_lanes_at LESS 0)
    message(FATAL_ERROR "exchange_info for '${_tag}' has no lane table.\nFile: ${_info}")
  endif()
  string(SUBSTRING "${_info_contents}" ${_lanes_at} -1 _lanes_section)
  string(REGEX MATCHALL "\n  [A-Za-z_0-9]+ +(state|flux) +(yes|no|n/a) +[^\n]*"
         _rows "${_lanes_section}")
  if(NOT _rows)
    message(FATAL_ERROR
      "exchange_info for '${_tag}' has a lane table with no lane rows in it.\n"
      "File: ${_info}")
  endif()

  set(_recorded "")
  foreach(_row IN LISTS _rows)
    if(NOT _row MATCHES "\n  ([A-Za-z_0-9]+) +(state|flux) +(yes|no|n/a) +(.*)$")
      message(FATAL_ERROR "Unparsable lane row in ${_info}:\n${_row}")
    endif()
    set(_name "${CMAKE_MATCH_1}")
    set(_character "${CMAKE_MATCH_2}")
    set(_averaged "${CMAKE_MATCH_3}")
    string(STRIP "${CMAKE_MATCH_4}" _represents)
    list(APPEND _recorded "${_name}")

    # Character is the lane's identity. A lane whose index failed to resolve is
    # reported as state, so a flux lane appearing as state means the row is a
    # default, not a description.
    if("${_name}" IN_LIST _flux_character)
      set(_expect_character "flux")
      set(_expect_averaged "${_flux_averaged}")
    else()
      set(_expect_character "state")
      set(_expect_averaged "${_state_averaged}")
    endif()
    if(NOT _character STREQUAL _expect_character)
      message(FATAL_ERROR
        "Run '${_tag}': lane '${_name}' is recorded as ${_character}-character, "
        "but the contract makes it ${_expect_character}-character. Character "
        "decides whether a lane is averaged, so a lane carrying the wrong one is "
        "forced at the wrong instant. A lane whose index did not resolve reports "
        "'state' by default, which is what this looks like.\nFile: ${_info}")
    endif()
    if(NOT _averaged STREQUAL _expect_averaged)
      message(FATAL_ERROR
        "Run '${_tag}': lane '${_name}' (${_character}-character) is recorded as "
        "averaged=${_averaged}, expected ${_expect_averaged} for this "
        "combination of contract and averaging flags.\nFile: ${_info}")
    endif()

    # 'represents' is the authoritative statement of the instant, so it has to
    # follow the averaged column rather than drift from it.
    if(_averaged STREQUAL "yes")
      if(NOT _represents MATCHES "^interval mean over \\[")
        message(FATAL_ERROR
          "Run '${_tag}': lane '${_name}' is averaged but represents "
          "'${_represents}', which does not describe an interval.\nFile: ${_info}")
      endif()
    else()
      if(NOT _represents MATCHES "^instantaneous at t=")
        message(FATAL_ERROR
          "Run '${_tag}': lane '${_name}' is not averaged but represents "
          "'${_represents}', which does not describe an instant.\nFile: ${_info}")
      endif()
    endif()
  endforeach()

  # ---- exactly one source per lane ----------------------------------------
  set(_deduped ${_recorded})
  list(REMOVE_DUPLICATES _deduped)
  list(LENGTH _recorded _n_recorded)
  list(LENGTH _deduped _n_deduped)
  if(NOT _n_recorded EQUAL _n_deduped)
    string(REPLACE ";" ", " _recorded_str "${_recorded}")
    message(FATAL_ERROR
      "Run '${_tag}': exchange_info lists a lane more than once (${_recorded_str}). "
      "A lane with two rows has two sources, and nothing says which one REMORA "
      "was given.\nFile: ${_info}")
  endif()

  # ---- the three surfaces have to agree -----------------------------------
  string(REPLACE ";" ", " _announced_str "${_announced}")
  string(REPLACE ";" ", " _written_str "${_written}")
  string(REPLACE ";" ", " _recorded_str "${_recorded}")

  if(NOT _announced_str STREQUAL _written_str)
    message(FATAL_ERROR
      "Run '${_tag}': the driver announced the active view as\n"
      "  ${_announced_str}\n"
      "but wrote\n"
      "  ${_written_str}\n"
      "to ${_plot}. REMORA received a different set of lanes than the driver "
      "reported selecting.")
  endif()
  if(NOT _recorded_str STREQUAL _written_str)
    message(FATAL_ERROR
      "Run '${_tag}': the exchange plotfile carries\n"
      "  ${_written_str}\n"
      "but its exchange_info records\n"
      "  ${_recorded_str}\n"
      "A written lane with no row is forcing applied with no provenance; a row "
      "with no written lane is provenance for a lane ERF never wrote.\n"
      "File: ${_info}")
  endif()

  # ---- and every lane has to belong to the contract, in its order ---------
  set(_expected_order "")
  foreach(_name IN LISTS _vocab)
    if("${_name}" IN_LIST _recorded)
      list(APPEND _expected_order "${_name}")
    endif()
  endforeach()
  string(REPLACE ";" ", " _expected_str "${_expected_order}")
  string(REPLACE ";" ", " _vocab_str "${_vocab}")
  if(NOT _recorded_str STREQUAL _expected_str)
    message(FATAL_ERROR
      "Run '${_tag}': the ${_contract} contract declares\n"
      "  ${_vocab_str}\n"
      "but the active view is\n"
      "  ${_recorded_str}\n"
      "Expected the contract's lanes in the contract's order, so this view "
      "either names a lane outside the contract or reorders it. Lane order is "
      "the plotfile component order REMORA reads.")
  endif()

  message(STATUS "lane provenance '${_tag}' (${_contract}): ${_recorded_str}")
endfunction()

# ---- flux contract: every lane is flux-character ---------------------------
# Nothing is averaged, so every lane must report an instant. This is also the
# case that would hide a resolution failure: a defaulted lane reports "state"
# here, and there is no state lane in this contract for it to hide among.
_check_case("flux_plain" "flux" "" "no" "no")

# ---- state contract, flux lanes averaged -----------------------------------
# Six state-character lanes and three flux-character ones under a single
# plotfile time. The table has to split its own view: rain, SWrad and LWrad
# report an interval mean while the other six report an instant.
_check_case("state_avg_flux" "state"
  "driver.time_average_atm2ocn_flux=1" "yes" "no")

# ---- state contract, state lanes averaged ----------------------------------
# The mirror image, so neither column can be a constant that happens to fit.
_check_case("state_avg_state" "state"
  "driver.time_average_atm2ocn_state=1" "no" "yes")

# ---- one-way means one way -------------------------------------------------
# The deck couples ERF->REMORA only. The provenance record has to say so, and no
# ocean->atmos exchange may have been written for REMORA to have fed back.
file(GLOB _o2a_plots "${_work_dir}/*_o2a?????")
file(READ "${_work_dir}/flux_plain_a2o00000/exchange_info" _first_info)
if(_first_info MATCHES "two-way coupling +: +yes")
  if(NOT _o2a_plots)
    message(FATAL_ERROR
      "exchange_info records two-way coupling on a deck configured "
      "ERFtoREMORA, and no ocean->atmos exchange was written. The provenance "
      "record contradicts the run.")
  endif()
elseif(_o2a_plots)
  string(REPLACE ";" ", " _o2a_str "${_o2a_plots}")
  message(FATAL_ERROR
    "The deck is one-way (ERFtoREMORA) and exchange_info agrees, but an "
    "ocean->atmos exchange was written anyway: ${_o2a_str}")
endif()

message(STATUS
  "===== LANE PROVENANCE SUMMARY =====\n"
  "  - announced view, written plotfile and exchange_info agree in all 3 runs\n"
  "  - every lane belongs to its contract, once, in contract order\n"
  "  - character and averaging agree with the contract and the flags\n"
  "  - one-way deck wrote no ocean->atmos exchange")
