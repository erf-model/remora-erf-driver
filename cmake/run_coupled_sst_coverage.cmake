# Regression check for the coupled-SST withholding contract (driver issue #13).
#
# Coupled ocean SST must reach ERF only on the cells the ocean grid actually
# covers. Where it does not reach, ERF keeps its own lower-boundary value. The
# failure this guards against is silent: REMORA's remap zero-fills the
# destination, so an uncovered ERF cell used to arrive as a literal 0 K that was
# indistinguishable from a computed zero, and was applied as though it were one.
#
# Three runs over the same deck:
#   full    - ocean and atmosphere share extents, so every surface cell is covered
#   partial - ocean shrunk to a quarter of the atmosphere footprint, so most are not
#   sliver  - the partial ocean shifted so its edge clips one ERF column by 5e-4
#
# The partial run is the one that matters. It asserts both halves of the
# contract: that coverage is reported as partial rather than total, and that
# t_surf over the whole domain never drops to the remap's zero fill.

if(NOT DEFINED ERF_REMORA_DRIVER_EXE)
  message(FATAL_ERROR "ERF_REMORA_DRIVER_EXE is required")
endif()
if(NOT DEFINED ERF_REMORA_SOURCE_DIR)
  message(FATAL_ERROR "ERF_REMORA_SOURCE_DIR is required")
endif()
if(NOT DEFINED ERF_REMORA_INPUT_FILE)
  message(FATAL_ERROR "ERF_REMORA_INPUT_FILE is required")
endif()
if(NOT DEFINED FEXTREMA_EXE)
  message(FATAL_ERROR "FEXTREMA_EXE is required")
endif()
if(NOT DEFINED ERF_REMORA_MAX_STEP)
  set(ERF_REMORA_MAX_STEP "2")
endif()
# Any uncovered cell that took the remap's zero fill drops t_surf to ~0 K, so a
# floor well below every physical sea-surface temperature separates the two
# outcomes without pinning the test to a particular deck's values.
if(NOT DEFINED ERF_REMORA_TSURF_FLOOR_K)
  set(ERF_REMORA_TSURF_FLOOR_K "100.0")
endif()

if(NOT EXISTS "${ERF_REMORA_DRIVER_EXE}")
  message(FATAL_ERROR "Missing driver executable: ${ERF_REMORA_DRIVER_EXE}")
endif()
if(NOT EXISTS "${FEXTREMA_EXE}")
  message(FATAL_ERROR "Missing fextrema executable: ${FEXTREMA_EXE}")
endif()

set(_input_file "${ERF_REMORA_SOURCE_DIR}/${ERF_REMORA_INPUT_FILE}")
if(NOT EXISTS "${_input_file}")
  message(FATAL_ERROR "Missing coverage deck: ${_input_file}")
endif()

get_filename_component(_input_name "${ERF_REMORA_INPUT_FILE}" NAME)
get_filename_component(_deck_dir "${_input_file}" DIRECTORY)

set(_work_subdir "coupled_sst_coverage")
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
    message(FATAL_ERROR "Missing coverage deck dependency: ${_src}")
  endif()
  file(COPY "${_src}" DESTINATION "${_work_dir}")
endforeach()

set(_common_args
  "max_step=${ERF_REMORA_MAX_STEP}"
  "erf.plot_int_1=-1"
  "remora.plot_int=-1"
)
if(DEFINED ERF_REMORA_RRTMGP_DIR AND NOT "${ERF_REMORA_RRTMGP_DIR}" STREQUAL "")
  if(NOT IS_DIRECTORY "${ERF_REMORA_RRTMGP_DIR}")
    message(FATAL_ERROR "Missing RRTMGP directory: ${ERF_REMORA_RRTMGP_DIR}")
  endif()
  list(APPEND _common_args "erf.rrtmgp_file_path=${ERF_REMORA_RRTMGP_DIR}")
else()
  list(APPEND _common_args "erf.radiation_model=None")
endif()

# Run the deck and return the coverage counts the driver reports.
function(_run_case _tag _extra_args _out_covered _out_total)
  set(_log "${_work_dir}/${_tag}.log")
  set(_cmd "${ERF_REMORA_DRIVER_EXE}" "${_work_dir}/${_input_name}")
  list(APPEND _cmd ${_common_args} ${_extra_args})

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
      "Coverage run '${_tag}' failed with exit code ${_result}\n"
      "Log: ${_log}\n${_log_contents}")
  endif()

  file(READ "${_log}" _log_contents)
  # "Coupled SST apply at t=0 s: covered 250 of 1000 surface cells, ..."
  if(NOT _log_contents MATCHES "covered ([0-9]+) of ([0-9]+) surface cells")
    message(FATAL_ERROR
      "Coverage run '${_tag}' never reported a coupled SST apply. The coupled "
      "path did not run at all, so this test proves nothing.\nLog: ${_log}")
  endif()
  set(${_out_covered} "${CMAKE_MATCH_1}" PARENT_SCOPE)
  set(${_out_total}   "${CMAKE_MATCH_2}" PARENT_SCOPE)
endfunction()

# ---- full coverage: ocean and atmosphere share extents ----------------------
_run_case("full" "" _full_covered _full_total)

if(NOT _full_covered EQUAL _full_total)
  message(FATAL_ERROR
    "Full-coverage run reported ${_full_covered} of ${_full_total} cells covered. "
    "The deck's ocean and atmosphere share extents, so every surface cell should "
    "have an ocean donor; a shortfall means the coverage flag is under-reporting.")
endif()

# ---- partial coverage: ocean shrunk to a quarter of the footprint -----------
set(_partial_args
  "remora.prob_hi=25000. 10000. 0."
  "remora.n_cell=25 10 40"
  "erf.plot2d_vars_1=t_surf"
  "erf.plot2d_file_1=tsurf"
  "erf.plot2d_int_1=1"
)
_run_case("partial" "${_partial_args}" _partial_covered _partial_total)

if(_partial_covered EQUAL 0)
  message(FATAL_ERROR
    "Partial-coverage run reported no covered cells at all. The ocean still "
    "overlaps a quarter of the atmosphere footprint, so the coverage flag is "
    "empty when it should be partial.")
endif()
if(NOT _partial_covered LESS _partial_total)
  message(FATAL_ERROR
    "Partial-coverage run reported ${_partial_covered} of ${_partial_total} cells "
    "covered. The ocean was shrunk to a quarter of the atmosphere footprint, so "
    "coverage must be partial; reporting everything covered means the flag is not "
    "tracking the remap stencils.")
endif()

# ---- the withholding half: uncovered cells kept a real temperature ----------
file(GLOB _tsurf_plots "${_work_dir}/tsurf?????")
if(NOT _tsurf_plots)
  message(FATAL_ERROR "Partial-coverage run produced no t_surf plotfile in ${_work_dir}")
endif()
list(SORT _tsurf_plots)
list(GET _tsurf_plots -1 _tsurf_last)

execute_process(
  COMMAND "${FEXTREMA_EXE}" "${_tsurf_last}"
  OUTPUT_VARIABLE _extrema_out
  ERROR_VARIABLE _extrema_err
  RESULT_VARIABLE _extrema_result
)
if(NOT _extrema_result EQUAL 0)
  message(FATAL_ERROR "fextrema failed on ${_tsurf_last}\n${_extrema_out}${_extrema_err}")
endif()

if(NOT _extrema_out MATCHES "t_surf[ \t]+([0-9.eE+-]+)[ \t]+([0-9.eE+-]+)")
  message(FATAL_ERROR "Could not read t_surf extrema from:\n${_extrema_out}")
endif()
set(_tsurf_min "${CMAKE_MATCH_1}")
set(_tsurf_max "${CMAKE_MATCH_2}")

if(_tsurf_min LESS ERF_REMORA_TSURF_FLOOR_K)
  message(FATAL_ERROR
    "t_surf fell to ${_tsurf_min} K (floor ${ERF_REMORA_TSURF_FLOOR_K} K) in "
    "${_tsurf_last}.\n"
    "Only ${_partial_covered} of ${_partial_total} surface cells have an ocean "
    "donor. The uncovered cells took the remap's zero fill instead of keeping "
    "ERF's own lower-boundary temperature, which is the defect issue #13 fixed.")
endif()

# ---- sliver coverage: an ocean edge that barely clips an ERF column --------
#
# The ocean is offset by 999.5 m on a 1000 m ERF grid, so ERF's first column is
# overlapped across 0.5 m -- an area fraction of 5e-4. Such a cell must be
# reported uncovered so ERF keeps its own lower-boundary value: assigning it an
# ocean SST attributes a whole cell's temperature on a twentieth of a percent of
# geometric overlap.
#
# This is the case a coverage flag keyed on "does any overlap exist" gets wrong,
# and it is invisible to the aligned-grid cases above, where every cell is either
# wholly covered or wholly outside.
#
# The check counts covered cells, NOT an inspection of t_surf. Normalized weights
# give a sliver cell a physically plausible SST that no floor can distinguish from
# a correct value. The shifted ocean has the partial run's area but touches one
# more ERF column, so it must report the same number of covered cells.
#
# This geometry does not exercise driver.remap_min_coverage_frac. The driver marks
# ERF cells whose centres lie outside the ocean domain as land, so an edge column
# overlapped by less than half is land whatever the threshold. The threshold acts
# only on ERF water cells with a small wet overlap, which takes a REMORA land mask
# across an ERF water cell; this deck has none.
set(_sliver_args
  "remora.prob_lo=999.5 0. -300."
  "remora.prob_hi=25999.5 10000. 0."
  "remora.n_cell=25 10 40"
  "erf.plot2d_vars_1=t_surf"
  "erf.plot2d_file_1=tsurf_sliver"
  "erf.plot2d_int_1=1"
)
_run_case("sliver" "${_sliver_args}" _sliver_covered _sliver_total)

if(NOT _sliver_covered EQUAL _partial_covered)
  message(FATAL_ERROR
    "Sliver-coverage run reported ${_sliver_covered} of ${_sliver_total} cells "
    "covered; the unshifted ocean of the same size covered ${_partial_covered}. "
    "The shifted ocean clips ERF's first column across 5e-4 of its area, and that "
    "column must not be counted as covered.")
endif()

file(GLOB _sliver_plots "${_work_dir}/tsurf_sliver?????")
if(NOT _sliver_plots)
  message(FATAL_ERROR "Sliver-coverage run produced no t_surf plotfile in ${_work_dir}")
endif()
list(SORT _sliver_plots)
list(GET _sliver_plots -1 _sliver_last)

execute_process(
  COMMAND "${FEXTREMA_EXE}" "${_sliver_last}"
  OUTPUT_VARIABLE _sliver_extrema
  ERROR_VARIABLE  _sliver_extrema_err
  RESULT_VARIABLE _sliver_extrema_result
)
if(NOT _sliver_extrema_result EQUAL 0)
  message(FATAL_ERROR
    "fextrema failed on ${_sliver_last}\n${_sliver_extrema}${_sliver_extrema_err}")
endif()
if(NOT _sliver_extrema MATCHES "t_surf[ \t]+([0-9.eE+-]+)[ \t]+([0-9.eE+-]+)")
  message(FATAL_ERROR "Could not read t_surf extrema from:\n${_sliver_extrema}")
endif()
set(_sliver_tsurf_min "${CMAKE_MATCH_1}")

if(_sliver_tsurf_min LESS ERF_REMORA_TSURF_FLOOR_K)
  message(FATAL_ERROR
    "t_surf fell to ${_sliver_tsurf_min} K (floor ${ERF_REMORA_TSURF_FLOOR_K} K) in "
    "${_sliver_last}.\n"
    "The ocean clips ERF's first column across 5e-4 of its area. That cell was "
    "treated as covered and took a fraction of one ocean cell's SST instead of "
    "keeping ERF's own value.")
endif()

message(STATUS
  "===== COUPLED SST COVERAGE SUMMARY =====\n"
  "  - full coverage:    ${_full_covered} of ${_full_total} cells covered\n"
  "  - partial coverage: ${_partial_covered} of ${_partial_total} cells covered\n"
  "  - sliver coverage:  ${_sliver_covered} of ${_sliver_total} cells covered "
  "(must equal the partial run's ${_partial_covered})\n"
  "  - t_surf range, partial run: ${_tsurf_min} / ${_tsurf_max} K\n"
  "  - t_surf min, sliver run:    ${_sliver_tsurf_min} K\n"
  "  - uncovered and sliver-covered cells both retained ERF's own value "
  "(floor ${ERF_REMORA_TSURF_FLOOR_K} K)")
