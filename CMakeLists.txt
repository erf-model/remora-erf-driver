cmake_minimum_required(VERSION 3.24)

project(erf-remora-driver CXX C)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ERF_REMORA_ENABLE_MPI "Enable MPI across AMReX/ERF/REMORA" ON)
option(ERF_REMORA_ENABLE_OPENMP "Enable OpenMP across AMReX/ERF/REMORA" OFF)
option(ERF_REMORA_ENABLE_TESTS "Enable top-level tests" OFF)
set(_erf_remora_enable_fcompare_default OFF)
if(ERF_REMORA_ENABLE_TESTS)
  set(_erf_remora_enable_fcompare_default ON)
endif()
option(ERF_REMORA_ENABLE_FCOMPARE
       "Enable fcompare build for top-level regression tests"
       ${_erf_remora_enable_fcompare_default})
set(ERF_REMORA_TEST_FCOMPARE_RTOL "2.0e-10" CACHE STRING
    "Relative tolerance for top-level ERF fcompare-based regression checks")
set(ERF_REMORA_TEST_FCOMPARE_ATOL "2.0e-10" CACHE STRING
    "Absolute tolerance for top-level ERF fcompare-based regression checks")
set(ERF_REMORA_RRTMGP_DIR "" CACHE PATH
    "Path to RRTMGP data directory for top-level restart regression tests")
message(STATUS "ERF-REMORA symbol policy: prefixed (fixed)")

if(ERF_REMORA_ENABLE_OPENMP)
  set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "" FORCE)
  set(Kokkos_ENABLE_SERIAL ON CACHE BOOL "" FORCE)
endif()

set(AMReX_SPACEDIM 3 CACHE STRING "" FORCE)
set(AMReX_MPI ${ERF_REMORA_ENABLE_MPI} CACHE BOOL "" FORCE)
set(AMReX_OMP ${ERF_REMORA_ENABLE_OPENMP} CACHE BOOL "" FORCE)
set(AMReX_PRECISION DOUBLE CACHE STRING "" FORCE)
set(AMReX_EB ON CACHE BOOL "" FORCE)
set(AMReX_FORTRAN OFF CACHE BOOL "" FORCE)
set(AMReX_FORTRAN_INTERFACES OFF CACHE BOOL "" FORCE)
set(AMReX_AMRDATA OFF CACHE BOOL "" FORCE)
set(AMReX_PARTICLES ON CACHE BOOL "" FORCE)
set(AMReX_LINEAR_SOLVERS ON CACHE BOOL "" FORCE)
set(AMReX_PLOTFILE_TOOLS ${ERF_REMORA_ENABLE_FCOMPARE} CACHE BOOL "" FORCE)
set(AMReX_INSTALL OFF CACHE BOOL "" FORCE)
set(AMReX_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
set(AMReX_ASSERTIONS ON CACHE BOOL "" FORCE)
set(AMReX_FPE OFF CACHE BOOL "" FORCE)
set(AMReX_GPU_BACKEND NONE CACHE STRING "" FORCE)

option(ERF_ENABLE_RRTMGP "Enable RTE-RRTMGP Radiation" ON)
option(ERF_ENABLE_NETCDF "Enable NETCDF" ON)
option(REMORA_ENABLE_PNETCDF "Enable REMORA PNetCDF IO" ON)

add_subdirectory(Submodules/amrex amrex-build)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Submodules/amrex/Tools/CMake")
unset(AMReX_DIR CACHE)

set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

set(ERF_DIM 3 CACHE STRING "" FORCE)
set(ERF_USE_INTERNAL_AMREX OFF CACHE BOOL "" FORCE)
set(ERF_ENABLE_MPI ${ERF_REMORA_ENABLE_MPI} CACHE BOOL "" FORCE)
set(ERF_ENABLE_OPENMP ${ERF_REMORA_ENABLE_OPENMP} CACHE BOOL "" FORCE)
set(ERF_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
set(ERF_ENABLE_DOCUMENTATION OFF CACHE BOOL "" FORCE)
set(ERF_ENABLE_HDF5 OFF CACHE BOOL "" FORCE)
set(ERF_ENABLE_PARTICLES OFF CACHE BOOL "" FORCE)
set(ERF_ENABLE_FCOMPARE ${ERF_REMORA_ENABLE_FCOMPARE} CACHE BOOL "" FORCE)
set(ERF_PRECISION DOUBLE CACHE STRING "" FORCE)
set(ERF_BUILD_LIBRARY_ONLY ON CACHE BOOL "" FORCE)
set(ERF_ENABLE_MULTIBLOCK ON CACHE BOOL "" FORCE)
set(ERF_ENABLE_TINY_PROFILE ON CACHE BOOL "" FORCE)

set(REMORA_DIM 3 CACHE STRING "" FORCE)
set(REMORA_USE_INTERNAL_AMREX OFF CACHE BOOL "" FORCE)
set(REMORA_ENABLE_MPI ${ERF_REMORA_ENABLE_MPI} CACHE BOOL "" FORCE)
set(REMORA_ENABLE_OPENMP ${ERF_REMORA_ENABLE_OPENMP} CACHE BOOL "" FORCE)
set(REMORA_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
set(REMORA_ENABLE_DOCUMENTATION OFF CACHE BOOL "" FORCE)
set(REMORA_ENABLE_FCOMPARE OFF CACHE BOOL "" FORCE)
set(REMORA_PRECISION DOUBLE CACHE STRING "" FORCE)
set(REMORA_BUILD_EXECUTABLES OFF CACHE BOOL "" FORCE)
set(REMORA_BUILD_LIBRARY_ONLY ON CACHE BOOL "" FORCE)
set(REMORA_ENABLE_TINY_PROFILE ON CACHE BOOL "" FORCE)
set(REMORA_ENABLE_PNETCDF ${REMORA_ENABLE_PNETCDF} CACHE BOOL "" FORCE)

option(ERF_ENABLE_PNETCDF "Enable PNETCDF" ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Submodules/ERF/CMake")
add_subdirectory(Submodules/ERF erf-build)
add_subdirectory(Submodules/REMORA remora-build)

if(ERF_REMORA_ENABLE_FCOMPARE)
  if(WIN32)
    set(FCOMPARE_EXE
        "${CMAKE_BINARY_DIR}/amrex-build/Tools/Plotfile/*/amrex_fcompare.exe"
        CACHE STRING "Path to fcompare executable for regression tests" FORCE)
  else()
    set(FCOMPARE_EXE
        "${CMAKE_BINARY_DIR}/amrex-build/Tools/Plotfile/amrex_fcompare"
        CACHE STRING "Path to fcompare executable for regression tests" FORCE)
  endif()
  message(STATUS "Top-level fcompare executable: ${FCOMPARE_EXE}")

  # Built alongside fcompare by the same AMReX Plotfile target. Used where a
  # check needs a field's range rather than a comparison against a reference.
  if(WIN32)
    set(FEXTREMA_EXE
        "${CMAKE_BINARY_DIR}/amrex-build/Tools/Plotfile/*/amrex_fextrema.exe"
        CACHE STRING "Path to fextrema executable for regression tests" FORCE)
  else()
    set(FEXTREMA_EXE
        "${CMAKE_BINARY_DIR}/amrex-build/Tools/Plotfile/amrex_fextrema"
        CACHE STRING "Path to fextrema executable for regression tests" FORCE)
  endif()
endif()

# Bake the library paths this build linked against into the driver's RPATH.
#
# ERF already does this for its own targets — Submodules/ERF/CMake/SetRpath.cmake,
# included at Submodules/ERF/CMakeLists.txt:688 — but CMAKE_INSTALL_RPATH and
# add_link_options() are both directory-scoped. ERF's copy runs inside the scope
# opened by add_subdirectory(Submodules/ERF) above, which closes before
# erf-remora-driver is defined here, so none of it reaches this executable. It has
# to be included again in this scope, ahead of add_executable(), because
# BUILD_WITH_INSTALL_RPATH and INSTALL_RPATH are read when the target is created.
#
# Without it the driver carries no RPATH for the Cray PE, and the loader resolves
# libmpi/libmpi_gtl_cuda through /opt/cray/pe/lib64 — the *default* PE version —
# regardless of what `module load` selected. On Perlmutter that silently ignores
# the cray-mpich/9.0.1 and cudatoolkit/12.9 pins in
# Submodules/ERF/Build/machines/perlmutter_erf.profile: the run picks up the
# default cray-mpich, whose GTL wants libcudart.so.13, and dies at startup with
# "libcudart.so.13: cannot open shared object file" against the 12.9 that is
# actually loaded. The pins bind at compile time; this is what makes them bind at
# run time too.
#
# The file stands down on its own when CRAY_ADD_RPATH=yes, since the Cray wrappers
# inject RPATH themselves in that case.
include(${CMAKE_CURRENT_SOURCE_DIR}/Submodules/ERF/CMake/SetRpath.cmake)

add_executable(
  erf-remora-driver
  Source/main.cpp
  Source/CoupledRestartUtils.cpp
  Source/ERFRemoraMultiBlockContainer.cpp
  Source/RunERFStandalone.cpp
  Source/RunREMORAStandalone.cpp
)

target_include_directories(
  erf-remora-driver
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Source
    ${CMAKE_CURRENT_SOURCE_DIR}/Submodules/ERF/Source
    ${CMAKE_CURRENT_SOURCE_DIR}/Submodules/REMORA/Source
)

target_link_libraries(
  erf-remora-driver
  PRIVATE
    AMReX::amrex
    erf_api
    remora_api
)

if(TARGET Kokkos::kokkos)
  target_link_libraries(erf-remora-driver PRIVATE Kokkos::kokkos)
elseif(TARGET kokkoscore)
  target_link_libraries(erf-remora-driver PRIVATE kokkoscore)
endif()

install(
  TARGETS erf-remora-driver
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

if(ERF_REMORA_ENABLE_TESTS)
  enable_testing()

  function(erf_remora_add_restart_test _name _input _work_subdir)
    set(options)
    set(oneValueArgs
      ERF_RESTART_SUFFIX
      REMORA_RESTART_SUFFIX
      ERF_COMPARE_SUFFIX
      REMORA_COMPARE_SUFFIX
      EXCHANGE_COMPARE_SUFFIX
      ERF_PLOT_INT
      ERF_CHECK_INT
      REMORA_PLOT_INT
      REMORA_CHECK_INT
      EXCHANGE_PLOT_INT)
    set(multiValueArgs RUNTIME_ARGS)
    cmake_parse_arguments(ERF_REMORA "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    list(JOIN ERF_REMORA_RUNTIME_ARGS ";" _args_joined)
    add_test(
      NAME ${_name}
      COMMAND
        "${CMAKE_COMMAND}"
        -DERF_REMORA_DRIVER_EXE=${ERF_REMORA_DRIVER_EXE}
        -DERF_REMORA_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
        -DERF_REMORA_INPUT_FILE=${_input}
        -DFCOMPARE_EXE=${FCOMPARE_EXE}
        -DERF_REMORA_TEST_FCOMPARE_RTOL=${ERF_REMORA_TEST_FCOMPARE_RTOL}
        -DERF_REMORA_TEST_FCOMPARE_ATOL=${ERF_REMORA_TEST_FCOMPARE_ATOL}
        -DERF_REMORA_RRTMGP_DIR=${ERF_REMORA_RRTMGP_DIR}
        -DERF_REMORA_WORK_SUBDIR=${_work_subdir}
        -DERF_REMORA_ERF_RESTART_SUFFIX=${ERF_REMORA_ERF_RESTART_SUFFIX}
        -DERF_REMORA_REMORA_RESTART_SUFFIX=${ERF_REMORA_REMORA_RESTART_SUFFIX}
        -DERF_REMORA_ERF_COMPARE_SUFFIX=${ERF_REMORA_ERF_COMPARE_SUFFIX}
        -DERF_REMORA_REMORA_COMPARE_SUFFIX=${ERF_REMORA_REMORA_COMPARE_SUFFIX}
        -DERF_REMORA_EXCHANGE_COMPARE_SUFFIX=${ERF_REMORA_EXCHANGE_COMPARE_SUFFIX}
        -DERF_REMORA_ERF_PLOT_INT=${ERF_REMORA_ERF_PLOT_INT}
        -DERF_REMORA_ERF_CHECK_INT=${ERF_REMORA_ERF_CHECK_INT}
        -DERF_REMORA_REMORA_PLOT_INT=${ERF_REMORA_REMORA_PLOT_INT}
        -DERF_REMORA_REMORA_CHECK_INT=${ERF_REMORA_REMORA_CHECK_INT}
        -DERF_REMORA_EXCHANGE_PLOT_INT=${ERF_REMORA_EXCHANGE_PLOT_INT}
        "-DERF_REMORA_ARGS=${_args_joined}"
        -P "${ERF_REMORA_RESTART_TEST_SCRIPT}"
    )
    set_tests_properties(${_name}
      PROPERTIES
        LABELS "reproducer;restart"
        PROCESSORS 1
    )
  endfunction()

  function(erf_remora_add_parallel_tests _test_prefix _input _work_subdir_prefix)
    set(options)
    set(oneValueArgs
      ERF_COMPARE_SUFFIX
      REMORA_COMPARE_SUFFIX
      EXCHANGE_COMPARE_SUFFIX
      ERF_PLOT_INT
      REMORA_PLOT_INT
      EXCHANGE_PLOT_INT
      MPI_NPROCS)
    set(multiValueArgs RUNTIME_ARGS)
    cmake_parse_arguments(ERF_REMORA "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    list(JOIN ERF_REMORA_RUNTIME_ARGS ";" _args_joined)
    add_test(
      NAME ${_test_prefix}_np${ERF_REMORA_MPI_NPROCS}_smoke
      COMMAND
        "${CMAKE_COMMAND}"
        -DERF_REMORA_DRIVER_EXE=${ERF_REMORA_DRIVER_EXE}
        -DERF_REMORA_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
        -DERF_REMORA_INPUT_FILE=${_input}
        -DERF_REMORA_RRTMGP_DIR=${ERF_REMORA_RRTMGP_DIR}
        -DMPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}
        -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
        -DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}
        -DERF_REMORA_WORK_SUBDIR=${_work_subdir_prefix}_smoke
        -DERF_REMORA_MPI_NPROCS=${ERF_REMORA_MPI_NPROCS}
        -DERF_REMORA_ERF_COMPARE_SUFFIX=${ERF_REMORA_ERF_COMPARE_SUFFIX}
        -DERF_REMORA_REMORA_COMPARE_SUFFIX=${ERF_REMORA_REMORA_COMPARE_SUFFIX}
        -DERF_REMORA_EXCHANGE_COMPARE_SUFFIX=${ERF_REMORA_EXCHANGE_COMPARE_SUFFIX}
        -DERF_REMORA_ERF_PLOT_INT=${ERF_REMORA_ERF_PLOT_INT}
        -DERF_REMORA_REMORA_PLOT_INT=${ERF_REMORA_REMORA_PLOT_INT}
        -DERF_REMORA_EXCHANGE_PLOT_INT=${ERF_REMORA_EXCHANGE_PLOT_INT}
        "-DERF_REMORA_ARGS=${_args_joined}"
        -P "${ERF_REMORA_MPI_SMOKE_TEST_SCRIPT}"
    )
    set_tests_properties(${_test_prefix}_np${ERF_REMORA_MPI_NPROCS}_smoke
      PROPERTIES
        LABELS "parallel"
        PROCESSORS ${ERF_REMORA_MPI_NPROCS}
    )

    add_test(
      NAME ${_test_prefix}_np${ERF_REMORA_MPI_NPROCS}_compare_to_np1
      COMMAND
        "${CMAKE_COMMAND}"
        -DERF_REMORA_DRIVER_EXE=${ERF_REMORA_DRIVER_EXE}
        -DERF_REMORA_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
        -DERF_REMORA_INPUT_FILE=${_input}
        -DFCOMPARE_EXE=${FCOMPARE_EXE}
        -DERF_REMORA_TEST_FCOMPARE_RTOL=${ERF_REMORA_TEST_FCOMPARE_RTOL}
        -DERF_REMORA_TEST_FCOMPARE_ATOL=${ERF_REMORA_TEST_FCOMPARE_ATOL}
        -DERF_REMORA_RRTMGP_DIR=${ERF_REMORA_RRTMGP_DIR}
        -DMPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}
        -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
        -DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}
        -DERF_REMORA_WORK_SUBDIR=${_work_subdir_prefix}_compare
        -DERF_REMORA_MPI_NPROCS=${ERF_REMORA_MPI_NPROCS}
        -DERF_REMORA_ERF_COMPARE_SUFFIX=${ERF_REMORA_ERF_COMPARE_SUFFIX}
        -DERF_REMORA_REMORA_COMPARE_SUFFIX=${ERF_REMORA_REMORA_COMPARE_SUFFIX}
        -DERF_REMORA_EXCHANGE_COMPARE_SUFFIX=${ERF_REMORA_EXCHANGE_COMPARE_SUFFIX}
        -DERF_REMORA_ERF_PLOT_INT=${ERF_REMORA_ERF_PLOT_INT}
        -DERF_REMORA_REMORA_PLOT_INT=${ERF_REMORA_REMORA_PLOT_INT}
        -DERF_REMORA_EXCHANGE_PLOT_INT=${ERF_REMORA_EXCHANGE_PLOT_INT}
        "-DERF_REMORA_ARGS=${_args_joined}"
        -P "${ERF_REMORA_MPI_COMPARE_TEST_SCRIPT}"
    )
    set_tests_properties(${_test_prefix}_np${ERF_REMORA_MPI_NPROCS}_compare_to_np1
      PROPERTIES
        LABELS "parallel"
        PROCESSORS ${ERF_REMORA_MPI_NPROCS}
    )
  endfunction()

  set(ERF_REMORA_DRIVER_EXE "$<TARGET_FILE:erf-remora-driver>")
  set(ERF_REMORA_TEST_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/cmake/run_reproducer.cmake")
  set(ERF_REMORA_RESTART_TEST_SCRIPT
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/run_restart_reproducer.cmake")
  set(ERF_REMORA_MPI_SMOKE_TEST_SCRIPT
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/run_mpi_upwelling_smoke.cmake")
  set(ERF_REMORA_MPI_COMPARE_TEST_SCRIPT
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/run_mpi_upwelling_compare.cmake")
  set(ERF_REMORA_SST_COVERAGE_TEST_SCRIPT
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/run_coupled_sst_coverage.cmake")
  set(ERF_REMORA_LANE_PROVENANCE_TEST_SCRIPT
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/run_lane_provenance.cmake")

  add_test(
    NAME compile_erf_remora_driver
    COMMAND
      "${CMAKE_COMMAND}"
      --build "${CMAKE_BINARY_DIR}"
      --target erf-remora-driver
  )
  set_tests_properties(compile_erf_remora_driver
    PROPERTIES
      LABELS "compile"
  )

  # Coupled SST must reach ERF only where the ocean grid covers it, and ERF must
  # keep its own lower-boundary value everywhere else. See issue #13.
  add_test(
    NAME reproducer_coupled_sst_coverage
    COMMAND
      "${CMAKE_COMMAND}"
      -DERF_REMORA_DRIVER_EXE=${ERF_REMORA_DRIVER_EXE}
      -DERF_REMORA_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
      -DERF_REMORA_INPUT_FILE=Exec/Upwelling/inputs_twoway
      -DFEXTREMA_EXE=${FEXTREMA_EXE}
      -DERF_REMORA_RRTMGP_DIR=${ERF_REMORA_RRTMGP_DIR}
      -DERF_REMORA_MAX_STEP=2
      -P "${ERF_REMORA_SST_COVERAGE_TEST_SCRIPT}"
  )
  set_tests_properties(reproducer_coupled_sst_coverage
    PROPERTIES
      LABELS "reproducer;coupling"
      PROCESSORS 1
  )

  # Every atmos->ocean forcing lane must have exactly one source, and the driver
  # must not hand REMORA a lane ERF never wrote. See issue #4.
  add_test(
    NAME reproducer_lane_provenance
    COMMAND
      "${CMAKE_COMMAND}"
      -DERF_REMORA_DRIVER_EXE=${ERF_REMORA_DRIVER_EXE}
      -DERF_REMORA_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
      -DERF_REMORA_INPUT_FILE=Exec/Upwelling/inputs_oneway
      -DERF_REMORA_RRTMGP_DIR=${ERF_REMORA_RRTMGP_DIR}
      -DERF_REMORA_MAX_STEP=1
      -P "${ERF_REMORA_LANE_PROVENANCE_TEST_SCRIPT}"
  )
  set_tests_properties(reproducer_lane_provenance
    PROPERTIES
      LABELS "reproducer;coupling"
      PROCESSORS 1
  )

  # The reproducibility tests below ask whether the coupled run gives the same
  # answer twice: once across rank counts, once across a checkpoint/restart.
  # Both require every other source of run-to-run variation to be absent, so
  # they run the deck's own PBL scheme rather than overriding it, and the deck
  # scheme has to be one that holds still. erf.pbl_type=MRF does. Measured on
  # the Upwelling deck: cold vs restarted and np1 vs np4 both fcompare to
  # exactly zero absolute error on every ERF, REMORA and exchange field.
  #
  # These tests used to force erf.pbl_type=None, because the deck's previous
  # scheme, MYNN25, is neither rank- nor restart-reproducible:
  #
  #   - across ranks, np1 vs npN is identical at step 0, diverges to 7e-06 in
  #     Kmv by step 10 and 2.4e-02 (0.37% relative) by step 20, saturating at
  #     the same magnitude for np2, np4 and np8
  #   - across restart, cold vs restarted diverges to 1.3 in Kmv (19% relative)
  #
  # That is still an open ERF issue and nothing here covers it; it is simply no
  # longer in the path of these tests. The override also had to supply constant
  # molecular diffusion, because MOST asserts has_diff (ERF.cpp) and disabling
  # the PBL scheme removed the only diffusion the deck had. MRF satisfies
  # has_diff on its own, so the deck's erf.molec_diff_type=None stands.

  if(FCOMPARE_EXE)
    set(_erf_remora_upwelling_cases
      "oneway|Exec/Upwelling/inputs_oneway"
      "twoway|Exec/Upwelling/inputs_twoway")
    set(_erf_remora_atm2ocn_modes state flux)
    foreach(_case_entry IN LISTS _erf_remora_upwelling_cases)
      string(REPLACE "|" ";" _case_fields "${_case_entry}")
      list(GET _case_fields 0 _case_name)
      list(GET _case_fields 1 _case_input)
      foreach(_atm2ocn_mode IN LISTS _erf_remora_atm2ocn_modes)
        erf_remora_add_restart_test(
          reproducer_restart_${_case_name}_${_atm2ocn_mode}_matched
          ${_case_input}
          restart_${_case_name}_${_atm2ocn_mode}_matched
          ERF_RESTART_SUFFIX 00001
          REMORA_RESTART_SUFFIX 00001
          ERF_COMPARE_SUFFIX 00002
          REMORA_COMPARE_SUFFIX 00002
          EXCHANGE_COMPARE_SUFFIX 00001
          ERF_PLOT_INT 1
          ERF_CHECK_INT 1
          REMORA_PLOT_INT 1
          REMORA_CHECK_INT 1
          EXCHANGE_PLOT_INT 1
          RUNTIME_ARGS
            "driver.atm2ocn_mode=${_atm2ocn_mode}"
            "erf.fixed_dt=1"
            "remora.fixed_dt=1"
            "remora.ndtfast=1"
            "max_step=2"
        )

        erf_remora_add_restart_test(
          reproducer_restart_${_case_name}_${_atm2ocn_mode}_substep
          ${_case_input}
          restart_${_case_name}_${_atm2ocn_mode}_substep
          ERF_RESTART_SUFFIX 00010
          REMORA_RESTART_SUFFIX 00001
          ERF_COMPARE_SUFFIX 00020
          REMORA_COMPARE_SUFFIX 00002
          EXCHANGE_COMPARE_SUFFIX 00001
          ERF_PLOT_INT 10
          ERF_CHECK_INT 10
          REMORA_PLOT_INT 1
          REMORA_CHECK_INT 1
          EXCHANGE_PLOT_INT 1
          RUNTIME_ARGS
            "driver.atm2ocn_mode=${_atm2ocn_mode}"
            "erf.fixed_dt=1"
            "remora.fixed_dt=10"
            "remora.ndtfast=1"
            "max_step=2"
        )

        if(_atm2ocn_mode STREQUAL "state")
          set(_a2o_average_variant average_state_lanes)
          set(_a2o_average_arg "driver.time_average_atm2ocn_state=1")
        else()
          set(_a2o_average_variant average_flux_lanes)
          set(_a2o_average_arg "driver.time_average_atm2ocn_flux=1")
        endif()
        erf_remora_add_restart_test(
          reproducer_restart_${_case_name}_${_atm2ocn_mode}_${_a2o_average_variant}_substep
          ${_case_input}
          restart_${_case_name}_${_atm2ocn_mode}_${_a2o_average_variant}_substep
          ERF_RESTART_SUFFIX 00010
          REMORA_RESTART_SUFFIX 00001
          ERF_COMPARE_SUFFIX 00020
          REMORA_COMPARE_SUFFIX 00002
          EXCHANGE_COMPARE_SUFFIX 00001
          ERF_PLOT_INT 10
          ERF_CHECK_INT 10
          REMORA_PLOT_INT 1
          REMORA_CHECK_INT 1
          EXCHANGE_PLOT_INT 1
          RUNTIME_ARGS
            "driver.atm2ocn_mode=${_atm2ocn_mode}"
            "${_a2o_average_arg}"
            "erf.fixed_dt=1"
            "remora.fixed_dt=10"
            "remora.ndtfast=1"
            "max_step=2"
        )
      endforeach()
    endforeach()
  else()
    message(STATUS
      "Skipping top-level restart reproducer test registration because "
      "FCOMPARE_EXE is unset")
  endif()

  if(ERF_REMORA_ENABLE_MPI
     AND MPIEXEC_EXECUTABLE
     AND FCOMPARE_EXE)
    set(_erf_remora_upwelling_cases
      "oneway|Exec/Upwelling/inputs_oneway"
      "twoway|Exec/Upwelling/inputs_twoway")
    set(_erf_remora_atm2ocn_modes state flux)

    # See the comment above the restart tests for why these run the deck's
    # own PBL scheme rather than overriding it.
    foreach(_case_entry IN LISTS _erf_remora_upwelling_cases)
      string(REPLACE "|" ";" _case_fields "${_case_entry}")
      list(GET _case_fields 0 _case_name)
      list(GET _case_fields 1 _case_input)
      foreach(_atm2ocn_mode IN LISTS _erf_remora_atm2ocn_modes)
        foreach(_mpi_nprocs 2 4 6 8)
          erf_remora_add_parallel_tests(
            parallel_upwelling_${_case_name}_${_atm2ocn_mode}_matched
            ${_case_input}
            mpi_${_case_name}_${_atm2ocn_mode}_matched_np${_mpi_nprocs}
            MPI_NPROCS ${_mpi_nprocs}
            ERF_COMPARE_SUFFIX 00002
            REMORA_COMPARE_SUFFIX 00002
            EXCHANGE_COMPARE_SUFFIX 00001
            ERF_PLOT_INT 1
            REMORA_PLOT_INT 1
            EXCHANGE_PLOT_INT 1
            RUNTIME_ARGS
              "driver.atm2ocn_mode=${_atm2ocn_mode}"
              "erf.fixed_dt=1"
              "remora.fixed_dt=1"
              "remora.ndtfast=1"
              "max_step=2"
          )

          erf_remora_add_parallel_tests(
            parallel_upwelling_${_case_name}_${_atm2ocn_mode}_substep
            ${_case_input}
            mpi_${_case_name}_${_atm2ocn_mode}_substep_np${_mpi_nprocs}
            MPI_NPROCS ${_mpi_nprocs}
            ERF_COMPARE_SUFFIX 00020
            REMORA_COMPARE_SUFFIX 00002
            EXCHANGE_COMPARE_SUFFIX 00001
            ERF_PLOT_INT 10
            REMORA_PLOT_INT 1
            EXCHANGE_PLOT_INT 1
            RUNTIME_ARGS
              "driver.atm2ocn_mode=${_atm2ocn_mode}"
              "erf.fixed_dt=1"
              "remora.fixed_dt=10"
              "remora.ndtfast=1"
              "max_step=2"
          )

          if(_atm2ocn_mode STREQUAL "state")
            set(_a2o_average_variant average_state_lanes)
            set(_a2o_average_arg "driver.time_average_atm2ocn_state=1")
          else()
            set(_a2o_average_variant average_flux_lanes)
            set(_a2o_average_arg "driver.time_average_atm2ocn_flux=1")
          endif()
          erf_remora_add_parallel_tests(
            parallel_upwelling_${_case_name}_${_atm2ocn_mode}_${_a2o_average_variant}_substep
            ${_case_input}
            mpi_${_case_name}_${_atm2ocn_mode}_${_a2o_average_variant}_substep_np${_mpi_nprocs}
            MPI_NPROCS ${_mpi_nprocs}
            ERF_COMPARE_SUFFIX 00020
            REMORA_COMPARE_SUFFIX 00002
            EXCHANGE_COMPARE_SUFFIX 00001
            ERF_PLOT_INT 10
            REMORA_PLOT_INT 1
            EXCHANGE_PLOT_INT 1
            RUNTIME_ARGS
              "driver.atm2ocn_mode=${_atm2ocn_mode}"
              "${_a2o_average_arg}"
              "erf.fixed_dt=1"
              "remora.fixed_dt=10"
              "remora.ndtfast=1"
              "max_step=2"
          )
        endforeach()
      endforeach()
    endforeach()
  else()
    message(STATUS
      "Skipping top-level MPI upwelling test registration because "
      "ERF_REMORA_ENABLE_MPI is OFF or MPIEXEC_EXECUTABLE or FCOMPARE_EXE is unset")
  endif()
endif()
