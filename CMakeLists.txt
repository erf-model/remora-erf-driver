cmake_minimum_required(VERSION 3.24)

project(erf-remora-driver CXX C)

include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(ERF_REMORA_ENABLE_MPI "Enable MPI across AMReX/ERF/REMORA" ON)
option(ERF_REMORA_ENABLE_OPENMP "Enable OpenMP across AMReX/ERF/REMORA" OFF)
option(ERF_REMORA_ENABLE_TESTS "Enable top-level tests" OFF)
message(STATUS "ERF-REMORA symbol policy: prefixed (fixed)")

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
set(AMReX_PLOTFILE_TOOLS OFF CACHE BOOL "" FORCE)
set(AMReX_INSTALL OFF CACHE BOOL "" FORCE)
set(AMReX_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
set(AMReX_ASSERTIONS ON CACHE BOOL "" FORCE)
set(AMReX_FPE OFF CACHE BOOL "" FORCE)
set(AMReX_GPU_BACKEND NONE CACHE STRING "" FORCE)

option(ERF_ENABLE_RRTMGP "Enable RTE-RRTMGP Radiation" ON)
option(ERF_ENABLE_NETCDF "Enable NETCDF" ON)

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
set(ERF_ENABLE_FCOMPARE OFF CACHE BOOL "" FORCE)
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
set(REMORA_LIBRARY_PROBLEM "Upwelling_Coupling" CACHE STRING "" FORCE)
set(REMORA_ENABLE_TINY_PROFILE ON CACHE BOOL "" FORCE)

option(ERF_ENABLE_PNETCDF "Enable PNETCDF" ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Submodules/ERF/CMake")
add_subdirectory(Submodules/ERF erf-build)
add_subdirectory(Submodules/REMORA remora-build)

add_executable(
  erf-remora-driver
  Source/main.cpp
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

install(
  TARGETS erf-remora-driver
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

if(ERF_REMORA_ENABLE_TESTS)
  enable_testing()

  set(ERF_REMORA_DRIVER_EXE "$<TARGET_FILE:erf-remora-driver>")
  set(ERF_REMORA_TEST_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/cmake/run_reproducer.cmake")

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

  add_test(
    NAME reproducer_coupling_infrastructure_check
    COMMAND
      "${CMAKE_COMMAND}"
      -DERF_REMORA_DRIVER_EXE=${ERF_REMORA_DRIVER_EXE}
      -DERF_REMORA_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
      -DERF_REMORA_INPUT_FILE=Exec/CouplingReproducers/inputs_coupling_infrastructure_check
      -DERF_REMORA_ARGS=
      -P "${ERF_REMORA_TEST_SCRIPT}"
  )
  set_tests_properties(reproducer_coupling_infrastructure_check
    PROPERTIES
      LABELS "reproducer"
  )

  add_test(
    NAME reproducer_one_step_solver_initialization_check
    COMMAND
      "${CMAKE_COMMAND}"
      -DERF_REMORA_DRIVER_EXE=${ERF_REMORA_DRIVER_EXE}
      -DERF_REMORA_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
      -DERF_REMORA_INPUT_FILE=Exec/CouplingReproducers/inputs_one_step_solver_initialization_check
      -DERF_REMORA_ARGS=max_step=1
      -P "${ERF_REMORA_TEST_SCRIPT}"
  )
  set_tests_properties(reproducer_one_step_solver_initialization_check
    PROPERTIES
      LABELS "reproducer"
  )

  add_test(
    NAME reproducer_erftoremora_coupleToERF
    COMMAND
      "${CMAKE_COMMAND}"
      -DERF_REMORA_DRIVER_EXE=${ERF_REMORA_DRIVER_EXE}
      -DERF_REMORA_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
      -DERF_REMORA_INPUT_FILE=Exec/CouplingReproducers/inputs_ERFtoREMORA_coupleToERF
      -DERF_REMORA_ARGS=max_step=10
      -P "${ERF_REMORA_TEST_SCRIPT}"
  )
  set_tests_properties(reproducer_erftoremora_coupleToERF
    PROPERTIES
      LABELS "reproducer"
  )

  add_test(
    NAME reproducer_two_way_coupleToERF
    COMMAND
      "${CMAKE_COMMAND}"
      -DERF_REMORA_DRIVER_EXE=${ERF_REMORA_DRIVER_EXE}
      -DERF_REMORA_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
      -DERF_REMORA_INPUT_FILE=Exec/CouplingReproducers/inputs_two_way_coupleToERF
      "-DERF_REMORA_ARGS=max_step=600;remora.v=1;erf.v=1;amr.v=1;amrex.tiny_profiler=1"
      -P "${ERF_REMORA_TEST_SCRIPT}"
  )
  set_tests_properties(reproducer_two_way_coupleToERF
    PROPERTIES
      LABELS "reproducer"
  )

  add_test(
    NAME reproducer_two_way_ABL
    COMMAND
      "${CMAKE_COMMAND}"
      -DERF_REMORA_DRIVER_EXE=${ERF_REMORA_DRIVER_EXE}
      -DERF_REMORA_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
      -DERF_REMORA_INPUT_FILE=Exec/CouplingReproducers/inputs_two_way_ABL
      -DERF_REMORA_ARGS=
      -P "${ERF_REMORA_TEST_SCRIPT}"
  )
  set_tests_properties(reproducer_two_way_ABL
    PROPERTIES
      LABELS "reproducer"
  )

  add_test(
    NAME reproducer_two_way_couette
    COMMAND
      "${CMAKE_COMMAND}"
      -DERF_REMORA_DRIVER_EXE=${ERF_REMORA_DRIVER_EXE}
      -DERF_REMORA_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
      -DERF_REMORA_INPUT_FILE=Exec/CouplingReproducers/inputs_two_way_couette
      -DERF_REMORA_ARGS=
      -P "${ERF_REMORA_TEST_SCRIPT}"
  )
  set_tests_properties(reproducer_two_way_couette
    PROPERTIES
      LABELS "reproducer"
  )
endif()
