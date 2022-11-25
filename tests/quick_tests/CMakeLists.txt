## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2021 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#
# A minimalistic set of tests:
#
enable_testing()

include_directories(
  ${CMAKE_BINARY_DIR}/include/
  ${CMAKE_SOURCE_DIR}/include/
  ${DEAL_II_BUNDLED_INCLUDE_DIRS}
  ${DEAL_II_INCLUDE_DIRS}
  )

# Use the first available build type (this prefers debug mode if available):
list(GET DEAL_II_BUILD_TYPES 0 _mybuild)
message(STATUS "Setting up quick_tests in ${_mybuild} mode")

set(ALL_TESTS) # clean variable

# define a macro to set up a quick test:
macro(make_quicktest test_basename build_name mpi_run)
  string(TOLOWER ${build_name} _build_lowercase)
  set(_target ${test_basename}.${_build_lowercase})
  list(APPEND ALL_TESTS "${_target}")
  add_executable(${_target} EXCLUDE_FROM_ALL ${test_basename}.cc)
  deal_ii_insource_setup_target(${_target} ${build_name})

  if("${mpi_run}" STREQUAL "")
    set(_command ${_target})
  else()
    if(CMAKE_SYSTEM_NAME MATCHES "Windows")
      set(_command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${mpi_run} ${MPIEXEC_PREFLAGS} ${_target})
    else()
      set(_command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${mpi_run} ${MPIEXEC_PREFLAGS} ./${_target})
    endif()
  endif()
  add_custom_target(${_target}.run
    DEPENDS ${_target}
    COMMAND
      ${_command} > ${_target}-OK 2>&1
      ||(echo "${_target}: RUN failed. Output:"
         && cat ${_target}-OK
         && rm ${_target}-OK
         && exit 1)
    COMMAND echo "${_target}: PASSED."
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

  # this is a hack to make sure the -OK file is deleted
  # even if compilation fails.
  add_custom_target(reset-${_target}-OK
    COMMAND ${CMAKE_COMMAND} -E remove -f ${_target}-OK
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
  add_dependencies(${_target} reset-${_target}-OK)

  add_test(NAME ${_target}
    COMMAND ${CMAKE_COMMAND} -DTRGT=${_target}.run -DTEST=${_target}
      -DBINARY_DIR=${CMAKE_BINARY_DIR}
      -P ${CMAKE_SOURCE_DIR}/cmake/scripts/run_test.cmake
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
  set_tests_properties(${_target} PROPERTIES
    LABEL "sanity checks"
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  # ^^^ Allow oversubscription for MPI (needed for Openmpi@3.0)
endmacro()


# Simple assembly/solver test. This makes sure we can compile and link correctly
# in debug and release.
foreach(_build ${DEAL_II_BUILD_TYPES})
  make_quicktest("step" ${_build} "")
endforeach()

# Test whether thread affinity is well behaved
make_quicktest("affinity" ${_mybuild} "")

# Test if MPI is configured correctly
if (DEAL_II_WITH_MPI)
  make_quicktest("mpi" ${_mybuild} 2)
endif()

# Test if TBB works correctly
if (DEAL_II_WITH_TBB)
  make_quicktest("tbb" ${_mybuild} "")
endif()

# test taskflow
if (DEAL_II_WITH_TASKFLOW)
  make_quicktest("taskflow" ${_mybuild} "")
endif()

# Test p4est. This test exposes a bug in OpenMPI 1.3 and 1.4
# Update to OpenMPI 1.5 or newer.
if (DEAL_II_WITH_P4EST)
  make_quicktest("p4est" ${_mybuild} 10)
endif()

# Test petsc
if (DEAL_II_WITH_PETSC)
  make_quicktest("step-petsc" ${_mybuild} "")
endif()

# Test slepc
if (DEAL_II_WITH_PETSC AND DEAL_II_WITH_SLEPC)
  make_quicktest("step-slepc" ${_mybuild} "")
endif()

# Test trilinos
if (DEAL_II_WITH_TRILINOS)
  make_quicktest("step-trilinos" ${_mybuild} "")
endif()

# Test metis
if (DEAL_II_WITH_METIS AND DEAL_II_WITH_MPI)
  make_quicktest("step-metis" ${_mybuild} 2)
endif()

# Test LAPACK
if (DEAL_II_WITH_LAPACK)
  make_quicktest("lapack" ${_mybuild} "")
endif()

# Test Umfpack
if (DEAL_II_WITH_UMFPACK)
  make_quicktest("umfpack" ${_mybuild} "")
endif()

# Test GSL
if (DEAL_II_WITH_GSL)
  make_quicktest("gsl" ${_mybuild} "")
endif()

# Test HDF5
if (DEAL_II_WITH_HDF5)
  make_quicktest("hdf5" ${_mybuild} "")
endif()

# Test Arpack
if (DEAL_II_WITH_ARPACK AND DEAL_II_WITH_UMFPACK)
  make_quicktest("arpack" ${_mybuild} "")
endif()

# Test Adol-C
if (DEAL_II_WITH_ADOLC)
  make_quicktest("adolc" ${_mybuild} "")
endif()

# Test SUNDIALS
if (DEAL_II_WITH_SUNDIALS)
  make_quicktest("sundials-ida" ${_mybuild} "")
endif()

# Test Assimp
if (DEAL_II_WITH_ASSIMP)
  make_quicktest("assimp" ${_mybuild} "")
endif()

# Test CUDA
if (DEAL_II_WITH_CUDA)
  make_quicktest("cuda" ${_mybuild} "")
endif()

# Test ScaLAPACK
if (DEAL_II_WITH_SCALAPACK)
  make_quicktest("scalapack" ${_mybuild} 5)
endif()

# Test Zlib and Boost
if (DEAL_II_WITH_ZLIB)
  make_quicktest("boost_zlib" ${_mybuild} "")
endif()

# Test GMSH
if (DEAL_II_WITH_GMSH)
  make_quicktest("gmsh" ${_mybuild} "")
endif()

# Test SymEngine
if (DEAL_II_WITH_SYMENGINE)
  make_quicktest("symengine" ${_mybuild} "")
endif()

# Test Kokkos
make_quicktest("kokkos" ${_mybuild} "")

# A custom test target:
add_custom_target(test
  COMMAND ${CMAKE_COMMAND} -D ALL_TESTS="${ALL_TESTS}" -DCMAKE_BUILD_TYPE=${_mybuild} -P ${CMAKE_CURRENT_SOURCE_DIR}/run.cmake
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Running quicktests..."
  )

message(STATUS "Setting up quick_tests in ${_mybuild} mode - Done")
