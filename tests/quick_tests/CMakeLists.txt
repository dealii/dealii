## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2019 by the deal.II authors
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
ENABLE_TESTING()

INCLUDE_DIRECTORIES(
  ${CMAKE_BINARY_DIR}/include/
  ${CMAKE_SOURCE_DIR}/include/
  ${DEAL_II_BUNDLED_INCLUDE_DIRS}
  ${DEAL_II_INCLUDE_DIRS}
  )

# Use the first available build type (this prefers debug mode if available):
LIST(GET DEAL_II_BUILD_TYPES 0 _mybuild)
MESSAGE(STATUS "Setting up quick_tests in ${_mybuild} mode")

SET(ALL_TESTS) # clean variable

# define a macro to set up a quick test:
MACRO(make_quicktest test_basename build_name mpi_run)
  STRING(TOLOWER ${build_name} _build_lowercase)
  SET(_target ${test_basename}.${_build_lowercase})
  LIST(APPEND ALL_TESTS "${_target}")
  IF(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${test_basename}.cc")
    ADD_EXECUTABLE(${_target} EXCLUDE_FROM_ALL ${test_basename}.cc)
  ELSE()
    ADD_EXECUTABLE(${_target} EXCLUDE_FROM_ALL ${test_basename}.cu)
  ENDIF()
  DEAL_II_INSOURCE_SETUP_TARGET(${_target} ${build_name})

  IF("${mpi_run}" STREQUAL "")
    SET(_command ${_target})
  ELSE()
    IF(CMAKE_SYSTEM_NAME MATCHES "Windows")
      SET(_command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${mpi_run} ${MPIEXEC_PREFLAGS} ${_target})
    ELSE()
      SET(_command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${mpi_run} ${MPIEXEC_PREFLAGS} ./${_target})
    ENDIF()
  ENDIF()
  ADD_CUSTOM_TARGET(${_target}.run
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
  ADD_CUSTOM_TARGET(kill-${_target}-OK
        COMMAND rm -f ${_target}-OK
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
  ADD_DEPENDENCIES(${_target} kill-${_target}-OK)

  ADD_TEST(NAME ${_target}
    COMMAND ${CMAKE_COMMAND} -DTRGT=${_target}.run -DTEST=${_target}
      -DBINARY_DIR=${CMAKE_BINARY_DIR}
      -P ${CMAKE_SOURCE_DIR}/cmake/scripts/run_test.cmake
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
  SET_TESTS_PROPERTIES(${_target} PROPERTIES
    LABEL "sanity checks"
    ENVIRONMENT "OMPI_MCA_rmaps_base_oversubscribe=1")
  # ^^^ Allow oversubscription for MPI (needed for Openmpi@3.0)
ENDMACRO()


# Simple assembly/solver test. This makes sure we can compile and link correctly
# in debug and release.
FOREACH(_build ${DEAL_II_BUILD_TYPES})
  make_quicktest("step" ${_build} "")
ENDFOREACH()

# Test whether thread affinity is well behaved
IF (DEAL_II_WITH_THREADS)
make_quicktest("affinity" ${_mybuild} "")
ENDIF()

# Test if MPI is configured correctly
IF (DEAL_II_WITH_MPI)
  make_quicktest("mpi" ${_mybuild} 2)
ENDIF()

# Test if TBB works correctly
IF (DEAL_II_WITH_THREADS)
  make_quicktest("tbb" ${_mybuild} "")
ENDIF()

# Test p4est. This test exposes a bug in OpenMPI 1.3 and 1.4
# Update to OpenMPI 1.5 or newer.
IF (DEAL_II_WITH_P4EST)
  make_quicktest("p4est" ${_mybuild} 10)
ENDIF()

# Test petsc
IF (DEAL_II_WITH_PETSC)
  make_quicktest("step-petsc" ${_mybuild} "")
ENDIF()

# Test slepc
IF (DEAL_II_WITH_PETSC AND DEAL_II_WITH_SLEPC)
  make_quicktest("step-slepc" ${_mybuild} "")
ENDIF()

# Test trilinos
IF (DEAL_II_WITH_TRILINOS)
  make_quicktest("step-trilinos" ${_mybuild} "")
ENDIF()

# Test metis
IF (DEAL_II_WITH_METIS AND DEAL_II_WITH_MPI)
  make_quicktest("step-metis" ${_mybuild} 2)
ENDIF()

# Test LAPACK
IF (DEAL_II_WITH_LAPACK)
  make_quicktest("lapack" ${_mybuild} "")
ENDIF()

# Test Umfpack
IF (DEAL_II_WITH_UMFPACK)
  make_quicktest("umfpack" ${_mybuild} "")
ENDIF()

# Test GSL
IF (DEAL_II_WITH_GSL)
  make_quicktest("gsl" ${_mybuild} "")
ENDIF()

# Test HDF5
IF (DEAL_II_WITH_HDF5)
  make_quicktest("hdf5" ${_mybuild} "")
ENDIF()

# Test Arpack
IF (DEAL_II_WITH_ARPACK AND DEAL_II_WITH_UMFPACK)
  make_quicktest("arpack" ${_mybuild} "")
ENDIF()

# Test Adol-C
IF (DEAL_II_WITH_ADOLC)
  make_quicktest("adolc" ${_mybuild} "")
ENDIF()

# Test SUNDIALS
IF (DEAL_II_WITH_SUNDIALS)
  make_quicktest("sundials-ida" ${_mybuild} "")
ENDIF()

# Test Assimp
IF (DEAL_II_WITH_ASSIMP)
  make_quicktest("assimp" ${_mybuild} "")
ENDIF()

# Test CUDA
IF (DEAL_II_WITH_CUDA)
  make_quicktest("cuda" ${_mybuild} "")
ENDIF()

# Test ScaLAPACK
IF (DEAL_II_WITH_SCALAPACK)
  make_quicktest("scalapack" ${_mybuild} 5)
ENDIF()

# Test Zlib and Boost
IF (DEAL_II_WITH_ZLIB)
  make_quicktest("boost_zlib" ${_mybuild} "")
ENDIF()

# Test GMSH
IF (DEAL_II_WITH_GMSH)
  make_quicktest("gmsh" ${_mybuild} "")
ENDIF()

# Test SymEngine
IF (DEAL_II_WITH_SYMENGINE)
  make_quicktest("symengine" ${_mybuild} "")
ENDIF()

# A custom test target:
ADD_CUSTOM_TARGET(test
  COMMAND ${CMAKE_COMMAND} -D ALL_TESTS="${ALL_TESTS}" -DCMAKE_BUILD_TYPE=${_mybuild} -P ${CMAKE_CURRENT_SOURCE_DIR}/run.cmake
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Running quicktests..."
  )

MESSAGE(STATUS "Setting up quick_tests in ${_mybuild} mode - Done")
