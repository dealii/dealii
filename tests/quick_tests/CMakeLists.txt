## ---------------------------------------------------------------------
##
## Copyright (C) 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
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
  ADD_EXECUTABLE(${_target} EXCLUDE_FROM_ALL ${test_basename}.cc)
  DEAL_II_INSOURCE_SETUP_TARGET(${_target} ${build_name})

  IF("${mpi_run}" STREQUAL "")
    SET(_command ./${_target})
  ELSE()
    SET(_command ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${mpi_run} ${MPIEXEC_PREFLAGS} ./${_target})
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
      -DDEAL_II_BINARY_DIR=${CMAKE_BINARY_DIR}
      -P ${CMAKE_SOURCE_DIR}/tests/run_test.cmake
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
  SET_TESTS_PROPERTIES(${_target} PROPERTIES LABEL "sanity checks")
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

# A custom test target:
ADD_CUSTOM_TARGET(test
  COMMAND ${CMAKE_COMMAND} -D ALL_TESTS="${ALL_TESTS}" -P ${CMAKE_CURRENT_SOURCE_DIR}/run.cmake
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

MESSAGE(STATUS "Setting up quick_tests in ${_mybuild} mode - Done")
