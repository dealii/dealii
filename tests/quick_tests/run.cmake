## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2016 by the deal.II authors
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

# This file is run when "make test" is executed by the user and is
# responsible for running the tests and printing some helpful
# error messages.
include(ProcessorCount)
PROCESSORCOUNT(_n_processors)
IF(_n_processors EQUAL 0)
  SET(_n_processors "1")
ENDIF()


# Windows quick tests have a race condition, so disable compiling/running
# tests in parallel. This avoid errors like:
#
# error MSB3491: Could not write lines to file
# "obj_boost_system_debug.dir\Debug\obj_boos.4A356C5C.tlog\obj_boost_system_debug.lastbuildstate". The
# process cannot access the file '...' because it is being used by another
# process.
IF(CMAKE_HOST_SYSTEM_NAME MATCHES "Windows")
  SET(_n_processors "1")
ENDIF()

SEPARATE_ARGUMENTS(ALL_TESTS)

EXECUTE_PROCESS(COMMAND ${CMAKE_CTEST_COMMAND} -j${_n_processors}
  -C ${CMAKE_BUILD_TYPE}
  --force-new-ctest-process
  --output-on-failure
  -O quicktests.log
  RESULT_VARIABLE res_var)

IF(NOT "${res_var}" STREQUAL "0")
  MESSAGE("
***************************************************************************
**                 Error: Some of the quick tests failed.                **
***************************************************************************

Please scroll up or check the file tests/quick_tests/quicktests.log for the
error messages. If you are unable to fix the problems, see the FAQ or write
to the mailing list linked at http://www.dealii.org\n"
    )

  FOREACH(test ${ALL_TESTS})  
    IF (${test} MATCHES "^affinity" AND NOT EXISTS ${test}-OK)
      MESSAGE("
The affinity test can fail when you are linking in a library like BLAS
which uses OpenMP. Even without calling any BLAS functions, OpenMP messes
with the thread affinity which causes TBB to run single-threaded only. You
can fix this by exporting OMP_NUM_THREADS=1. Also see GOMP_CPU_AFFINITY 
and OMP_PROC_BIND.\n"
        )
    ENDIF()

    IF (${test} MATCHES "^step-petsc" AND NOT EXISTS ${test}-OK)
      MESSAGE("
Additional information about PETSc issues is available
at:\nhttp://www.dealii.org/developer/external-libs/petsc.html\n"
        )
    ENDIF()

    IF (${test} MATCHES "^p4est" AND NOT EXISTS ${test}-OK)
      MESSAGE("
The p4est test can fail if you are running an OpenMPI version before 1.5.
This is a known problem and the only work around is to update to a more
recent version or use a different MPI library like MPICH.\n"
        )
    ENDIF()

  ENDFOREACH()

  # The CMake command MESSAGE(SEND_ERROR ...) is, to the best of the authors'
  # knowledge, the only way to set the exit status of CMake to a nonzero value.
  # If we used MESSAGE(SEND_ERROR ...) at the top (with the actual error
  # message) then subsequent messages (i.e., the test specific help) would not
  # be printed. Hence, do it down here.
  MESSAGE(SEND_ERROR "")

ENDIF()
