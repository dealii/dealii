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

# This file is run when "make test" is executed by the user and is 
# responsible for running the tests and printing some helpful
# error messages.

SEPARATE_ARGUMENTS(ALL_TESTS)

EXECUTE_PROCESS(COMMAND ${CMAKE_CTEST_COMMAND} --force-new-ctest-process --output-on-failure -O quicktests.log RESULT_VARIABLE res_var)

if(NOT "${res_var}" STREQUAL "0")
  MESSAGE( "

*******************************     WARNING     *******************************

Some of the tests failed!

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
can fix this by exporting OMP_NUM_THREADS=1.\n"
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

ENDIF()
