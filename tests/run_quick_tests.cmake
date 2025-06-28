## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2013 - 2024 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

#
# This file is run when "make test" is executed by the user and is
# responsible for running the tests and printing some helpful error
# messages.
#

include(ProcessorCount)
PROCESSORCOUNT(_n_processors)

if(_n_processors EQUAL 0)
  set(_n_processors "1")
endif()

# Avoid race conditions with native Windows build tools:
if(CMAKE_HOST_SYSTEM_NAME MATCHES "Windows")
  set(_n_processors "1")
endif()

if("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE "Debug")
endif()
message(STATUS "Running quick_tests in ${CMAKE_BUILD_TYPE} mode with -j${_n_processors}:")

#
# Redirect quick tests to an output.log for CMake versions 3.18 onwards
# (when we can actually instruct CMake to also print to the standard
# terminal).
#

set(_redirect_output "")
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
  set(_redirect_output
    OUTPUT_VARIABLE _output ECHO_OUTPUT_VARIABLE
    ERROR_VARIABLE _output ECHO_ERROR_VARIABLE
    )
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" _build_type)
execute_process(COMMAND ${CMAKE_CTEST_COMMAND}
  -j${_n_processors} -C ${CMAKE_BUILD_TYPE} --force-new-ctest-process
  -R "(a-framework|quick_tests|examples)/"
  --output-on-failure
  ${_redirect_output}
  RESULT_VARIABLE _return_value
  )
file(WRITE quick_tests.log ${_output})

if("${_return_value}" STREQUAL "0")
  message("
**************************************************************************

                 🎉  🎉  🎉  All quick tests passed.  🎉  🎉  🎉

**************************************************************************"
    )
else()
  message("
**************************************************************************
**                                                                      **
**                Error: Some of the quick tests failed.                **
**                                                                      **
**************************************************************************

Check the terminal output above. You can run all failing quick tests again
via  $ ctest --rerun-failed --output-on-failure -R quick_tests/

If you need help with this problem, open a bug report on the issue tracker

               https://github.com/dealii/dealii/issues/new

or write to the mailing list linked at https://www.dealii.org/mail.html
Please also add the files quick_tests.log , Testing/Temporary/LastTest.log
to your bug report."
    )

  set(_affinity_issue FALSE)
  string(REPLACE "\n" ";" _output "${_output}")
  foreach(_line ${_output})
    if ("${_line}" MATCHES " Test .*quick_tests/affinity.*Failed.*sec")
      set(_affinity_issue TRUE)
    endif()
  endforeach()
  if(_affinity_issue)
    message("
The affinity test can fail when you are linking in a library like BLAS
which uses OpenMP. Even without calling any BLAS functions, OpenMP messes
with the thread affinity which causes TBB to run single-threaded only. You
can fix this by exporting OMP_NUM_THREADS=1. Also see GOMP_CPU_AFFINITY
and OMP_PROC_BIND."
      )
  endif()

  message("\n")

  # ensure that this script exits with a non-zero exit code
  message(FATAL_ERROR "quick tests failed")
endif()
