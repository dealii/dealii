## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2026 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

#
# A small wrapper around an executable to do stdin/stdout redirection. This
# works around an issue with the ctest launcher that gets enabled with
# CTEST_USE_LAUNCHER: "ctest --launch" fails to forward its own stdin to
# the launched process on many CMake/CTest releases (3.28 through 4.1+,
# backported inconsistently, so we cannot detect it).
#
# https://gitlab.kitware.com/cmake/cmake/-/issues/27610
#
# The following variables must be set:
#
#   EXE             - the expand_instantiations executable to run
#   CLASS_LIST_FILE - the template-arguments file (passed as a command
#                      line argument to EXE)
#   IN              - the .inst.in file to feed to EXE via stdin
#   OUT             - the file to capture EXE's stdout in
#

execute_process(
  COMMAND ${EXE} ${CLASS_LIST_FILE}
  INPUT_FILE ${IN}
  OUTPUT_FILE ${OUT}
  RESULT_VARIABLE _result
  ERROR_VARIABLE _error
  )

if(NOT _result EQUAL 0)
  message(FATAL_ERROR
    "expand_instantiations failed (exit code ${_result}) while processing "
    "'${IN}':\n${_error}")
endif()
