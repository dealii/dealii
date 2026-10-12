## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2023 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

#
# Glob all test summary files and remove them
#

if("${_summary_files_prefix}" STREQUAL "")
  message(FATAL_ERROR "The variable _summary_files_prefix is not set.")
endif()

file(GLOB _files "${_summary_files_prefix}_*")
if(NOT "${_files}" STREQUAL "")
  file(REMOVE ${_files})
endif()
