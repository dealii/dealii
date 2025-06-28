## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2013 - 2025 by the deal.II authors
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
# Configuration for thread support in deal.II with the help of the Taskflow
# library:
#

set(DEAL_II_WITH_TASKFLOW ON CACHE BOOL "")

macro(feature_taskflow_find_external var)
  find_package(DEAL_II_TASKFLOW)

  if(TASKFLOW_FOUND)
    set(${var} TRUE)
  endif()

  if(TASKFLOW_VERSION VERSION_LESS "3.10")
    # Clear the previously determined version numbers to avoid confusion
    set(TASKFLOW_VERSION "bundled")
    set(TASKFLOW_VERSION_MAJOR "")
    set(TASKFLOW_VERSION_MINOR "")

    message(STATUS
      "The externally provided Taskflow library is older than version 3.10, "
      "which cannot be used with deal.II."
      )
    set(TASKFLOW_ADDITIONAL_ERROR_STRING
      "The externally provided Taskflow library is older than version\n"
      "3.10, which is the oldest version compatible with deal.II."
      )
    set(${var} FALSE)
  endif()
endmacro()


configure_feature(TASKFLOW)
