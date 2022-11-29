## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2021 by the deal.II authors
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
# Configuration for thread support in deal.II with the help of the Taskflow
# library:
#

#
# Disallow the default detection of taskflow for the time being to avoid
# unpleasant surprises on user side.
#
set(DEAL_II_WITH_TASKFLOW OFF CACHE BOOL "")

macro(feature_taskflow_find_external var)
  find_package(DEAL_II_TASKFLOW)

  if(TASKFLOW_FOUND)
    set(${var} TRUE)
  endif()

  if(TASKFLOW_VERSION VERSION_LESS "2.4")
    # Clear the previously determined version numbers to avoid confusion
    set(TASKFLOW_VERSION "bundled")
    set(TASKFLOW_VERSION_MAJOR "")
    set(TASKFLOW_VERSION_MINOR "")

    message(STATUS
      "The externally provided Taskflow library is older than version 2.4, "
      "which cannot be used with deal.II."
      )
    set(TASKFLOW_ADDITIONAL_ERROR_STRING
      "The externally provided Taskflow library is older than version\n"
      "2.4, which is the oldest version compatible with deal.II."
      )
    set(${var} FALSE)
  endif()


  if(NOT TASKFLOW_VERSION VERSION_LESS "3.0" AND NOT DEAL_II_HAVE_CXX17)
    # Clear the previously determined version numbers to avoid confusion
    set(TASKFLOW_VERSION "bundled")
    set(TASKFLOW_VERSION_MAJOR "")
    set(TASKFLOW_VERSION_MINOR "")

    message(STATUS
      "The externally provided Taskflow library (version 3.0 onwards)
      requires C++17 support, which has not been configured."
      )
    set(TASKFLOW_ADDITIONAL_ERROR_STRING
      "The externally provided Taskflow library (version 3.0 onwards) "
      "requires C++17 support, but no C++17 support had been detected "
      "during configuration.\n"
      "Try to set -DDEAL_II_CXX_FLAGS=\"-std=c++17\" by hand.\n"
      )
    set(${var} FALSE)
  endif()
endmacro()


configure_feature(TASKFLOW)

