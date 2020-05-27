## ---------------------------------------------------------------------
##
## Copyright (C) 2020 by the deal.II authors
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
# Try to find the Threading Building Blocks library
#
# This module exports
#
#   CPP_TASKFLOW_INCLUDE_DIRS
#   CPP_TASKFLOW_VERSION
#

SET(CPP_TASKFLOW_DIR "" CACHE PATH "An optional hint to a Cpp Taskflow installation")
SET_IF_EMPTY(CPP_TASKFLOW_DIR "$ENV{CPP_TASKFLOW_DIR}")

FIND_PACKAGE(CPP_TASKFLOW_CONFIG
  CONFIG QUIET
  NAMES Cpp-Taskflow
  HINTS
    ${CPP_TASKFLOW_DIR}/lib/cmake/Cpp-Taskflow
    ${CPP_TASKFLOW_DIR}
  PATH_SUFFIXES
    lib64/cmake/Cpp-Taskflow
    lib/cmake/Cpp-Taskflow
    lib${LIB_SUFFIX}/cmake/Cpp-Taskflow
  NO_SYSTEM_ENVIRONMENT_PATH
  )

SET(CPP_TASKFLOW_INCLUDE_DIR ${Cpp-Taskflow_INCLUDE_DIR})

#
# Extract version numbers:
#
SET(CPP_TASKFLOW_VERSION "${CPP_TASKFLOW_CONFIG_VERSION}")
STRING(REGEX REPLACE
  "^([0-9]+).*$" "\\1"
  CPP_TASKFLOW_VERSION_MAJOR "${CPP_TASKFLOW_CONFIG_VERSION}")
STRING(REGEX REPLACE
  "^[0-9]+\\.([0-9]+).*$" "\\1"
  CPP_TASKFLOW_VERSION_MINOR "${CPP_TASKFLOW_CONFIG_VERSION}")

DEAL_II_PACKAGE_HANDLE(CPP_TASKFLOW
  INCLUDE_DIRS REQUIRED CPP_TASKFLOW_INCLUDE_DIR
  USER_INCLUDE_DIRS REQUIRED CPP_TASKFLOW_INCLUDE_DIR
  CLEAR CPP_TASKFLOW_CONFIG_DIR
  )
