## ---------------------------------------------------------------------
##
## Copyright (C) 2020 - 2021 by the deal.II authors
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
# Try to find the Taskflow library
#
# This module exports
#
#   TASKFLOW_INCLUDE_DIRS
#   TASKFLOW_VERSION
#

SET(TASKFLOW_DIR "" CACHE PATH "An optional hint to a Taskflow installation")
SET_IF_EMPTY(TASKFLOW_DIR "$ENV{TASKFLOW_DIR}")

FIND_PACKAGE(TASKFLOW_CONFIG
  CONFIG QUIET
  NAMES Taskflow
  HINTS
    ${TASKFLOW_DIR}/lib/cmake/Taskflow
    ${TASKFLOW_DIR}
  PATH_SUFFIXES
    lib64/cmake/Taskflow
    lib/cmake/Taskflow
    lib${LIB_SUFFIX}/cmake/Taskflow
  NO_SYSTEM_ENVIRONMENT_PATH
  )

SET(TASKFLOW_INCLUDE_DIR ${Taskflow_INCLUDE_DIR})

#
# Extract version numbers:
#
SET(TASKFLOW_VERSION "${TASKFLOW_CONFIG_VERSION}")
STRING(REGEX REPLACE
  "^([0-9]+).*$" "\\1"
  TASKFLOW_VERSION_MAJOR "${TASKFLOW_CONFIG_VERSION}")
STRING(REGEX REPLACE
  "^[0-9]+\\.([0-9]+).*$" "\\1"
  TASKFLOW_VERSION_MINOR "${TASKFLOW_CONFIG_VERSION}")

DEAL_II_PACKAGE_HANDLE(TASKFLOW
  INCLUDE_DIRS REQUIRED TASKFLOW_INCLUDE_DIR
  USER_INCLUDE_DIRS REQUIRED TASKFLOW_INCLUDE_DIR
  CLEAR TASKFLOW_CONFIG_DIR
  )
