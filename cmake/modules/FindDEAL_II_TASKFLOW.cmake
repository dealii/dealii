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

set(TASKFLOW_DIR "" CACHE PATH "An optional hint to a Taskflow installation")
set_if_empty(TASKFLOW_DIR "$ENV{TASKFLOW_DIR}")

find_package(TASKFLOW_CONFIG
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

set(TASKFLOW_INCLUDE_DIR ${Taskflow_INCLUDE_DIR})

#
# Extract version numbers:
#
set(TASKFLOW_VERSION "${TASKFLOW_CONFIG_VERSION}")
string(REGEX REPLACE
  "^([0-9]+).*$" "\\1"
  TASKFLOW_VERSION_MAJOR "${TASKFLOW_CONFIG_VERSION}")
string(REGEX REPLACE
  "^[0-9]+\\.([0-9]+).*$" "\\1"
  TASKFLOW_VERSION_MINOR "${TASKFLOW_CONFIG_VERSION}")

process_feature(TASKFLOW
  INCLUDE_DIRS REQUIRED TASKFLOW_INCLUDE_DIR
  CLEAR TASKFLOW_CONFIG_DIR
  )
