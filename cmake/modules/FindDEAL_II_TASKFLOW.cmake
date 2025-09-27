## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2020 - 2022 by the deal.II authors
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
