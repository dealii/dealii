## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2017 - 2025 by the deal.II authors
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
# Try to find the SUNDIALS libraries
#
# This module exports
#
#   SUNDIALS_LIBRARIES
#   SUNDIALS_INCLUDE_DIR
#   SUNDIALS_WITH_IDAS
#   SUNDIALS_WITH_MPI
#   SUNDIALS_VERSION
#   SUNDIALS_VERSION_MAJOR
#   SUNDIALS_VERSION_MINOR
#   SUNDIALS_VERSION_PATCH
#
# Note that sundials headers are typically installed in several directories,
# e.g.,
#
# /prefix/include/ida/ida.h
# /prefix/include/kinsol/kinsol.h
# /prefix/include/sundials/sundials_nvector.h
#
# Here SUNDIALS_INCLUDE_DIR is just '/prefix/include/', not
# '/prefix/include/sundials/'.

set(SUNDIALS_DIR "" CACHE PATH "An optional hint to a SUNDIALS_DIR installation")
set_if_empty(SUNDIALS_DIR "$ENV{SUNDIALS_DIR}")

deal_ii_find_library(SUNDIALS_LIB_CORE NAMES sundials_core
  HINTS ${SUNDIALS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

deal_ii_find_library(SUNDIALS_LIB_IDAS NAMES sundials_idas
  HINTS ${SUNDIALS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

deal_ii_find_library(SUNDIALS_LIB_IDA NAMES sundials_ida
  HINTS ${SUNDIALS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

deal_ii_find_library(SUNDIALS_LIB_ARKODE NAMES sundials_arkode
  HINTS ${SUNDIALS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

deal_ii_find_library(SUNDIALS_LIB_KINSOL NAMES sundials_kinsol
  HINTS ${SUNDIALS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

deal_ii_find_library(SUNDIALS_LIB_SER NAMES sundials_nvecserial
  HINTS ${SUNDIALS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

deal_ii_find_path(SUNDIALS_INCLUDE_DIR sundials/sundials_version.h
  HINTS ${SUNDIALS_DIR}
  PATH_SUFFIXES include
)

set(_sundials_lib_par)
if(DEAL_II_WITH_MPI)
  deal_ii_find_library(SUNDIALS_LIB_PAR NAMES sundials_nvecparallel
    HINTS ${SUNDIALS_DIR}
    PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
    )
  set(_sundials_lib_par "SUNDIALS_LIB_PAR")
endif()

#
# If IDAS is available, we prefer said library over IDA. We have to make
# sure to only link one of the two library variants:
#
if(NOT SUNDIALS_LIB_IDAS MATCHES "-NOTFOUND")
  set(_sundials_lib_ida "SUNDIALS_LIB_IDAS")
  set(SUNDIALS_WITH_IDAS TRUE)
else()
  set(_sundials_lib_ida "SUNDIALS_LIB_IDA")
  set(SUNDIALS_WITH_IDAS FALSE)
endif()

#
# Extract SUNDIALS version.
#
deal_ii_find_file(SUNDIALS_CONFIG_H
  NAMES sundials_config.h
  HINTS ${SUNDIALS_INCLUDE_DIR}/sundials
  )
if(NOT SUNDIALS_CONFIG_H MATCHES "-NOTFOUND")
  file(STRINGS "${SUNDIALS_CONFIG_H}" SUNDIALS_VERSION_MAJOR_STRING
    REGEX "#define.*SUNDIALS_VERSION_MAJOR"
    )
  string(REGEX REPLACE "^.*SUNDIALS_VERSION_MAJOR.*([0-9]+).*" "\\1"
    SUNDIALS_VERSION_MAJOR "${SUNDIALS_VERSION_MAJOR_STRING}"
    )
  file(STRINGS "${SUNDIALS_CONFIG_H}" SUNDIALS_VERSION_MINOR_STRING
    REGEX "#define.*SUNDIALS_VERSION_MINOR"
    )
  string(REGEX REPLACE "^.*SUNDIALS_VERSION_MINOR.*([0-9]+).*" "\\1"
    SUNDIALS_VERSION_MINOR "${SUNDIALS_VERSION_MINOR_STRING}"
    )
  file(STRINGS "${SUNDIALS_CONFIG_H}" SUNDIALS_VERSION_PATCH_STRING
    REGEX "#define.*SUNDIALS_VERSION_PATCH"
    )
  string(REGEX REPLACE "^.*SUNDIALS_VERSION_PATCH.*([0-9]+).*" "\\1"
    SUNDIALS_VERSION_PATCH "${SUNDIALS_VERSION_PATCH_STRING}"
    )

  set(SUNDIALS_VERSION
    "${SUNDIALS_VERSION_MAJOR}.${SUNDIALS_VERSION_MINOR}.${SUNDIALS_VERSION_PATCH}"
    )

  #
  # Determine whether Sundials was configured with MPI:
  #
  file(STRINGS "${SUNDIALS_CONFIG_H}" SUNDIALS_MPI_STRING
    REGEX "^[ \t]*#[ \t]*define[ \t]+SUNDIALS_MPI_ENABLED[ \t]1")
  if("${SUNDIALS_MPI_STRING}" STREQUAL "")
    set(SUNDIALS_WITH_MPI FALSE)
  else()
    set(SUNDIALS_WITH_MPI TRUE)
  endif()
endif()

#
# sundials_core is only required as of version 7.0.0.
#
set(_sundials_lib_core)
if(SUNDIALS_VERSION VERSION_GREATER_EQUAL 7)
  set(_sundials_lib_core "SUNDIALS_LIB_CORE")
endif()

process_feature(SUNDIALS
  LIBRARIES REQUIRED
    ${_sundials_lib_core}
    ${_sundials_lib_ida}
    SUNDIALS_LIB_ARKODE
    SUNDIALS_LIB_KINSOL
    SUNDIALS_LIB_SER
    ${_sundials_lib_par}
  INCLUDE_DIRS REQUIRED
    SUNDIALS_INCLUDE_DIR
  CLEAR
    SUNDIALS_LIB_CORE
    SUNDIALS_LIB_IDA
    SUNDIALS_LIB_IDAS
    SUNDIALS_LIB_ARKODE
    SUNDIALS_LIB_KINSOL
    SUNDIALS_LIB_SER
    ${_sundials_lib_par}
    SUNDIALS_INCLUDE_DIR
    SUNDIALS_CONFIG_H
  )
