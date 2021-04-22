## ---------------------------------------------------------------------
##
## Copyright (C) 2017 - 2019 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal2lkit library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#
# Try to find the SUNDIALS libraries
#
# This module exports
#
#   SUNDIALS_LIBRARIES
#   SUNDIALS_INCLUDE_DIR
#   SUNDIALS_WITH_IDAS
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

SET(SUNDIALS_DIR "" CACHE PATH "An optional hint to a SUNDIALS_DIR installation")
SET_IF_EMPTY(SUNDIALS_DIR "$ENV{SUNDIALS_DIR}")

DEAL_II_FIND_LIBRARY(SUNDIALS_LIB_IDAS NAMES sundials_idas
  HINTS ${SUNDIALS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

DEAL_II_FIND_LIBRARY(SUNDIALS_LIB_IDA NAMES sundials_ida
  HINTS ${SUNDIALS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

DEAL_II_FIND_LIBRARY(SUNDIALS_LIB_ARKODE NAMES sundials_arkode
  HINTS ${SUNDIALS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

DEAL_II_FIND_LIBRARY(SUNDIALS_LIB_KINSOL NAMES sundials_kinsol
  HINTS ${SUNDIALS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

DEAL_II_FIND_LIBRARY(SUNDIALS_LIB_SER NAMES sundials_nvecserial
  HINTS ${SUNDIALS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

DEAL_II_FIND_PATH(SUNDIALS_INCLUDE_DIR sundials/sundials_version.h
  HINTS ${SUNDIALS_DIR}
  PATH_SUFFIXES include
)

SET(_sundials_lib_par)
IF(DEAL_II_WITH_MPI)
  DEAL_II_FIND_LIBRARY(SUNDIALS_LIB_PAR NAMES sundials_nvecparallel
    HINTS ${SUNDIALS_DIR}
    PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
    )
  SET(_sundials_lib_par "SUNDIALS_LIB_PAR")
ENDIF()

#
# If IDAS is available, we prefer said library over IDA. We have to make
# sure to only link one of the two library variants:
#
IF(NOT SUNDIALS_LIB_IDAS MATCHES "-NOTFOUND")
  SET(_sundials_lib_ida "SUNDIALS_LIB_IDAS")
  SET(SUNDIALS_WITH_IDAS TRUE)
ELSE()
  SET(_sundials_lib_ida "SUNDIALS_LIB_IDA")
  SET(SUNDIALS_WITH_IDAS FALSE)
ENDIF()

#
# Extract SUNDIALS version.
#
DEAL_II_FIND_FILE(SUNDIALS_CONFIG_H
  NAMES sundials_config.h
  HINTS ${SUNDIALS_INCLUDE_DIR}/sundials
  )
IF(NOT SUNDIALS_CONFIG_H MATCHES "-NOTFOUND")
  FILE(STRINGS "${SUNDIALS_CONFIG_H}" SUNDIALS_VERSION_MAJOR_STRING
    REGEX "#define.*SUNDIALS_VERSION_MAJOR"
    )
  STRING(REGEX REPLACE "^.*SUNDIALS_VERSION_MAJOR.*([0-9]+).*" "\\1"
    SUNDIALS_VERSION_MAJOR "${SUNDIALS_VERSION_MAJOR_STRING}"
    )
  FILE(STRINGS "${SUNDIALS_CONFIG_H}" SUNDIALS_VERSION_MINOR_STRING
    REGEX "#define.*SUNDIALS_VERSION_MINOR"
    )
  STRING(REGEX REPLACE "^.*SUNDIALS_VERSION_MINOR.*([0-9]+).*" "\\1"
    SUNDIALS_VERSION_MINOR "${SUNDIALS_VERSION_MINOR_STRING}"
    )
  FILE(STRINGS "${SUNDIALS_CONFIG_H}" SUNDIALS_VERSION_PATCH_STRING
    REGEX "#define.*SUNDIALS_VERSION_PATCH"
    )
  STRING(REGEX REPLACE "^.*SUNDIALS_VERSION_PATCH.*([0-9]+).*" "\\1"
    SUNDIALS_VERSION_PATCH "${SUNDIALS_VERSION_PATCH_STRING}"
    )
ENDIF()
IF(NOT "${SUNDIALS_VERSION_MAJOR}")
  SET(SUNDIALS_VERSION_MAJOR "0")
ENDIF()
IF(NOT "${SUNDIALS_VERSION_MINOR}")
  SET(SUNDIALS_VERSION_MINOR "0")
ENDIF()
IF(NOT "${SUNDIALS_VERSION_PATCH}")
  SET(SUNDIALS_VERSION_PATCH "0")
ENDIF()
SET(SUNDIALS_VERSION
    "${SUNDIALS_VERSION_MAJOR}.${SUNDIALS_VERSION_MINOR}.${SUNDIALS_VERSION_PATCH}"
    )

DEAL_II_PACKAGE_HANDLE(SUNDIALS
  LIBRARIES REQUIRED
    ${_sundials_lib_ida}
    SUNDIALS_LIB_ARKODE
    SUNDIALS_LIB_KINSOL
    SUNDIALS_LIB_SER
    ${_sundials_lib_par}
  INCLUDE_DIRS REQUIRED
    SUNDIALS_INCLUDE_DIR
  USER_INCLUDE_DIRS REQUIRED
    SUNDIALS_INCLUDE_DIR
  CLEAR
    SUNDIALS_LIB_IDA
    SUNDIALS_LIB_IDAS
    SUNDIALS_LIB_ARKODE
    SUNDIALS_LIB_KINSOL
    SUNDIALS_LIB_SER
    ${_sundials_lib_par}
    SUNDIALS_INCLUDE_DIR
    SUNDIALS_CONFIG_H
  )
