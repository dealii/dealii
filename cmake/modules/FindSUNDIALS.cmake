## ---------------------------------------------------------------------
##
## Copyright (C) 2017 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal2lkit library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# Try to find the SUNDIALS libraries
#
# This module exports
#
#   SUNDIALS_LIB_IDA
#   SUNDIALS_LIB_KINSOL
#   SUNDIALS_LIB_PAR
#   SUNDIALS_LIB_SER
#   SUNDIALS_INCLUDE_DIR
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

DEAL_II_FIND_LIBRARY(SUNDIALS_LIB_IDA NAMES sundials_ida
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

#
# only define the extra hint if ${SUNDIALS_DIR} is nonempty so that we do not
# try to search through the top-level '/include/' directory, should it happen to
# exist
#
STRING(COMPARE EQUAL "${SUNDIALS_DIR}" "" _sundials_dir_is_empty)
IF(NOT ${_sundials_dir_is_empty})
  SET(_sundials_include_hint_dir "${SUNDIALS_DIR}/include/")
ENDIF()

DEAL_II_FIND_PATH(SUNDIALS_INCLUDE_DIR sundials/sundials_nvector.h
  HINTS ${SUNDIALS_DIR} "${_sundials_include_hint_dir}"
)

SET(_sundials_additional_libs)
IF(DEAL_II_WITH_MPI)
  DEAL_II_FIND_LIBRARY(SUNDIALS_LIB_PAR NAMES sundials_nvecparallel
    HINTS ${SUNDIALS_DIR}
    PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
    )
  SET(_sundials_additional_libs SUNDIALS_LIB_PAR)
ENDIF()

DEAL_II_PACKAGE_HANDLE(SUNDIALS
  LIBRARIES REQUIRED
    SUNDIALS_LIB_IDA SUNDIALS_LIB_KINSOL SUNDIALS_LIB_SER SUNDIALS_LIB_PAR
    ${_sundials_additional_libs}
  INCLUDE_DIRS REQUIRED SUNDIALS_INCLUDE_DIR
  USER_INCLUDE_DIRS REQUIRED SUNDIALS_INCLUDE_DIR
  CLEAR
    SUNDIALS_LIB_IDA SUNDIALS_LIB_KINSOL SUNDIALS_LIB_SER SUNDIALS_LIB_PAR
    SUN_INC
  )
