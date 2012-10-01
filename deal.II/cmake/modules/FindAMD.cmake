#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

#
# Try to find the AMD library
#
# This is a helper module for FindUMFPACK.cmake
#
# This module exports
#
#   AMD_LIBRARY
#   AMD_INCLUDE_DIR
#

INCLUDE(FindPackageHandleStandardArgs)

SET_IF_EMPTY(AMD_DIR "$ENV{AMD_DIR}")
SET_IF_EMPTY(UMFPACK_DIR "$ENV{UMFPACK_DIR}")

FIND_PATH(AMD_INCLUDE_DIR amd.h
  HINTS
    ${AMD_DIR}
    ${UMFPACK_DIR}
    ${UMFPACK_DIR}/../AMD/
  PATH_SUFFIXES
    amd include/amd include Include AMD/Include
)

FIND_LIBRARY(AMD_LIBRARY
  NAMES amd
  HINTS
    ${AMD_DIR}
    ${UMFPACK_DIR}
    ${UMFPACK_DIR}/../AMD
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib Lib AMD/Lib
  )

FIND_PACKAGE_HANDLE_STANDARD_ARGS(AMD DEFAULT_MSG AMD_LIBRARY AMD_INCLUDE_DIR)

IF(AMD_FOUND)
  MARK_AS_ADVANCED(
    AMD_LIBRARY
    AMD_INCLUDE_DIR
    AMD_DIR
  )
ELSE()
  SET(AMD_DIR "" CACHE STRING
    "An optional hint to an AMD directory"
    )
ENDIF()

