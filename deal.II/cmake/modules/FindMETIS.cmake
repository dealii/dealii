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
# Try to find the (serial) METIS library
#
# This module exports
#
#   METIS_LIBRARIES
#   METIS_INCLUDE_DIRS
#

INCLUDE(FindPackageHandleStandardArgs)

SET_IF_EMPTY(METIS_DIR "$ENV{METIS_DIR}")

#
# TODO: Metis is usually pretty self contained. So no external dependencies
# so far... But there could be dependencies on pcre and mpi...
#

FIND_PATH(METIS_INCLUDE_DIR metis.h
  HINTS
    ${METIS_DIR}
  PATH_SUFFIXES
    metis include/metis include
  )

FIND_LIBRARY(METIS_LIBRARY
  NAMES metis
  HINTS
    ${METIS_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib
    # This is a hint, isn't it?
    build/${CMAKE_CXX_PLATFORM_ID}-${CMAKE_SYSTEM_PROCESSOR}/libmetis
  )

SET(METIS_LIBRARIES ${METIS_LIBRARY})
SET(METIS_INCLUDE_DIRS ${METIS_INCLUDE_DIRS})

FIND_PACKAGE_HANDLE_STANDARD_ARGS(METIS DEFAULT_MSG METIS_LIBRARIES METIS_INCLUDE_DIRS)

IF(METIS_FOUND)
  #
  # Extract the version number out of metis.h
  #
  FILE(STRINGS "${METIS_INCLUDE_DIR}/metis.h" METIS_MAJOR_STRING
    REGEX "METIS_VER_MAJOR")
  STRING(REGEX REPLACE "^.*METIS_VER_MAJOR.* ([0-9]+).*" "\\1" METIS_MAJOR "${METIS_MAJOR_STRING}")

  FILE(STRINGS "${METIS_INCLUDE_DIR}/metis.h" METIS_MINOR_STRING
    REGEX "METIS_VER_MINOR")
  STRING(REGEX REPLACE "^.*METIS_VER_MINOR.* ([0-9]+).*" "\\1" METIS_MINOR "${METIS_MINOR_STRING}")

  MARK_AS_ADVANCED(
    METIS_LIBRARY
    METIS_INCLUDE_DIR
    METIS_DIR
  )
ELSE()
  SET(METIS_DIR "" CACHE STRING
    "An optional hint to a metis directory"
    )
ENDIF()

