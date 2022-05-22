## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2022 by the deal.II authors
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
# Try to find the MUPARSER library
#
# This module exports
#
#   MUPARSER_LIBRARIES
#   MUPARSER_INCLUDE_DIRS
#   MUPARSER_VERSION
#   MUPARSER_VERSION_MAJOR
#   MUPARSER_VERSION_MINOR
#   MUPARSER_VERSION_SUBMINOR
#

SET(MUPARSER_DIR "" CACHE PATH "An optional hint to a MUPARSER installation")
SET_IF_EMPTY(MUPARSER_DIR "$ENV{MUPARSER_DIR}")

DEAL_II_FIND_LIBRARY(MUPARSER_LIBRARY
  NAMES muparser muparserd
  HINTS ${MUPARSER_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

DEAL_II_FIND_PATH(MUPARSER_INCLUDE_DIR muParserDef.h
  HINTS ${MUPARSER_DIR}
  PATH_SUFFIXES include
  )

IF(EXISTS ${MUPARSER_INCLUDE_DIR}/muParserDef.h)
  FILE(STRINGS "${MUPARSER_INCLUDE_DIR}/muParserDef.h" MUPARSER_VERSION_STRING_LINE
    # Try to match the line
    #
    #     #define MUP_VERSION _T("2.2.4")
    REGEX "^[ \t]*#[ \t]*define[ \t]+MUP_VERSION _T"
    )

  IF("${MUPARSER_VERSION_STRING_LINE}" STREQUAL "")
    # try again with the newer version format (starting in at least 2.3.2),
    # which matches the line
    #
    #     static const string_type ParserVersion = string_type(_T("2.3.2"));
    FILE(STRINGS "${MUPARSER_INCLUDE_DIR}/muParserDef.h" MUPARSER_VERSION_STRING_LINE
      REGEX "string_type ParserVersion = string_type"
      )
  ENDIF()

  STRING(REGEX REPLACE ".*\"(.*)\".*" "\\1"
    _VERSION_STRING "${MUPARSER_VERSION_STRING_LINE}"
    )

  STRING(REPLACE "." ";" _VERSION_LIST ${_VERSION_STRING})
  LIST(GET _VERSION_LIST 0 MUPARSER_VERSION_MAJOR)
  LIST(GET _VERSION_LIST 1 MUPARSER_VERSION_MINOR)
  LIST( LENGTH _VERSION_LIST _LISTLEN ) 
  IF (${_LISTLEN} GREATER 2)
    LIST(GET _VERSION_LIST 2 MUPARSER_VERSION_SUBMINOR)
  ELSE()
    SET(MUPARSER_VERSION_SUBMINOR "0")
  ENDIF()

  SET(MUPARSER_VERSION
    "${MUPARSER_VERSION_MAJOR}.${MUPARSER_VERSION_MINOR}.${MUPARSER_VERSION_SUBMINOR}"
    )
ENDIF()

DEAL_II_PACKAGE_HANDLE(MUPARSER
  LIBRARIES REQUIRED MUPARSER_LIBRARY
  INCLUDE_DIRS REQUIRED MUPARSER_INCLUDE_DIR
  CLEAR MUPARSER_LIBRARY MUPARSER_INCLUDE_DIR
  )
