## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2013 - 2022 by the deal.II authors
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

set(MUPARSER_DIR "" CACHE PATH "An optional hint to a MUPARSER installation")
set_if_empty(MUPARSER_DIR "$ENV{MUPARSER_DIR}")

deal_ii_find_library(MUPARSER_LIBRARY
  NAMES muparser muparserd
  HINTS ${MUPARSER_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

deal_ii_find_path(MUPARSER_INCLUDE_DIR muParserDef.h
  HINTS ${MUPARSER_DIR}
  PATH_SUFFIXES include
  )

if(EXISTS ${MUPARSER_INCLUDE_DIR}/muParserDef.h)
  file(STRINGS "${MUPARSER_INCLUDE_DIR}/muParserDef.h" MUPARSER_VERSION_STRING_LINE
    # Try to match the line
    #
    #     #define MUP_VERSION _T("2.2.4")
    REGEX "^[ \t]*#[ \t]*define[ \t]+MUP_VERSION _T"
    )

  if("${MUPARSER_VERSION_STRING_LINE}" STREQUAL "")
    # try again with the newer version format (starting in at least 2.3.2),
    # which matches the line
    #
    #     static const string_type ParserVersion = string_type(_T("2.3.2"));
    file(STRINGS "${MUPARSER_INCLUDE_DIR}/muParserDef.h" MUPARSER_VERSION_STRING_LINE
      REGEX "string_type ParserVersion = string_type"
      )
  endif()

  string(REGEX REPLACE ".*\"(.*)\".*" "\\1"
    _VERSION_STRING "${MUPARSER_VERSION_STRING_LINE}"
    )

  string(REPLACE "." ";" _VERSION_LIST ${_VERSION_STRING})
  list(GET _VERSION_LIST 0 MUPARSER_VERSION_MAJOR)
  list(GET _VERSION_LIST 1 MUPARSER_VERSION_MINOR)
  list( LENGTH _VERSION_LIST _LISTLEN ) 
  if (${_LISTLEN} GREATER 2)
    list(GET _VERSION_LIST 2 MUPARSER_VERSION_SUBMINOR)
  else()
    set(MUPARSER_VERSION_SUBMINOR "0")
  endif()

  set(MUPARSER_VERSION
    "${MUPARSER_VERSION_MAJOR}.${MUPARSER_VERSION_MINOR}.${MUPARSER_VERSION_SUBMINOR}"
    )
endif()

process_feature(MUPARSER
  LIBRARIES REQUIRED MUPARSER_LIBRARY
  INCLUDE_DIRS REQUIRED MUPARSER_INCLUDE_DIR
  CLEAR MUPARSER_LIBRARY MUPARSER_INCLUDE_DIR
  )
