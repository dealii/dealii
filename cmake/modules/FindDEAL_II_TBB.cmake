## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2022 by the deal.II authors
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
# Try to find the Threading Building Blocks library
#
# This module exports
#
#   TBB_LIBRARIES
#   TBB_INCLUDE_DIRS
#   TBB_WITH_DEBUGLIB
#   TBB_WITH_ONEAPI
#   TBB_VERSION
#   TBB_VERSION_MAJOR
#   TBB_VERSION_MINOR
#

set(TBB_DIR "" CACHE PATH "An optional hint to a TBB installation")
set_if_empty(TBB_DIR "$ENV{TBB_DIR}")

file(GLOB _path ${TBB_DIR}/build/*_release ${TBB_DIR}/lib/intel64/gcc*)
deal_ii_find_library(TBB_LIBRARY
  NAMES tbb
  HINTS
    ${_path}
    ${TBB_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

#
# Also search for the debug library:
#
file(GLOB _path ${TBB_DIR}/build/*_debug ${TBB_DIR}/lib/intel64/gcc*)
deal_ii_find_library(TBB_DEBUG_LIBRARY
  NAMES tbb_debug
  HINTS
    ${_path}
    ${TBB_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )
if(NOT TBB_DEBUG_LIBRARY MATCHES "-NOTFOUND")
  set(TBB_WITH_DEBUGLIB TRUE)
  set(_libraries debug TBB_DEBUG_LIBRARY optimized TBB_LIBRARY)
else()
  set(_libraries TBB_LIBRARY)
endif()

#
# Check for old TBB header layout:
#

deal_ii_find_path(TBB_INCLUDE_DIR tbb/tbb_stddef.h
  HINTS
    ${TBB_DIR}
  PATH_SUFFIXES include include/tbb tbb
  )

if(EXISTS ${TBB_INCLUDE_DIR}/tbb/tbb_stddef.h)
  file(STRINGS "${TBB_INCLUDE_DIR}/tbb/tbb_stddef.h" TBB_VERSION_MAJOR_STRING
    REGEX "#define.*TBB_VERSION_MAJOR")
  string(REGEX REPLACE "^.*TBB_VERSION_MAJOR +([0-9]+).*" "\\1"
    TBB_VERSION_MAJOR "${TBB_VERSION_MAJOR_STRING}"
    )
  file(STRINGS "${TBB_INCLUDE_DIR}/tbb/tbb_stddef.h" TBB_VERSION_MINOR_STRING
    REGEX "#define.*TBB_VERSION_MINOR")
  string(REGEX REPLACE "^.*TBB_VERSION_MINOR +([0-9]+).*" "\\1"
    TBB_VERSION_MINOR "${TBB_VERSION_MINOR_STRING}"
    )
  set(TBB_VERSION
    "${TBB_VERSION_MAJOR}.${TBB_VERSION_MINOR}"
    )

  set(TBB_WITH_ONEAPI FALSE)
else()

  #
  # Check for new oneAPI TBB header layout:
  #

  deal_ii_find_path(TBB_INCLUDE_DIR oneapi/tbb/version.h
    HINTS
      ${TBB_DIR}
    PATH_SUFFIXES include include/tbb tbb
    )

  if(EXISTS ${TBB_INCLUDE_DIR}/oneapi/tbb/version.h)
    file(STRINGS "${TBB_INCLUDE_DIR}/oneapi/tbb/version.h" TBB_VERSION_MAJOR_STRING
      REGEX "#define.*TBB_VERSION_MAJOR")
    string(REGEX REPLACE "^.*TBB_VERSION_MAJOR +([0-9]+).*" "\\1"
      TBB_VERSION_MAJOR "${TBB_VERSION_MAJOR_STRING}"
      )
    file(STRINGS "${TBB_INCLUDE_DIR}/oneapi/tbb/version.h" TBB_VERSION_MINOR_STRING
      REGEX "#define.*TBB_VERSION_MINOR")
    string(REGEX REPLACE "^.*TBB_VERSION_MINOR +([0-9]+).*" "\\1"
      TBB_VERSION_MINOR "${TBB_VERSION_MINOR_STRING}"
      )
    set(TBB_VERSION
      "${TBB_VERSION_MAJOR}.${TBB_VERSION_MINOR}"
      )
  endif()

  set(TBB_WITH_ONEAPI TRUE)
endif()

process_feature(TBB
  LIBRARIES REQUIRED ${_libraries}
  INCLUDE_DIRS REQUIRED TBB_INCLUDE_DIR
  CLEAR TBB_DEBUG_LIBRARY TBB_LIBRARY TBB_INCLUDE_DIR
  )
