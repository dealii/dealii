## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2023 by the deal.II authors
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
# Try to find the Threading Building Blocks library
#
# This module exports
#
#   TBB_LIBRARIES
#   TBB_INCLUDE_DIRS
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
  set(_libraries
    LIBRARIES_RELEASE REQUIRED TBB_LIBRARY
    LIBRARIES_DEBUG REQUIRED TBB_DEBUG_LIBRARY
    )
else()
  set(_libraries LIBRARIES REQUIRED TBB_LIBRARY)
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
  ${_libraries}
  INCLUDE_DIRS REQUIRED TBB_INCLUDE_DIR
  CLEAR TBB_DEBUG_LIBRARY TBB_LIBRARY TBB_INCLUDE_DIR
  )
