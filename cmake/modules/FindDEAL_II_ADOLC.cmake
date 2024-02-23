## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2017 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

# Note: This file is adapted from https://raw.githubusercontent.com/joaoleal/CppADCodeGen/master/cmake/FindADOLC.cmake

#
# Try to find Adolc
#
# This module exports
#
#   ADOLC_INCLUDE_DIR
#   ADOLC_LIBRARY
#   ADOLC_WITH_ADVANCED_BRANCHING
#   ADOLC_WITH_ATRIG_ERF
#   ADOLC_WITH_BOOST_ALLOCATOR
#

set(ADOLC_DIR "" CACHE PATH "An optional hint to an ADOL-C installation")
set_if_empty(ADOLC_DIR "$ENV{ADOLC_DIR}")

deal_ii_find_path(ADOLC_INCLUDE_DIR
  NAMES adolc/adolc.h
  HINTS ${ADOLC_DIR}
  PATH_SUFFIXES include
  )

deal_ii_find_library(ADOLC_LIBRARY
  NAMES adolc
  HINTS ${ADOLC_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

#
# Look for adolc_settings.h - we'll query it to determine supported features:
#

deal_ii_find_file(ADOLC_SETTINGS_H adolc_settings.h
  HINTS ${ADOLC_INCLUDE_DIR} "${ADOLC_INCLUDE_DIR}/adolc/internal"
  NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH NO_CMAKE_FIND_ROOT_PATH
  )

if(EXISTS ${ADOLC_SETTINGS_H})
  #
  # Check whether ADOL-C is configured with extra trig functions
  #
  file(STRINGS "${ADOLC_SETTINGS_H}" ADOLC_ATRIG_ERF_STRING
    REGEX "^[ \t]*#[ \t]*define[ \t]+ATRIG_ERF"
    )
  if(NOT "${ADOLC_ATRIG_ERF_STRING}" STREQUAL "")
    set(ADOLC_WITH_ATRIG_ERF TRUE)
  else()
    set(ADOLC_WITH_ATRIG_ERF FALSE)
  endif()

  #
  # Check whether ADOL-C is configured with advanced branching
  #
  file(STRINGS "${ADOLC_SETTINGS_H}" ADOLC_ADVANCED_BRANCHING_STRING
    REGEX "^[ \t]*#[ \t]*define[ \t]+ADOLC_ADVANCED_BRANCHING"
    )
  if(NOT "${ADOLC_ADVANCED_BRANCHING_STRING}" STREQUAL "")
    set(ADOLC_WITH_ADVANCED_BRANCHING TRUE)
  else()
    set(ADOLC_WITH_ADVANCED_BRANCHING FALSE)
  endif()

  #
  # Check whether ADOL-C is configured with tapeless number reference counting
  #
  file(STRINGS "${ADOLC_SETTINGS_H}" ADOLC_WITH_TAPELESS_REFCOUNTING_STRING
    REGEX "^[ \t]*#[ \t]*define[ \t]+USE_ADTL_REFCOUNTING 1"
    )
  if(NOT "${ADOLC_WITH_TAPELESS_REFCOUNTING_STRING}" STREQUAL "")
    set(ADOLC_WITH_TAPELESS_REFCOUNTING TRUE)
  else()
    set(ADOLC_WITH_TAPELESS_REFCOUNTING FALSE)
  endif()

  #
  # Check whether ADOL-C is configured to use the Boost pool allocator
  #
  file(STRINGS "${ADOLC_SETTINGS_H}" ADOLC_BOOST_POOL_STRING
    REGEX "^[ \t]*#[ \t]*define[ \t]+USE_BOOST_POOL 1"
    )
  if(NOT "${ADOLC_BOOST_POOL_STRING}" STREQUAL "")
    set(ADOLC_WITH_BOOST_ALLOCATOR TRUE)
    set(_additional_include_dirs OPTIONAL BOOST_INCLUDE_DIRS)
    set(_additional_library OPTIONAL BOOST_LIBRARIES)
  else()
    set(ADOLC_WITH_BOOST_ALLOCATOR FALSE)
    set(_additional_include_dirs)
    set(_additional_library)
  endif()
endif()


process_feature(ADOLC
  LIBRARIES
    REQUIRED ADOLC_LIBRARY
    ${_additional_library}
  INCLUDE_DIRS
    REQUIRED ADOLC_INCLUDE_DIR
    ${_additional_include_dirs}
  CLEAR ADOLC_INCLUDE_DIR ADOLC_LIBRARY ADOLC_SETTINGS_H
    ADOLC_DOUBLE_CAST_CHECK ADOLC_ADOUBLE_OSTREAM_CHECK # clean up checks in configure_adolc.cmake
  )
