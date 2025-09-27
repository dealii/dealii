## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2022 by the deal.II authors
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
# This file implements the DEAL_II_INITIALIZE_VARIABLES macro, which is
# part of the deal.II library.
#
# Usage:
#       deal_ii_initialize_cached_variables()
#
# This sets some cached variables to the values used for compiling the
# deal.II library.
#
# This macro has to be called before project()!
#

macro(deal_ii_initialize_cached_variables)

  if(NOT DEAL_II_PROJECT_CONFIG_INCLUDED)
    message(FATAL_ERROR
      "\nDEAL_II_INITIALIZE_CACHED_VARIABLES can only be called in external "
      "projects after the inclusion of deal.IIConfig.cmake. It is not "
      "intended for internal use.\n\n"
      )
  endif()

  #
  # Set build type according to available libraries
  #
  if(DEAL_II_BUILD_TYPE MATCHES "Debug")
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
      "Choose the type of build, options are: Debug, Release"
      )
  else()
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug, Release"
      )
  endif()

  #
  # Reset build type if unsupported, i.e. if it is not Debug, Release, or
  # DebugRelease, or if the library doesn't support it
  #
  if( NOT "${CMAKE_BUILD_TYPE}" MATCHES "^(Debug|Release|DebugRelease)$"
      OR NOT "${DEAL_II_BUILD_TYPE}" MATCHES "${CMAKE_BUILD_TYPE}" )

    if("${DEAL_II_BUILD_TYPE}" STREQUAL "DebugRelease")
      set(_new_build_type "Debug")
    else()
      set(_new_build_type "${DEAL_II_BUILD_TYPE}")
    endif()

    message(
"###
#
#  WARNING:
#
#  CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\" unsupported by current installation!
#  deal.II was configured with \"${DEAL_II_BUILD_TYPE}\".
#
#  CMAKE_BUILD_TYPE was forced to \"${_new_build_type}\".
#
###"
      )
    set(CMAKE_BUILD_TYPE "${_new_build_type}" CACHE STRING
      "Choose the type of build, options are: Debug, Release"
      FORCE
      )

  endif()


  set(CMAKE_CXX_COMPILER ${DEAL_II_CXX_COMPILER} CACHE STRING
    "CXX Compiler.")
  set(CMAKE_CXX_FLAGS "" CACHE STRING
    "Flags used by the compiler during all build types."
    )
  set(CMAKE_CXX_FLAGS_DEBUG "" CACHE STRING
    "Flags used by the compiler during debug builds."
    )
  set(CMAKE_CXX_FLAGS_RELEASE "" CACHE STRING
    "Flags used by the compiler during release builds."
    )


  mark_as_advanced(CMAKE_INSTALL_PREFIX)

endmacro()
