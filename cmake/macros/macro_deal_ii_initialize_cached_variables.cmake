## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2015 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# This file implements the DEAL_II_INITIALIZE_VARIABLES macro, which is
# part of the deal.II library.
#
# Usage:
#       DEAL_II_INITIALIZE_CACHED_VARIABLES()
#
# This sets some cached variables to the values used for compiling the
# deal.II library.
#
# This macro has to be called before PROJECT()!
#

MACRO(DEAL_II_INITIALIZE_CACHED_VARIABLES)

  IF(NOT DEAL_II_PROJECT_CONFIG_INCLUDED)
    MESSAGE(FATAL_ERROR
      "\nDEAL_II_INITIALIZE_CACHED_VARIABLES can only be called in external "
      "projects after the inclusion of deal.IIConfig.cmake. It is not "
      "intended for internal use.\n\n"
      )
  ENDIF()

  #
  # Set build type according to build type of deal.II
  #
  SET(CMAKE_BUILD_TYPE ${DEAL_II_BUILD_TYPE} CACHE STRING
    "Choose the type of build, options are: Debug, Release, DebugRelease")

  #
  # Reset build type if unsupported, i.e. if it is not (case insensitively
  # equal to Debug or Release or unsupported by the current build type:
  #
  STRING(TOLOWER "${CMAKE_BUILD_TYPE}" _cmake_build_type)

  IF(NOT "${_cmake_build_type}" MATCHES "^(debug|release|debugrelease)$")


    MESSAGE(
"###
#
#  WARNING:
#
#  CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\" unsupported by current installation!
#  deal.II was built with CMAKE_BUILD_TYPE \"${DEAL_II_BUILD_TYPE}\".
#
#  CMAKE_BUILD_TYPE is forced to \"${DEAL_II_BUILD_TYPE}\".
#
###"
      )
    SET(CMAKE_BUILD_TYPE ${DEAL_II_BUILD_TYPE} CACHE STRING
      "Choose the type of build, options are: Debug, Release, DebugRelease"
      FORCE
      )

  ENDIF()


  SET(CMAKE_CXX_COMPILER ${DEAL_II_CXX_COMPILER} CACHE STRING
    "CXX Compiler.")

  SET(CMAKE_C_COMPILER ${DEAL_II_C_COMPILER} CACHE STRING
    "C Compiler.")

  SET(CMAKE_CXX_FLAGS "" CACHE STRING
    "Flags used by the compiler during all build types."
    )

  SET(CMAKE_CXX_FLAGS_DEBUG "" CACHE STRING
    "Flags used by the compiler during debug builds."
    )

  SET(CMAKE_CXX_FLAGS_RELEASE "" CACHE STRING
    "Flags used by the compiler during release builds."
    )

  MARK_AS_ADVANCED(CMAKE_INSTALL_PREFIX)

ENDMACRO()

