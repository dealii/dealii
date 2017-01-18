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
# This file implements the DEAL_II_SETUP_TARGET macro, which is
# part of the deal.II library.
#
# Usage:
#       DEAL_II_SETUP_TARGET(target)
#       DEAL_II_SETUP_TARGET(target DEBUG|RELEASE)
#
# This appends necessary include directories, linker flags, compile flags
# and compile definitions and the deal.II library link interface to the
# given target. In particular:
#
# INCLUDE_DIRECTORIES is appended with
#   "${DEAL_II_INCLUDE_DIRS}"
#
# COMPILE_FLAGS is appended with
#   "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_<build type>}"
#
# LINK_FLAGS is appended with
#   "${DEAL_II_LINKER_FLAGS ${DEAL_II_LINKER_FLAGS_<build type>}"
#
# COMPILE_DEFINITIONS is appended with
#   "${DEAL_II_USER_DEFINITIONS};${DEAL_II_USER_DEFINITIONS_<build type>}"
#
# If no "DEBUG" or "RELEASE" keyword is specified after the target, the
# current CMAKE_BUILD_TYPE is used instead. A CMAKE_BUILD_TYPE "Debug" is
# equivalent to the DEBUG keyword, a CMAKE_BUILD_TYPE "Release" is
# equivalent to the RELEASE keyword.
#
# This macro throws a FATAL_ERROR in case no DEBUG/RELEASE keyword is set
# and the build type is different from "Debug", or "Release".
#
# If the requested build type is not available (e.g. DEBUG request but
# deal.II was compiled with release mode only), the macro throws a
# FATAL_ERROR.
#

MACRO(DEAL_II_SETUP_TARGET _target)

  IF(NOT DEAL_II_PROJECT_CONFIG_INCLUDED)
    MESSAGE(FATAL_ERROR
      "\nDEAL_II_SETUP_TARGET can only be called in external projects after "
      "the inclusion of deal.IIConfig.cmake. It is not intended for "
      "internal use.\n\n"
      )
  ENDIF()

  IF(NOT DEAL_II_TARGET_CONFIG_INCLUDED)
    INCLUDE(${DEAL_II_TARGET_CONFIG})
    SET(DEAL_II_TARGET_CONFIG_INCLUDED TRUE)
  ENDIF()

  # Necessary for setting INCLUDE_DIRECTORIES via SET_PROPERTY
  CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

  #
  # Set build type with the help of the specified keyword, or
  # CMAKE_BUILD_TYPE:
  #
  IF("${ARGN}" MATCHES "^(DEBUG|RELEASE)$")
    SET(_build "${ARGN}")
  ELSEIF("${ARGN}" STREQUAL "")
    IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
      SET(_build "DEBUG")
    ELSEIF(CMAKE_BUILD_TYPE STREQUAL "Release")
      SET(_build "RELEASE")
    ELSE()
      MESSAGE(FATAL_ERROR
        "\nDEAL_II_SETUP_TARGET cannot determine DEBUG, or RELEASE flavor "
        "for target. CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\" is neither "
        "equal to \"Debug\", nor \"Release\"\n"
        "Set CMAKE_BUILD_TYPE accordingly, or use an explicit annotation: "
        "  DEAL_II_SETUP_TARGET(<target> DEBUG|RELEASE)\n\n"
        )
    ENDIF()
  ELSE()
    MESSAGE(FATAL_ERROR
      "\nDEAL_II_SETUP_TARGET called with invalid second argument. "
      "Valid arguments are (empty), DEBUG, or RELEASE\n\n"
      )
  ENDIF()


  #
  # We can only append DEBUG link flags and compile definitions if deal.II
  # was built with the Debug or DebugRelease build type. So test for this:
  #
  IF("${_build}" STREQUAL "DEBUG" AND NOT DEAL_II_BUILD_TYPE MATCHES "Debug")
    SET(_build "RELEASE")
  ENDIF()

  SET_PROPERTY(TARGET ${_target} APPEND PROPERTY
    INCLUDE_DIRECTORIES "${DEAL_II_INCLUDE_DIRS}"
    )
  SET_PROPERTY(TARGET ${_target} APPEND_STRING PROPERTY
    COMPILE_FLAGS "${DEAL_II_CXX_FLAGS} ${DEAL_II_CXX_FLAGS_${_build}}"
    )
  SET_PROPERTY(TARGET ${_target} APPEND_STRING PROPERTY
    LINK_FLAGS " ${DEAL_II_LINKER_FLAGS} ${DEAL_II_LINKER_FLAGS_${_build}}"
    )
  SET_PROPERTY(TARGET ${_target} APPEND PROPERTY
    COMPILE_DEFINITIONS "${DEAL_II_USER_DEFINITIONS};${DEAL_II_USER_DEFINITIONS_${_build}}"
    )

  #
  # Set up the link interface:
  #
  GET_PROPERTY(_type TARGET ${_target} PROPERTY TYPE)
  IF(NOT "${_type}" STREQUAL "OBJECT_LIBRARY")
    TARGET_LINK_LIBRARIES(${_target} ${DEAL_II_TARGET_${_build}})
  ENDIF()

  #
  # If DEAL_II_STATIC_EXECUTABLE is set, switch the final link type to
  # static:
  #
  IF(DEAL_II_STATIC_EXECUTABLE)
    SET_PROPERTY(TARGET ${_target} PROPERTY
      LINK_SEARCH_END_STATIC TRUE
      )
  ENDIF()

ENDMACRO()
