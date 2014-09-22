## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2013 by the deal.II authors
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
#
# This appends necessary include directories, linker flags, compile
# definitions and the deal.II library link interface to the given target.
#
# The current CMAKE_BUILD_TYPE is respected by setting the appropriate
# debug or release variants (if available).
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
  # Append include directories, and build-type independent linker flags and
  # compile definitions
  #
  SET_PROPERTY(TARGET ${_target} APPEND PROPERTY
    INCLUDE_DIRECTORIES "${DEAL_II_INCLUDE_DIRS}"
    )
  SET_PROPERTY(TARGET ${_target} APPEND_STRING PROPERTY
    LINK_FLAGS " ${DEAL_II_LINKER_FLAGS}"
    )
  SET_PROPERTY(TARGET ${_target} APPEND PROPERTY
    COMPILE_DEFINITIONS "${DEAL_II_USER_DEFINITIONS}"
    )

  #
  # Append build type dependent flags and definitions.
  #

  #
  # For this we obey the behaviour of the "optimized" and "debug"
  # keywords and this is a bit tricky:
  #
  # If the global property DEBUG_CONFIGURATIONS is set all build types that
  # (case insensitively) match one of the listed build types is considered
  # a "debug" build. The rest is "optimized".
  #
  # Otherwise every build type that (case insensitively) matches "debug" is
  # considered a debug build.
  #
  GET_PROPERTY(_debug_configurations_set
    GLOBAL PROPERTY DEBUG_CONFIGURATIONS SET
    )
  IF(_debug_configurations_set)
    STRING(TOLOWER "${CMAKE_BUILD_TYPE}" _cmake_build_type)
    GET_PROPERTY(_debug_configurations
      GLOBAL PROPERTY DEBUG_CONFIGURATIONS
      )
    FOREACH(_debug_type ${_debug_configurations})
      STRING(TOLOWER "${_debug_type}" _debug_type)
      IF("${_cmake_build_type}" STREQUAL "${_debug_type}")
        SET(_on_debug_build TRUE)
        BREAK()
      ENDIF()
    ENDFOREACH()
  ELSE()
    STRING(TOLOWER "${CMAKE_BUILD_TYPE}" _cmake_build_type)
    IF("${_cmake_build_type}" STREQUAL "debug")
      SET(_on_debug_build TRUE)
    ENDIF()
  ENDIF()

  #
  # We can only append DEBUG link flags and compile definitions if deal.II
  # was build with Debug or DebugRelease build type. So test for this:
  #
  IF( _on_debug_build
      AND DEAL_II_BUILD_TYPE MATCHES "Debug" )
    SET_PROPERTY(TARGET ${_target} APPEND_STRING PROPERTY
      LINK_FLAGS " ${DEAL_II_LINKER_FLAGS_DEBUG}"
      )
    SET_PROPERTY(TARGET ${_target} APPEND PROPERTY
      COMPILE_DEFINITIONS "${DEAL_II_USER_DEFINITIONS_DEBUG}"
      )
  ELSE()
    SET_PROPERTY(TARGET ${_target} APPEND_STRING PROPERTY
      LINK_FLAGS " ${DEAL_II_LINKER_FLAGS_RELEASE}"
      )
    SET_PROPERTY(TARGET ${_target} APPEND PROPERTY
      COMPILE_DEFINITIONS "${DEAL_II_USER_DEFINITIONS_RELEASE}"
      )
  ENDIF()

  #
  # Set up the link interface:
  #
  GET_PROPERTY(_type TARGET ${_target} PROPERTY TYPE)
  IF(NOT "${_type}" STREQUAL "OBJECT_LIBRARY")
    TARGET_LINK_LIBRARIES(${_target} ${DEAL_II_TARGET})
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
