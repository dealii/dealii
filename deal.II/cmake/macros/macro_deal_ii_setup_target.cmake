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
# This file implements the DEAL_II_SETUP_TARGET macro, which is
# part of the deal.II library.
#
# Usage:
#       DEAL_II_SETUP_TARGET(target)
#
# This sets the necessary include directories, linker flags, compile
# definitions and the deal.II library the target will be linked against.
#
#

MACRO(DEAL_II_SETUP_TARGET target)

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

  # Necessary for setting INCLUDE_DIRECTORIES via SET_TARGET_PROPERTIES
  CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

  GET_DIRECTORY_PROPERTY(inc_dirs_so_far INCLUDE_DIRECTORIES)
  SET_TARGET_PROPERTIES(${target} PROPERTIES
    INCLUDE_DIRECTORIES
      "${inc_dirs_so_far};${DEAL_II_EXTERNAL_INCLUDE_DIRS};${DEAL_II_INCLUDE_DIRS}"
    LINK_FLAGS
      "${DEAL_II_LINKER_FLAGS}"
    COMPILE_DEFINITIONS
      "${DEAL_II_USER_DEFINITIONS}"
    )

  #
  # Set build type dependend flags and definitions:
  #
  FOREACH(build ${DEAL_II_BUILD_TYPES})
    SET_TARGET_PROPERTIES(${target} PROPERTIES
      LINK_FLAGS_${build}
        "${DEAL_II_LINKER_FLAGS_${build}}"
      COMPILE_DEFINITIONS_${build}
        "${DEAL_II_USER_DEFINITIONS_${build}}"
      )

  ENDFOREACH()

  #
  # Link againgst the correct deal.II library target for the current
  # CMAKE_BUILD_TYPE:
  #
  STRING(TOUPPER "${CMAKE_BUILD_TYPE}" build)
  IF(NOT "${DEAL_II_TARGET_${build}}" STREQUAL "")
    TARGET_LINK_LIBRARIES(${target}
      ${DEAL_II_TARGET_${build}}
      )
  ELSE()
    MESSAGE(FATAL_ERROR
      "\nNo matching deal.II library target for current build type: "
      "\"${CMAKE_BUILD_TYPE}\"\n"
      "Candidates are (case insensitive): ${DEAL_II_BUILD_TYPES}\n\n"
      )
  ENDIF()

ENDMACRO()

