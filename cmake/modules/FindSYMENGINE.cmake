## ---------------------------------------------------------------------
##
## Copyright (C) 2019 - 2022 by the deal.II authors
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
# - Try to find SymEngine
#
# This module exports
#
#   SYMENGINE_INCLUDE_DIR
#   SYMENGINE_LIBRARY
#   SYMENGINE_WITH_LLVM
#

SET(SYMENGINE_DIR "" CACHE PATH "An optional hint to a SymEngine installation")
SET_IF_EMPTY(SYMENGINE_DIR "$ENV{SYMENGINE_DIR}")

#
# SymEngine overwrites the CMake module path, so we save
# and restore it after this library is found and configured.
#
SET (DEAL_II_CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH})

#
# Include the SymEngine:
#
FIND_PACKAGE(SymEngine
  CONFIG QUIET
  HINTS ${SYMENGINE_DIR}
  PATH_SUFFIXES lib/cmake/symengine
  NO_SYSTEM_ENVIRONMENT_PATH
  )

#
# Reset the CMake module path
#
SET (CMAKE_MODULE_PATH ${DEAL_II_CMAKE_MODULE_PATH})


#
# Look for symengine_config.h - we'll query it to determine supported features:
#
IF(SymEngine_FOUND)
  DEAL_II_FIND_FILE(SYMENGINE_SETTINGS_H symengine_config.h
    HINTS ${SYMENGINE_INCLUDE_DIRS}
    PATH_SUFFIXES symengine
    NO_DEFAULT_PATH
    NO_CMAKE_ENVIRONMENT_PATH
    NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH
    NO_CMAKE_FIND_ROOT_PATH
    )
ENDIF()

#
# Version check
#
IF(EXISTS ${SYMENGINE_SETTINGS_H})

  FILE(STRINGS "${SYMENGINE_SETTINGS_H}" SYMENGINE_VERSION_MAJOR_STRING
    REGEX "#define.*SYMENGINE_MAJOR_VERSION")
  STRING(REGEX REPLACE "^.*SYMENGINE_MAJOR_VERSION.*([0-9]+).*" "\\1"
    SYMENGINE_VERSION_MAJOR "${SYMENGINE_VERSION_MAJOR_STRING}"
    )
  FILE(STRINGS "${SYMENGINE_SETTINGS_H}" SYMENGINE_VERSION_MINOR_STRING
    REGEX "#define.*SYMENGINE_MINOR_VERSION")
  STRING(REGEX REPLACE "^.*SYMENGINE_MINOR_VERSION.*([0-9]+).*" "\\1"
    SYMENGINE_VERSION_MINOR "${SYMENGINE_VERSION_MINOR_STRING}"
    )
  FILE(STRINGS "${SYMENGINE_SETTINGS_H}" SYMENGINE_VERSION_PATCH_STRING
    REGEX "#define.*SYMENGINE_PATCH_VERSION")
  STRING(REGEX REPLACE "^.*SYMENGINE_PATCH_VERSION.*([0-9]+).*" "\\1"
    SYMENGINE_VERSION_PATCH "${SYMENGINE_VERSION_PATCH_STRING}"
    )
    
  SET(SYMENGINE_VERSION ${SymEngine_VERSION})
ENDIF()

#
# Feature checks
#

MACRO(_symengine_feature_check _var _regex)
  IF(EXISTS ${SYMENGINE_SETTINGS_H})
    FILE(STRINGS "${SYMENGINE_SETTINGS_H}" SYMENGINE_${_var}_STRING
      REGEX "${_regex}")
    IF("${SYMENGINE_${_var}_STRING}" STREQUAL "")
      SET(SYMENGINE_WITH_${_var} FALSE)
    ELSE()
      SET(SYMENGINE_WITH_${_var} TRUE)
    ENDIF()
  ENDIF()
ENDMACRO()

# Other possible features of interest: BOOST, GMP
_symengine_feature_check(LLVM "#define.*HAVE_SYMENGINE_LLVM")

#
# Sanitize include dirs:
#

STRING(REGEX REPLACE
  "(lib64|lib)\\/cmake\\/symengine\\/\\.\\.\\/\\.\\.\\/\\.\\.\\/" ""
  SYMENGINE_INCLUDE_DIRS  "${SYMENGINE_INCLUDE_DIRS}"
  )
REMOVE_DUPLICATES(SYMENGINE_INCLUDE_DIRS)

#
# Get the full path for the SYMENGINE_LIBRARIES. Some of these libraries are
# CMake targets, so we can query them directly for this information.
#
FOREACH(SYMENGINE_LIBRARY_NAME ${SYMENGINE_LIBRARIES})
   IF (TARGET ${SYMENGINE_LIBRARY_NAME})
       GET_PROPERTY(SYMENGINE_LIBRARY TARGET ${SYMENGINE_LIBRARY_NAME} PROPERTY LOCATION)
   ELSE ()
       SET(SYMENGINE_LIBRARY ${SYMENGINE_LIBRARY_NAME})
   ENDIF()

  SET(_symengine_libraries ${_symengine_libraries} ${SYMENGINE_LIBRARY})
ENDFOREACH()
SET(SYMENGINE_LIBRARIES ${_symengine_libraries})

DEAL_II_PACKAGE_HANDLE(SYMENGINE
  LIBRARIES REQUIRED SYMENGINE_LIBRARIES
  INCLUDE_DIRS REQUIRED SYMENGINE_INCLUDE_DIRS
  USER_INCLUDE_DIRS REQUIRED SYMENGINE_INCLUDE_DIRS
  CLEAR SYMENGINE_SETTINGS_H SYMENGINE_SKIP_DEPENDENCIES SymEngine_DIR
)
