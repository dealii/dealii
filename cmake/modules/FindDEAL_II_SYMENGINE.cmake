## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2019 - 2022 by the deal.II authors
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
# - Try to find SymEngine
#
# This module exports
#
#   SYMENGINE_INCLUDE_DIR
#   SYMENGINE_LIBRARY
#   SYMENGINE_WITH_LLVM
#

set(SYMENGINE_DIR "" CACHE PATH "An optional hint to a SymEngine installation")
set_if_empty(SYMENGINE_DIR "$ENV{SYMENGINE_DIR}")

#
# SymEngine overwrites the CMake module path, so we save
# and restore it after this library is found and configured.
#
set(_cmake_module_path ${CMAKE_MODULE_PATH})

#
# Include the SymEngine:
#
find_package(SymEngine
  CONFIG QUIET
  HINTS ${SYMENGINE_DIR}
  PATH_SUFFIXES lib/cmake/symengine
  NO_SYSTEM_ENVIRONMENT_PATH
  )

#
# Reset the CMake module path
#
set(CMAKE_MODULE_PATH ${_cmake_module_path})


#
# Look for symengine_config.h - we'll query it to determine supported features:
#
if(SymEngine_FOUND)
  deal_ii_find_file(SYMENGINE_SETTINGS_H symengine_config.h
    HINTS ${SYMENGINE_INCLUDE_DIRS}
    PATH_SUFFIXES symengine
    NO_DEFAULT_PATH
    NO_CMAKE_ENVIRONMENT_PATH
    NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH
    NO_CMAKE_FIND_ROOT_PATH
    )
endif()

#
# Version check
#
if(EXISTS ${SYMENGINE_SETTINGS_H})

  file(STRINGS "${SYMENGINE_SETTINGS_H}" SYMENGINE_VERSION_MAJOR_STRING
    REGEX "#define.*SYMENGINE_MAJOR_VERSION")
  string(REGEX REPLACE "^.*SYMENGINE_MAJOR_VERSION.*([0-9]+).*" "\\1"
    SYMENGINE_VERSION_MAJOR "${SYMENGINE_VERSION_MAJOR_STRING}"
    )
  file(STRINGS "${SYMENGINE_SETTINGS_H}" SYMENGINE_VERSION_MINOR_STRING
    REGEX "#define.*SYMENGINE_MINOR_VERSION")
  string(REGEX REPLACE "^.*SYMENGINE_MINOR_VERSION.*([0-9]+).*" "\\1"
    SYMENGINE_VERSION_MINOR "${SYMENGINE_VERSION_MINOR_STRING}"
    )
  file(STRINGS "${SYMENGINE_SETTINGS_H}" SYMENGINE_VERSION_PATCH_STRING
    REGEX "#define.*SYMENGINE_PATCH_VERSION")
  string(REGEX REPLACE "^.*SYMENGINE_PATCH_VERSION.*([0-9]+).*" "\\1"
    SYMENGINE_VERSION_PATCH "${SYMENGINE_VERSION_PATCH_STRING}"
    )
    
  set(SYMENGINE_VERSION ${SymEngine_VERSION})
endif()

#
# Feature checks
#

macro(_symengine_feature_check _var _regex)
  if(EXISTS ${SYMENGINE_SETTINGS_H})
    file(STRINGS "${SYMENGINE_SETTINGS_H}" SYMENGINE_${_var}_STRING
      REGEX "${_regex}")
    if("${SYMENGINE_${_var}_STRING}" STREQUAL "")
      set(SYMENGINE_WITH_${_var} FALSE)
    else()
      set(SYMENGINE_WITH_${_var} TRUE)
    endif()
  endif()
endmacro()

# Other possible features of interest: BOOST, GMP
_symengine_feature_check(LLVM "#define.*HAVE_SYMENGINE_LLVM")

#
# Sanitize include dirs:
#

string(REGEX REPLACE
  "(lib64|lib)\\/cmake\\/symengine\\/\\.\\.\\/\\.\\.\\/\\.\\.\\/" ""
  SYMENGINE_INCLUDE_DIRS  "${SYMENGINE_INCLUDE_DIRS}"
  )
remove_duplicates(SYMENGINE_INCLUDE_DIRS)

#
# Get the full path for the SYMENGINE_LIBRARIES. Some of these libraries are
# CMake targets, so we can query them directly for this information.
#
foreach(SYMENGINE_LIBRARY_NAME ${SYMENGINE_LIBRARIES})
   if (TARGET ${SYMENGINE_LIBRARY_NAME})
       get_property(SYMENGINE_LIBRARY TARGET ${SYMENGINE_LIBRARY_NAME} PROPERTY LOCATION)
   else ()
       set(SYMENGINE_LIBRARY ${SYMENGINE_LIBRARY_NAME})
   endif()

  set(_symengine_libraries ${_symengine_libraries} ${SYMENGINE_LIBRARY})
endforeach()
set(SYMENGINE_LIBRARIES ${_symengine_libraries})

process_feature(SYMENGINE
  LIBRARIES REQUIRED SYMENGINE_LIBRARIES
  INCLUDE_DIRS REQUIRED SYMENGINE_INCLUDE_DIRS
  CLEAR SYMENGINE_SETTINGS_H SYMENGINE_SKIP_DEPENDENCIES SymEngine_DIR
)
