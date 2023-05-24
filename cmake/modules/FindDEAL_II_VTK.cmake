## ---------------------------------------------------------------------
##
## Copyright (C) 2022 by the deal.II authors
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
# Try to find the VTK libraries
#
# This module exports:
#
# VTK_INCLUDE_DIR
# VTK_VERSION_MAJOR
# VTK_VERSION_MINOR
# VTK_LIBRARIES
#

set(VTK_DIR "" CACHE PATH "An optional hint to a VTK installation")
set_if_empty(VTK_DIR "$ENV{VTK_DIR}")

if(NOT "${VTK_DIR}" STREQUAL "")
  set(VTK_DIR ${VTK_DIR})
endif()

find_package(VTK 9.0.0 QUIET HINTS ${VTK_DIR})

if(VTK_FOUND)
  set(VTK_VERSION "${VTK_VERSION}")
  set(VTK_MAJOR_VERSION "${VTK_MAJOR_VERSION}")
  set(VTK_MINOR_VERSION "${VTK_MINOR_VERSION}")

  set(VTK_INCLUDE_DIR
      ${VTK_PREFIX_PATH}/include/vtk-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION})

  #
  # We'd like to have the full library names but the VTK package only
  # exports a list with namespace-qualified names, e.g. VTK::CommonCore.
  # So we check again for every lib and store the full path:
  #
  set(_libraries "")

  foreach(_component ${VTK_AVAILABLE_COMPONENTS})
    deal_ii_find_library(VTK_LIBRARY_${_component}
      NAMES libvtk${_component}-${VTK_VERSION_MAJOR}.${VTK_VERSION_MINOR}.so
      HINTS ${VTK_PREFIX_PATH}/lib
      NO_DEFAULT_PATH
      NO_CMAKE_ENVIRONMENT_PATH
      NO_CMAKE_PATH
      NO_SYSTEM_ENVIRONMENT_PATH
      NO_CMAKE_SYSTEM_PATH
      NO_CMAKE_FIND_ROOT_PATH
    )

    if (NOT "${VTK_LIBRARY_${_component}}" MATCHES "-NOTFOUND")
      list(APPEND _libraries VTK_LIBRARY_${_component})
    else()
      # VTK_AVAILABLE_COMPONENTS contains also header-only modules.
      # If the library has not been found, check if corresponding
      # headers exist.
      if (NOT EXISTS "${VTK_INCLUDE_DIR}/${_component}" AND
          NOT EXISTS "${VTK_INCLUDE_DIR}/vtk_${_component}.h")
        message(FATAL_ERROR "VTK: component ${_component} not found.")
      endif()
    endif()
  endforeach()
endif()

process_feature(VTK
  LIBRARIES REQUIRED ${_libraries}
  INCLUDE_DIRS REQUIRED VTK_INCLUDE_DIR
  CLEAR VTK_INCLUDE_DIR VTK_LIBRARIES ${_libraries}
)
