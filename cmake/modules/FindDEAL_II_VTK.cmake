## ---------------------------------------------------------------------
##
## Copyright (C) 2023 by the deal.II authors
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

  # Try to find full paths from targets contained in VTK_LIBRARIES.
  set(_libraries)
  foreach(_library ${VTK_LIBRARIES})
    if(NOT ${_library} MATCHES "Python"
	AND NOT ${_library} MATCHES "MPI4Py")
      get_target_property(_configurations ${_library} IMPORTED_CONFIGURATIONS)

      if(_configurations)
	foreach(_configuration ${_configurations})
          get_target_property(_imported_location ${_library} IMPORTED_LOCATION_${_configuration})
          list(APPEND _libraries ${_imported_location})
	endforeach()
      endif()
    endif()
  endforeach()
endif()

process_feature(VTK
  LIBRARIES REQUIRED _libraries
  INCLUDE_DIRS REQUIRED VTK_INCLUDE_DIR
  CLEAR VTK_INCLUDE_DIR VTK_LIBRARIES _libraries
)
