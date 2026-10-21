## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2023 - 2026 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

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

# Call find_package(VTK). In an unfriendly twist, the VTK package
# configuration *unsets* the MPI version variables that we have
# previously found and set. So we have to save them around the call to
# find_package().
set (_mpi_version ${MPI_VERSION})
set (_mpi_version_major ${MPI_VERSION_MAJOR})
set (_mpi_version_minor ${MPI_VERSION_MINOR})

# Only request the VTK modules we actually use.
#
# Without an explicit list of components VTK_LIBRARIES contains every
# module of the VTK installation, including the Python wrapping libraries
# that leave the Python symbols unresolved and thus lead to underlinkage.
#
find_package(VTK 9.0.0 QUIET HINTS ${VTK_DIR}
  COMPONENTS
    CommonCore CommonDataModel CommonExecutionModel FiltersCore IOCore IOLegacy IOXML
  OPTIONAL_COMPONENTS
    FiltersGeneral
  )

set (MPI_VERSION ${_mpi_version})
set (MPI_VERSION_MAJOR ${_mpi_version_major})
set (MPI_VERSION_MINOR ${_mpi_version_minor})


if(VTK_FOUND)
  set(VTK_VERSION "${VTK_VERSION}")
  set(VTK_MAJOR_VERSION "${VTK_MAJOR_VERSION}")
  set(VTK_MINOR_VERSION "${VTK_MINOR_VERSION}")

  set(VTK_INCLUDE_DIR
    ${VTK_PREFIX_PATH}/include/vtk-${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION})

  #
  # Collect the requested modules together with the VTK modules they depend on:
  #
  set(_targets)
  set(_queue ${VTK_LIBRARIES})
  while(_queue)
    list(POP_FRONT _queue _target)
    if(NOT TARGET ${_target} OR "${_target}" IN_LIST _targets)
      continue()
    endif()
    list(APPEND _targets ${_target})
    get_target_property(_dependencies ${_target} INTERFACE_LINK_LIBRARIES)
    if(_dependencies)
      foreach(_dependency ${_dependencies})
        if(_dependency MATCHES "^VTK::")
          list(APPEND _queue ${_dependency})
        endif()
      endforeach()
    endif()
  endwhile()

  # Try to find full paths from the collected targets.
  set(_libraries)
  foreach(_library ${_targets})
    get_target_property(_configurations ${_library} IMPORTED_CONFIGURATIONS)

    if(_configurations)
      foreach(_configuration ${_configurations})
        get_target_property(_imported_location ${_library} IMPORTED_LOCATION_${_configuration})
        list(APPEND _libraries ${_imported_location})
      endforeach()
    endif()
  endforeach()
endif()

process_feature(VTK
  LIBRARIES REQUIRED _libraries
  INCLUDE_DIRS REQUIRED VTK_INCLUDE_DIR
  CLEAR VTK_INCLUDE_DIR VTK_LIBRARIES _libraries _targets _queue
)
