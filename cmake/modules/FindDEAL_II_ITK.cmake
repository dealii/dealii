## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2026 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

set(ITK_DIR "" CACHE PATH "An optional hint to an ITK installation")
set_if_empty(ITK_DIR "$ENV{ITK_DIR}")

if(NOT "${ITK_DIR}" STREQUAL "")
  set(ITK_DIR ${ITK_DIR})
endif()

find_package(ITK 5.3.0 QUIET HINTS ${ITK_DIR})
if(ITK_FOUND)
  set(ITK_VERSION "${ITK_VERSION}")
  set(ITK_MAJOR_VERSION "${ITK_MAJOR_VERSION}")
  set(ITK_MINOR_VERSION "${ITK_MINOR_VERSION}")

  set(_libraries)
  foreach(_library ${ITK_LIBRARIES})
    # Ignore potential Python wrapping libraries
    if(NOT ${_library} MATCHES "Python")
      
      if(TARGET ${_library})
        get_target_property(_configurations ${_library} IMPORTED_CONFIGURATIONS)
        
        if(_configurations)
          foreach(_configuration ${_configurations})
            get_target_property(_imported_location ${_library} IMPORTED_LOCATION_${_configuration})
            if(_imported_location)
              list(APPEND _libraries ${_imported_location})
            endif()
          endforeach()
        else()
          # Fallback if there are no specific imported configurations
          get_target_property(_imported_location ${_library} IMPORTED_LOCATION)
          if(_imported_location)
            list(APPEND _libraries ${_imported_location})
          endif()
        endif()
      else()
        list(APPEND _libraries ${_library})
      endif()
    endif()
  endforeach()

  # Clean up any potential duplicates
  if(_libraries)
    list(REMOVE_DUPLICATES _libraries)
  endif()
endif()

process_feature(ITK
  LIBRARIES REQUIRED _libraries
  INCLUDE_DIRS REQUIRED ITK_INCLUDE_DIRS
  CLEAR ITK_INCLUDE_DIRS ITK_LIBRARIES _libraries
)
