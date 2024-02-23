## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2023 by the deal.II authors
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
# Try to find the HDF5 library
#
# This module exports
#
#   HDF5_LIBRARIES
#   HDF5_INCLUDE_DIRS
#   HDF5_IS_PARALLEL
#

set(HDF5_DIR "" CACHE PATH "An optional hint to an hdf5 directory")
set_if_empty(HDF5_DIR "$ENV{HDF5_DIR}")

if( NOT "${HDF5_DIR}" STREQUAL "" AND
    NOT "${HDF5_DIR}" STREQUAL "HDF5_DIR-NOTFOUND" )
  set(HDF5_ROOT "${HDF5_DIR}")
endif()

set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5)

set(_include_dirs "${HDF5_INCLUDE_DIRS}")
set(_libraries_tmp "${HDF5_LIBRARIES};${HDF5_HL_LIBRARIES}")

# HDF5_LIBRARIES and HDF5_HL_LIBRARIES might contain targets or full paths to libraries
# try to find full paths in the former case
set(_libraries)
foreach(_library ${_libraries_tmp})
  if(TARGET ${_library})
    get_target_property(_configurations ${_library} IMPORTED_CONFIGURATIONS)
    if(_configurations)
      foreach(_configuration ${_configurations})
        get_target_property(_imported_location ${_library} IMPORTED_LOCATION_${_configuration})
        list(APPEND _libraries ${_imported_location})
      endforeach()
    else()
      get_target_property(_imported_location ${_library} IMPORTED_LOCATION)
      list(APPEND _libraries ${_imported_location})
    endif()
  else()
    list(APPEND _libraries ${_library})
  endif()
endforeach()

process_feature(HDF5
  LIBRARIES
    REQUIRED _libraries
    OPTIONAL MPI_C_LIBRARIES
  INCLUDE_DIRS
    REQUIRED _include_dirs
  CLEAR
    HDF5_C_COMPILER_EXECUTABLE HDF5_C_LIBRARY_dl HDF5_C_LIBRARY_hdf5
    HDF5_C_LIBRARY_m HDF5_C_LIBRARY_mpi HDF5_C_LIBRARY_sz HDF5_C_LIBRARY_z
    HDF5_DIFF_EXECUTABLE HDF5_PUBCONF
  )
