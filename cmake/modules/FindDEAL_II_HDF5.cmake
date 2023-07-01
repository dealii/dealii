## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2023 by the deal.II authors
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
