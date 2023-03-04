## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2022 by the deal.II authors
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
set(_libraries "hdf5") #${HDF5_LIBRARIES};${HDF5_HL_LIBRARIES}")

#
# We'd like to have the full library names but the HDF5 package only
# exports a list with short names. So check again for every lib and store
# the full path:
#
#set(_libraries "")
#foreach(_library hdf5 ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
#  list(APPEND _libraries HDF5_LIBRARY_${_library})
#  deal_ii_find_library(HDF5_LIBRARY_${_library}
#    NAMES ${_library}
#    HINTS ${HDF5_INSTALL_LIBRARY_DIR}
#    NO_DEFAULT_PATH
#    NO_CMAKE_ENVIRONMENT_PATH
#    NO_CMAKE_PATH
#    NO_SYSTEM_ENVIRONMENT_PATH
#    NO_CMAKE_SYSTEM_PATH
#    NO_CMAKE_FIND_ROOT_PATH
#    )
#endforeach()
bla
MESSAGE(FATAL_ERROR "${_libraries}")

process_feature(HDF5
  LIBRARIES
    REQUIRED ${_libraries}
    OPTIONAL MPI_C_LIBRARIES
  INCLUDE_DIRS
    REQUIRED _include_dirs
  CLEAR
    HDF5_C_COMPILER_EXECUTABLE HDF5_C_LIBRARY_dl HDF5_C_LIBRARY_hdf5
    HDF5_C_LIBRARY_m HDF5_C_LIBRARY_mpi HDF5_C_LIBRARY_sz HDF5_C_LIBRARY_z
    HDF5_DIFF_EXECUTABLE HDF5_PUBCONF
  )
