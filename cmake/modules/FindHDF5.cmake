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

SET(HDF5_DIR "" CACHE PATH "An optional hint to an hdf5 directory")
SET_IF_EMPTY(HDF5_DIR "$ENV{HDF5_DIR}")

# temporarily disable ${CMAKE_SOURCE_DIR}/cmake/modules for module lookup
LIST(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)

IF(NOT "${HDF5_DIR}" STREQUAL "")
  SET(HDF5_ROOT "${HDF5_DIR}")
ENDIF()

SET(HDF5_PREFER_PARALLEL TRUE)
FIND_PACKAGE(HDF5)

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules/)

SET(_include_dirs "${HDF5_INCLUDE_DIRS}")
SET(_libraries "${HDF5_LIBRARIES};${HDF5_HL_LIBRARIES}")

DEAL_II_PACKAGE_HANDLE(HDF5
  LIBRARIES
    REQUIRED _libraries
    OPTIONAL MPI_C_LIBRARIES
  INCLUDE_DIRS
    REQUIRED _include_dirs
  USER_INCLUDE_DIRS
    REQUIRED _include_dirs
  CLEAR
    HDF5_C_COMPILER_EXECUTABLE HDF5_C_LIBRARY_dl HDF5_C_LIBRARY_hdf5
    HDF5_C_LIBRARY_m HDF5_C_LIBRARY_mpi HDF5_C_LIBRARY_sz HDF5_C_LIBRARY_z
    HDF5_DIFF_EXECUTABLE HDF5_PUBCONF
  )
