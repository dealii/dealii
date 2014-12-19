## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2014 by the deal.II authors
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
# Try to find the HDF5 library
#
# This module exports
#
#   HDF5_LIBRARIES
#   HDF5_INCLUDE_DIRS
#   HDF5_WITH_MPI
#

SET(HDF5_DIR "" CACHE PATH "An optional hint to an hdf5 directory")
SET_IF_EMPTY(HDF5_DIR "$ENV{HDF5_DIR}")

DEAL_II_FIND_PATH(HDF5_INCLUDE_DIR hdf5.h
  HINTS ${HDF5_DIR}
  PATH_SUFFIXES hdf5 hdf5/include include/hdf5 include
  )

DEAL_II_FIND_LIBRARY(HDF5_LIBRARY NAMES hdf5
  HINTS ${HDF5_DIR}
  PATH_SUFFIXES hdf5/lib lib${LIB_SUFFIX} lib64 lib
  )

DEAL_II_FIND_LIBRARY(HDF5_HL_LIBRARY NAMES hdf5_hl
  HINTS ${HDF5_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

DEAL_II_FIND_FILE(HDF5_PUBCONF NAMES H5pubconf.h H5pubconf-64.h
  HINTS ${HDF5_INCLUDE_DIR} ${HDF5_DIR}
  PATH_SUFFIXES hdf5 hdf5/include include/hdf5 include
  )

IF(EXISTS ${HDF5_PUBCONF})
  #
  # Is hdf5 compiled with support for mpi?
  #
  FILE(STRINGS ${HDF5_PUBCONF} HDF5_MPI_STRING
    REGEX "#define.*H5_HAVE_PARALLEL 1"
    )
  IF("${HDF5_MPI_STRING}" STREQUAL "")
    SET(HDF5_WITH_MPI FALSE)
  ELSE()
    SET(HDF5_WITH_MPI TRUE)
  ENDIF()
ENDIF()

DEAL_II_PACKAGE_HANDLE(HDF5
  LIBRARIES
    REQUIRED HDF5_HL_LIBRARY HDF5_LIBRARY
    OPTIONAL MPI_C_LIBRARIES
  INCLUDE_DIRS
    REQUIRED HDF5_INCLUDE_DIR
  USER_INCLUDE_DIRS
    REQUIRED HDF5_INCLUDE_DIR
  CLEAR HDF5_HL_LIBRARY HDF5_LIBRARY HDF5_INCLUDE_DIR HDF5_PUBCONF
  )
