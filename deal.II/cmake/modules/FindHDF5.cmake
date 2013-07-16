#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

#
# Try to find the HDF5 library
#
# This module exports
#
#   HDF5_LIBRARIES
#   HDF5_INCLUDE_DIRS
#   HDF5_WITH_MPI
#

INCLUDE(FindPackageHandleStandardArgs)

SET_IF_EMPTY(HDF5_DIR "$ENV{HDF5_DIR}")

FIND_PATH(HDF5_INCLUDE_DIR hdf5.h
  HINTS
    ${HDF5_DIR}
  PATH_SUFFIXES
    hdf5 hdf5/include include/hdf5 include
  )

FIND_LIBRARY(HDF5_LIBRARY NAMES hdf5
  HINTS
    ${HDF5_DIR}
  PATH_SUFFIXES
    hdf5/lib lib${LIB_SUFFIX} lib64 lib
  )

FIND_LIBRARY(HDF5_HL_LIBRARY NAMES hdf5_hl
  HINTS
    ${HDF5_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib
  )

SET(_output ${HDF5_HL_LIBRARY} ${HDF5_LIBRARY})
FIND_PACKAGE_HANDLE_STANDARD_ARGS(HDF5 DEFAULT_MSG
  _output # Cosmetic: Gives nice output
  HDF5_HL_LIBRARY
  HDF5_LIBRARY
  HDF5_INCLUDE_DIR
  )

MARK_AS_ADVANCED(
  HDF5_LIBRARY
  HDF5_HL_LIBRARY
  HDF5_INCLUDE_DIR
  )

IF(HDF5_FOUND)
  SET(HDF5_INCLUDE_DIRS
    ${HDF5_INCLUDE_DIR}
    )
  SET(HDF5_LIBRARIES
    ${HDF5_HL_LIBRARY}
    ${HDF5_LIBRARY}
    ${MPI_C_LIBRARIES} # for good measure
    )

  #
  # Is hdf5 compiled with support for mpi?
  #
  FILE(STRINGS "${HDF5_INCLUDE_DIR}/H5pubconf.h" HDF5_MPI_STRING
    REGEX "#define.*H5_HAVE_PARALLEL 1")
  IF("${HDF5_MPI_STRING}" STREQUAL "")
    SET(HDF5_WITH_MPI FALSE)
  ELSE()
    SET(HDF5_WITH_MPI TRUE)
  ENDIF()

  MARK_AS_ADVANCED(HDF5_DIR)
ELSE()
  SET(HDF5_DIR "" CACHE PATH
    "An optional hint to an hdf5 directory"
    )
ENDIF()

