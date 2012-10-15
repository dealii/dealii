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
#

INCLUDE(FindPackageHandleStandardArgs)

FIND_PATH(HDF5_INCLUDE_DIR hdf5.h
  HINTS
    ${HDF5_DIR}
  PATH_SUFFIXES
    hdf5 include/hdf5 include
  )

FIND_LIBRARY(HDF5_LIBRARY
  NAMES hdf5
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
)

FIND_LIBRARY(HDF5_HL_LIBRARY
  NAMES hdf5_hl
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(HDF5 DEFAULT_MSG
  HDF5_INCLUDE_DIR
  HDF5_LIBRARY
  HDF5_HL_LIBRARY
  )

SET(HDF5_INCLUDE_DIRS
  ${HDF5_INCLUDE_DIR}
  )

SET(HDF5_LIBRARIES
  ${HDF5_LIBRARY}
  ${HDF5_HL_LIBRARY}
  )

IF(HDF5_FOUND)
  MARK_AS_ADVANCED(
  HDF5_INCLUDE_DIR
  HDF5_LIBRARY
  HDF5_HL_LIBRARY
  )
ENDIF()

