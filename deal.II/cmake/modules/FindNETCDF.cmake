#####
##
## Copyright (C) 2012, 2013 by the deal.II authors
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
# Try to find the NETCDF C and C++ libraries
#
# This module exports
#
#   NETCDF_LIBRARIES
#   NETCDF_INCLUDE_DIRS
#

INCLUDE(FindPackageHandleStandardArgs)

SET_IF_EMPTY(NETCDF_DIR "$ENV{NETCDF_DIR}")

FIND_PATH(NETCDF_INCLUDE_DIR netcdf.hh
  HINTS
    ${NETCDF_DIR}
  PATH_SUFFIXES
    netcdf include
  )

#
# TODO:
#
# - netcdf might externally depend on hdf5. Check and fix this.
# - separate C++ and C library search
#

FIND_LIBRARY(NETCDF_CPLUSPLUS_LIBRARY NAMES netcdf_c++ netcdf_cpp
  HINTS
    ${NETCDF_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib
  )

FIND_LIBRARY(NETCDF_C_LIBRARY NAMES netcdf
  HINTS
    ${NETCDF_DIR}
  PATH_SUFFIXES
    lib${LIB_SUFFIX} lib64 lib
  )

SET(_output ${NETCDF_CPLUSPLUS_LIBRARY} ${NETCDF_C_LIBRARY})
FIND_PACKAGE_HANDLE_STANDARD_ARGS(NETCDF DEFAULT_MSG
  _output # Cosmetic: Gives nice output
  NETCDF_CPLUSPLUS_LIBRARY
  NETCDF_C_LIBRARY
  NETCDF_INCLUDE_DIR
  )

IF(NETCDF_FOUND)
  SET(NETCDF_INCLUDE_DIRS
    ${NETCDF_INCLUDE_DIR}
    )
  SET(NETCDF_LIBRARIES
    ${NETCDF_CPLUSPLUS_LIBRARY}
    ${NETCDF_C_LIBRARY}
    )

  MARK_AS_ADVANCED(
    NETCDF_CPLUSPLUS_LIBRARY
    NETCDF_C_LIBRARY
    NETCDF_INCLUDE_DIR
    NETCDF_DIR
  )
ELSE()
  SET(NETCDF_DIR "" CACHE PATH
    "An optional hint to a NETCDF installation"
    )
ENDIF()
