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
# Try to find the NETCDF C and C++ libraries
#
# This module exports
#
#   NETCDF_LIBRARIES
#   NETCDF_INCLUDE_DIRS
#

SET(NETCDF_DIR "" CACHE PATH "An optional hint to a NETCDF installation")
SET_IF_EMPTY(NETCDF_DIR "$ENV{NETCDF_DIR}")

DEAL_II_FIND_PATH(NETCDF_INCLUDE_DIR netcdfcpp.h
  HINTS ${NETCDF_DIR}
  PATH_SUFFIXES netcdf include
  )

#
# TODO:
#
# - netcdf might externally depend on hdf5. Check and fix this.
# - separate C++ and C library search
#

DEAL_II_FIND_LIBRARY(NETCDF_CPLUSPLUS_LIBRARY NAMES netcdf_c++ netcdf_cpp
  HINTS ${NETCDF_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

DEAL_II_FIND_LIBRARY(NETCDF_C_LIBRARY NAMES netcdf
  HINTS ${NETCDF_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

DEAL_II_PACKAGE_HANDLE(NETCDF
  LIBRARIES REQUIRED NETCDF_CPLUSPLUS_LIBRARY NETCDF_C_LIBRARY
  INCLUDE_DIRS REQUIRED NETCDF_INCLUDE_DIR
  CLEAR NETCDF_CPLUSPLUS_LIBRARY NETCDF_C_LIBRARY NETCDF_INCLUDE_DIR
  )
