## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2025 by the deal.II authors
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
# Try to find the NetCDF library.
#
# This module exports 
#     target netCDF::netcdf
#     NETCDF_FOUND
#     NETCDF_HAS_PARALLEL4
#

# not all distributions provide cmake config, and some (ubuntu 22) are broken
# so we just search for headers the old school way
find_path(_nc_include_dir
  NAMES netcdf_meta.h
  HINTS ${netCDF_DIR} ${NETCDF_DIR} ENV NETCDF_DIR ENV netCDF_DIR
  PATH_SUFFIXES netcdf include
  )
deal_ii_find_library(_nc_library
  NAMES netcdf
  HINTS ${netCDF_DIR} ${NETCDF_DIR} ENV NETCDF_DIR ENV netCDF_DIR
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

if(_nc_include_dir AND _nc_library)
  # netcdf found, define interface target
  set(netCDF_FOUND True)
  add_library(_netcdf INTERFACE)
  target_link_libraries(_netcdf INTERFACE ${_nc_library})
  target_include_directories(_netcdf INTERFACE ${_nc_include_dir})
  add_library(netCDF::netcdf ALIAS _netcdf)
  set(NETCDF_TARGETS "netCDF::netcdf")

  # scan features
  file(STRINGS "${_nc_include_dir}/netcdf_meta.h" 
    _nc_def_has_par 
    REGEX "#define[ \t]+NC_HAS_PARALLEL4[ \t]+1"
    )
  if(_nc_def_has_par)
    set(NETCDF_HAS_PARALLEL4 TRUE) # same variable as netcdf cmake config would define 
  else()
    set(NETCDF_HAS_PARALLEL4 FALSE)
  endif()
endif()

process_feature(NETCDF
  TARGETS
    REQUIRED NETCDF_TARGETS
  )
