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
#     NETCDF_HAS_PARALLEL
#

find_package(netCDF 4.0.0 CONFIG QUIET)

if(NOT TARGET netCDF::netcdf)
  # not all distributions provide cmake config, so we fall back to searching files  
  find_path(_nc_include_dir
    NAMES netcdf_meta.h
    HINTS ${NETCDF_DIR} ENV NETCDF_DIR
    PATH_SUFFIXES netcdf include
    )
  deal_ii_find_library(_nc_library
    NAMES netcdf
    HINTS ${NETCDF_DIR} ENV NETCDF_DIR
    PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
    )

  if(_nc_include_dir AND _nc_library)
    # netcdf found, define interface target
    set(netCDF_FOUND True)
    add_library(_netcdf INTERFACE)
    target_link_libraries(_netcdf INTERFACE ${_nc_library})
    target_include_directories(_netcdf INTERFACE ${_nc_include_dir})
    add_library(netCDF::netcdf ALIAS _netcdf)

    # scan features
    file(STRINGS "${_nc_include_dir}/netcdf_meta.h" 
      _nc_def_has_par 
      REGEX "#define[ \t]+NC_HAS_PARALLEL[ \t]+1"
      )
    if(_nc_def_has_par)
      set(netCDF_HAS_PARALLEL TRUE) # same variable as netcdf cmake config would define 
    else()
      set(netCDF_HAS_PARALLEL FALSE)
    endif()
  endif()
endif()

if (TARGET netCDF::netcdf)
  set(netCDF_TARGETS "netCDF::netcdf")
else()
  unset(netCDF_TARGETS)
endif()

process_feature(NETCDF
  TARGETS
    REQUIRED netCDF_TARGETS
  )
set(NETCDF_HAS_PARALLEL ${netCDF_HAS_PARALLEL}) # uppercase for consistency with other features
