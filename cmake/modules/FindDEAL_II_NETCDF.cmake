## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2023 by the deal.II authors
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
# Try to find the NetCDF library using nc-config executable.
#
# This module exports
#
#   NETCDF_LIBRARIES
#   NETCDF_INCLUDE_DIRS
#   NETCDF_HAS_PARALLEL
#

find_program(DEAL_II_NETCDF_CONFIG_PROGRAM nc-config)
if (DEAL_II_NETCDF_CONFIG_PROGRAM)
  execute_process(
    COMMAND ${DEAL_II_NETCDF_CONFIG_PROGRAM} "--has-nc4" 
    OUTPUT_VARIABLE _has_nc4)
  if(NOT _has_nc4)
    message(WARNING "Ignoring the discovered NetCDF library because it does not support format version 4.")
  else()
    execute_process(
      COMMAND ${DEAL_II_NETCDF_CONFIG_PROGRAM} "--includedir" 
      OUTPUT_VARIABLE _netcdf_include_dirs)
    string(STRIP "${_netcdf_include_dirs}" _netcdf_include_dirs)
    execute_process(
      COMMAND ${DEAL_II_NETCDF_CONFIG_PROGRAM} "--libs" 
      OUTPUT_VARIABLE _netcdf_libraries)
    string(STRIP "${_netcdf_libraries}" _netcdf_libraries)
    execute_process(
      COMMAND ${DEAL_II_NETCDF_CONFIG_PROGRAM} "--has-parallel" 
      OUTPUT_VARIABLE _netcdf_has_parallel)
  endif()
endif()

process_feature(NETCDF
  LIBRARIES
    REQUIRED _netcdf_libraries
  INCLUDE_DIRS
    REQUIRED _netcdf_include_dirs
  )
set(NETCDF_HAS_PARALLEL ${_netcdf_has_parallel})
