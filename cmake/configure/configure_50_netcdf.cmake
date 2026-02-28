## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2022 by the deal.II authors
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
# Configuration for the netcdf library:
#

set(FEATURE_NETCDF_DEPENDS MPI)

macro(feature_netcdf_find_external var)
  find_package(DEAL_II_NETCDF)
  
  if(NETCDF_FOUND)
    set(${var} TRUE)
    
    if(NOT NETCDF_HAS_PARALLEL)
      set(${var} FALSE)
      message(STATUS "Insufficient netCDF installation found: "
        "netCDF has to be configured with MPI support."
        )
      set(NETCDF_ADDITIONAL_ERROR_STRING
        "Insufficient netCDF installation found!\n"
        "netCDF has to be configured with MPI support.\n"
        )
    endif()
  
    check_mpi_interface(NETCDF ${var})
  endif()
endmacro()

configure_feature(NETCDF)
