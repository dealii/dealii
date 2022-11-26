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
# Configuration for the hdf5 library:
#

set(FEATURE_HDF5_DEPENDS MPI)

macro(feature_hdf5_find_external var)
  find_package(DEAL_II_HDF5)

  if(HDF5_FOUND)
    set(${var} TRUE)

    if(NOT HDF5_IS_PARALLEL)
      message(STATUS "Insufficient hdf5 installation found: "
        "hdf5 has to be configured with MPI support."
        )
      set(HDF5_ADDITIONAL_ERROR_STRING
        "Insufficient hdf5 installation found!\n"
        "hdf5 has to be configured with MPI support.\n"
        )
      set(${var} FALSE)
    endif()

    check_mpi_interface(HDF5 ${var})
  endif()
endmacro()


configure_feature(HDF5)
