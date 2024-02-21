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
