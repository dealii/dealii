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

SET(FEATURE_HDF5_DEPENDS MPI)

MACRO(FEATURE_HDF5_FIND_EXTERNAL var)
  FIND_PACKAGE(HDF5)

  IF(HDF5_FOUND)
    SET(${var} TRUE)

    IF(NOT HDF5_IS_PARALLEL)
      MESSAGE(STATUS "Insufficient hdf5 installation found: "
        "hdf5 has to be configured with MPI support."
        )
      SET(HDF5_ADDITIONAL_ERROR_STRING
        "Insufficient hdf5 installation found!\n"
        "hdf5 has to be configured with MPI support.\n"
        )
      SET(${var} FALSE)
    ENDIF()

    CHECK_MPI_INTERFACE(HDF5 ${var})
  ENDIF()
ENDMACRO()


CONFIGURE_FEATURE(HDF5)
