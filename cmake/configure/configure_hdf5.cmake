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
# Configuration for the hdf5 library:
#

SET(FEATURE_HDF5_AFTER MPI)


MACRO(FEATURE_HDF5_FIND_EXTERNAL var)
  FIND_PACKAGE(HDF5)

  IF(HDF5_FOUND)
    SET(${var} TRUE)

    IF( (HDF5_WITH_MPI AND NOT DEAL_II_WITH_MPI) OR
        (NOT HDF5_WITH_MPI AND DEAL_II_WITH_MPI) )
      MESSAGE(STATUS "Insufficient hdf5 installation found: "
        "hdf5 has to be configured with the same MPI configuration as deal.II."
        )
      SET(HDF5_ADDITIONAL_ERROR_STRING
        "Insufficient hdf5 installation found!\n"
        "hdf5 has to be configured with the same MPI configuration as deal.II, but found:\n"
        "  DEAL_II_WITH_MPI = ${DEAL_II_WITH_MPI}\n"
        "  HDF5_WITH_MPI    = ${HDF5_WITH_MPI}\n"
        )
      SET(${var} FALSE)
    ENDIF()

    CHECK_MPI_INTERFACE(HDF5 ${var})
  ENDIF()
ENDMACRO()


CONFIGURE_FEATURE(HDF5)
