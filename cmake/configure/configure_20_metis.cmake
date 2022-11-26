## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2017 by the deal.II authors
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

set(FEATURE_METIS_AFTER MPI)

#
# Configuration for the metis library:
#

macro(feature_metis_find_external var)
  find_package(DEAL_II_METIS)

  if(METIS_FOUND)
    set(${var} TRUE)

    if(NOT METIS_VERSION_MAJOR GREATER 4)
      message(STATUS "Insufficient metis installation found: "
        "Version 5.x required!"
        )
      set(METIS_ADDITIONAL_ERROR_STRING
        "Could not find a sufficiently modern metis installation: "
        "Version 5.x required!\n"
        )
      set(${var} FALSE)
    endif()

    check_mpi_interface(METIS ${var})
  endif()
endmacro()

configure_feature(METIS)
