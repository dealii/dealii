## ---------------------------------------------------------------------
##
## Copyright (C) 2012 - 2021 by the deal.II authors
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
# Configuration for the p4est and sc libraries:
#

set(FEATURE_P4EST_DEPENDS MPI)


macro(feature_p4est_find_external var)
  find_package(DEAL_II_P4EST)

  if(P4EST_FOUND)
    set(${var} TRUE)

    #
    # We require at least version 2.0
    #
    set(_version_required 2.0)
    if(P4EST_VERSION VERSION_LESS ${_version_required})
      message(STATUS "Insufficient p4est installation found: "
        "At least version ${_version_required} is required."
        )
      set(P4EST_ADDITIONAL_ERROR_STRING
        "Insufficient p4est installation found!\n"
        "At least version ${_version_required} is required.\n"
        )
      set(${var} FALSE)
    endif()

    #
    # Check whether p4est supports mpi:
    #
    if(NOT P4EST_WITH_MPI)
      message(STATUS "Insufficient p4est installation found: "
        "p4est has to be configured with MPI enabled."
        )
      set(P4EST_ADDITIONAL_ERROR_STRING
        ${P4EST_ADDITIONAL_ERROR_STRING}
        "Insufficient p4est installation found!\n"
        "p4est has to be configured with MPI enabled.\n"
        )
      set(${var} FALSE)
    endif()

    #
    # Check whether p4est is built against zlib:
    #
    if(NOT P4EST_WITH_ZLIB)
      message(STATUS "Insufficient p4est installation found: "
        "p4est has to be configured with enabled zlib support."
        )
      set(P4EST_ADDITIONAL_ERROR_STRING
        ${P4EST_ADDITIONAL_ERROR_STRING}
        "Insufficient p4est installation found!\n"
        "p4est has to be configured with enabled zlib support.\n"
        )
      set(${var} FALSE)
    endif()

    check_mpi_interface(P4EST ${var})
  endif()
endmacro()

macro(feature_p4est_configure_external)
  set(DEAL_II_P4EST_WITH_VTK_BINARY ${P4EST_WITH_VTK_BINARY})
  set(DEAL_II_P4EST_WITH_SEARCH_LOCAL ${P4EST_WITH_SEARCH_LOCAL})
endmacro()


configure_feature(P4EST)
