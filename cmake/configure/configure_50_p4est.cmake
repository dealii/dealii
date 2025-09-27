## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2025 by the deal.II authors
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
# Configuration for the p4est and sc libraries:
#

set(FEATURE_P4EST_DEPENDS MPI)


macro(feature_p4est_find_external var)
  find_package(DEAL_II_P4EST)

  if(P4EST_FOUND)
    set(${var} TRUE)

    #
    # We require at least version 2.2
    #
    set(_version_required 2.2)
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
endmacro()


configure_feature(P4EST)
