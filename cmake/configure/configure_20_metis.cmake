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
