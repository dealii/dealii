## ---------------------------------------------------------------------
##
## Copyright (C) 2021 by the deal.II authors
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
# Configuration for ArborX support in deal.II:
#

set(FEATURE_ARBORX_AFTER MPI)
set(FEATURE_ARBORX_DEPENDS KOKKOS)

macro(feature_arborx_find_external var)
  find_package(DEAL_II_ARBORX)

  if(ARBORX_FOUND)
    #
    # So, we have a library. Let's see whether we can use it:
    #
    set(${var} TRUE)

    #
    # ArborX has to be configured with the same MPI configuration as
    # deal.II.
    #
    if((NOT ARBORX_WITH_MPI AND DEAL_II_WITH_MPI) OR (ARBORX_WITH_MPI AND NOT DEAL_II_WITH_MPI))
      message(STATUS "Could not find a sufficient ArborX installation: "
        "ArborX has to be configured with the same MPI configuration as deal.II."
        )
      set(ARBORX_ADDITIONAL_ERROR_STRING
        ${ARBORX_ADDITIONAL_ERROR_STRING}
        "Could not find a sufficient ArborX installation:\n"
        "ArborX has to be configured with the same MPI configuration as deal.II, but found:\n"
        "  DEAL_II_WITH_MPI = ${DEAL_II_WITH_MPI}\n"
        "  ARBORX_WITH_MPI  = ${ARBORX_WITH_MPI}\n"
        )
      set(${var} FALSE)
    endif()
  endif()

  set(DEAL_II_ARBORX_WITH_MPI ${ARBORX_WITH_MPI})
endmacro()

configure_feature(ARBORX)

