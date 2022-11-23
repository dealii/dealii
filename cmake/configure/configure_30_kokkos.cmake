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
# Configuration for Kokkos support in deal.II:
#

SET(DEAL_II_WITH_KOKKOS ON # Always true. We need it :-]
  CACHE BOOL "Build deal.II with support for Kokkos." FORCE
  )

MACRO(FEATURE_KOKKOS_CONFIGURE_BUNDLED)
  SET(KOKKOS_BUNDLED_INCLUDE_DIRS
    ${KOKKOS_FOLDER}/algorithms/src
    ${KOKKOS_FOLDER}/containers/src
    ${KOKKOS_FOLDER}/core/src
    ${KOKKOS_FOLDER}/simd/src
    ${KOKKOS_FOLDER}/tpls/desul/include
    )
ENDMACRO()

CONFIGURE_FEATURE(KOKKOS)

#
# DEAL_II_WITH_KOKKOS is always required.
#
IF(NOT DEAL_II_WITH_KOKKOS)
  IF(DEAL_II_FEATURE_AUTODETECTION)
    FEATURE_ERROR_MESSAGE("KOKKOS")
  ELSE()
    MESSAGE(FATAL_ERROR "\n"
      "Unmet configuration requirements: "
      "DEAL_II_WITH_KOKKOS required, but set to OFF!\n\n"
      )
  ENDIF()
ENDIF()
