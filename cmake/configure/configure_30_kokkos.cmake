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

set(DEAL_II_WITH_KOKKOS ON # Always true. We need it :-]
  CACHE BOOL "Build deal.II with support for Kokkos." FORCE
  )

macro(feature_kokkos_configure_bundled)
  set(KOKKOS_BUNDLED_INCLUDE_DIRS
    ${KOKKOS_FOLDER}/algorithms/src
    ${KOKKOS_FOLDER}/containers/src
    ${KOKKOS_FOLDER}/core/src
    ${KOKKOS_FOLDER}/simd/src
    ${KOKKOS_FOLDER}/tpls/desul/include
    )

  enable_if_supported(DEAL_II_CXX_FLAGS -Wno-float-conversion)
  enable_if_supported(DEAL_II_CXX_FLAGS -Wno-missing-field-initializers)
  enable_if_supported(DEAL_II_CXX_FLAGS -Wno-suggest-override)
  enable_if_supported(DEAL_II_CXX_FLAGS -Wno-unused-but-set-parameter)
endmacro()

configure_feature(KOKKOS)

#
# DEAL_II_WITH_KOKKOS is always required.
#
if(NOT DEAL_II_WITH_KOKKOS)
  if(DEAL_II_FEATURE_AUTODETECTION)
    feature_error_message("KOKKOS")
  else()
    message(FATAL_ERROR "\n"
      "Unmet configuration requirements: "
      "DEAL_II_WITH_KOKKOS required, but set to OFF!\n\n"
      )
  endif()
endif()
