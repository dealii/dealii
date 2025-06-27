## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2021 - 2025 by the deal.II authors
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
# Configuration for Kokkos support in deal.II:
#

set(DEAL_II_WITH_KOKKOS ON # Always true. We need it :-]
  CACHE BOOL "Build deal.II with support for Kokkos." FORCE
  )


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

#
# Sanity check: Make sure we do not accidentally end up using bundled
# Kokkos with Trilinos and PETSc built against an external Kokkos library.
# We should have never gotten into this situation. So all we can do now is
# to issue a FATAL_ERROR.
#

if(DEAL_II_FEATURE_KOKKOS_BUNDLED_CONFIGURED AND DEAL_II_WITH_TRILINOS AND TRILINOS_WITH_KOKKOS)
  message(FATAL_ERROR "\n"
    "Internal build system error: We have selected Trilinos shipping with (or "
    "built against) an external Kokkos library, but ended up selecting our "
    "bundled Kokkos library.\n\n"
    )
endif()

if(DEAL_II_FEATURE_KOKKOS_BUNDLED_CONFIGURED AND DEAL_II_WITH_PETSC AND PETSC_WITH_KOKKOS)
  message(FATAL_ERROR "\n"
    "Internal build system error: We have selected PETSc built against an "
    "external Kokkos library, but ended up selecting our bundled Kokkos "
    "library.\n\n"
    )
endif()
