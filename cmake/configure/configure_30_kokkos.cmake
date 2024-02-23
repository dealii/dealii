## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2021 - 2022 by the deal.II authors
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

