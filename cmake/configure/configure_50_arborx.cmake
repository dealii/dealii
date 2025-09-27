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

    list(APPEND CMAKE_REQUIRED_LIBRARIES
      ArborX::ArborX
    )

    if(ArborX_VERSION VERSION_LESS 2.0.0)
      check_cxx_compiler_bug(
        "
        #include <ArborX.hpp>
        int main() {
          Kokkos::View<ArborX::Point*, Kokkos::HostSpace> points(\"points\", 0);
          [[maybe_unused]] ArborX::BVH<Kokkos::HostSpace> bvh(Kokkos::DefaultHostExecutionSpace{}, points);
        }
        "
        DEAL_II_ARBORX_CXX20_BUG)
      reset_cmake_required()

      if(DEAL_II_ARBORX_CXX20_BUG)
        message(STATUS "Could not find a sufficient ArborX installation: "
          "The ArborX version doesn't work with C++20 or higher."
          )
        set(ARBORX_ADDITIONAL_ERROR_STRING
          ${ARBORX_ADDITIONAL_ERROR_STRING}
          "Could not find a sufficient ArborX installation:\n"
          "The ArborX version doesn't work with C++20 or higher. "
          "Try using a later ArborX release or try specifying a lower C++ standard.\n"
          )
        set(${var} FALSE)
      endif()
    endif()

    if(ArborX_VERSION VERSION_GREATER_EQUAL 2.0.0)
      if(NOT DEAL_II_HAVE_CXX20)
        message(STATUS "Could not find a sufficient ArborX installation: "
          "The ArborX version ${ArborX_VERSION} requires C++20 or higher."
          )
        set(ARBORX_ADDITIONAL_ERROR_STRING
          ${ARBORX_ADDITIONAL_ERROR_STRING}
          "Could not find a sufficient ArborX installation:\n"
          "The ArborX version ${ArborX_VERSION} requires C++20 or higher. "
          "Try using an earlier ArborX release or try specifying a higher C++ standard.\n"
          )
        set(${var} FALSE)
      endif()
    endif()
  endif()

  set(DEAL_II_ARBORX_WITH_MPI ${ARBORX_WITH_MPI})
endmacro()

configure_feature(ARBORX)
