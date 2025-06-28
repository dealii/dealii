## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2017 - 2025 by the deal.II authors
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
# Configuration for the SUNDIALS library:
#

macro(feature_sundials_find_external var)
  find_package(DEAL_II_SUNDIALS)

  if(SUNDIALS_FOUND)
    set(${var} TRUE)

    #
    # We require at least sundials 5.4.0
    #
    set(_version_required 5.4.0)
    if(SUNDIALS_VERSION VERSION_LESS ${_version_required})
      message(STATUS "Could not find a sufficient Sundials installation: "
        "deal.II requires at least version ${_version_required}, "
        "but version ${SUNDIALS_VERSION} was found."
        )
      set(SUNDIALS_ADDITIONAL_ERROR_STRING
        ${SUNDIALS_ADDITIONAL_ERROR_STRING}
        "The SUNDIALS installation (found at \"${SUNDIALS_DIR}\")\n"
        "with version ${SUNDIALS_VERSION} is too old.\n"
        "deal.II requires at least version ${_version_required}.\n\n"
        )
      set(${var} FALSE)
    endif()

    #
    # Sundials has to be configured with the same MPI configuration as
    # deal.II.
    #
    if((SUNDIALS_WITH_MPI AND NOT DEAL_II_WITH_MPI) OR (NOT SUNDIALS_WITH_MPI AND DEAL_II_WITH_MPI))
      message(STATUS "Could not find a sufficient Sundials installation: "
        "deal.II and Sundials must have MPI support either both enabled or both disabled."
      )
      set(SUNDIALS_ADDITIONAL_ERROR_STRING
        ${SUNDIALS_ADDITIONAL_ERROR_STRING}
        "The Sundials installation (found at \"${SUNDIALS_DIR}\") has:\n"
        "  SUNDIALS_WITH_MPI = ${SUNDIALS_WITH_MPI}\n"
        " While deal.II has:\n"
        "  DEAL_II_WITH_MPI = ${DEAL_II_WITH_MPI}\n"
      )
      set(${var} FALSE)
    endif()
  endif()
endmacro()

macro(feature_sundials_configure_external)
  set(DEAL_II_SUNDIALS_WITH_IDAS ${SUNDIALS_WITH_IDAS})
endmacro()

configure_feature(SUNDIALS)
