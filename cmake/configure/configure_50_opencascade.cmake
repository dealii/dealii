## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2012 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

#
# Configuration for the OpenCASCADE library:
#

macro(feature_opencascade_find_external var)
    find_package(DEAL_II_OPENCASCADE)

    if(OPENCASCADE_FOUND)
        set(${var} TRUE)

        #
        # We require at least OpenCASCADE 7.0.0
        #
        set(_version_required 7.0.0)
        if(OPENCASCADE_VERSION VERSION_LESS ${_version_required})
            message(STATUS "Could not find a sufficient OpenCASCADE installation: "
                    "deal.II requires at least version ${_version_required}, "
                    "but version ${OPENCASCADE_VERSION} was found."
            )
            set(OPENCASCADE_ADDITIONAL_ERROR_STRING
                    ${OPENCASCADE_ADDITIONAL_ERROR_STRING}
                    "The OpenCASCADE installation (found at \"${OPENCASCADE_DIR}\")\n"
                    "with version ${OPENCASCADE_VERSION} is too old.\n"
                    "deal.II requires at least version ${_version_required}.\n\n"
            )
            set(${var} FALSE)
        endif()

    endif()
endmacro()

configure_feature(OPENCASCADE)
