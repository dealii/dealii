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
# Configuration for the PSBLAS library:
#

set(FEATURE_PSBLAS_DEPENDS MPI)

macro(feature_psblas_find_external var)
    find_package(DEAL_II_PSBLAS)

    if(PSBLAS_FOUND)
        set(${var} TRUE)

        set(_version_required 3.9.0)
        if(PSBLAS_VERSION VERSION_LESS ${_version_required})
            message(STATUS "Insufficient PSBLAS installation found: "
                "At least version ${_version_required} is required."
                )
            set(PSBLAS_ADDITIONAL_ERROR_STRING
                "Insufficient PSBLAS installation found!\n"
                "At least version ${_version_required} is required.\n"
                )
            set(${var} FALSE)

        endif()
    endif()
endmacro()
configure_feature(PSBLAS)
