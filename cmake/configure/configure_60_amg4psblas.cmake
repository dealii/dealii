## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2017 - 2022 by the deal.II authors
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

set(FEATURE_AMG4PSBLAS_DEPENDS PSBLAS)

macro(feature_amg4psblas_find_external var)
    find_package(DEAL_II_AMG4PSBLAS)

    if(AMG4PSBLAS_FOUND)
        set(${var} TRUE)

            set(_version_required 1.2.0)
            if(AMG4PSBLAS_VERSION VERSION_LESS ${_version_required})
                message(STATUS "Insufficient AMG4PSBLAS installation found: "
                "At least version ${_version_required} is required. "
                "Detected version: ${AMG4PSBLAS_VERSION}."
                )
                set(AMG4PSBLAS_ADDITIONAL_ERROR_STRING
                "Insufficient AMG4PSBLAS installation found!\n"
                "At least version ${_version_required} is required.\n"
                "Detected version: ${AMG4PSBLAS_VERSION}.\n"
                )
                set(${var} FALSE)
            endif()
    endif()
endmacro()

macro(feature_amg4psblas_configure_external)
  set(DEAL_II_EXPAND_AMG4PSBLAS_PRECONDITIONER "PSCToolkitWrappers::PreconditionAMG")
endmacro()

configure_feature(AMG4PSBLAS)