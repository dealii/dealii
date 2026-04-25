## -----------------------------------------------------------------------------
##
## SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
## Copyright (C) 2017 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Detailed license information governing the source code and contributions
## can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
##
## -----------------------------------------------------------------------------

#
# Configuration for the PSBLAS library:
#

set(FEATURE_PSBLAS_DEPENDS MPI LAPACK)

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

macro(feature_psblas_configure_external)
  set(DEAL_II_EXPAND_PSBLAS_VECTOR "PSCToolkitWrappers::Vector")
  set(DEAL_II_EXPAND_PSBLAS_SPARSE_MATRICES "PSCToolkitWrappers::SparseMatrix")
  set(DEAL_II_EXPAND_PSBLAS_SPARSITY_PATTERN "PSCToolkitWrappers::SparsityPattern")
endmacro()

configure_feature(PSBLAS)
