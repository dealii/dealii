## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2019 - 2022 by the deal.II authors
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
# Configuration for the SymEngine library:
#


macro(feature_symengine_find_external var)
  find_package(DEAL_II_SYMENGINE)

  if(SYMENGINE_FOUND)
    set(${var} TRUE)

    #
    # We require at least version 0.6 of the symengine library:
    #
    set(_version_required "0.6")

    if(SYMENGINE_VERSION VERSION_LESS ${_version_required})
      message(STATUS "Insufficient SymEngine installation found: "
              "At least version ${_version_required} is required "
              "but version ${SYMENGINE_VERSION} was found."
             )
      set(SYMENGINE_ADDITIONAL_ERROR_STRING
          "Insufficient SymEngine installation found!\n"
          "At least version ${_version_required} is required "
          "but version ${SYMENGINE_VERSION} was found.\n"
         )
      set(${var} FALSE)
    endif()
  endif()
endmacro()

macro(feature_symengine_configure_external)
  set(DEAL_II_SYMENGINE_WITH_LLVM ${SYMENGINE_WITH_LLVM})

  if(DEAL_II_SYMENGINE_WITH_LLVM)
    message(STATUS "Configured with SymEngine LLVM capabilities.")
  endif()

  #
  # Overwrite the compiler flags imported from SymEngine
  #
  set(SYMENGINE_CXX_FLAGS)
  set(SYMENGINE_CXX_FLAGS_DEBUG)
  set(SYMENGINE_CXX_FLAGS_RELEASE)
endmacro()


configure_feature(SYMENGINE)
