## ---------------------------------------------------------------------
##
## Copyright (C) 2019 - 2020 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

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
