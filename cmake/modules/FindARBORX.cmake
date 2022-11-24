## ---------------------------------------------------------------------
##
## Copyright (C) 2021 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#
# Try to find the ArborX library
#
# This module exports
#
#   ARBORX_INCLUDE_DIRS
#   ARBORX_INTERFACE_LINK_FLAGS
#

set(ARBORX_DIR "" CACHE PATH "An optional hint to an ArborX installation")
set_if_empty(ARBORX_DIR "$ENV{ARBORX_DIR}")


find_package(ArborX
  HINTS ${ARBORX_DIR} ${ArborX_DIR} $ENV{ArborX_DIR}
  )

if(ArborX_FOUND)
  get_property(ARBORX_INSTALL_INCLUDE_DIR TARGET ArborX::ArborX PROPERTY INTERFACE_INCLUDE_DIRECTORIES)

  #
  # Check whether ArborX was compiled with MPI support
  #
  message(STATUS
    "Checking whether the found ArborX has MPI support:"
    )

  #
  # Look for ArborX_Config.hpp - we'll query it to determine MPI support:
  #
  deal_ii_find_file(ARBORX_CONFIG_HPP ArborX_Config.hpp
    HINTS ${ARBORX_INSTALL_INCLUDE_DIR}
    NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH NO_CMAKE_FIND_ROOT_PATH
    )

  if(EXISTS ${ARBORX_CONFIG_HPP})
    #
    # Determine whether ArborX was configured with MPI:
    #
    file(STRINGS "${ARBORX_CONFIG_HPP}" ARBORX_MPI_STRING
      REGEX "#define ARBORX_ENABLE_MPI")
    if("${ARBORX_MPI_STRING}" STREQUAL "")
      message(STATUS "ArborX has no MPI support")
    else()
      set(DEAL_II_ARBORX_WITH_MPI TRUE)
      message(STATUS "ArborX has MPI support")
    endif()
  endif()
endif()

deal_ii_package_handle(ARBORX
  # ArborX is a header-only library
  INCLUDE_DIRS REQUIRED ARBORX_INSTALL_INCLUDE_DIR
  USER_INCLUDE_DIRS REQUIRED ARBORX_INSTALL_INCLUDE_DIR
  CLEAR ARBORX_DIR ArborX_DIR
  )
