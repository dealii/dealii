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

SET(ARBORX_DIR "" CACHE PATH "An optional hint to an ArborX installation")
SET_IF_EMPTY(ARBORX_DIR "$ENV{ARBORX_DIR}")


FIND_PACKAGE(ArborX
  HINTS ${ARBORX_DIR} ${ArborX_DIR} $ENV{ArborX_DIR}
  )

IF(ArborX_FOUND)
  GET_PROPERTY(ARBORX_INSTALL_INCLUDE_DIR TARGET ArborX::ArborX PROPERTY INTERFACE_INCLUDE_DIRECTORIES)

  #
  # Check whether ArborX was compiled with MPI support
  #
  MESSAGE(STATUS
    "Checking whether the found ArborX has MPI support:"
    )

  #
  # Look for ArborX_Config.hpp - we'll query it to determine MPI support:
  #
  DEAL_II_FIND_FILE(ARBORX_CONFIG_HPP ArborX_Config.hpp
    HINTS ${ARBORX_INSTALL_INCLUDE_DIR}
    NO_DEFAULT_PATH NO_CMAKE_ENVIRONMENT_PATH NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH NO_CMAKE_FIND_ROOT_PATH
    )

  IF(EXISTS ${ARBORX_CONFIG_HPP})
    #
    # Determine whether ArborX was configured with MPI:
    #
    FILE(STRINGS "${ARBORX_CONFIG_HPP}" ARBORX_MPI_STRING
      REGEX "#define ARBORX_ENABLE_MPI")
    IF("${ARBORX_MPI_STRING}" STREQUAL "")
      MESSAGE(STATUS "ArborX has no MPI support")
    ELSE()
      SET(DEAL_II_ARBORX_WITH_MPI TRUE)
      MESSAGE(STATUS "ArborX has MPI support")
    ENDIF()
  ENDIF()
ENDIF()

DEAL_II_PACKAGE_HANDLE(ARBORX
  # ArborX is a header-only library
  INCLUDE_DIRS REQUIRED ARBORX_INSTALL_INCLUDE_DIR
  USER_INCLUDE_DIRS REQUIRED ARBORX_INSTALL_INCLUDE_DIR
  CLEAR ARBORX_DIR ArborX_DIR
  )
