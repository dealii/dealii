## ---------------------------------------------------------------------
##
## Copyright (C) 2018 - 2020 by the deal.II authors
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
# Try to find the GINKGO library
#
# This module exports
#
#   GINKGO_INCLUDE_DIRS
#   GINKGO_INTERFACE_LINK_FLAGS
#

SET(GINKGO_DIR "" CACHE PATH "An optional hint to a GINKGO installation")
SET_IF_EMPTY(GINKGO_DIR "$ENV{GINKGO_DIR}")

#
# Save and restore the ${CMAKE_MODULE_PATH} variable. The Ginkgo project
# configuration unfortunately overrides the variable which causes
# subsequent configuration to fail.
#
SET(_cmake_module_path ${CMAKE_MODULE_PATH})
FIND_PACKAGE(Ginkgo
  HINTS ${GINKGO_DIR} ${Ginkgo_DIR} $ENV{Ginkgo_DIR}
  )
SET(CMAKE_MODULE_PATH ${_cmake_module_path})

#
# Cosmetic clean up: Let's remove all variables beginning with "GINKGO_"
# that are actually not used during configuration but show up in
# detailed.log
#
unset(GINKGO_CXX_COMPILER)

DEAL_II_PACKAGE_HANDLE(GINKGO
  LIBRARIES
    REQUIRED GINKGO_INTERFACE_LINK_LIBRARIES
  INCLUDE_DIRS
    REQUIRED GINKGO_INSTALL_INCLUDE_DIR
  USER_INCLUDE_DIRS
    REQUIRED GINKGO_INSTALL_INCLUDE_DIR
  CLEAR
    Ginkgo_DIR
  )
