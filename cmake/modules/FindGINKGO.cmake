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

find_package(Ginkgo
  HINTS ${GINKGO_DIR} ${Ginkgo_DIR} $ENV{Ginkgo_DIR}
  )

DEAL_II_PACKAGE_HANDLE(GINKGO
  LIBRARIES
    REQUIRED GINKGO_INTERFACE_LINK_FLAGS
  INCLUDE_DIRS
    REQUIRED GINKGO_INSTALL_INCLUDE_DIR
  USER_INCLUDE_DIRS
    REQUIRED GINKGO_INSTALL_INCLUDE_DIR
  CLEAR
    GINKGO_INSTALL_INCLUDE_DIR GINKGO_INTERFACE_LINK_FLAGS
  )
