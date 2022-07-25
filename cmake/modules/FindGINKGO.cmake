## ---------------------------------------------------------------------
##
## Copyright (C) 2018 - 2021 by the deal.II authors
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
#   GINKGO_VERSION
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

#
# We'd like to have the full library names but the Ginkgo package only
# exports a list with short names. So check again for every lib and store
# the full path:
#
SET(_libraries "")
FOREACH(_library ginkgo ${GINKGO_INTERFACE_LINK_LIBRARIES})
  LIST(APPEND _libraries GINKGO_LIBRARY_${_library})
  DEAL_II_FIND_LIBRARY(GINKGO_LIBRARY_${_library}
    NAMES ${_library}
    HINTS ${GINKGO_INSTALL_LIBRARY_DIR}
    NO_DEFAULT_PATH
    NO_CMAKE_ENVIRONMENT_PATH
    NO_CMAKE_PATH
    NO_SYSTEM_ENVIRONMENT_PATH
    NO_CMAKE_SYSTEM_PATH
    NO_CMAKE_FIND_ROOT_PATH
    )
ENDFOREACH()

#
# Get ginkgo version number
#
IF(Ginkgo_FOUND)
  SET(GINKGO_VERSION "${GINKGO_PROJECT_VERSION}")
ENDIF()

DEAL_II_PACKAGE_HANDLE(GINKGO
  LIBRARIES REQUIRED ${_libraries}
  INCLUDE_DIRS REQUIRED GINKGO_INSTALL_INCLUDE_DIR
  USER_INCLUDE_DIRS REQUIRED GINKGO_INSTALL_INCLUDE_DIR
  CLEAR Ginkgo_DIR ${_libraries}
  )
