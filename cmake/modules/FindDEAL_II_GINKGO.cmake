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
# Try to find the GINKGO library
#
# This module exports
#
#   GINKGO_INCLUDE_DIRS
#   GINKGO_INTERFACE_LINK_FLAGS
#   GINKGO_VERSION
#

set(GINKGO_DIR "" CACHE PATH "An optional hint to a GINKGO installation")
set_if_empty(GINKGO_DIR "$ENV{GINKGO_DIR}")

#
# Save and restore the ${CMAKE_MODULE_PATH} variable. The Ginkgo project
# configuration unfortunately overrides the variable which causes
# subsequent configuration to fail.
#
set(_cmake_module_path ${CMAKE_MODULE_PATH})
find_package(Ginkgo QUIET
  HINTS ${GINKGO_DIR} ${Ginkgo_DIR} $ENV{Ginkgo_DIR}
  )
set(CMAKE_MODULE_PATH ${_cmake_module_path})

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
set(_libraries "")
foreach(_library ginkgo ${GINKGO_INTERFACE_LINK_LIBRARIES})
  # Make sure to only pick up Ginkgo's own libraries here, skipping
  # the MPI libraries that are listed here as of Ginkgo 1.5.0.
  if(_library MATCHES "ginkgo.*")
    list(APPEND _libraries GINKGO_LIBRARY_${_library})
    deal_ii_find_library(GINKGO_LIBRARY_${_library}
      NAMES ${_library}
      HINTS ${GINKGO_INSTALL_LIBRARY_DIR}
      NO_DEFAULT_PATH
      NO_CMAKE_ENVIRONMENT_PATH
      NO_CMAKE_PATH
      NO_SYSTEM_ENVIRONMENT_PATH
      NO_CMAKE_SYSTEM_PATH
      NO_CMAKE_FIND_ROOT_PATH
      )
  endif()
endforeach()

#
# Get ginkgo version number
#
if(Ginkgo_FOUND)
  set(GINKGO_VERSION "${GINKGO_PROJECT_VERSION}")
endif()

process_feature(GINKGO
  LIBRARIES REQUIRED ${_libraries}
  INCLUDE_DIRS REQUIRED GINKGO_INSTALL_INCLUDE_DIR
  CLEAR Ginkgo_DIR ${_libraries}
  )
