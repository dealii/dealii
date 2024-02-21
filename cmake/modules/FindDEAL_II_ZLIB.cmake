## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2014 - 2022 by the deal.II authors
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
# Try to find the ZLIB library
#
# This module exports
#
#   ZLIB_LIBRARIES
#   ZLIB_INCLUDE_DIRS
#   ZLIB_VERSION
#

set(ZLIB_DIR "" CACHE PATH "An optional hint to a ZLIB installation")
set_if_empty(ZLIB_DIR "$ENV{ZLIB_DIR}")

if(NOT "${ZLIB_DIR}" STREQUAL "")
  set(ZLIB_ROOT ${ZLIB_DIR})
endif()
find_package(ZLIB)

set(ZLIB_VERSION ${ZLIB_VERSION_STRING})

process_feature(ZLIB
  LIBRARIES REQUIRED ZLIB_LIBRARY
  INCLUDE_DIRS REQUIRED ZLIB_INCLUDE_DIR
  CLEAR ZLIB_INCLUDE_DIR ZLIB_LIBRARY
  )
