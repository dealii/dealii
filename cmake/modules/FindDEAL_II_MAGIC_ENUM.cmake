## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2024 - 2025 by the deal.II authors
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
# Try to find the MAGIC_ENUM library
#
# This module exports
#
#   MAGIC_ENUM_FOUND
#   MAGIC_ENUM_INCLUDE_DIRS
#

set(MAGIC_ENUM_DIR "" CACHE PATH "An optional hint to a MAGIC_ENUM installation")
set_if_empty(MAGIC_ENUM_DIR "$ENV{MAGIC_ENUM_DIR}")

deal_ii_find_path(MAGIC_ENUM_INCLUDE_DIR magic_enum.hpp
  HINTS ${MAGIC_ENUM_DIR}
  PATH_SUFFIXES include include/magic_enum
  )

process_feature(MAGIC_ENUM
  INCLUDE_DIRS
    REQUIRED MAGIC_ENUM_INCLUDE_DIR
  CLEAR MAGIC_ENUM_INCLUDE_DIR
  )
