## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2017 - 2022 by the deal.II authors
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
# Try to find the ASSIMP libraries
#
# This module exports
#
#   ASSIMP_LIBRARIES
#   ASSIMP_INCLUDE_DIRS
#

set(ASSIMP_DIR "" CACHE PATH "An optional hint to a Assimp installation")
set_if_empty(ASSIMP_DIR "$ENV{ASSIMP_DIR}")

deal_ii_find_library(ASSIMP_LIB NAMES assimp
  HINTS ${ASSIMP_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

deal_ii_find_path(ASSIMP_INC assimp/defs.h
  HINTS ${ASSIMP_DIR}
  PATH_SUFFIXES include
  )

process_feature(ASSIMP
  LIBRARIES REQUIRED ASSIMP_LIB
  INCLUDE_DIRS REQUIRED ASSIMP_INC
  CLEAR ASSIMP_LIB ASSIMP_INC
  )
