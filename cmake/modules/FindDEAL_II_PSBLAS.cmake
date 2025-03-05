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
# Try to find the PSBLAS library
#
# This module exports
#
#   PSBLAS_LIBRARY
#   PSBLAS_INCLUDE_DIR
#

set(PSBLAS_DIR "" CACHE PATH "An optional hint to a PSBLAS installation containing the PSBLAS include directory and libraries")
set_if_empty(PSBLAS_DIR "$ENV{PSBLAS_DIR}")

deal_ii_find_library(PSBLAS_LIBRARY
  NAMES psb_base
  HINTS ${PSBLAS_DIR}
  PATH_SUFFIXES lib${LIB_SUFFIX} lib64 lib
  )

deal_ii_find_path(PSBLAS_INCLUDE_DIR psb_c_base.h
  HINTS ${PSBLAS_DIR}
  PATH_SUFFIXES include
  )

process_feature(PSBLAS
  LIBRARIES 
    REQUIRED PSBLAS_LIBRARY
  INCLUDE_DIRS 
    REQUIRED PSBLAS_INCLUDE_DIR
  CLEAR
    PSBLAS_LIBRARY PSBLAS_INCLUDE_DIR
  )
