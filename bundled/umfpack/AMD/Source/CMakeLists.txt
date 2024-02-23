## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2023 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

strip_flag(DEAL_II_WARNING_FLAGS "-Wfloat-conversion")

set(src_amd
  amd_1.cc
  amd_2.cc
  amd_aat.cc
  amd_control.cc
  amd_defaults.cc
  amd_dump.cc
  amd_info.cc
  amd_order.cc
  amd_postorder.cc
  amd_post_tree.cc
  amd_preprocess.cc
  amd_valid.cc
  )

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../Include)
define_object_library(bundled_amd_int OBJECT ${src_amd})
deal_ii_add_definitions(bundled_amd_int "DINT")

define_object_library(bundled_amd_long OBJECT ${src_amd})
deal_ii_add_definitions(bundled_amd_long "DLONG")

define_object_library(bundled_amd_global OBJECT "amd_global.cc")

